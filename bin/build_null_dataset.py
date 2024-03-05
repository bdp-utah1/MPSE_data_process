#!/usr/bin/env python3

from os import mkdir
import os.path as path
import sys
import argparse
import random
import json
import csv 
import re

from pyhpo.ontology import Ontology
from pyhpo.set import BasicHPOSet

from collections import Counter, defaultdict
from joblib import dump, load

import numpy as np
import pandas as pd

csv.field_size_limit(sys.maxsize)


def argue():
    """Parses command line arguments.

    Returns:
        ArgumentParser: An ArgumentParser object that can be used to define and parse command line arguments.
    """
    parser = argparse.ArgumentParser("""
    This script takes an MPSE input dataset and permutes the HPO terms for each sample in the cohort.
    This is done to calculate the ability of MPSE to classify cases/controls using randomized junk data.
    There are two ways to permute: 1) sample from the entire HPO or 2) sample from the original cohort's represented terms (a subset of the HPO).
    The desired outcome is poor performance with permuted data.
    Note: the number of terms in each sample is maintained, only the term identities are permuted.
    """)
    parser.add_argument("-t", "--training",
            default="data/test/fake_training_data.tsv", 
            help="Case/control training data in standard format.")
    parser.add_argument("-p", "--prospective",
            required=False,
            help="Prospective data in standard format.")
    return parser


def ready(ftag, delim="\t", drop_header=False):
    """Reads a delimited file and returns its contents as a list of rows.

    Args:
        ftag (str): The path or filename of the delimited file to be read.
        delim (str, optional): The delimiter used in the delimited file. Defaults to "\t".
        drop_header (bool, optional): Whether to drop the header row from the delimited file. Defaults to False.

    Returns:
        list: A list containing the rows of the delimited file, where each row is a list of its fields.
    """
    with open(ftag, newline="") as f:
        reader = csv.reader(f, delimiter=delim)
        if drop_header:
            out = [row for row in reader][1:]
        else:
            out = [row for row in reader]
    return out 


def get_column_positions(data, col_names):
    """Returns a dictionary mapping column names to their corresponding indices in the data.

    Args:
        data (list): A list of lists representing tabular data.
        col_names (list): A list of column names.

    Returns:
        dict: A dictionary mapping column names to their corresponding indices in the data.
    """
    idx_dic = {name: data[0].index(name) for name in col_names}
    return idx_dic


def parse_hpo(hpo_str):
    """Parses and formats an HPO code from a string.

    Args:
        hpo_str (str): A string containing an HPO code with variable formatting.

    Returns:
        str: The parsed and formatted HPO code in the format 'HP:#######'.
    """
    pattern = r'(?i)([hp]{2}):?(\d{7})'
    replacement = r'HP:\2'
    return re.sub(pattern, replacement, hpo_str)


def remove_parent_terms(hpo_lst):
    """Retrieves the most specific HPO term of each subtree.

    Args:
        hpo_lst (list): A list of HPO terms.

    Returns:
        list: A sorted list of child terms' IDs.
    """
    # BasicHPOSet does the following:
    # -removes parent terms
    # -removes modifier terms
    # -replaces obsolete terms
    hpo_set = BasicHPOSet.from_queries(hpo_lst)
    #hpo_subset = hpo_set.child_nodes()
    hpo_dic = hpo_set.toJSON()
    hpo_str = sorted([x["id"] for x in hpo_dic])
    return hpo_str


def clean_codes(codes, keep_others=False):
    """Basic pre-processing of input codes in preparation for modeling.

    Args:
        codes (iterable): An iterable of codes.
        keep_others (bool, optional): Whether to keep non-HPO/ICD-10-CM codes. Defaults to False.

    Returns:
        list: A filtered list of codes.
    """
    valid_hpo = Ontology.to_dataframe().index.tolist()
    hpo = []
    icd= []
    other = []
    for cde in codes:
        if cde in valid_hpo:
            hpo.append(cde)
        elif cm.is_valid_item(cde):
            icd.append(cm.add_dot(cde))
        else:
            other.append(cde)
    hpo_clean = remove_parent_terms(hpo)
    icd_clean = sorted(icd)
    other_clean = sorted(other) if keep_others else []
    return hpo_clean + icd_clean + other_clean


def make_compliant(data, dataset_name, col_idx, check_cols=None, keep_all_codes=False):
    """Performs data compliance checks and modifications on the provided dataset.

    Args:
        data (list): A list of lists representing tabular data.
        dataset_name (str): The name of the dataset.
        col_idx (dict): A dictionary mapping column names to their corresponding indices.
        check_cols (list, optional): A list of column names to perform value set checks on. Defaults to None.
        keep_all_codes (bool, optional): Whether to keep all codes or remove unrecognized codes. Defaults to False.

    Returns:
        list: The modified dataset after performing compliance checks and modifications.
    """
    # remove rows with no codes
    data = [row for row in data if row[col_idx["codes"]] != ""]

    # check identifiers are unique
    ids = [x[col_idx["pid"]] for x in data[1:]]
    if len(ids) != len(set(ids)):
        msg = "Warning: the dataset '{0}' PIDs are not unique. Please check this is expected..."
        print(msg.format(dataset_name), file=sys.stderr)

    # check value sets for seq_status, diagnostic, incidental
    # fill "" with 0
    if check_cols:
        for row in data[1:]:
            for col in check_cols:
                if row[col_idx[col]] == "":
                    row[col_idx[col]] = "0"
                elif row[col_idx[col]] not in ["0","1"]:
                    msg = "{0}: non-compliant value for column '{1}'; must be '0' or '1'"
                    raise ValueError(msg.format(dataset_name, col))
                    #sys.exit(msg.format(dataset_name, col))
    
    # order code list and remove duplicate codes
    # clean codes
    data[0].append("codes_clean")
    for row in data[1:]:
        dirty = sorted(set(row[col_idx["codes"]].split(";")))
        clean = clean_codes(dirty, keep_others=keep_all_codes)
        row[col_idx["codes"]] = ";".join(dirty)
        row.append(";".join(clean))
    return data


def build_null_dataset(data, col_idx, keep_terms, weights):
    term_cnts = [len(x[col_idx["codes_clean"]].split(";")) for x in data[1:]]
    #null_samp = [np.random.choice(keep_terms, size=n, replace=False, p=weights) for n in term_cnts]
    null_samp = [np.random.choice(keep_terms, size=n, replace=False) for n in term_cnts]
    null_samp = [["codes_permuted"]] + [[";".join(x)] for x in null_samp]
    null_out = [x+y for x,y in zip(data, null_samp)]
    return null_out


def main():
    parser = argue()
    args = parser.parse_args()
    _ = Ontology()

    train = ready(args.training)
    col_pos_names = ["pid","seq_status","diagnostic","codes"]
    train_col_idx = get_column_positions(train, col_pos_names)
    train = make_compliant(train, 
            "train_data", 
            train_col_idx, 
            ["seq_status","diagnostic"],
            False)
    train_col_idx = get_column_positions(train, col_pos_names + ["codes_clean"])

    term_col = [x[-1] for x in train[1:]]
    term_lst = ";".join(term_col).split(";")
    term_cnt = len(term_lst)
    term_dic = Counter(term_lst)
    a = term_dic.keys()
    p = [x/term_cnt for x in term_dic.values()]

    null = build_null_dataset(train, train_col_idx, keep_terms=list(a), weights=list(p))
    null_col_idx = get_column_positions(null, ["pid","seq_status","diagnostic","codes_permuted"])

    print("\t".join(["pid","seq_status","diagnostic","codes"]))
    for row in null[1:]:
        out_row = [row[null_col_idx["pid"]],
                   row[null_col_idx["seq_status"]],
                   row[null_col_idx["diagnostic"]],
                   row[null_col_idx["codes_permuted"]]]
        print("\t".join(out_row))



if __name__ == "__main__":
    main()
