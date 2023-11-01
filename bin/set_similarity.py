#!/usr/bin/env python3

import os
import sys
import csv
import random
import argparse
from pyhpo.ontology import Ontology
from pyhpo.set import BasicHPOSet


def argue():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--datasets", nargs=2, required=True,
                        help="Two datasets for computing Jaccard similarity.")
    args = parser.parse_args()
    return args


def ready(ftag, delim="\t", drop_header=False):
    with open(ftag, newline="") as f:
        reader = csv.reader(f, delimiter=delim)
        if drop_header:
            out = [row for row in reader][1:]
        else:
            out = [row for row in reader]
    return out 


def get_column_positions(data, col_names):
    idx_dic = {name: data[0].index(name) for name in col_names}
    return idx_dic


def make_id_lookup(data, col_pos):
    lookup = {}
    for idx,row in enumerate(data):
        lookup[row[col_pos]] = idx
    return lookup


def add_parent_terms(s):
    s_parents = set()
    for term in s:
        for parent in term.parents:
            s_parents.add(parent.id)
    return s_parents


def find_parents_recursive(node, visited):
    if node in visited:
        return set()

    parents = set()
    for parent in Ontology.get_hpo_object(node).parents:
        parents.add(parent.id)
        if parent.id != "HP:0000001":
            parents.update(find_parents_recursive(parent.id, visited))
    return parents


def get_random_hpo_sample(cnt, opts):
    samp = random.sample(opts, cnt)
    return samp


def calculate_jaccard(list1, list2):
    s1 = set(list1)
    s2 = set(list2)
    intersection = s1.intersection(s2)
    union = s1.union(s2)
    return len(intersection) / len(union)


def main():
    args = argue()
    _ = Ontology()

    df1_name = os.path.basename(args.datasets[0])
    df2_name = os.path.basename(args.datasets[1])

    df1 = ready(args.datasets[0])
    df2 = ready(args.datasets[1])

    df1_col_pos = get_column_positions(df1, ["pid", "codes"])
    df2_col_pos = get_column_positions(df2, ["pid", "codes"])

    df1_lookup = make_id_lookup(df1, df1_col_pos["pid"])
    df2_lookup = make_id_lookup(df2, df2_col_pos["pid"])

    id1 = list(df1_lookup.keys())
    id2 = list(df2_lookup.keys())
    ids = set(id1).intersection(set(id2))

    valid_hpo = _.to_dataframe().index.tolist()

    rows = []
    for id in ids:
        if id != "pid":
            raw1 = df1[df1_lookup[id]][df1_col_pos["codes"]].split(";")
            raw2 = df2[df2_lookup[id]][df2_col_pos["codes"]].split(";")

            real1 = [x for x in raw1 if x in valid_hpo]
            real2 = [x for x in raw2 if x in valid_hpo]

            set1 = BasicHPOSet.from_queries(real1)
            set2 = BasicHPOSet.from_queries(real2)

            lst1 = [x["id"] for x in set1.toJSON()]
            lst2 = [x["id"] for x in set2.toJSON()]

            len1 = len(lst1)
            len2 = len(lst2)

            set1_p = set()
            visited1 = set()
            for term in lst1:
                set1_p.add(term)
                set1_p.update(find_parents_recursive(term, visited1))
                visited1.update(set1_p)

            set2_p = set()
            visited2 = set()
            for term in lst2:
                set2_p.add(term)
                set2_p.update(find_parents_recursive(term, visited2))


            rand_samp1 = get_random_hpo_sample(len1, valid_hpo)
            rand_samp2 = get_random_hpo_sample(len2, valid_hpo)

            rand_set1 = BasicHPOSet.from_queries(rand_samp1)
            rand_set2 = BasicHPOSet.from_queries(rand_samp2)

            rand_samp1_p = set()
            rand_visited1 = set()
            for term in rand_samp1:
                rand_samp1_p.add(term)
                rand_samp1_p.update(find_parents_recursive(term, rand_visited1))

            rand_samp2_p = set()
            rand_visited2 = set()
            for term in rand_samp2:
                rand_samp2_p.add(term)
                rand_samp2_p.update(find_parents_recursive(term, rand_visited2))

            jac = calculate_jaccard(lst1, lst2)
            jac_p = calculate_jaccard(set1_p, set2_p)
            sim = set1.similarity(set2)

            rand_jac = calculate_jaccard(rand_samp1, rand_samp2)
            rand_jac_p = calculate_jaccard(rand_samp1_p, rand_samp2_p)
            rand_sim = rand_set1.similarity(rand_set2)

            rows.append([id, len1, len2, round(jac, 3), round(rand_jac, 3), round(jac_p, 3), round(rand_jac_p, 3), round(sim, 3), round(rand_sim, 3)])


    print(f"# set 1 -- {args.datasets[0]}")
    print(f"# set 2 -- {args.datasets[1]}")
    print()
    header = ["pid","cnt1","cnt2","jac","rand_jac","jac_p","rand_jac_p","onto_sim","rand_onto_sim"]
    print("\t".join(header))
    for x in rows:
        print("\t".join([str(y) for y in x]))



if __name__=="__main__":
    main()
