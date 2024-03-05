#!/usr/bin/env python3

import os
import sys
import csv
import argparse

from pyhpo.ontology import Ontology
from pyhpo.set import HPOSet, BasicHPOSet
from pyhpo.stats import EnrichmentModel
import simple_icd_10_cm as cm


def argue():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dataset", required=True,
                        help="Dataset for computing HPO set statistics.")
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


def main():
    args = argue()
    _ = Ontology()

    df_name = os.path.basename(args.dataset)
    df = ready(args.dataset)

    df_col_pos = get_column_positions(df, ["pid", "codes"])
    df_lookup = make_id_lookup(df, df_col_pos["pid"])
    ids = list(df_lookup.keys())

    valid_hpo = _.to_dataframe().index.tolist()

    rows = []
    for id in ids:
        if id != "pid":
            codes = df[df_lookup[id]][df_col_pos["codes"]].split(";")
            hpo = []
            icd = []
            other = []

            for cde in codes:
                if cde in valid_hpo:
                    hpo.append(cde)
                elif cm.is_valid_item(cde):
                    icd.append(cm.add_dot(cde))
                else:
                    other.append(cde)

            all_hpo_cnt = len(set(hpo))
            icd_cnt = len(set(icd))
            other_cnt = len(set(other))
            total_cnt = len(set(codes))

            if all_hpo_cnt == 0:
                rows.append((id, 0, 0, icd_cnt, other_cnt, total_cnt, 0, 0, None, None))
            else:
                s = BasicHPOSet.from_queries(hpo)

                terminal_hpo_cnt = len(s.toJSON())

                gene_model = EnrichmentModel("gene")
                genes = gene_model.enrichment(method="hypergeom", hposet=s)
                top10_genes = [x["item"] for x in genes[:20]]

                omim_model = EnrichmentModel("omim")
                omims = omim_model.enrichment(method="hypergeom", hposet=s)
                top10_omims = [x["item"] for x in omims[:20]]

                info_con = s.information_content(kind="omim")

                rows.append((id, 
                             all_hpo_cnt, 
                             terminal_hpo_cnt,
                             icd_cnt,
                             other_cnt,
                             total_cnt,
                             round(info_con["mean"], 3),
                             round(info_con["total"], 3),
                             ";".join([str(x.toJSON()["id"]) for x in top10_genes]),
                             ";".join([str(x.toJSON()["id"]) for x in top10_omims])
                             ))

    header = ["pid","all_hpo_cnt","terminal_hpo_cnt","icd_cnt","other_cnt","total_cnt",
              "mean_info_cont","total_info_cont","top10genes","top10omims"]
    print("\t".join(header))
    for x in rows:
        print("\t".join([str(y) for y in x]))



if __name__=="__main__":
    main()
