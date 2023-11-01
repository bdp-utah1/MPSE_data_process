#!/usr/bin/env python3

import os
import sys
import csv
import argparse
from pyhpo.ontology import Ontology
from pyhpo.set import HPOSet, BasicHPOSet
from pyhpo.stats import EnrichmentModel


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
            raw = df[df_lookup[id]][df_col_pos["codes"]].split(";")
            real = [x for x in raw if x in valid_hpo]
            # s = HPOSet.from_queries(real)
            s = BasicHPOSet.from_queries(real)

            term_cnt = len(s.toJSON())

            #omim_dis = s.omim_diseases()
            #omim_str = [str(x.toJSON()["id"]) for x in omim_dis]

            gene_model = EnrichmentModel("gene")
            genes = gene_model.enrichment(method="hypergeom", hposet=s)
            top10_genes = [x["item"] for x in genes[:20]]

            omim_model = EnrichmentModel("omim")
            omims = omim_model.enrichment(method="hypergeom", hposet=s)
            top10_omims = [x["item"] for x in omims[:20]]

            info_con = s.information_content(kind="omim")
            #variance = s.variance()

            rows.append((id, 
                         term_cnt, 
                         round(info_con["mean"], 3),
                         round(info_con["total"], 3),
                         # variance, 
                         # ";".join(omim_str),
                         ";".join([str(x.toJSON()["id"]) for x in top10_genes]),
                         ";".join([str(x.toJSON()["id"]) for x in top10_omims])
                         ))

    header = ["pid","term_count","mean_info_cont","total_info_cont","top10genes","top10omims"]
    print("\t".join(header))
    for x in rows:
        print("\t".join([str(y) for y in x]))



if __name__=="__main__":
    main()
