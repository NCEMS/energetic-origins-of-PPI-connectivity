import os, sys
import numpy as np
import pandas as pd
import pint
import pint_pandas
import argparse
import typing

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--nodes", required=True, help="Input nodes file"
    )
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where results will be saved",
    )
    parser.add_argument("--output_prefix", default="0", help="Prefix for output files")
    parser.add_argument(
        "--organism_tag",
        help="Organism label for this run")
    parser.add_argument("--trans_speed_db")
    args = parser.parse_args()

    #nodes_df = pd.read_pickle(args.nodes)
    nodes_df = pd.read_csv(args.nodes)

    trans_speed_df = pd.read_csv(args.trans_speed_db, sep="\t", header=None, names=["gene", "num_ncs", "raw_counts"])

    #print (type(trans_speed_df["raw_counts"].iloc[0]))
    trans_speed_df["raw_counts"] = trans_speed_df["raw_counts"].apply(lambda s: [float(x) for x in s.split(",")])
    trans_speed_df["canonical_frame_counts"] = trans_speed_df["float_list"].apply(lambda x: x[::3])

    #trans_speed_df_cols = ["ENSG", "Degradation rates (min-1)", "R2 (quality of curve fitting)", "t1/2 (min)"]
    #merge_col1 = "node"
    #merge_col2 = "ENSG"
    #nodes_df = add_trans_speed(nodes_df, halflife_df, merge_col1, merge_col2, halflife_df_cols)

    #nodes_df.to_pickle(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife-speed.pkl")

if __name__ == "__main__":
    main()
