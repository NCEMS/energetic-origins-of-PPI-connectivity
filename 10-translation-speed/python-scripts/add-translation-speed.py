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
    parser.add_argument("--trans_speed_file")
    args = parser.parse_args()

    # read in the nodes file
    nodes_df = pd.read_pickle(args.nodes)
    #nodes_df = pd.read_csv(args.nodes)

    # read in the translation speed information
    trans_speed_df = pd.read_pickle(args.trans_speed_file)

    # merge the two pd.DataFrame objects to insert translation speed information into the nodes_df
    nodes_df = nodes_df.merge(trans_speed_df[["gene", "translation_speed_score"]], how="left", left_on="node", right_on="gene")

    nodes_df.to_pickle(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife-expr-speed.pkl")

if __name__ == "__main__":
    main()
