import os, sys
import argparse
import pandas as pd
from functools import reduce
import pint
import pint_pandas

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--input_nodes")
    parser.add_argument("--input_lipms")
    args = parser.parse_args()

    # load data into memory
    lipms = pd.read_csv(args.input_lipms)
    nodes_df = pd.read_pickle(args.input_nodes)

    # create the merged pd.DataFrame
    nodes_df = nodes_df.merge(lipms, how="left", right_on="Protein ID", left_on="node")

    print(nodes_df.columns)

    nodes_df = nodes_df.drop("Protein ID", axis=1)

    # save the result to file
    nodes_df.to_pickle(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife-expr-speed-PTMGPT2-LiPMS.pkl")

if __name__ == "__main__":
    main()
