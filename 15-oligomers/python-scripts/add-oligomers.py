import os, sys
import argparse
import pandas as pd
import numpy as np
import pint
import pint_pandas


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--input_nodes")
    parser.add_argument("--olig_data")
    args = parser.parse_args()

    # read in the pre-processed oligomer information from Complex Portal
    olig_df = pd.read_csv(args.olig_data)

    # read in the nodes_df from the previous pipeline step
    nodes_df = pd.read_pickle(args.input_nodes)

    # merge the nodes_df and olig_df on node
    nodes_df = nodes_df.merge(olig_df, on="node", how="left")

    # save the updated nodes_df to file
    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife-expr-speed-PTMGPT2-LiPMS-entanglement-chaperones-oligomers.pkl"
    )


if __name__ == "__main__":
    main()
