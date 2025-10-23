import os
import sys
import argparse
import pandas as pd


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_suffix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--input_nodes", required=True)
    parser.add_argument("--domain_data", required=True)
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.input_nodes)
    domain_df = pd.read_pickle(args.domain_data)

    nodes_df = nodes_df.merge(domain_df, on="UniProtKB-AC", how="left")

    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )


if __name__ == "__main__":
    main()
