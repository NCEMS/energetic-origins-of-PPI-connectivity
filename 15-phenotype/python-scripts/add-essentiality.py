import os, sys
import argparse
import pandas as pd


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_suffix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--input_data", required=True)
    parser.add_argument("--input_nodes", required=True)
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.input_nodes)
    ess_df = pd.read_csv(args.input_data)

    essential_nodes = set(ess_df["node"])

    nodes_df["essential"] = nodes_df["node"].isin(essential_nodes).astype(int)

    details_dict = ess_df.set_index("node")["Details"].to_dict()

    nodes_df["essentiality_details"] = nodes_df["node"].map(details_dict)

    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )


if __name__ == "__main__":
    main()
