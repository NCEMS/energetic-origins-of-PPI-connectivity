import os, sys
import numpy as np
import pandas as pd
import pint
import pint_pandas
import argparse
import typing


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", required=True, help="Input nodes file")
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where results will be saved",
    )
    parser.add_argument("--output_prefix", default="0", help="Prefix for output files")
    parser.add_argument(
        "--output_suffix", default="step10", help="Suffix for output files"
    )
    parser.add_argument("--organism_tag", help="Organism label for this run")
    parser.add_argument("--input_file")
    args = parser.parse_args()

    # read in the nodes file
    nodes_df = pd.read_pickle(args.nodes)

    # read in the translation efficiency information
    TE_df = pd.read_csv(args.input_file, sep="\t")

    # merge the two pd.DataFrame objects to insert translation speed information into the nodes_df
    nodes_df = nodes_df.merge(
        TE_df[["gene", "log2_TE"]], how="left", left_on="node", right_on="gene"
    )

    # save the resulting pd.DataFrame containing log2_TE information
    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )


if __name__ == "__main__":
    main()
