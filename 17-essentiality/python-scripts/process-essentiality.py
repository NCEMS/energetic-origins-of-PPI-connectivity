import os
import sys
import argparse
import pandas as pd

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--input_data", required=True)
    parser.add_argument("--strain_filter", required=True)
    args = parser.parse_args()

    # load file
    df = pd.read_csv(args.input_data, sep="\t", on_bad_lines="warn")

    # filter to get strain of interest
    df = df[df["Strain Background"] == args.strain_filter]
    print (df.info())

    # rename for consistency with previous pipeline steps
    df.rename(columns={"Gene Systematic Name":"node"}, inplace=True)

    # extract columns of interest
    df = df[["node", "Strain Background", "Details"]]

    #collapsed_df = df.groupby("node", as_index=False)["Details"].agg(" | ".join)
    collapsed_df = df.groupby("node", as_index=False)["Details"].agg(
        lambda x: " | ".join(x.dropna().astype(str))
    )

    # save to output csv
    collapsed_df.to_csv(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-essentiality-processed.csv", index=False)

if __name__ == "__main__":
    main()
