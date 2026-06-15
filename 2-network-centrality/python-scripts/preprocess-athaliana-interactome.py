#!/usr/bin/env python3

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="Preprocess Arabidopsis PPI table into edge and node files."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input tab-delimited Arabidopsis PPI file",
    )

    parser.add_argument(
        "--edges",
        required=True,
        help="Output edge CSV file",
    )

    parser.add_argument(
        "--nodes",
        required=True,
        help="Output node CSV file",
    )

    args = parser.parse_args()

    # read input PPI table
    df = pd.read_csv(args.input, sep="\t")

    # build edge table
    edges_df = df[["arabidopsis_locus_a", "arabidopsis_locus_b", "confidence_tier"]].copy()
    edges_df = edges_df.rename(
        columns={
            "arabidopsis_locus_a": "source",
            "arabidopsis_locus_b": "target",
        }
    )

    # build node table from unique edge endpoints
    nodes = pd.unique(edges_df[["source", "target"]].values.ravel())

    nodes_df = pd.DataFrame({"name": nodes})

    # put in a dummy column for _wkshell, which is expected by the next program in line
    nodes_df["_wkshell"] = 0.0

    # write output files
    edges_df.to_csv(args.edges, index=False)
    nodes_df.to_csv(args.nodes, index=False)

    print(f"Wrote {len(edges_df):,} edges to {args.edges}")
    print(f"Wrote {len(nodes_df):,} nodes to {args.nodes}")


if __name__ == "__main__":
    main()
