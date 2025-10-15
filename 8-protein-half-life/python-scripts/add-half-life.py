import pandas as pd
import numpy as np
import pint
import pint_pandas
import argparse
import typing
from typing import List


def add_halflife_db1(
    nodes_df: pd.DataFrame,
    halflife_df: pd.DataFrame,
    merge_col1,
    merge_col2,
    halflife_df_cols,
) -> pd.DataFrame:
    """
    Merge protein half life information into the nodes_df

    Args:
        nodes_df (pd.DataFrame): nodes_df to be merged with halflife information
        halflife_df (pd.DataFrame): dataframe with halflife information to be merged with nodes_df
        merge_col1 (str): column within nodes_df to use as ID when merging
        merge_col2 (str): column within halflife_df to use as ID when merging
        halflife_df_cols (List[str]): list of the columns within halflife_df to be merged into nodes_df

    Returns:
        Updated nodes_df (pd.DataFrame) with half life information
    """

    halflife_df = halflife_df[halflife_df_cols]

    # drop duplicates; the half-life data used for s288c have two entries for YPR033C & YER168C for different isoforms
    # this code keeps the first instance, which corresponds to the first isoform in both cases
    halflife_df = halflife_df.drop_duplicates(subset="ENSG", keep="first")

    nodes_df = nodes_df.merge(
        halflife_df, how="left", left_on=merge_col1, right_on=merge_col2
    )

    return nodes_df


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
        "--output_suffix", default="step8", help="Suffix for output files"
    )
    parser.add_argument("--organism_tag", help="Organism label for this run")
    parser.add_argument("--halflife_db")
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)

    # ADD FIRST HALF-LIFE DATABASE HERE
    halflife_df1 = pd.read_csv(args.halflife_db)
    halflife_df1_cols = [
        "ENSG",
        "Degradation rates (min-1)",
        "R2 (quality of curve fitting)",
        "t1/2 (min)",
    ]
    merge_col1 = "node"
    merge_col2 = "ENSG"
    nodes_df = add_halflife_db1(
        nodes_df, halflife_df1, merge_col1, merge_col2, halflife_df1_cols
    )

    # ADD SECOND HALF-LIFE DATABASE HERE

    # output the results to file
    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )


if __name__ == "__main__":
    main()
