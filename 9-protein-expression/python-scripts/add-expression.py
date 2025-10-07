import pandas as pd
import numpy as np
import pint
import pint_pandas
import argparse
import typing
from typing import List

def add_expression(nodes_df: pd.DataFrame, expression_df: pd.DataFrame, merge_col1, merge_col2, expression_df_cols) -> pd.DataFrame:
    """
    Merge protein expression information into the nodes_df

    Args:
        nodes_df (pd.DataFrame): nodes_df to be merged with expression information
        expression_df (pd.DataFrame): dataframe with expression information to be merged with nodes_df
        merge_col1 (str): column within nodes_df to use as ID when merging
        merge_col2 (str): column within expression_df to use as ID when merging
        expression_df_cols (List[str]): list of the columns within expression_df to be merged into nodes_df

    Returns:
        Updated nodes_df (pd.DataFrame) with half life information
    """

    expression_df = expression_df[expression_df_cols]
    nodes_df = nodes_df.merge(expression_df, how="left", left_on=merge_col1, right_on=merge_col2)

    return nodes_df

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
    parser.add_argument("--output_prefix", default="0", help="Prefix for output file")
    parser.add_argument("--output_suffix", default="0", help="Suffix for output file")
    parser.add_argument(
        "--organism_tag",
        help="Organism label for this run")
    parser.add_argument("--expression_db")
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)

    expression_df = pd.read_csv(args.expression_db)

    expression_df_cols = ["Systematic Name", "Mean molecules per cell", "Median molecules per cell"]

    merge_col1 = "node"

    merge_col2 = "Systematic Name"

    nodes_df = add_expression(nodes_df, expression_df, merge_col1, merge_col2, expression_df_cols)

    nodes_df.to_pickle(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl")

if __name__ == "__main__":
    main()
