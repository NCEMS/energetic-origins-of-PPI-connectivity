import pandas as pd
import typing
from typing import List
import argparse
import pint
from pint import UnitRegistry
import pint_pandas
import numpy as np


def flatten_nodes(nodes_df: pd.DataFrame, nested_columns: List[str]) -> pd.DataFrame:
    """
    Takes a nodes_df pd.DataFrame and flattens it such that each row corresponds to one IDR

    Args:
        nodes_df (pd.DataFrame): input nodes_df from previous pipeline step
        nested_columns (List[str]): list of the column names that are nested dictionaries to be flattened

    Returns:
        The flattened pd.DataFrame to be output as a .csv file
    """

    flattened_rows = []

    for _, row in nodes_df.iterrows():

        if not isinstance(row[nested_columns[0]], dict):
            continue

        idr_keys = row[nested_columns[0]].keys()

        for idr_index in idr_keys:
            new_row = {}

            for col in nodes_df.columns:
                if col not in nested_columns:
                    new_row[col] = row[col]

            for col in nested_columns:
                value = row[col].get(idr_index)
                if isinstance(value, dict):
                    for subkey, subval in value.items():
                        new_row[f"{col}_{subkey}"] = subval
                else:
                    new_row[col] = value

            new_row["IDR_index"] = idr_index
            flattened_rows.append(new_row)

    idr_df = pd.DataFrame(flattened_rows)

    return idr_df


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
        "--output_suffix", default="final", help="Suffix for output files"
    )
    parser.add_argument("--organism_tag", help="Organism label for this run")
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)

    # print the columns in the pd.DataFrame
    junk = nodes_df.columns
    for j in junk:
        print(j)

    # drop some unneeded columns
    columns_to_drop = [
        "dH",
        "dCp",
        "dS",
    ]

    nodes_df = nodes_df.drop(columns=columns_to_drop)

    # remove pint units from Ghosh-Dill-dG
    ureg = UnitRegistry()
    nodes_df["Ghosh-Dill-dG"] = nodes_df["Ghosh-Dill-dG"].apply(
        lambda x: (
            x.to("kilocalorie / mole").magnitude if isinstance(x, pint.Quantity) else x
        )
    )

    # save as a pickle file with reprocessed columns
    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )

    # save to a .csv file
    nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.csv",
        index=False,
    )

    # reformat to expand dictionaries across N IDRs
    nested_columns = ["albatross", "cider", "IDR_sequences", "IDR_ranges"]

    flat_nodes_df = flatten_nodes(nodes_df, nested_columns)

    flat_nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-per-IDR-{args.output_suffix}.csv",
        index=False,
    )

    flat_nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-per-IDR-{args.output_suffix}.pkl"
    )


if __name__ == "__main__":

    main()
