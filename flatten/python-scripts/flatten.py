import pandas as pd
import typing
from typing import List
import argparse
import pint
from pint import UnitRegistry
import pint_pandas

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

            new_row['IDR_index'] = idr_index
            flattened_rows.append(new_row)

    idr_df = pd.DataFrame(flattened_rows)

    return idr_df


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
    parser.add_argument("--output_prefix", default="0", help="Prefix for output files")
    parser.add_argument(
        "--organism_tag",
        help="Organism label for this run")
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)

    columns_to_drop = ["ENSG", "Systematic Name", "dH", "dCp", "dS"]

    ureg = UnitRegistry()

    #nodes_df["Ghosh-Dill-dG"] = nodes_df['Ghosh-Dill-dG'].apply(lambda x: ureg(x).to_base_units().magnitude)
    nodes_df["Ghosh-Dill-dG"] = nodes_df["Ghosh-Dill-dG"].apply(
        lambda x: x.to("kilocalorie / mole").magnitude if isinstance(x, pint.Quantity) else x
    )

    if nodes_df["signalP_trimmed_sequence_y"].equals(nodes_df["signalP_trimmed_sequence_x"]):
        nodes_df = nodes_df.rename(columns={"signalP_trimmed_sequence_x":"signalP_trimmed_sequence"})
        columns_to_drop.append("signalP_trimmed_sequence_y")

    nodes_df = nodes_df.drop(columns=columns_to_drop)

    nodes_df.to_csv(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-final-per-node.csv", index=False)

    nested_columns = ["albatross", "cider", "IDR_sequences"]

    flat_nodes_df = flatten_nodes(nodes_df, nested_columns)

    flat_nodes_df.to_csv(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-final-per-IDR.csv", index=False)


if __name__ == "__main__":

    main()
