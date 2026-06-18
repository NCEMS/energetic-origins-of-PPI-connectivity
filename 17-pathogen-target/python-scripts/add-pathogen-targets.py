#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Add pathogen effector target annotations to an annotated "
            "A. thaliana node dataframe."
        )
    )

    parser.add_argument(
        "--input_nodes",
        required=True,
        help="Input annotated nodes dataframe as a pandas pickle file."
    )

    parser.add_argument(
        "--pathogen_data",
        required=True,
        help=(
            "CSV file containing pathogen effector-target data. "
            "Required columns: Effectors, Target_Proteins."
        )
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory where the updated nodes pickle will be written."
    )

    parser.add_argument(
        "--output_prefix",
        required=True,
        help="Output filename prefix."
    )

    parser.add_argument(
        "--organism_tag",
        required=True,
        help="Organism tag used in the output filename."
    )

    parser.add_argument(
        "--output_suffix",
        required=True,
        help="Output filename suffix."
    )

    return parser.parse_args()


def infer_node_column(nodes_df: pd.DataFrame) -> str:
    """
    Infer the column in nodes_df containing the AGI/node identifier.

    This is written defensively because different pipeline steps may use
    slightly different node identifier column names.
    """
    candidate_columns = [
        "node",
        "Node",
        "name",
        "Name",
        "protein",
        "Protein",
        "gene",
        "Gene",
        "agi",
        "AGI",
        "locus",
        "Locus",
        "arabidopsis_locus",
        "Arabidopsis_Locus",
    ]

    for col in candidate_columns:
        if col in nodes_df.columns:
            return col

    raise ValueError(
        "Could not infer the node identifier column in the input nodes dataframe. "
        f"Available columns are: {list(nodes_df.columns)}"
    )


def read_nodes(input_nodes: str) -> pd.DataFrame:
    """
    Read the input nodes dataframe.
    """
    input_nodes = Path(input_nodes)

    if not input_nodes.exists():
        raise FileNotFoundError(f"Input nodes file does not exist: {input_nodes}")

    if input_nodes.suffix == ".pkl":
        return pd.read_pickle(input_nodes)

    if input_nodes.suffix == ".csv":
        return pd.read_csv(input_nodes)

    raise ValueError(
        "Unsupported input_nodes file type. Expected .pkl or .csv, "
        f"but got: {input_nodes}"
    )


def read_pathogen_data(pathogen_data: str) -> pd.DataFrame:
    """
    Read and validate the pathogen effector-target table.
    """
    pathogen_data = Path(pathogen_data)

    if not pathogen_data.exists():
        raise FileNotFoundError(f"Pathogen data file does not exist: {pathogen_data}")

    patho_df = pd.read_csv(pathogen_data)

    required_columns = {"Effectors", "Target_Proteins"}
    missing = required_columns - set(patho_df.columns)

    if missing:
        raise ValueError(
            "Pathogen data file is missing required columns: "
            f"{sorted(missing)}. Available columns are: {list(patho_df.columns)}"
        )

    patho_df = patho_df[["Effectors", "Target_Proteins"]].copy()

    # Normalize identifiers and effector names.
    patho_df["Effectors"] = patho_df["Effectors"].astype("string").str.strip()
    patho_df["Target_Proteins"] = (
        patho_df["Target_Proteins"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    # Remove incomplete rows.
    patho_df = patho_df.dropna(subset=["Effectors", "Target_Proteins"])
    patho_df = patho_df[
        (patho_df["Effectors"] != "") &
        (patho_df["Target_Proteins"] != "")
    ]

    # Remove exact duplicate effector-target pairs.
    patho_df = patho_df.drop_duplicates()

    return patho_df


def aggregate_effectors(patho_df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse the pathogen table to one row per target AGI.

    Output columns:
        node
        pathogen_target
        pathogen_effectors
        pathogen_effector_count
    """
    aggregated = (
        patho_df
        .groupby("Target_Proteins", as_index=False)
        .agg(
            pathogen_effectors=(
                "Effectors",
                lambda x: ";".join(sorted(set(x.dropna())))
            ),
            pathogen_effector_count=(
                "Effectors",
                lambda x: len(set(x.dropna()))
            )
        )
        .rename(columns={"Target_Proteins": "node"})
    )

    aggregated["pathogen_target"] = 1

    return aggregated[
        [
            "node",
            "pathogen_target",
            "pathogen_effectors",
            "pathogen_effector_count",
        ]
    ]


def add_pathogen_annotations(
    nodes_df: pd.DataFrame,
    aggregated_pathogen_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Left-join pathogen target annotations onto the node dataframe.
    """
    node_col = infer_node_column(nodes_df)

    annotated_df = nodes_df.copy()

    # Create a normalized temporary join key.
    annotated_df["_pathogen_join_key"] = (
        annotated_df[node_col]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    merge_df = aggregated_pathogen_df.copy()
    merge_df["_pathogen_join_key"] = (
        merge_df["node"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    merge_df = merge_df.drop(columns=["node"])

    annotated_df = annotated_df.merge(
        merge_df,
        on="_pathogen_join_key",
        how="left",
        validate="m:1",
    )

    annotated_df = annotated_df.drop(columns=["_pathogen_join_key"])

    # Fill non-targets.
    annotated_df["pathogen_target"] = (
        annotated_df["pathogen_target"]
        .fillna(0)
        .astype(int)
    )

    annotated_df["pathogen_effectors"] = (
        annotated_df["pathogen_effectors"]
        .fillna("")
        .astype(str)
    )

    annotated_df["pathogen_effector_count"] = (
        annotated_df["pathogen_effector_count"]
        .fillna(0)
        .astype(int)
    )

    return annotated_df


def main():
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = (
        output_dir /
        f"{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )

    nodes_df = read_nodes(args.input_nodes)
    patho_df = read_pathogen_data(args.pathogen_data)

    aggregated_pathogen_df = aggregate_effectors(patho_df)

    annotated_df = add_pathogen_annotations(
        nodes_df=nodes_df,
        aggregated_pathogen_df=aggregated_pathogen_df,
    )

    annotated_df.to_pickle(output_file)

    n_nodes = len(annotated_df)
    n_targets = int(annotated_df["pathogen_target"].sum())
    n_unique_target_agis = aggregated_pathogen_df.shape[0]
    n_unique_effectors = patho_df["Effectors"].nunique()

    print(f"Wrote annotated nodes dataframe to: {output_file}")
    print(f"Total nodes: {n_nodes}")
    print(f"Nodes annotated as pathogen targets: {n_targets}")
    print(f"Unique AGIs in pathogen target file: {n_unique_target_agis}")
    print(f"Unique effectors in pathogen target file: {n_unique_effectors}")


if __name__ == "__main__":
    main()
