#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import List

import pandas as pd
import pint
import pint_pandas


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Add preprocessed protein half-life datasets to a node dataframe. "
            "Each half-life CSV must contain an AGI column, which is merged "
            "against the node column of the input nodes dataframe."
        )
    )

    parser.add_argument("--nodes", required=True, help="Input nodes pickle file")
    parser.add_argument("--halflife_db1", required=True, help="Processed half-life DB1 CSV")
    parser.add_argument("--halflife_db2", required=True, help="Processed half-life DB2 CSV")
    parser.add_argument("--halflife_db3", required=True, help="Processed half-life DB3 CSV")

    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where results will be saved",
    )
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files")
    parser.add_argument("--output_suffix", required=True, help="Suffix for output files")
    parser.add_argument("--organism_tag", required=True, help="Organism label for this run")

    return parser.parse_args()


def normalize_agi(series: pd.Series) -> pd.Series:
    """
    Normalize AGI-like identifiers for merging.

    This preserves the original node column in the output, but makes the merge
    robust to lowercase identifiers, surrounding whitespace, and accidental
    isoform suffixes such as AT1G01010.1.
    """
    return (
        series.astype("string")
        .str.strip()
        .str.upper()
        .str.replace(r"\.\d+$", "", regex=True)
    )


def read_halflife_csv(path: str, label: str) -> pd.DataFrame:
    """
    Read one preprocessed half-life CSV and validate that it is merge-ready.
    """
    df = pd.read_csv(path)

    # Drop accidental index columns from CSV writing, if present.
    df = df.loc[:, ~df.columns.str.match(r"^Unnamed")]

    if "AGI" not in df.columns:
        raise ValueError(
            f"{label}: expected column 'AGI' in {path}, but found columns: "
            f"{list(df.columns)}"
        )

    df = df.copy()
    df["AGI"] = normalize_agi(df["AGI"])

    # Remove rows with missing or empty AGI identifiers.
    before = len(df)
    df = df[df["AGI"].notna() & (df["AGI"] != "")]
    removed = before - len(df)
    if removed > 0:
        print(f"{label}: removed {removed:,} rows with missing/empty AGI values")

    # The preprocessing scripts should already have collapsed duplicates.
    # If duplicates remain, fail loudly rather than expanding the node table.
    duplicated = df["AGI"].duplicated(keep=False)
    if duplicated.any():
        examples = sorted(df.loc[duplicated, "AGI"].dropna().unique())[:20]
        raise ValueError(
            f"{label}: found duplicated AGI values after preprocessing. "
            f"The final merge expects one row per AGI. Example duplicates: {examples}"
        )

    return df


def add_halflife_table(
    nodes_df: pd.DataFrame,
    halflife_df: pd.DataFrame,
    label: str,
    node_col: str = "node",
    agi_col: str = "AGI",
) -> pd.DataFrame:
    """
    Left merge one preprocessed half-life table into nodes_df.

    The merge is many-to-one:
        nodes_df[node_col] -> halflife_df[agi_col]

    This preserves all rows in nodes_df and adds all non-AGI columns from the
    half-life table.
    """
    if node_col not in nodes_df.columns:
        raise ValueError(
            f"nodes_df is missing required merge column '{node_col}'. "
            f"Found columns: {list(nodes_df.columns)}"
        )

    if agi_col not in halflife_df.columns:
        raise ValueError(
            f"{label}: halflife_df is missing required merge column '{agi_col}'."
        )

    feature_cols: List[str] = [c for c in halflife_df.columns if c != agi_col]

    overlapping_cols = sorted(set(feature_cols).intersection(nodes_df.columns))
    if overlapping_cols:
        raise ValueError(
            f"{label}: these half-life columns already exist in nodes_df and would "
            f"cause ambiguous duplicate columns: {overlapping_cols}"
        )

    left_key = f"__{label}_node_merge_key"
    right_key = f"__{label}_agi_merge_key"
    indicator_col = f"__{label}_merge_status"

    reserved = {left_key, right_key, indicator_col}
    collisions = reserved.intersection(nodes_df.columns).union(
        reserved.intersection(halflife_df.columns)
    )
    if collisions:
        raise ValueError(f"{label}: temporary column name collision: {collisions}")

    left_df = nodes_df.copy()
    right_df = halflife_df.copy()

    left_df[left_key] = normalize_agi(left_df[node_col])
    right_df[right_key] = normalize_agi(right_df[agi_col])

    # Drop AGI from the right table so the final annotated node file does not
    # accumulate AGI, AGI_x, AGI_y, etc.
    right_df = right_df.drop(columns=[agi_col])

    before_rows = len(left_df)

    merged_df = left_df.merge(
        right_df,
        how="left",
        left_on=left_key,
        right_on=right_key,
        validate="m:1",
        indicator=indicator_col,
    )

    matched_nodes = int((merged_df[indicator_col] == "both").sum())

    merged_df = merged_df.drop(columns=[left_key, right_key, indicator_col])

    after_rows = len(merged_df)
    if after_rows != before_rows:
        raise RuntimeError(
            f"{label}: merge changed row count from {before_rows:,} to {after_rows:,}. "
            "This should not happen for a validated left merge."
        )

    print(
        f"{label}: added {len(feature_cols):,} columns; "
        f"matched {matched_nodes:,}/{before_rows:,} node rows"
    )

    return merged_df


def main() -> None:
    args = parse_args()

    nodes_df = pd.read_pickle(args.nodes)

    db_paths = [
        ("halflife_db1", args.halflife_db1),
        ("halflife_db2", args.halflife_db2),
        ("halflife_db3", args.halflife_db3),
    ]

    for label, path in db_paths:
        halflife_df = read_halflife_csv(path, label)
        nodes_df = add_halflife_table(nodes_df, halflife_df, label)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = (
        output_dir
        / f"{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )

    nodes_df.to_pickle(output_file)
    print(f"Wrote annotated node dataframe to: {output_file}")


if __name__ == "__main__":
    main()
