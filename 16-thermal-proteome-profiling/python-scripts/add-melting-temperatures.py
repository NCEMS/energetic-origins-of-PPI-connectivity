#!/usr/bin/env python3

"""
Integrate preprocessed protein melting-temperature datasets into a node table.

Each preprocessed dataset must contain:
    gene

The node annotation table must contain:
    node

All non-key columns from each preprocessed dataset are retained. The datasets
are left-joined independently so that every row from the node table remains in
the final output.

Expected inputs:
    - Pickled pandas node table
    - Lyu 2023 preprocessed CSV
    - Meltome Atlas preprocessed CSV
    - Volkening 2019 preprocessed CSV

Output:
    - Pickled pandas DataFrame containing the original node annotations plus
      all columns from the three preprocessed melting-temperature datasets.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Left-join preprocessed melting-temperature datasets into a "
            "pickled node annotation table."
        )
    )

    parser.add_argument(
        "--nodes",
        required=True,
        type=Path,
        help="Input pickled pandas node annotation table.",
    )
    parser.add_argument(
        "--lyu_data",
        required=True,
        type=Path,
        help="Preprocessed Lyu 2023 CSV file.",
    )
    parser.add_argument(
        "--meltome_data",
        required=True,
        type=Path,
        help="Preprocessed Meltome Atlas CSV file.",
    )
    parser.add_argument(
        "--volkening_data",
        required=True,
        type=Path,
        help="Preprocessed Volkening 2019 CSV file.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output path for the integrated pickled node table.",
    )
    parser.add_argument(
        "--node_col",
        default="node",
        help="Merge column in the node annotation table. Default: node",
    )
    parser.add_argument(
        "--gene_col",
        default="gene",
        help="Merge column in each preprocessed CSV. Default: gene",
    )

    return parser.parse_args()


def require_file(path: Path, label: str) -> None:
    """Raise a clear error when an expected input file is missing."""

    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist or is not a file: {path}")


def normalize_agi(series: pd.Series) -> pd.Series:
    """
    Normalize Arabidopsis AGI identifiers to uppercase locus-level IDs.

    Examples:
        AT1G01010   -> AT1G01010
        at1g01010   -> AT1G01010
        AT1G01010.1 -> AT1G01010

    Missing values remain pandas missing values.
    """

    normalized = series.astype("string").str.strip().str.upper()

    # Remove a terminal transcript or isoform suffix if one is present.
    normalized = normalized.str.replace(r"\.\d+$", "", regex=True)

    # Convert empty strings back to missing values.
    normalized = normalized.replace("", pd.NA)

    return normalized


def read_preprocessed_dataset(
    path: Path,
    dataset_name: str,
    gene_col: str,
) -> pd.DataFrame:
    """
    Read and validate one preprocessed melting-temperature dataset.

    The returned table contains one row per normalized AGI. Duplicate AGIs are
    treated as an error because they could multiply node rows during merging.
    """

    require_file(path, dataset_name)

    df = pd.read_csv(path, low_memory=False)

    if gene_col not in df.columns:
        raise KeyError(
            f"{dataset_name} is missing required merge column "
            f"{gene_col!r}: {path}"
        )

    df = df.copy()
    df[gene_col] = normalize_agi(df[gene_col])

    missing_gene_count = int(df[gene_col].isna().sum())
    if missing_gene_count:
        raise ValueError(
            f"{dataset_name} contains {missing_gene_count:,} row(s) with a "
            f"missing or empty {gene_col!r} value."
        )

    duplicate_mask = df.duplicated(subset=[gene_col], keep=False)

    if duplicate_mask.any():
        duplicate_examples = (
            df.loc[duplicate_mask, gene_col]
            .value_counts()
            .head(20)
            .to_dict()
        )

        raise ValueError(
            f"{dataset_name} contains "
            f"{df.loc[duplicate_mask, gene_col].nunique():,} duplicated AGI(s). "
            "The integration requires one row per AGI. Example duplicate "
            f"counts: {duplicate_examples}"
        )

    if df.columns.duplicated().any():
        duplicate_columns = df.columns[df.columns.duplicated()].tolist()
        raise ValueError(
            f"{dataset_name} contains duplicate column names: "
            f"{duplicate_columns}"
        )

    return df


def check_column_collisions(
    nodes_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    gene_col: str,
    dataset_name: str,
) -> None:
    """
    Ensure annotation columns will not overwrite existing node columns.

    The annotation merge key is excluded because it is removed after merging.
    """

    annotation_columns = set(annotation_df.columns) - {gene_col}
    collisions = sorted(annotation_columns.intersection(nodes_df.columns))

    if collisions:
        raise ValueError(
            f"{dataset_name} contains annotation column(s) already present in "
            f"the node table: {collisions}. Refusing to create merge suffixes "
            "or overwrite existing annotations."
        )


def add_dataset(
    nodes_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    dataset_name: str,
    node_col: str,
    gene_col: str,
) -> pd.DataFrame:
    """
    Left-join one preprocessed dataset into the node table.

    The merge is validated as many-to-one:
        - the node table may theoretically contain repeated node identifiers;
        - each preprocessed table must contain at most one row per AGI.

    In the expected pipeline, the node table itself should have one row per node.
    """

    check_column_collisions(
        nodes_df=nodes_df,
        annotation_df=annotation_df,
        gene_col=gene_col,
        dataset_name=dataset_name,
    )

    source_row_count = len(nodes_df)
    annotation_columns = [
        column for column in annotation_df.columns if column != gene_col
    ]

    merged = nodes_df.merge(
        annotation_df,
        how="left",
        left_on=node_col,
        right_on=gene_col,
        validate="many_to_one",
        sort=False,
    )

    if len(merged) != source_row_count:
        raise RuntimeError(
            f"{dataset_name} merge changed the node-row count from "
            f"{source_row_count:,} to {len(merged):,}."
        )

    merged = merged.drop(columns=[gene_col])

    matched_rows = int(merged[annotation_columns].notna().any(axis=1).sum())

    print(
        f"{dataset_name}: added {len(annotation_columns):,} column(s); "
        f"{matched_rows:,} of {len(merged):,} node row(s) received at least "
        "one annotation."
    )

    return merged


def report_added_columns(
    original_columns: Sequence[str],
    final_df: pd.DataFrame,
) -> None:
    """Print a concise report of columns added during integration."""

    original_set = set(original_columns)
    added_columns = [
        column for column in final_df.columns if column not in original_set
    ]

    print(f"Original node columns: {len(original_columns):,}")
    print(f"Added melting-temperature columns: {len(added_columns):,}")
    print(f"Final node columns: {len(final_df.columns):,}")

    if added_columns:
        print("Added columns:")
        for column in added_columns:
            nonmissing = int(final_df[column].notna().sum())
            print(f"  {column}: {nonmissing:,} nonmissing value(s)")


def main() -> None:
    """Run the integration workflow."""

    args = parse_args()

    require_file(args.nodes, "Node table")

    print(f"Reading node table: {args.nodes}")
    nodes_df = pd.read_pickle(args.nodes)

    if args.node_col not in nodes_df.columns:
        raise KeyError(
            f"Node table is missing required merge column "
            f"{args.node_col!r}: {args.nodes}"
        )

    if nodes_df.columns.duplicated().any():
        duplicate_columns = nodes_df.columns[
            nodes_df.columns.duplicated()
        ].tolist()
        raise ValueError(
            f"Node table contains duplicate column names: {duplicate_columns}"
        )

    nodes_df = nodes_df.copy()
    nodes_df[args.node_col] = normalize_agi(nodes_df[args.node_col])

    missing_node_count = int(nodes_df[args.node_col].isna().sum())
    if missing_node_count:
        raise ValueError(
            f"Node table contains {missing_node_count:,} row(s) with missing "
            f"or empty {args.node_col!r} identifiers."
        )

    duplicate_node_count = int(
        nodes_df.duplicated(subset=[args.node_col], keep=False).sum()
    )
    if duplicate_node_count:
        print(
            "WARNING: node table contains "
            f"{duplicate_node_count:,} row(s) belonging to duplicated node "
            "identifiers. The merges remain safe because the annotation "
            "tables are constrained to one row per AGI."
        )

    original_columns = nodes_df.columns.tolist()
    original_index = nodes_df.index.copy()

    datasets = [
        (
            "Lyu 2023",
            args.lyu_data,
        ),
        (
            "Meltome Atlas",
            args.meltome_data,
        ),
        (
            "Volkening 2019",
            args.volkening_data,
        ),
    ]

    integrated_df = nodes_df

    for dataset_name, path in datasets:
        print(f"Reading {dataset_name}: {path}")

        annotation_df = read_preprocessed_dataset(
            path=path,
            dataset_name=dataset_name,
            gene_col=args.gene_col,
        )

        integrated_df = add_dataset(
            nodes_df=integrated_df,
            annotation_df=annotation_df,
            dataset_name=dataset_name,
            node_col=args.node_col,
            gene_col=args.gene_col,
        )

    # pandas merge creates a new RangeIndex. Restore the original node-table
    # index so downstream code sees the same row/index relationship.
    integrated_df.index = original_index

    report_added_columns(
        original_columns=original_columns,
        final_df=integrated_df,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    integrated_df.to_pickle(args.output)

    print(f"Integrated node table written to: {args.output}")
    print(f"Final dimensions: {integrated_df.shape[0]:,} rows × "
          f"{integrated_df.shape[1]:,} columns")


if __name__ == "__main__":
    main()
