#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
import numpy as np


def infer_node_column(df: pd.DataFrame, file_label: str) -> str:
    """
    Infer the node identifier column.

    Accepts common node column names used across the pipeline and Cytoscape.
    """
    candidates = ["node", "name", "AGI", "agi"]

    for col in candidates:
        if col in df.columns:
            return col

    raise ValueError(
        f"Could not infer node column for {file_label}. "
        f"Tried {candidates}. Available columns: {list(df.columns)}"
    )


def clean_node_ids(series: pd.Series) -> pd.Series:
    """
    Normalize node IDs for joining.
    """
    return series.astype(str).str.strip()


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Left-join precomputed Cytoscape wk-shell results onto a canonical "
            "node file."
        )
    )

    parser.add_argument(
        "--nodes",
        required=True,
        help="Canonical input node CSV.",
    )

    parser.add_argument(
        "--wkshell",
        required=True,
        help=(
            "Precomputed Cytoscape wk-shell node table. Expected to contain "
            "a node/name column and _wkshell."
        ),
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output nodes CSV with _wkshell added.",
    )

    parser.add_argument(
        "--fill-missing",
        type=float,
        default=np.nan,
        help=(
            "Value to use for nodes without a precomputed wk-shell value. "
            "Default: NaN. Existing network-centrality.py fills _wkshell NaN "
            "with 0.0 internally."
        ),
    )

    args = parser.parse_args()

    nodes_df = pd.read_csv(args.nodes)
    wkshell_df = pd.read_csv(args.wkshell)

    nodes_node_col = infer_node_column(nodes_df, "canonical nodes file")
    wkshell_node_col = infer_node_column(wkshell_df, "wk-shell file")

    if "_wkshell" not in wkshell_df.columns:
        raise ValueError(
            "wk-shell file must contain a column named '_wkshell'. "
            f"Available columns: {list(wkshell_df.columns)}"
        )

    nodes_df = nodes_df.copy()
    wkshell_df = wkshell_df.copy()

    nodes_df["_join_node"] = clean_node_ids(nodes_df[nodes_node_col])
    wkshell_df["_join_node"] = clean_node_ids(wkshell_df[wkshell_node_col])

    if nodes_df["_join_node"].duplicated().any():
        examples = (
            nodes_df.loc[nodes_df["_join_node"].duplicated(), "_join_node"]
            .head(10)
            .tolist()
        )
        raise ValueError(
            "Canonical nodes file contains duplicated node IDs. "
            f"Examples: {examples}"
        )

    if wkshell_df["_join_node"].duplicated().any():
        examples = (
            wkshell_df.loc[wkshell_df["_join_node"].duplicated(), "_join_node"]
            .head(10)
            .tolist()
        )
        raise ValueError(
            "wk-shell file contains duplicated node IDs. "
            f"Examples: {examples}"
        )

    # Keep only columns needed for this preprocessing step.
    wkshell_keep = wkshell_df[["_join_node", "_wkshell"]].copy()

    # If the raw node file already has a dummy _wkshell column, replace it.
    if "_wkshell" in nodes_df.columns:
        nodes_df = nodes_df.drop(columns=["_wkshell"])

    merged = nodes_df.merge(
        wkshell_keep,
        on="_join_node",
        how="left",
        validate="one_to_one",
    )

    if not pd.isna(args.fill_missing):
        merged["_wkshell"] = merged["_wkshell"].fillna(args.fill_missing)

    n_nodes = len(merged)
    n_with_wkshell = merged["_wkshell"].notna().sum()
    n_missing_wkshell = merged["_wkshell"].isna().sum()

    print(f"Input nodes:              {n_nodes:,}")
    print(f"Nodes with _wkshell:      {n_with_wkshell:,}")
    print(f"Nodes missing _wkshell:   {n_missing_wkshell:,}")

    # Drop internal join column before writing.
    merged = merged.drop(columns=["_join_node"])

    output_file = Path(args.output)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    merged.to_csv(output_file, index=False, na_rep=np.nan)

    print(f"Wrote nodes with wk-shell to: {output_file}")


if __name__ == "__main__":
    main()
