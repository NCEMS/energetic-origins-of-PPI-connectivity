#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

# These imports are retained because your node pickle files require them.
import pint
import pint_pandas


GENE_ID_COL = "Gene_ID"
NODE_COL = "node"


def clean_gene_id(value) -> str | None:
    """
    Clean Gene_ID values for merging to node-level AGI IDs.

    Behavior:
        At1g01010  -> AT1G01010
        AT1G01010  -> AT1G01010
        AT1G01010.1 -> AT1G01010

    The isoform-suffix removal is included as a guard, even though these TE files
    appear to use base AGI identifiers.
    """
    if value is None or pd.isna(value):
        return None

    gene_id = str(value).strip().upper()

    if not gene_id:
        return None

    if gene_id.lower() in {"nan", "none", "na", "null"}:
        return None

    # Remove optional isoform suffix, e.g. AT1G01010.1 -> AT1G01010
    gene_id = gene_id.split(".", 1)[0]

    return gene_id


def make_dataset_col_name(raw_col: str, existing_cols: set[str]) -> str:
    """
    Convert one SRX column name into a clean TE column name.

    Usually:
        SRX5164421 -> TE_SRX5164421

    If a duplicate column name somehow appears across files, append a numeric suffix.
    """
    base = f"TE_{str(raw_col).strip()}"
    col = base

    i = 2
    while col in existing_cols:
        col = f"{base}_{i}"
        i += 1

    existing_cols.add(col)
    return col


def read_te_file(
    path: str,
    existing_dataset_cols: set[str],
) -> pd.DataFrame:
    """
    Read one translation-efficiency file.

    Input format:
        Gene_ID SRX... SRX... SRX...

    Each SRX column is treated as an independent TE dataset.

    Returns:
        DataFrame with one row per Gene_ID and TE_<SRX> columns.
    """
    path = str(path)

    # Files are whitespace/tab-delimited. "NA" is converted to NaN.
    df = pd.read_csv(
        path,
        sep=r"\s+",
        engine="python",
        na_values=["NA", "NaN", "nan", "N/A", ""],
        keep_default_na=True,
    )

    df.columns = [str(col).strip() for col in df.columns]

    if GENE_ID_COL not in df.columns:
        raise ValueError(
            f"{path} is missing required column '{GENE_ID_COL}'. "
            f"Observed columns: {list(df.columns)}"
        )

    raw_dataset_cols = [col for col in df.columns if col != GENE_ID_COL]

    if not raw_dataset_cols:
        raise ValueError(f"{path} contains no TE dataset columns after '{GENE_ID_COL}'.")

    n_raw_rows = len(df)

    df[GENE_ID_COL] = df[GENE_ID_COL].apply(clean_gene_id)
    df = df[df[GENE_ID_COL].notna()].copy()

    rename_map = {}

    for raw_col in raw_dataset_cols:
        clean_col = make_dataset_col_name(raw_col, existing_dataset_cols)
        rename_map[raw_col] = clean_col

    df = df[[GENE_ID_COL] + raw_dataset_cols].rename(columns=rename_map)

    te_cols = list(rename_map.values())

    # Convert values:
    #   0.0 stays 0.0
    #   NA becomes NaN
    #   nonnumeric strings become NaN
    for col in te_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # If a Gene_ID appears more than once within a file, collapse per TE column
    # using the median to avoid later merge expansion.
    n_before_collapse = len(df)

    df = (
        df.groupby(GENE_ID_COL, as_index=False)[te_cols]
        .median(numeric_only=True)
        .reset_index(drop=True)
    )

    n_after_collapse = len(df)

    n_duplicate_gene_rows = n_before_collapse - n_after_collapse

    n_values = df[te_cols].notna().sum().sum()
    n_genes_with_any_value = df[te_cols].notna().any(axis=1).sum()

    print(f"\n{path}")
    print(f"  Raw rows: {n_raw_rows:,}")
    print(f"  TE dataset columns: {len(te_cols):,}")
    print(f"  Duplicate Gene_ID rows collapsed within file: {n_duplicate_gene_rows:,}")
    print(f"  Unique Gene_ID values after cleaning: {df[GENE_ID_COL].nunique():,}")
    print(f"  Gene_ID values with at least one TE value: {n_genes_with_any_value:,}")
    print(f"  Total non-NA TE values: {n_values:,}")

    return df


def build_harmonized_translation_efficiency(
    input_files: List[str],
) -> pd.DataFrame:
    """
    Combine all translation-efficiency files into one Gene_ID-level table.

    For each Gene_ID:
        - each TE_<SRX> column is an independent dataset value
        - 0.0 values are retained as true zeroes
        - NA values are represented as NaN
        - if one non-NA value exists, harmonized value is that value
        - if multiple non-NA values exist, harmonized value is the median
        - standard deviation is reported only where 2 or more non-NA values exist
    """
    if len(input_files) == 0:
        raise ValueError("No input TE files were provided.")

    existing_dataset_cols = set()
    te_dfs = []

    for path in input_files:
        te_dfs.append(
            read_te_file(
                path=path,
                existing_dataset_cols=existing_dataset_cols,
            )
        )

    harmonized_df = te_dfs[0]

    for next_df in te_dfs[1:]:
        harmonized_df = harmonized_df.merge(
            next_df,
            on=GENE_ID_COL,
            how="outer",
            validate="one_to_one",
        )

    te_cols = [
        col for col in harmonized_df.columns
        if col.startswith("TE_")
    ]

    if not te_cols:
        raise ValueError("No TE columns were found after reading input files.")

    harmonized_df["n_translation_efficiency_values"] = (
        harmonized_df[te_cols].notna().sum(axis=1)
    )

    harmonized_df["harmonized_translation_efficiency"] = (
        harmonized_df[te_cols].median(axis=1, skipna=True)
    )

    has_multiple_values = harmonized_df["n_translation_efficiency_values"].ge(2)

    harmonized_df["translation_efficiency_std"] = np.where(
        has_multiple_values,
        harmonized_df[te_cols].std(axis=1, skipna=True, ddof=1),
        np.nan,
    )

    harmonized_df["translation_efficiency_range"] = np.where(
        has_multiple_values,
        harmonized_df[te_cols].max(axis=1, skipna=True)
        - harmonized_df[te_cols].min(axis=1, skipna=True),
        np.nan,
    )

    harmonized_df = harmonized_df.sort_values(GENE_ID_COL).reset_index(drop=True)

    print("\nHarmonized translation-efficiency summary:")
    print(f"  Total Gene_ID values across TE files: {len(harmonized_df):,}")
    print(f"  Independent TE dataset columns: {len(te_cols):,}")
    print(
        "  Gene_ID values with at least one TE value: "
        f"{harmonized_df['n_translation_efficiency_values'].gt(0).sum():,}"
    )
    print(
        "  Gene_ID values with two or more TE values: "
        f"{harmonized_df['n_translation_efficiency_values'].ge(2).sum():,}"
    )
    print(
        "  Gene_ID values with TE standard deviation: "
        f"{harmonized_df['translation_efficiency_std'].notna().sum():,}"
    )

    return harmonized_df


def add_translation_efficiency(
    nodes_df: pd.DataFrame,
    harmonized_te_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge harmonized translation-efficiency information into nodes_df without
    changing row count.
    """
    if NODE_COL not in nodes_df.columns:
        raise ValueError(f"nodes_df must contain a '{NODE_COL}' column.")

    n_before = len(nodes_df)

    nodes_df = nodes_df.copy()
    nodes_df[NODE_COL] = nodes_df[NODE_COL].astype(str).str.strip().str.upper()

    harmonized_te_df = harmonized_te_df.copy()
    harmonized_te_df[GENE_ID_COL] = (
        harmonized_te_df[GENE_ID_COL].astype(str).str.strip().str.upper()
    )

    nodes_df = nodes_df.merge(
        harmonized_te_df,
        how="left",
        left_on=NODE_COL,
        right_on=GENE_ID_COL,
        validate="one_to_one",
    )

    if GENE_ID_COL in nodes_df.columns:
        nodes_df = nodes_df.drop(columns=[GENE_ID_COL])

    n_after = len(nodes_df)

    if n_after != n_before:
        raise ValueError(
            f"Translation-efficiency merge changed row count: "
            f"{n_before:,} -> {n_after:,}. "
            "This indicates duplicate node IDs or duplicate Gene_ID values."
        )

    print("\nFinal node annotation summary:")
    print(f"  Input/output node rows: {n_after:,}")
    print(
        "  Nodes with harmonized translation efficiency: "
        f"{nodes_df['harmonized_translation_efficiency'].notna().sum():,}"
    )
    print(
        "  Nodes with TE values from two or more datasets: "
        f"{nodes_df['n_translation_efficiency_values'].ge(2).sum():,}"
    )

    return nodes_df


def main():
    parser = argparse.ArgumentParser(
        description="Add harmonized translation-efficiency information to network nodes."
    )

    parser.add_argument("--nodes", required=True, help="Input nodes pickle file")

    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where results will be saved",
    )

    parser.add_argument("--output_prefix", default="0", help="Prefix for output file")
    parser.add_argument("--output_suffix", default="step10", help="Suffix for output file")
    parser.add_argument("--organism_tag", required=True, help="Organism label for this run")

    parser.add_argument(
        "--input_files",
        nargs="+",
        required=True,
        help="Translation-efficiency files to harmonize.",
    )

    parser.add_argument(
        "--harmonized_translation_efficiency_out",
        required=True,
        help="Path to write harmonized translation-efficiency CSV for inspection.",
    )

    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)

    if nodes_df[NODE_COL].duplicated().any():
        duplicated = nodes_df.loc[
            nodes_df[NODE_COL].duplicated(),
            NODE_COL,
        ].head(10).tolist()

        raise ValueError(
            "Input nodes file contains duplicate node IDs. "
            f"Examples: {duplicated}"
        )

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    Path(args.harmonized_translation_efficiency_out).parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    nodes_df = nodes_df.copy()
    nodes_df[NODE_COL] = nodes_df[NODE_COL].astype(str).str.strip().str.upper()

    harmonized_te_df = build_harmonized_translation_efficiency(
        input_files=args.input_files,
    )

    harmonized_te_df.to_csv(
        args.harmonized_translation_efficiency_out,
        index=False,
    )

    nodes_df = add_translation_efficiency(
        nodes_df=nodes_df,
        harmonized_te_df=harmonized_te_df,
    )

    output_pkl = (
        f"{args.output_dir}/"
        f"{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )

    nodes_df.to_pickle(output_pkl)

    print(f"\nWrote annotated nodes file: {output_pkl}")
    print(
        "Wrote harmonized translation-efficiency table: "
        f"{args.harmonized_translation_efficiency_out}"
    )


if __name__ == "__main__":
    main()
