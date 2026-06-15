import argparse
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd


ABUNDANCE_ID_COL = "string_external_id"
ABUNDANCE_VALUE_COL = "abundance"
NODE_COL = "node"
NODE_UNIPROT_COL = "UniProtKB-AC"


def parse_uniprot_from_string_external_id(value) -> Optional[str]:
    """
    Parse UniProt accession from STRING-style external ID.

    Example:
        3702.O03042 -> O03042
    """
    if value is None or pd.isna(value):
        return None

    value = str(value).strip()

    if not value:
        return None

    if "." not in value:
        return value.upper()

    return value.rsplit(".", 1)[1].strip().upper()


def parse_uniprot_accessions(value) -> list[str]:
    """
    Parse UniProtKB-AC field from nodes_df.

    Handles:
        Q9SP32
        F4HQG6;Q56Y50;Q9SP32
        NaN
    """
    if value is None or pd.isna(value):
        return []

    accessions = []

    for accession in str(value).split(";"):
        accession = accession.strip().upper()

        if not accession:
            continue

        if accession.lower() in {"nan", "none", "na", "null"}:
            continue

        accessions.append(accession)

    return accessions


def build_uniprot_to_node_map(nodes_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a long-form mapping table from UniProt accession to TAIR node.

    Returns:
        DataFrame with columns:
            node
            UniProtKB_AC_single
    """

    required_cols = {NODE_COL, NODE_UNIPROT_COL}
    missing_cols = required_cols - set(nodes_df.columns)

    if missing_cols:
        raise ValueError(
            f"nodes_df is missing required column(s): {', '.join(sorted(missing_cols))}"
        )

    mapping_rows = []

    for _, row in nodes_df[[NODE_COL, NODE_UNIPROT_COL]].iterrows():
        node = str(row[NODE_COL]).strip().upper()
        accessions = parse_uniprot_accessions(row[NODE_UNIPROT_COL])

        for accession in accessions:
            mapping_rows.append(
                {
                    NODE_COL: node,
                    "UniProtKB_AC_single": accession,
                }
            )

    mapping_df = pd.DataFrame(mapping_rows)

    if mapping_df.empty:
        raise ValueError(
            "No usable UniProtKB-AC mappings were found in nodes_df. "
            "Check the UniProtKB-AC column."
        )

    mapping_df = mapping_df.drop_duplicates()

    n_nodes_with_uniprot = mapping_df[NODE_COL].nunique()
    n_total_nodes = nodes_df[NODE_COL].nunique()
    n_nodes_multi_uniprot = (
        mapping_df.groupby(NODE_COL)["UniProtKB_AC_single"]
        .nunique()
        .gt(1)
        .sum()
    )

    print(f"Total unique nodes in nodes_df: {n_total_nodes:,}")
    print(f"Nodes with at least one UniProtKB-AC mapping: {n_nodes_with_uniprot:,}")
    print(f"Nodes with multiple UniProtKB-AC mappings: {n_nodes_multi_uniprot:,}")
    print(f"Unique UniProt accessions in node mapping: {mapping_df['UniProtKB_AC_single'].nunique():,}")

    return mapping_df


def read_abundance_file(
    path: str,
    file_index: int,
    uniprot_to_node_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Read one abundance file, map UniProt IDs to TAIR nodes, and collapse to one
    abundance value per TAIR node for that file.

    If multiple UniProt IDs map to the same TAIR node and have abundance values
    in the same file, the median is used.
    """

    value_col = f"abundance_file_{file_index}"

    # New files are whitespace-delimited.
    df = pd.read_csv(path, sep=r"\s+", engine="python")

    required_cols = {ABUNDANCE_ID_COL, ABUNDANCE_VALUE_COL}
    missing_cols = required_cols - set(df.columns)

    if missing_cols:
        raise ValueError(
            f"{path} is missing required column(s): "
            f"{', '.join(sorted(missing_cols))}. "
            f"Observed columns: {list(df.columns)}"
        )

    n_raw_rows = len(df)

    df = df[[ABUNDANCE_ID_COL, ABUNDANCE_VALUE_COL]].copy()

    df["UniProtKB_AC_single"] = df[ABUNDANCE_ID_COL].apply(
        parse_uniprot_from_string_external_id
    )

    df[value_col] = pd.to_numeric(df[ABUNDANCE_VALUE_COL], errors="coerce")

    df = df.dropna(subset=["UniProtKB_AC_single", value_col])

    # If the same UniProt appears more than once in one abundance file,
    # collapse first to avoid merge expansion.
    n_before_uniprot_collapse = len(df)

    df = (
        df.groupby("UniProtKB_AC_single", as_index=False)[value_col]
        .median()
    )

    n_after_uniprot_collapse = len(df)

    n_duplicate_uniprot_rows = n_before_uniprot_collapse - n_after_uniprot_collapse

    # Map UniProt accession to TAIR node.
    mapped = df.merge(
        uniprot_to_node_df,
        on="UniProtKB_AC_single",
        how="inner",
        validate="one_to_many",
    )

    n_mapped_uniprot = mapped["UniProtKB_AC_single"].nunique()
    n_mapped_nodes = mapped[NODE_COL].nunique()

    # Count nodes receiving abundance values through multiple UniProt IDs.
    n_nodes_multi_uniprot_abundance = (
        mapped.groupby(NODE_COL)["UniProtKB_AC_single"]
        .nunique()
        .gt(1)
        .sum()
    )

    # Collapse to one abundance value per TAIR node for this file.
    node_abundance = (
        mapped.groupby(NODE_COL, as_index=False)
        .agg(
            **{
                value_col: (value_col, "median"),
                f"n_uniprot_abundance_values_file_{file_index}": (
                    "UniProtKB_AC_single",
                    "nunique",
                ),
                f"matched_uniprot_accessions_file_{file_index}": (
                    "UniProtKB_AC_single",
                    lambda x: ";".join(sorted(set(x))),
                ),
            }
        )
    )

    print(f"\n{path}")
    print(f"  Raw abundance rows: {n_raw_rows:,}")
    print(f"  Duplicate UniProt rows collapsed within file: {n_duplicate_uniprot_rows:,}")
    print(f"  Unique UniProt accessions with numeric abundance: {len(df):,}")
    print(f"  Unique UniProt accessions matched to nodes: {n_mapped_uniprot:,}")
    print(f"  Unique TAIR nodes matched to abundance: {n_mapped_nodes:,}")
    print(
        "  TAIR nodes with abundance from multiple UniProt IDs in this file: "
        f"{n_nodes_multi_uniprot_abundance:,}"
    )

    return node_abundance


def build_harmonized_abundance(
    abundance_files: List[str],
    uniprot_to_node_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Combine abundance values from multiple abundance files.

    For each TAIR node:
        - abundance_file_1 ... abundance_file_7 are per-file TAIR-level values
        - if one file has a value, harmonized_abundance is that value
        - if 2 or more files have values, harmonized_abundance is the median
    """

    if len(abundance_files) == 0:
        raise ValueError("No abundance files were provided.")

    abundance_dfs = []

    for i, path in enumerate(abundance_files, start=1):
        abundance_dfs.append(
            read_abundance_file(
                path=path,
                file_index=i,
                uniprot_to_node_df=uniprot_to_node_df,
            )
        )

    harmonized_df = abundance_dfs[0]

    for next_df in abundance_dfs[1:]:
        harmonized_df = harmonized_df.merge(
            next_df,
            on=NODE_COL,
            how="outer",
            validate="one_to_one",
        )

    abundance_cols = [
        col for col in harmonized_df.columns
        if col.startswith("abundance_file_")
    ]

    harmonized_df["n_abundance_files_with_values"] = harmonized_df[
        abundance_cols
    ].notna().sum(axis=1)

    harmonized_df["harmonized_abundance"] = harmonized_df[abundance_cols].median(
        axis=1,
        skipna=True,
    )

    harmonized_df = harmonized_df.sort_values(NODE_COL).reset_index(drop=True)

    print("\nHarmonized abundance summary:")
    print(f"  Total TAIR nodes across abundance files: {len(harmonized_df):,}")
    print(
        "  TAIR nodes with at least one abundance value: "
        f"{harmonized_df['n_abundance_files_with_values'].gt(0).sum():,}"
    )
    print(
        "  TAIR nodes with abundance values from two or more files: "
        f"{harmonized_df['n_abundance_files_with_values'].ge(2).sum():,}"
    )

    return harmonized_df


def add_expression(
    nodes_df: pd.DataFrame,
    harmonized_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge harmonized abundance information into nodes_df without changing row count.
    """

    if NODE_COL not in nodes_df.columns:
        raise ValueError(f"nodes_df must contain a '{NODE_COL}' column.")

    n_before = len(nodes_df)

    nodes_df = nodes_df.copy()
    nodes_df[NODE_COL] = nodes_df[NODE_COL].astype(str).str.strip().str.upper()

    nodes_df = nodes_df.merge(
        harmonized_df,
        on=NODE_COL,
        how="left",
        validate="one_to_one",
    )

    n_after = len(nodes_df)

    if n_after != n_before:
        raise ValueError(
            f"Abundance merge changed row count: {n_before:,} -> {n_after:,}. "
            "This indicates duplicate node IDs in either nodes_df or harmonized_df."
        )

    print("\nFinal node annotation summary:")
    print(f"  Input/output node rows: {n_after:,}")
    print(
        "  Nodes with harmonized abundance: "
        f"{nodes_df['harmonized_abundance'].notna().sum():,}"
    )

    return nodes_df


def main():
    parser = argparse.ArgumentParser(
        description="Add harmonized protein abundance information to network nodes."
    )

    parser.add_argument("--nodes", required=True, help="Input nodes pickle file")

    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where results will be saved",
    )

    parser.add_argument("--output_prefix", default="0", help="Prefix for output file")
    parser.add_argument("--output_suffix", default="0", help="Suffix for output file")
    parser.add_argument("--organism_tag", required=True, help="Organism label for this run")

    parser.add_argument(
        "--abundance_files",
        nargs="+",
        required=True,
        help="Individual abundance files, e.g. 1-abundance.txt ... 7-abundance.txt",
    )

    parser.add_argument(
        "--harmonized_abundance_out",
        required=True,
        help="Path to write harmonized abundance CSV for inspection",
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
    Path(args.harmonized_abundance_out).parent.mkdir(parents=True, exist_ok=True)

    nodes_df = nodes_df.copy()
    nodes_df[NODE_COL] = nodes_df[NODE_COL].astype(str).str.strip().str.upper()

    uniprot_to_node_df = build_uniprot_to_node_map(nodes_df)

    harmonized_df = build_harmonized_abundance(
        abundance_files=args.abundance_files,
        uniprot_to_node_df=uniprot_to_node_df,
    )

    harmonized_df.to_csv(args.harmonized_abundance_out, index=False)

    nodes_df = add_expression(nodes_df, harmonized_df)

    output_pkl = (
        f"{args.output_dir}/"
        f"{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )

    nodes_df.to_pickle(output_pkl)

    print(f"\nWrote annotated nodes file: {output_pkl}")
    print(f"Wrote harmonized abundance table: {args.harmonized_abundance_out}")


if __name__ == "__main__":
    main()
