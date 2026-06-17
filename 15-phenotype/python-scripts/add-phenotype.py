#!/usr/bin/env python3

import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd

# Required because upstream node pickle files may contain pint/pint-pandas objects.
import pint
import pint_pandas


NODE_COL = "node"

PHENO_AGI_COL = "AGI"
PHENO_GROUP_COL = "Phenotypic functional group"

LETHAL_AGI_COL = "AGI"
LETHAL_COL = "is_lethal"


def normalize_agi(value) -> str | None:
    """
    Normalize Arabidopsis AGI/locus identifiers to node-level format.

    Examples:
        At1g01010   -> AT1G01010
        AT1G01010.1 -> AT1G01010
    """
    if value is None or pd.isna(value):
        return None

    value = str(value).strip().upper()

    if not value or value.lower() in {"nan", "none", "na", "null"}:
        return None

    # Remove optional isoform suffix.
    value = re.sub(r"\.\d+$", "", value)

    return value


def validate_required_columns(
    df: pd.DataFrame,
    required_cols: list[str],
    file_label: str,
) -> None:
    """
    Check that an input DataFrame contains expected columns.
    """
    missing = [col for col in required_cols if col not in df.columns]

    if missing:
        raise ValueError(
            f"{file_label} is missing required column(s): {missing}\n"
            f"Observed columns: {list(df.columns)}"
        )


def read_phenotype_data(path: str) -> pd.DataFrame:
    """
    Read phenotype functional group table.

    Expected columns:
        AGI
        Phenotypic functional group

    Output columns:
        node
        phenotypic_functional_group
        n_phenotypic_functional_groups

    If an AGI maps to multiple functional groups, groups are joined with ';'.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [col.strip() for col in df.columns]

    validate_required_columns(
        df=df,
        required_cols=[PHENO_AGI_COL, PHENO_GROUP_COL],
        file_label=path,
    )

    df = df[[PHENO_AGI_COL, PHENO_GROUP_COL]].copy()

    df[NODE_COL] = df[PHENO_AGI_COL].apply(normalize_agi)

    df["phenotypic_functional_group"] = (
        df[PHENO_GROUP_COL]
        .astype(str)
        .str.strip()
        .replace({"": np.nan, "nan": np.nan, "None": np.nan, "NA": np.nan})
    )

    df = df[df[NODE_COL].notna()].copy()
    df = df.dropna(subset=["phenotypic_functional_group"])
    df = df[[NODE_COL, "phenotypic_functional_group"]].drop_duplicates()

    if df.empty:
        print(f"\nPhenotype table {path}: no usable phenotype rows after cleaning.")
        return pd.DataFrame(
            columns=[
                NODE_COL,
                "phenotypic_functional_group",
                "n_phenotypic_functional_groups",
            ]
        )

    collapsed = (
        df.groupby(NODE_COL, as_index=False)
        .agg(
            phenotypic_functional_group=(
                "phenotypic_functional_group",
                lambda x: ";".join(sorted(set(x))),
            ),
            n_phenotypic_functional_groups=(
                "phenotypic_functional_group",
                lambda x: len(set(x)),
            ),
        )
    )

    print("\nPhenotype data summary:")
    print(f"  Raw phenotype rows: {len(df):,}")
    print(f"  Unique AGIs with phenotype annotation: {collapsed[NODE_COL].nunique():,}")
    print(
        "  AGIs with multiple phenotype functional groups: "
        f"{collapsed['n_phenotypic_functional_groups'].gt(1).sum():,}"
    )

    return collapsed


def normalize_lethal_value(value) -> str | float:
    """
    Normalize lethal-status values while preserving a user-readable Yes/No column.

    Returns:
        'Yes'
        'No'
        np.nan
    """
    if value is None or pd.isna(value):
        return np.nan

    value = str(value).strip()

    if not value:
        return np.nan

    value_lower = value.lower()

    if value_lower in {"yes", "y", "true", "1", "lethal"}:
        return "Yes"

    if value_lower in {"no", "n", "false", "0", "nonlethal", "non-lethal"}:
        return "No"

    return value


def read_lethal_data(path: str) -> pd.DataFrame:
    """
    Read lethal-status table.

    Expected columns:
        AGI
        is_lethal

    Output columns:
        node
        is_lethal
        is_lethal_bool

    Conflicting duplicate AGI statuses are treated as an error.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [col.strip() for col in df.columns]

    validate_required_columns(
        df=df,
        required_cols=[LETHAL_AGI_COL, LETHAL_COL],
        file_label=path,
    )

    df = df[[LETHAL_AGI_COL, LETHAL_COL]].copy()

    df[NODE_COL] = df[LETHAL_AGI_COL].apply(normalize_agi)
    df["is_lethal"] = df[LETHAL_COL].apply(normalize_lethal_value)

    df = df[df[NODE_COL].notna()].copy()
    df = df.dropna(subset=["is_lethal"])
    df = df[[NODE_COL, "is_lethal"]].drop_duplicates()

    if df.empty:
        print(f"\nLethal-status table {path}: no usable lethal-status rows after cleaning.")
        return pd.DataFrame(columns=[NODE_COL, "is_lethal", "is_lethal_bool"])

    status_counts = df.groupby(NODE_COL)["is_lethal"].nunique()
    conflicting_nodes = status_counts[status_counts.gt(1)].index.tolist()

    if conflicting_nodes:
        examples = df[df[NODE_COL].isin(conflicting_nodes)].sort_values(NODE_COL).head(30)
        raise ValueError(
            "Conflicting lethal-status annotations detected for one or more AGIs.\n"
            "This should be resolved before merging.\n\n"
            f"Example conflicting rows:\n{examples.to_string(index=False)}"
        )

    collapsed = df.drop_duplicates(subset=NODE_COL, keep="first").reset_index(drop=True)

    collapsed["is_lethal_bool"] = collapsed["is_lethal"].map(
        {
            "Yes": True,
            "No": False,
        }
    )

    print("\nLethal-status data summary:")
    print(f"  Unique AGIs with lethal-status annotation: {collapsed[NODE_COL].nunique():,}")
    print(f"  Lethal genes: {collapsed['is_lethal_bool'].eq(True).sum():,}")
    print(f"  Non-lethal genes: {collapsed['is_lethal_bool'].eq(False).sum():,}")

    return collapsed


def add_phenotype_annotations(
    nodes_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    lethal_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge phenotype and lethal-status annotations onto nodes_df without changing row count.
    """
    if NODE_COL not in nodes_df.columns:
        raise ValueError(f"Input nodes file must contain column '{NODE_COL}'.")

    nodes_df = nodes_df.copy()
    nodes_df[NODE_COL] = nodes_df[NODE_COL].apply(normalize_agi)

    if nodes_df[NODE_COL].isna().any():
        bad = nodes_df[nodes_df[NODE_COL].isna()].head(10)
        raise ValueError(
            "Input nodes file contains rows with unusable node identifiers.\n"
            f"Examples:\n{bad.to_string(index=False)}"
        )

    if nodes_df[NODE_COL].duplicated().any():
        examples = nodes_df.loc[
            nodes_df[NODE_COL].duplicated(),
            NODE_COL,
        ].head(10).tolist()

        raise ValueError(
            "Input nodes file contains duplicate node IDs after normalization. "
            f"Examples: {examples}"
        )

    n_before = len(nodes_df)

    nodes_df = nodes_df.merge(
        phenotype_df,
        on=NODE_COL,
        how="left",
        validate="one_to_one",
    )

    nodes_df = nodes_df.merge(
        lethal_df,
        on=NODE_COL,
        how="left",
        validate="one_to_one",
    )

    n_after = len(nodes_df)

    if n_after != n_before:
        raise ValueError(
            f"Phenotype merge changed row count: {n_before:,} -> {n_after:,}. "
            "This indicates duplicate keys in one of the merged tables."
        )

    print("\nFinal node annotation summary:")
    print(f"  Input/output node rows: {n_after:,}")
    print(
        "  Nodes with phenotype functional group: "
        f"{nodes_df['phenotypic_functional_group'].notna().sum():,}"
    )
    print(
        "  Nodes with lethal-status annotation: "
        f"{nodes_df['is_lethal'].notna().sum():,}"
    )
    print(
        "  Nodes annotated as lethal: "
        f"{nodes_df['is_lethal_bool'].eq(True).sum():,}"
    )

    return nodes_df


def main():
    parser = argparse.ArgumentParser(
        description="Add phenotype functional group and lethal-status annotations to nodes."
    )

    parser.add_argument("--input_nodes", required=True, help="Input nodes pickle file")
    parser.add_argument("--phenotype_data", required=True, help="Phenotype TSV file")
    parser.add_argument("--lethal_data", required=True, help="Lethal-status TSV file")

    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where results will be saved",
    )

    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_suffix", required=True)
    parser.add_argument("--organism_tag", required=True)

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    nodes_df = pd.read_pickle(args.input_nodes)

    phenotype_df = read_phenotype_data(args.phenotype_data)
    lethal_df = read_lethal_data(args.lethal_data)

    nodes_df = add_phenotype_annotations(
        nodes_df=nodes_df,
        phenotype_df=phenotype_df,
        lethal_df=lethal_df,
    )

    output_pkl = (
        output_dir
        / f"{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )

    nodes_df.to_pickle(output_pkl)

    print(f"\nWrote annotated nodes file: {output_pkl}")


if __name__ == "__main__":
    main()
