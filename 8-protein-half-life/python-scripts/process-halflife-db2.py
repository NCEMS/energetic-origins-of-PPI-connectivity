#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = [
    "AGI",
    "Tissue",
    "Fraction",
    "Diff. log2k",
    "Fold change in k",
    "ctrl.mean",
    "30C.mean",
    "p.value"
]


VALUE_COLUMNS = [
    "Diff. log2k",
    "Fold change in k",
    "ctrl.mean",
    "30C.mean",
]


RAW_TO_CLEAN_METRIC = {
    "Diff. log2k": "diff_log2k",
    "Fold change in k": "k_fold_change",
}


def clean_agi_string(value) -> list[str]:
    """
    Convert a raw AGI field into a list of clean uppercase base AGI IDs.

    Behavior:
        "At1g56070; At1g56075"
            -> ["AT1G56070", "AT1G56075"]

        "At1g01010.1; At1g01010.2"
            -> ["AT1G01010"]

        "At1g01010.1; At1g02020.2"
            -> ["AT1G01010", "AT1G02020"]

    This function:
        - splits distinct AGIs into separate IDs,
        - uppercases all IDs,
        - removes isoform suffixes such as .1, .2, .3,
        - removes duplicate AGIs within a single source row.
    """
    if pd.isna(value):
        return []

    value = str(value).strip()
    if not value:
        return []

    # Tolerate semicolon- or comma-delimited AGI lists.
    parts = re.split(r"[;,]", value)

    clean_parts = []
    seen = set()

    for part in parts:
        agi = part.strip().upper()

        if not agi:
            continue

        # Remove isoform suffixes such as .1, .2, .3.
        agi = re.sub(r"\.\d+$", "", agi)

        if agi and agi not in seen:
            clean_parts.append(agi)
            seen.add(agi)

    return clean_parts


def clean_label(value, fallback: str) -> str:
    """
    Convert tissue/fraction labels into safe column-name components.

    Examples:
        "Root" -> "root"
        "30 °C" -> "30_c"
        "100kg" -> "100kg"
        "membrane fraction" -> "membrane_fraction"
    """
    if pd.isna(value):
        return fallback

    label = str(value).strip()

    if not label:
        return fallback

    label = label.lower()

    # Make common symbols readable.
    label = label.replace("°", "")
    label = label.replace("%", "pct")

    # Replace whitespace with underscores.
    label = re.sub(r"\s+", "_", label)

    # Replace all remaining non-column-friendly characters.
    label = re.sub(r"[^a-z0-9._+-]+", "_", label)

    # Collapse repeated underscores and trim.
    label = re.sub(r"_+", "_", label).strip("_")

    if not label:
        return fallback

    return label


def validate_columns(df: pd.DataFrame, input_file: str) -> None:
    """
    Ensure all required columns are present in the input table.
    """
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]

    if missing:
        raise ValueError(
            f"Missing required columns in {input_file}: {missing}\n"
            f"Observed columns: {list(df.columns)}"
        )


def handle_duplicate_agi_condition_rows(
    df: pd.DataFrame,
    duplicate_policy: str,
) -> pd.DataFrame:
    """
    Handle duplicated AGI + Tissue + Fraction rows after AGI splitting.

    Different Tissue/Fraction combinations for the same AGI are expected and
    are preserved. This function only handles cases where the same AGI appears
    multiple times for the same Tissue/Fraction condition.

    duplicate_policy:
        error:
            Fail on duplicated AGI + Tissue + Fraction rows.

        first:
            Keep the first row for each AGI + Tissue + Fraction.

        mean:
            Average numeric values across duplicated AGI + Tissue + Fraction rows.
    """
    duplicate_key = ["AGI", "Tissue_label", "Fraction_label"]

    if not df.duplicated(subset=duplicate_key).any():
        return df.reset_index(drop=True)

    duplicate_rows = df[df.duplicated(subset=duplicate_key, keep=False)].sort_values(
        duplicate_key
    )

    if duplicate_policy == "error":
        example = duplicate_rows.head(30).to_string(index=False)
        raise ValueError(
            "Duplicate AGI + Tissue + Fraction values detected after splitting "
            "multi-AGI entries.\n"
            "Different Tissue/Fraction combinations for the same AGI are allowed, "
            "but repeated rows for the same AGI and same Tissue/Fraction require "
            "a policy decision.\n"
            "Inspect the source table or rerun with --duplicate_policy first or mean.\n\n"
            f"Example duplicate rows:\n{example}"
        )

    if duplicate_policy == "first":
        return df.drop_duplicates(subset=duplicate_key, keep="first").reset_index(
            drop=True
        )

    if duplicate_policy == "mean":
        return (
            df.groupby(duplicate_key, as_index=False)[VALUE_COLUMNS]
            .mean(numeric_only=True)
            .reset_index(drop=True)
        )

    if duplicate_policy == "median":
        grouped = (
            df.groupby(duplicate_key, as_index=False)
            .agg(
                {
                    "Diff. log2k": "median",
                }
            )
            .reset_index(drop=True)
        )

    grouped["Fold change in k"] = 2 ** grouped["Diff. log2k"]

    return grouped

    raise ValueError(f"Unsupported duplicate_policy: {duplicate_policy}")


def pivot_tissue_fraction_specific_columns(
    df: pd.DataFrame,
    source_prefix: str,
) -> pd.DataFrame:
    """
    Convert long AGI/Tissue/Fraction rows into one row per AGI with
    condition-specific columns.

    Example output columns:
        AGI
        Fan2024_root_100kg_diff_log2k
        Fan2024_root_100kg_k_fold_change
        Fan2024_shoot_100kg_diff_log2k
        Fan2024_shoot_100kg_k_fold_change
    """
    agis = (
        df[["AGI"]]
        .drop_duplicates()
        .sort_values("AGI")
        .reset_index(drop=True)
    )

    wide_df = agis.copy()

    condition_cols = ["Tissue_label", "Fraction_label"]
    conditions = (
        df[condition_cols]
        .drop_duplicates()
        .sort_values(condition_cols)
        .itertuples(index=False, name=None)
    )

    for tissue_label, fraction_label in conditions:
        condition_df = df[
            (df["Tissue_label"] == tissue_label)
            & (df["Fraction_label"] == fraction_label)
        ].copy()

        keep_cols = ["AGI"] + VALUE_COLUMNS
        condition_df = condition_df[keep_cols]

        condition_prefix = f"{source_prefix}_{tissue_label}_{fraction_label}"

        rename_map = {
            col: f"{condition_prefix}_{RAW_TO_CLEAN_METRIC[col]}"
            for col in VALUE_COLUMNS
        }

        condition_df = condition_df.rename(columns=rename_map)

        wide_df = wide_df.merge(condition_df, how="left", on="AGI")

    return wide_df


def process_db2(
    input_file: str,
    output_file: str,
    source_prefix: str,
    duplicate_policy: str,
) -> None:
    """
    Preprocess Fan2024 turnover-response data.

    Input columns expected include:
        AGI
        Tissue
        Fraction
        Diff. log2k
        Fold change in k

    Output:
        One row per uppercase base AGI.

    Output columns:
        AGI
        <source_prefix>_<tissue>_<fraction>_diff_log2k
        <source_prefix>_<tissue>_<fraction>_k_fold_change

    Multi-AGI rows are exploded first.
    Different Tissue/Fraction rows for the same AGI are pivoted into separate columns.
    """
    input_path = Path(input_file)
    output_path = Path(output_file)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # The provided example is tab-delimited.
    df = pd.read_csv(input_path, sep="\t", dtype=str)

    # Clean accidental whitespace from column names.
    # Important because the header can contain names such as "Diff. log2k ".
    df.columns = [col.strip() for col in df.columns]

    validate_columns(df, input_file)

    # Retain only the AGI key, condition labels, and requested measurement columns.
    df = df[REQUIRED_COLUMNS].copy()

    # Clean condition labels for use in output column names.
    df["Tissue"] = df["Tissue"].astype(str).str.strip()
    df["Fraction"] = df["Fraction"].astype(str).str.strip()

    df["Tissue_label"] = df["Tissue"].apply(lambda x: clean_label(x, "unknown_tissue"))
    df["Fraction_label"] = df["Fraction"].apply(
        lambda x: clean_label(x, "unknown_fraction")
    )

    # Split multi-AGI fields into one AGI per row.
    df["AGI"] = df["AGI"].apply(clean_agi_string)
    df = df.explode("AGI", ignore_index=True)

    # Remove rows with no usable AGI.
    df = df[df["AGI"].notna() & (df["AGI"].astype(str).str.len() > 0)].copy()

    # Convert retained measurement columns to numeric.
    for col in VALUE_COLUMNS:
        df[col] = pd.to_numeric(df[col].astype(str).str.strip(), errors="coerce")

    # Drop exact duplicate rows before checking true condition-level duplicates.
    df = df.drop_duplicates().reset_index(drop=True)

    # Different Tissue/Fraction combinations for the same AGI are allowed.
    # Only repeated AGI + same Tissue + same Fraction is considered a duplicate problem.
    df = handle_duplicate_agi_condition_rows(
        df=df,
        duplicate_policy=duplicate_policy,
    )

    # Pivot to one row per AGI with tissue/fraction-specific columns.
    wide_df = pivot_tissue_fraction_specific_columns(
        df=df,
        source_prefix=source_prefix,
    )

    wide_df.to_csv(output_path, index=False)

    print(f"Wrote processed DB2 table: {output_path}")
    print(f"Rows written: {len(wide_df):,}")
    print(f"Unique AGIs: {wide_df['AGI'].nunique():,}")

    observed_conditions = (
        df[["Tissue_label", "Fraction_label"]]
        .drop_duplicates()
        .sort_values(["Tissue_label", "Fraction_label"])
    )
    print("Tissue/Fraction conditions observed:")
    for row in observed_conditions.itertuples(index=False):
        print(f"  {row.Tissue_label} / {row.Fraction_label}")


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess Fan2024 protein turnover-response data."
    )

    parser.add_argument(
        "--input_file",
        required=True,
        help="Input Fan2024 table.",
    )

    parser.add_argument(
        "--output_file",
        required=True,
        help="Output processed CSV file.",
    )

    parser.add_argument(
        "--source_prefix",
        default="Fan2024",
        help="Prefix to add to retained source columns. Default: Fan2024",
    )

    parser.add_argument(
        "--duplicate_policy",
        choices=["error", "first", "mean", "median"],
        default="error",
        help=(
            "How to handle repeated rows for the same AGI, Tissue, and Fraction. "
            "Different Tissue/Fraction combinations for the same AGI are always "
            "retained as separate columns. Default: error."
        ),
    )

    args = parser.parse_args()

    process_db2(
        input_file=args.input_file,
        output_file=args.output_file,
        source_prefix=args.source_prefix,
        duplicate_policy=args.duplicate_policy,
    )


if __name__ == "__main__":
    main()
