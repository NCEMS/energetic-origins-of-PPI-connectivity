#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = [
    "AGI",
    "Average KD (d-1)",
]


VALUE_COLUMNS = [
    "Average KD (d-1)",
]


RAW_TO_CLEAN_METRIC = {
    "Average KD (d-1)": "k_deg_day_inverse",
}


def clean_agi_string(value) -> list[str]:
    """
    Convert a raw Li2017 AGI field into a list of clean uppercase base AGI IDs.

    Behavior:
        "AT1G01080.1,AT1G01080.2"
            -> ["AT1G01080"]

        "AT1G01080.1,AT1G01090.2"
            -> ["AT1G01080", "AT1G01090"]

        "At1g56070; At1g56075"
            -> ["AT1G56070", "AT1G56075"]

    This function:
        - splits comma- or semicolon-delimited identifiers,
        - uppercases all IDs,
        - removes isoform suffixes such as .1, .2, .3,
        - removes duplicate base AGIs within a single source row.
    """
    if pd.isna(value):
        return []

    value = str(value).strip()
    if not value:
        return []

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


def make_output_column_name(source_prefix: str, clean_metric_name: str) -> str:
    """
    Build clean output column names.

    Examples:
        source_prefix="Li2017_" -> Li2017_k_deg_day_inverse
        source_prefix="Li2017"  -> Li2017_k_deg_day_inverse
    """
    source_prefix = source_prefix.strip()

    if source_prefix.endswith("_"):
        return f"{source_prefix}{clean_metric_name}"

    return f"{source_prefix}_{clean_metric_name}"


def validate_columns(df: pd.DataFrame, input_file: str) -> None:
    """
    Ensure all required columns are present after column-name cleanup.
    """
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]

    if missing:
        raise ValueError(
            f"Missing required columns in {input_file}: {missing}\n"
            f"Observed columns: {list(df.columns)}"
        )


def handle_duplicate_agis(
    df: pd.DataFrame,
    duplicate_policy: str,
    value_columns: list[str],
) -> pd.DataFrame:
    """
    Ensure one row per AGI for clean downstream merging.

    duplicate_policy:
        error:
            Fail if the same AGI appears in multiple source rows.

        first:
            Keep the first row for each AGI.

        mean:
            Average numeric values across duplicate AGI rows.

        median:
            Take the median numeric value across duplicate AGI rows.
    """
    if not df["AGI"].duplicated().any():
        return df.reset_index(drop=True)

    duplicate_rows = df[df["AGI"].duplicated(keep=False)].sort_values("AGI")

    if duplicate_policy == "error":
        example = duplicate_rows.head(30).to_string(index=False)
        raise ValueError(
            "Duplicate AGI values detected after splitting multi-AGI entries.\n"
            "This would cause row multiplication during later merges.\n"
            "Inspect the source table or rerun with --duplicate_policy first, mean, or median.\n\n"
            f"Example duplicate rows:\n{example}"
        )

    if duplicate_policy == "first":
        return df.drop_duplicates(subset="AGI", keep="first").reset_index(drop=True)

    if duplicate_policy == "mean":
        return (
            df.groupby("AGI", as_index=False)[value_columns]
            .mean(numeric_only=True)
            .reset_index(drop=True)
        )

    if duplicate_policy == "median":
        return (
            df.groupby("AGI", as_index=False)[value_columns]
            .median(numeric_only=True)
            .reset_index(drop=True)
        )

    raise ValueError(f"Unsupported duplicate_policy: {duplicate_policy}")


def process_db3(
    input_file: str,
    output_file: str,
    source_prefix: str,
    duplicate_policy: str,
) -> None:
    """
    Preprocess Li2017 degradation-rate data.

    Input columns expected include:
        AGI
        Average KD (d-1)

    Output columns:
        AGI
        Li2017_k_deg_day_inverse

    Interpretation:
        Average KD (d-1) is treated as the source-reported degradation-rate
        constant in inverse days.
    """
    input_path = Path(input_file)
    output_path = Path(output_file)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Input example is tab-delimited.
    df = pd.read_csv(input_path, sep="\t", dtype=str)

    # Clean leading/trailing whitespace from column names.
    df.columns = [col.strip() for col in df.columns]

    validate_columns(df, input_file)

    df = df[REQUIRED_COLUMNS].copy()

    # Split multi-AGI fields and remove isoform suffixes.
    df["AGI"] = df["AGI"].apply(clean_agi_string)
    df = df.explode("AGI", ignore_index=True)

    # Remove rows with no usable AGI.
    df = df[df["AGI"].notna() & (df["AGI"].astype(str).str.len() > 0)].copy()

    # Convert retained measurement column to numeric.
    for col in VALUE_COLUMNS:
        df[col] = pd.to_numeric(df[col].astype(str).str.strip(), errors="coerce")

    # Add clean source-specific column name.
    rename_map = {
        col: make_output_column_name(source_prefix, RAW_TO_CLEAN_METRIC[col])
        for col in VALUE_COLUMNS
    }

    df = df.rename(columns=rename_map)

    prefixed_value_columns = [rename_map[col] for col in VALUE_COLUMNS]
    output_columns = ["AGI"] + prefixed_value_columns
    df = df[output_columns]

    # Drop exact duplicate rows first.
    df = df.drop_duplicates().reset_index(drop=True)

    # Ensure one row per AGI unless otherwise requested.
    df = handle_duplicate_agis(
        df=df,
        duplicate_policy=duplicate_policy,
        value_columns=prefixed_value_columns,
    )

    df.to_csv(output_path, index=False)

    print(f"Wrote processed DB3 table: {output_path}")
    print(f"Rows written: {len(df):,}")
    print(f"Unique AGIs: {df['AGI'].nunique():,}")
    print(f"Duplicate policy used: {duplicate_policy}")
    print(f"Output columns: {list(df.columns)}")


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess Li2017 protein degradation-rate data."
    )

    parser.add_argument(
        "--input_file",
        required=True,
        help="Input Li2017 table.",
    )

    parser.add_argument(
        "--output_file",
        required=True,
        help="Output processed CSV file.",
    )

    parser.add_argument(
        "--source_prefix",
        default="Li2017",
        help="Prefix to add to retained source columns. Default: Li2017",
    )

    parser.add_argument(
        "--duplicate_policy",
        choices=["error", "first", "mean", "median"],
        default="error",
        help=(
            "How to handle AGIs that appear in multiple source rows after splitting. "
            "Options: error, first, mean, median. Default: error."
        ),
    )

    args = parser.parse_args()

    process_db3(
        input_file=args.input_file,
        output_file=args.output_file,
        source_prefix=args.source_prefix,
        duplicate_policy=args.duplicate_policy,
    )


if __name__ == "__main__":
    main()
