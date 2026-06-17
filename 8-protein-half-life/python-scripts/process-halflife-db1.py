#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = [
    "Fraction",
    "AGI",
    "Mean of turnover rate (log2k)",
    "SD",
    "CV",
    "Half-life (hr-1)",
]


VALUE_COLUMNS = [
    "Mean of turnover rate (log2k)",
    "SD",
    "CV",
    "Half-life (hr-1)",
]


RAW_TO_CLEAN_METRIC = {
    "Mean of turnover rate (log2k)": "log2k_mean",
    "SD": "log2k_sd",
    "CV": "log2k_cv",
    "Half-life (hr-1)": "halflife_hr",
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

    parts = re.split(r"[;,]", value)

    clean_parts = []
    seen = set()

    for part in parts:
        agi = part.strip().upper()

        if not agi:
            continue

        # Remove isoform suffixes such as .1, .2, .3
        agi = re.sub(r"\.\d+$", "", agi)

        if agi and agi not in seen:
            clean_parts.append(agi)
            seen.add(agi)

    return clean_parts


def clean_fraction_label(value) -> str:
    """
    Convert a raw fraction label into a safe column-name component.

    Examples:
        "soluble" -> "soluble"
        "membrane fraction" -> "membrane_fraction"
        "100kg" -> "100kg"
    """
    if pd.isna(value):
        return "unknown_fraction"

    label = str(value).strip()

    if not label:
        return "unknown_fraction"

    # Make fraction labels safe for column names.
    label = re.sub(r"\s+", "_", label)
    label = re.sub(r"[^A-Za-z0-9._+-]+", "_", label)
    label = label.strip("_")

    if not label:
        return "unknown_fraction"

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


def handle_duplicate_agi_fraction_rows(
    df: pd.DataFrame,
    duplicate_policy: str,
) -> pd.DataFrame:
    """
    Handle duplicated AGI + Fraction rows after AGI splitting.

    Different fractions for the same AGI are expected and are preserved.
    This function only handles cases where the same AGI appears multiple
    times for the same fraction.

    duplicate_policy:
        error:
            Fail on duplicated AGI + Fraction rows.

        first:
            Keep the first row for each AGI + Fraction.

        mean:
            Average numeric values across duplicated AGI + Fraction rows.
    """
    duplicate_key = ["AGI", "Fraction_label"]

    if not df.duplicated(subset=duplicate_key).any():
        return df.reset_index(drop=True)

    duplicate_rows = df[df.duplicated(subset=duplicate_key, keep=False)].sort_values(
        duplicate_key
    )

    if duplicate_policy == "error":
        example = duplicate_rows.head(30).to_string(index=False)
        raise ValueError(
            "Duplicate AGI + Fraction values detected after splitting multi-AGI entries.\n"
            "Different fractions for the same AGI are allowed, but repeated rows for the\n"
            "same AGI and same fraction require a policy decision.\n"
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
        return (
            df.groupby("AGI", as_index=False)[value_columns]
            .median(numeric_only=True)
            .reset_index(drop=True)
        )

    raise ValueError(f"Unsupported duplicate_policy: {duplicate_policy}")


def pivot_fraction_specific_columns(
    df: pd.DataFrame,
    source_prefix: str,
) -> pd.DataFrame:
    """
    Convert long AGI/Fraction rows into one row per AGI with fraction-specific columns.

    Example output columns:
        AGI
        Fan2016_soluble_log2k_mean
        Fan2016_soluble_log2k_sd
        Fan2016_soluble_log2k_cv
        Fan2016_soluble_halflife_hr
    """
    agis = (
        df[["AGI"]]
        .drop_duplicates()
        .sort_values("AGI")
        .reset_index(drop=True)
    )

    wide_df = agis.copy()

    fraction_labels = sorted(df["Fraction_label"].dropna().unique())

    for fraction_label in fraction_labels:
        fraction_df = df[df["Fraction_label"] == fraction_label].copy()

        keep_cols = ["AGI"] + VALUE_COLUMNS
        fraction_df = fraction_df[keep_cols]

        rename_map = {
            col: f"{source_prefix}_{fraction_label}_{RAW_TO_CLEAN_METRIC[col]}"
            for col in VALUE_COLUMNS
        }

        fraction_df = fraction_df.rename(columns=rename_map)

        wide_df = wide_df.merge(fraction_df, how="left", on="AGI")

    return wide_df


def process_db1(
    input_file: str,
    output_file: str,
    source_prefix: str,
    duplicate_policy: str,
) -> None:
    """
    Preprocess Fan2016 protein half-life data.

    Input columns expected:
        Fraction
        UniProt ID
        Protein
        AGI
        # peptide
        Mean of turnover rate (log2k)
        SD
        CV
        Half-life (hr-1)

    Output:
        One row per uppercase base AGI.

    Output columns:
        AGI
        <source_prefix>_<fraction>_Mean of turnover rate (log2k)
        <source_prefix>_<fraction>_SD
        <source_prefix>_<fraction>_CV
        <source_prefix>_<fraction>_Half-life (hr-1)

    Multi-AGI rows are exploded first.
    Different fractions for the same AGI are pivoted into separate columns.
    """
    input_path = Path(input_file)
    output_path = Path(output_file)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # The provided example is tab-delimited.
    df = pd.read_csv(input_path, sep="\t", dtype=str)

    # Clean accidental whitespace from column names.
    df.columns = [col.strip() for col in df.columns]

    validate_columns(df, input_file)

    # Retain only the fraction, AGI key, and requested measurement columns.
    df = df[REQUIRED_COLUMNS].copy()

    # Clean fraction labels for use in output column names.
    df["Fraction"] = df["Fraction"].astype(str).str.strip()
    df["Fraction_label"] = df["Fraction"].apply(clean_fraction_label)

    # Split multi-AGI fields into one AGI per row.
    df["AGI"] = df["AGI"].apply(clean_agi_string)
    df = df.explode("AGI", ignore_index=True)

    # Remove rows with no usable AGI.
    df = df[df["AGI"].notna() & (df["AGI"].astype(str).str.len() > 0)].copy()

    # Convert retained measurement columns to numeric.
    for col in VALUE_COLUMNS:
        df[col] = pd.to_numeric(df[col].astype(str).str.strip(), errors="coerce")

    # Drop exact duplicate rows before checking true AGI + Fraction duplicates.
    df = df.drop_duplicates().reset_index(drop=True)

    # Different fractions for the same AGI are allowed.
    # Only repeated AGI + same Fraction is considered a duplicate problem.
    df = handle_duplicate_agi_fraction_rows(
        df=df,
        duplicate_policy=duplicate_policy,
    )

    # Pivot to one row per AGI with fraction-specific columns.
    wide_df = pivot_fraction_specific_columns(
        df=df,
        source_prefix=source_prefix,
    )

    wide_df.to_csv(output_path, index=False)

    print(f"Wrote processed DB1 table: {output_path}")
    print(f"Rows written: {len(wide_df):,}")
    print(f"Unique AGIs: {wide_df['AGI'].nunique():,}")
    print(f"Fractions observed: {sorted(df['Fraction_label'].dropna().unique())}")


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess Fan2016 protein half-life data."
    )

    parser.add_argument(
        "--input_file",
        required=True,
        help="Input Fan2016 half-life table.",
    )

    parser.add_argument(
        "--output_file",
        required=True,
        help="Output processed CSV file.",
    )

    parser.add_argument(
        "--source_prefix",
        default="Fan2016",
        help="Prefix to add to retained source columns. Default: Fan2016",
    )

    parser.add_argument(
        "--duplicate_policy",
        choices=["error", "first", "mean"],
        default="error",
        help=(
            "How to handle repeated rows for the same AGI and same Fraction. "
            "Different fractions for the same AGI are always retained as separate columns. "
            "Default: error."
        ),
    )

    args = parser.parse_args()

    process_db1(
        input_file=args.input_file,
        output_file=args.output_file,
        source_prefix=args.source_prefix,
        duplicate_policy=args.duplicate_policy,
    )


if __name__ == "__main__":
    main()
