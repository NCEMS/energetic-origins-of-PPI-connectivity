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
    "p.value",
]


MEASUREMENT_COLUMNS = [
    "Diff. log2k",
    "Fold change in k",
    "ctrl.mean",
    "30C.mean",
]


RAW_TO_CLEAN_METRIC = {
    "Diff. log2k": "diff_log2k",
    "Fold change in k": "k_fold_change",
    "ctrl.mean": "ctrl_log2k_mean",
    "30C.mean": "temp30_log2k_mean",
}


AGI_RE = re.compile(r"^AT[1-5CM]G\d{5}$")


FUSED_PVALUE_CTRL_RE = re.compile(
    r"^\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"([+-](?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"\s*$"
)


def clean_agi_string(value) -> list[str]:
    """
    Convert a raw AGI field into clean uppercase base AGI IDs.

    Examples:
        At4g33360 -> AT4G33360
        AT1G01080.1,AT1G01080.2 -> AT1G01080
        At1g56070; At1g56075 -> AT1G56070, AT1G56075
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

        agi = re.sub(r"\.\d+$", "", agi)

        if AGI_RE.match(agi) and agi not in seen:
            clean_parts.append(agi)
            seen.add(agi)

    return clean_parts


def clean_label(value, fallback: str) -> str:
    """
    Convert tissue/fraction labels into safe column-name components.
    """
    if pd.isna(value):
        return fallback

    label = str(value).strip()

    if not label:
        return fallback

    label = label.lower()
    label = label.replace("°", "")
    label = label.replace("%", "pct")
    label = re.sub(r"\s+", "_", label)
    label = re.sub(r"[^a-z0-9._+-]+", "_", label)
    label = re.sub(r"_+", "_", label).strip("_")

    if not label:
        return fallback

    return label


def validate_columns(df: pd.DataFrame, input_file: str) -> None:
    """
    Ensure all required columns are present after header cleanup.
    """
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]

    if missing:
        raise ValueError(
            f"Missing required columns in {input_file}: {missing}\n"
            f"Observed columns: {list(df.columns)}"
        )


def repair_fused_pvalue_ctrlmean_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Repair rows where the source file appears to have fused p.value and ctrl.mean.

    Example problematic field:
        p.value = '0.0181-5.159'

    Expected interpretation:
        p.value   = 0.0181
        ctrl.mean = -5.159

    Because the missing delimiter shifts later fields left, this also shifts:
        old ctrl.mean -> 30C.mean
        old 30C.mean  -> Tissue
        old Tissue    -> Fraction
    """
    df = df.copy()

    repaired = 0

    for idx, value in df["p.value"].items():
        match = FUSED_PVALUE_CTRL_RE.match(str(value))

        if not match:
            continue

        fraction_missing = (
            pd.isna(df.at[idx, "Fraction"])
            or str(df.at[idx, "Fraction"]).strip() == ""
        )

        if not fraction_missing:
            continue

        old_ctrl_mean = df.at[idx, "ctrl.mean"]
        old_30c_mean = df.at[idx, "30C.mean"]
        old_tissue = df.at[idx, "Tissue"]

        df.at[idx, "p.value"] = match.group(1)
        df.at[idx, "ctrl.mean"] = match.group(2)
        df.at[idx, "30C.mean"] = old_ctrl_mean
        df.at[idx, "Tissue"] = old_30c_mean
        df.at[idx, "Fraction"] = old_tissue

        repaired += 1

    if repaired:
        print(f"Repaired rows with fused p.value/ctrl.mean fields: {repaired:,}")

    return df


def coerce_numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert numeric source columns to floats.
    """
    df = df.copy()

    numeric_columns = MEASUREMENT_COLUMNS + ["p.value"]

    for col in numeric_columns:
        df[col] = pd.to_numeric(df[col].astype(str).str.strip(), errors="coerce")

    return df


def handle_duplicate_agi_condition_rows(
    df: pd.DataFrame,
    duplicate_policy: str,
) -> pd.DataFrame:
    """
    Handle duplicated AGI + Tissue + Fraction rows after AGI splitting.

    Different Tissue/Fraction combinations for the same AGI are expected and
    are preserved. This function only resolves repeated rows for the same
    AGI, Tissue, and Fraction condition.
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
            "Inspect the source table or rerun with --duplicate_policy first, "
            "mean, or median.\n\n"
            f"Example duplicate rows:\n{example}"
        )

    if duplicate_policy == "first":
        return (
            df.drop_duplicates(subset=duplicate_key, keep="first")
            .reset_index(drop=True)
        )

    if duplicate_policy not in {"mean", "median"}:
        raise ValueError(f"Unsupported duplicate_policy: {duplicate_policy}")

    agg_func = "mean" if duplicate_policy == "mean" else "median"

    grouped = (
        df.groupby(duplicate_key, as_index=False)
        .agg(
            {
                "Diff. log2k": agg_func,
                "ctrl.mean": agg_func,
                "30C.mean": agg_func,
            }
        )
        .reset_index(drop=True)
    )

    # Keep fold-change internally consistent with the aggregated log2 difference.
    grouped["Fold change in k"] = 2 ** grouped["Diff. log2k"]

    grouped = grouped[duplicate_key + MEASUREMENT_COLUMNS]

    return grouped


def pivot_tissue_fraction_specific_columns(
    df: pd.DataFrame,
    source_prefix: str,
) -> pd.DataFrame:
    """
    Convert long AGI/Tissue/Fraction rows into one row per AGI with
    condition-specific columns.
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

        keep_cols = ["AGI"] + MEASUREMENT_COLUMNS
        condition_df = condition_df[keep_cols]

        condition_prefix = f"{source_prefix}_{tissue_label}_{fraction_label}"

        rename_map = {
            col: f"{condition_prefix}_{RAW_TO_CLEAN_METRIC[col]}"
            for col in MEASUREMENT_COLUMNS
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
    Preprocess Fan2024 protein turnover-response data.

    Output:
        One row per uppercase base AGI.

    Output columns:
        AGI
        <source_prefix>_<tissue>_<fraction>_diff_log2k
        <source_prefix>_<tissue>_<fraction>_k_fold_change
        <source_prefix>_<tissue>_<fraction>_ctrl_log2k_mean
        <source_prefix>_<tissue>_<fraction>_temp30_log2k_mean
    """
    input_path = Path(input_file)
    output_path = Path(output_file)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t", dtype=str, engine="python")

    # Important because the header can contain names such as "Diff. log2k ".
    df.columns = [col.strip() for col in df.columns]

    validate_columns(df, input_file)

    # Repair apparent source-formatting issue before selecting/cleaning columns.
    df = repair_fused_pvalue_ctrlmean_rows(df)

    df = df[REQUIRED_COLUMNS].copy()

    df["Tissue"] = df["Tissue"].astype(str).str.strip()
    df["Fraction"] = df["Fraction"].astype(str).str.strip()

    df["Tissue_label"] = df["Tissue"].apply(lambda x: clean_label(x, "unknown_tissue"))
    df["Fraction_label"] = df["Fraction"].apply(
        lambda x: clean_label(x, "unknown_fraction")
    )

    df["AGI"] = df["AGI"].apply(clean_agi_string)
    df = df.explode("AGI", ignore_index=True)

    df = df[df["AGI"].notna() & (df["AGI"].astype(str).str.len() > 0)].copy()

    df = coerce_numeric_columns(df)

    # Drop exact duplicate rows before handling true condition-level duplicates.
    df = df.drop_duplicates().reset_index(drop=True)

    df = handle_duplicate_agi_condition_rows(
        df=df,
        duplicate_policy=duplicate_policy,
    )

    wide_df = pivot_tissue_fraction_specific_columns(
        df=df,
        source_prefix=source_prefix,
    )

    wide_df.to_csv(output_path, index=False)

    print(f"Wrote processed DB2 table: {output_path}")
    print(f"Rows written: {len(wide_df):,}")
    print(f"Unique AGIs: {wide_df['AGI'].nunique():,}")
    print(f"Duplicate policy used: {duplicate_policy}")

    observed_conditions = (
        df[["Tissue_label", "Fraction_label"]]
        .drop_duplicates()
        .sort_values(["Tissue_label", "Fraction_label"])
    )

    print("Tissue/Fraction conditions observed:")
    for row in observed_conditions.itertuples(index=False):
        print(f"  {row.Tissue_label} / {row.Fraction_label}")

    print(f"Output columns: {list(wide_df.columns)}")


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
