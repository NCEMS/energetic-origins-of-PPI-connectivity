#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


AGI_PATTERN = re.compile(r"AT[1-5CM]G\d{5}", flags=re.IGNORECASE)
REQUIRED_COLUMNS = {"AGI", "Tm_median", "Tm_sd", "n_modeled"}
REPLICATE_R2_PATTERN = re.compile(r"^C\d+\.r2$", flags=re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Preprocess the Volkening et al. 2019 Arabidopsis melting-temperature "
            "dataset into a compact, one-row-per-AGI table."
        )
    )
    parser.add_argument("--input", required=True, help="Input Volkening CSV file.")
    parser.add_argument("--output", required=True, help="Processed gene-level CSV.")
    parser.add_argument(
        "--audit_output",
        required=True,
        help="Audit CSV containing retained source values and processing status.",
    )
    parser.add_argument(
        "--duplicate_output",
        required=True,
        help="CSV containing duplicate AGI rows requiring review.",
    )
    parser.add_argument(
        "--summary_output",
        required=True,
        help="TSV containing preprocessing counts and selected columns.",
    )
    return parser.parse_args()


def ensure_parent_dirs(paths: Iterable[str | Path]) -> None:
    for path in paths:
        Path(path).parent.mkdir(parents=True, exist_ok=True)


def extract_agis(value: object) -> list[str]:
    """Extract unique locus-level AGIs from a possibly messy source field."""
    if pd.isna(value):
        return []

    matches = AGI_PATTERN.findall(str(value).upper())
    return sorted(set(match.upper() for match in matches))


def main() -> None:
    args = parse_args()
    output_paths = [
        args.output,
        args.audit_output,
        args.duplicate_output,
        args.summary_output,
    ]
    ensure_parent_dirs(output_paths)

    # Do not leave a stale principal output after a failed validation run.
    Path(args.output).unlink(missing_ok=True)

    input_path = Path(args.input)
    header = pd.read_csv(input_path, nrows=0)
    header_columns = list(header.columns)

    missing = sorted(REQUIRED_COLUMNS.difference(header_columns))
    if missing:
        raise KeyError(
            f"{input_path} is missing required column(s): {missing}. "
            f"Available columns: {header_columns}"
        )

    replicate_r2_cols = [
        column for column in header_columns if REPLICATE_R2_PATTERN.match(column)
    ]
    if not replicate_r2_cols:
        raise KeyError(
            f"No replicate R2 columns matching C1.r2, C2.r2, ... were found in "
            f"{input_path}. Available columns: {header_columns}"
        )

    usecols = ["AGI", "Tm_median", "Tm_sd", "n_modeled", *replicate_r2_cols]
    source = pd.read_csv(input_path, usecols=usecols)
    input_row_count = len(source)
    source.insert(0, "source_row", np.arange(1, len(source) + 1, dtype=int))
    source["source_AGI"] = source["AGI"]
    source["gene"] = source["AGI"].map(extract_agis)
    source["n_agis_in_source"] = source["gene"].map(len)
    source = source.explode("gene", ignore_index=True)

    for column in ["Tm_median", "Tm_sd", "n_modeled", *replicate_r2_cols]:
        source[column] = pd.to_numeric(source[column], errors="coerce")

    source["R2_median"] = source[replicate_r2_cols].median(axis=1, skipna=True)

    source["processing_status"] = np.select(
        [
            source["gene"].isna(),
            source["Tm_median"].isna(),
        ],
        [
            "missing_or_invalid_agi",
            "missing_tm",
        ],
        default="usable",
    )

    audit_columns = [
        "source_row",
        "source_AGI",
        "gene",
        "n_agis_in_source",
        "processing_status",
        "Tm_median",
        "Tm_sd",
        "n_modeled",
        "R2_median",
        *replicate_r2_cols,
    ]
    audit = source[audit_columns].copy()
    audit.to_csv(args.audit_output, index=False)

    usable = source.loc[
        source["processing_status"].eq("usable"),
        ["source_row", "gene", "Tm_median", "Tm_sd", "n_modeled", "R2_median"],
    ].copy()

    duplicate_mask = usable.duplicated(subset=["gene"], keep=False)
    duplicates = usable.loc[duplicate_mask].sort_values(
        ["gene", "source_row"], kind="stable"
    )
    duplicates.to_csv(args.duplicate_output, index=False)

    summary = pd.DataFrame(
        [
            {
                "source_file": input_path.name,
                "input_rows": input_row_count,
                "expanded_rows": len(source),
                "rows_with_valid_agi": int(source["gene"].notna().sum()),
                "rows_with_numeric_tm": int(source["Tm_median"].notna().sum()),
                "usable_rows": len(usable),
                "unique_usable_genes": int(usable["gene"].nunique()),
                "duplicate_rows": len(duplicates),
                "duplicate_genes": int(duplicates["gene"].nunique()),
                "replicate_r2_columns": ";".join(replicate_r2_cols),
                "retained_output_columns": (
                    "Volkening2019_Tm_median;Volkening2019_Tm_sd;"
                    "Volkening2019_n_modeled;Volkening2019_R2_median"
                ),
            }
        ]
    )
    summary.to_csv(args.summary_output, sep="\t", index=False)

    if not duplicates.empty:
        raise ValueError(
            f"Volkening preprocessing found {len(duplicates)} rows representing "
            f"{duplicates['gene'].nunique()} duplicated AGI(s). Review "
            f"{args.duplicate_output}. No duplicate values were averaged or selected."
        )

    processed = usable.drop(columns="source_row").rename(
        columns={
            "Tm_median": "Volkening2019_Tm_median",
            "Tm_sd": "Volkening2019_Tm_sd",
            "n_modeled": "Volkening2019_n_modeled",
            "R2_median": "Volkening2019_R2_median",
        }
    )

    processed["Volkening2019_n_modeled"] = (
        processed["Volkening2019_n_modeled"].round().astype("Int64")
    )
    processed = processed.sort_values("gene", kind="stable").reset_index(drop=True)
    processed.to_csv(args.output, index=False)

    print("Volkening et al. 2019 preprocessing complete")
    print(f"  Input rows: {summary.loc[0, 'input_rows']:,}")
    print(f"  Output genes: {len(processed):,}")
    print(f"  Replicate R2 columns summarized: {len(replicate_r2_cols)}")
    print(f"  Processed output: {args.output}")
    print(f"  Audit output: {args.audit_output}")
    print(f"  Duplicate output: {args.duplicate_output}")
    print(f"  Summary output: {args.summary_output}")


if __name__ == "__main__":
    main()
