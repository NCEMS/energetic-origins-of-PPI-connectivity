#!/usr/bin/env python3

"""Run an isolated discrepancy test for Arabidopsis PPI evidence curation.

The script reruns ../curate-athaliana-ppi-evidence.py into this directory's
outputs/ folder, then compares the regenerated in-vivo/in-vitro tables against
the existing manually curated files in ../data-files.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


OUTPUT_COLUMNS = [
    "protein_a",
    "protein_b",
    "label",
    "evidence_count",
    "source_files",
    "sources",
    "methods",
]

KEY_COLUMNS = ["protein_a", "protein_b", "label"]
METADATA_COLUMNS = ["evidence_count", "source_files", "sources", "methods"]


def read_edges(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str).fillna("")
    missing = [column for column in OUTPUT_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing columns: {missing}")

    df = df[OUTPUT_COLUMNS].copy()
    for column in ["protein_a", "protein_b", "label"]:
        df[column] = df[column].str.strip()
    for column in METADATA_COLUMNS:
        df[column] = df[column].astype(str).str.strip()
    df["evidence_count"] = df["evidence_count"].replace("", "0").astype(int).astype(str)
    return df.sort_values(OUTPUT_COLUMNS).reset_index(drop=True)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def rows_to_frame(rows: set[tuple[str, ...]]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    return pd.DataFrame(sorted(rows), columns=OUTPUT_COLUMNS)


def keys_to_frame(keys: set[tuple[str, ...]]) -> pd.DataFrame:
    if not keys:
        return pd.DataFrame(columns=KEY_COLUMNS)
    return pd.DataFrame(sorted(keys), columns=KEY_COLUMNS)


def compare_tables(label: str, manual_path: Path, generated_path: Path, out_dir: Path) -> list[str]:
    manual = read_edges(manual_path)
    generated = read_edges(generated_path)

    manual_rows = set(map(tuple, manual[OUTPUT_COLUMNS].itertuples(index=False, name=None)))
    generated_rows = set(map(tuple, generated[OUTPUT_COLUMNS].itertuples(index=False, name=None)))

    only_manual_rows = manual_rows - generated_rows
    only_generated_rows = generated_rows - manual_rows

    manual_keys = set(map(tuple, manual[KEY_COLUMNS].itertuples(index=False, name=None)))
    generated_keys = set(map(tuple, generated[KEY_COLUMNS].itertuples(index=False, name=None)))

    only_manual_keys = manual_keys - generated_keys
    only_generated_keys = generated_keys - manual_keys
    shared_keys = manual_keys & generated_keys

    manual_by_key = manual.set_index(KEY_COLUMNS)
    generated_by_key = generated.set_index(KEY_COLUMNS)
    mismatch_rows: list[dict[str, str]] = []
    for key in sorted(shared_keys):
        manual_row = manual_by_key.loc[key]
        generated_row = generated_by_key.loc[key]
        if isinstance(manual_row, pd.DataFrame) or isinstance(generated_row, pd.DataFrame):
            mismatch_rows.append(
                {
                    "protein_a": key[0],
                    "protein_b": key[1],
                    "label": key[2],
                    "field": "__duplicate_key__",
                    "manual": str(len(manual_row)) if isinstance(manual_row, pd.DataFrame) else "1",
                    "generated": str(len(generated_row)) if isinstance(generated_row, pd.DataFrame) else "1",
                }
            )
            continue

        for column in METADATA_COLUMNS:
            manual_value = str(manual_row[column])
            generated_value = str(generated_row[column])
            if manual_value != generated_value:
                mismatch_rows.append(
                    {
                        "protein_a": key[0],
                        "protein_b": key[1],
                        "label": key[2],
                        "field": column,
                        "manual": manual_value,
                        "generated": generated_value,
                    }
                )

    mismatch = pd.DataFrame(
        mismatch_rows,
        columns=["protein_a", "protein_b", "label", "field", "manual", "generated"],
    )

    write_csv(rows_to_frame(only_manual_rows), out_dir / f"{label}_only_manual_exact.csv")
    write_csv(rows_to_frame(only_generated_rows), out_dir / f"{label}_only_generated_exact.csv")
    write_csv(keys_to_frame(only_manual_keys), out_dir / f"{label}_keys_only_manual.csv")
    write_csv(keys_to_frame(only_generated_keys), out_dir / f"{label}_keys_only_generated.csv")
    write_csv(mismatch, out_dir / f"{label}_metadata_mismatch.csv")

    return [
        f"[{label}] manual rows:                 {len(manual):,}",
        f"[{label}] generated rows:              {len(generated):,}",
        f"[{label}] exact rows shared:           {len(manual_rows & generated_rows):,}",
        f"[{label}] exact rows only in manual:   {len(only_manual_rows):,}",
        f"[{label}] exact rows only in generated:{len(only_generated_rows):,}",
        f"[{label}] keys shared:                 {len(shared_keys):,}",
        f"[{label}] keys only in manual:         {len(only_manual_keys):,}",
        f"[{label}] keys only in generated:      {len(only_generated_keys):,}",
        f"[{label}] metadata mismatch cells:     {len(mismatch):,}",
    ]


def write_missing_inputs_report(missing: list[Path], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    report = out_dir / "MISSING_RAW_INPUTS.txt"
    lines = [
        "The discrepancy test did not run because raw source inputs are missing.",
        "",
        "Missing:",
    ]
    lines.extend(f"- {path}" for path in missing)
    lines.extend(
        [
            "",
            "Put those files in discrepancy-test/inputs/ or rerun with:",
            "",
            "python run-discrepancy-test.py --biogrid path/to/biogrid --table-s1 path/to/table_s1",
        ]
    )
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(report.read_text(encoding="utf-8"))


def run_curation(args: argparse.Namespace, generated_vivo: Path, generated_vitro: Path) -> None:
    command = [
        sys.executable,
        str(args.curation_script),
        "--biogrid",
        str(args.biogrid),
        "--table-s1",
        str(args.table_s1),
        "--in-vivo-output",
        str(generated_vivo),
        "--in-vitro-output",
        str(generated_vitro),
        "--unknown-methods-output",
        str(args.output_dir / "unknown_methods.csv"),
    ]
    command.extend(args.curation_args)
    subprocess.run(command, check=True)


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    download_dir = here.parent

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--biogrid",
        type=Path,
        default=here / "inputs" / "BIOGRID-ARABIDOPSIS-THALIANA-5.0.258.tab3.xlsx",
        help="Raw BioGRID table path.",
    )
    parser.add_argument(
        "--table-s1",
        type=Path,
        default=here / "inputs" / "Table S1.xlsx",
        help="Raw Table S1 path.",
    )
    parser.add_argument(
        "--curation-script",
        type=Path,
        default=download_dir / "curate-athaliana-ppi-evidence.py",
    )
    parser.add_argument(
        "--manual-in-vivo",
        type=Path,
        default=download_dir / "data-files" / "ppi_in_vivo_edges.csv",
    )
    parser.add_argument(
        "--manual-in-vitro",
        type=Path,
        default=download_dir / "data-files" / "ppi_in_vitro_edges.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=here / "outputs",
        help="Directory for regenerated files and discrepancy reports.",
    )
    parser.add_argument(
        "--skip-curation",
        action="store_true",
        help="Compare existing generated outputs without rerunning curation.",
    )
    parser.add_argument(
        "curation_args",
        nargs=argparse.REMAINDER,
        help="Optional extra args passed to curate-athaliana-ppi-evidence.py after --.",
    )

    args = parser.parse_args()
    if args.curation_args and args.curation_args[0] == "--":
        args.curation_args = args.curation_args[1:]
    return args


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    generated_vivo = args.output_dir / "generated_ppi_in_vivo_edges.csv"
    generated_vitro = args.output_dir / "generated_ppi_in_vitro_edges.csv"

    required = [args.curation_script, args.manual_in_vivo, args.manual_in_vitro]
    if not args.skip_curation:
        required.extend([args.biogrid, args.table_s1])
    else:
        required.extend([generated_vivo, generated_vitro])

    missing = [path for path in required if not path.exists()]
    if missing:
        write_missing_inputs_report(missing, args.output_dir)
        raise SystemExit(2)

    if not args.skip_curation:
        run_curation(args, generated_vivo, generated_vitro)

    lines = []
    lines.extend(compare_tables("in_vivo", args.manual_in_vivo, generated_vivo, args.output_dir))
    lines.append("")
    lines.extend(compare_tables("in_vitro", args.manual_in_vitro, generated_vitro, args.output_dir))

    summary = "\n".join(lines) + "\n"
    summary_path = args.output_dir / "discrepancy_summary.txt"
    summary_path.write_text(summary, encoding="utf-8")
    print(summary)
    print(f"Wrote discrepancy reports to: {args.output_dir}")


if __name__ == "__main__":
    main()
