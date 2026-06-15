#!/usr/bin/env python3

"""Curate A. thaliana PPI evidence from raw BioGRID and Table S1 tables.

The raw source workbooks are intentionally kept out of the GitHub repository.
This script recreates the two manually curated evidence tables consumed by
build-athaliana-ppi-networks.py:

* data-files/ppi_in_vivo_edges.csv
* data-files/ppi_in_vitro_edges.csv

Each output row is one undirected A. thaliana locus pair for one evidence class.
Duplicate raw evidence rows are aggregated into evidence_count, source_files,
sources, and methods columns.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


TAIR_LOCUS_RE = re.compile(r"\bAT(?:[1-5]G|MG|CG)\d{5}\b", re.IGNORECASE)

OUTPUT_COLUMNS = [
    "protein_a",
    "protein_b",
    "label",
    "evidence_count",
    "source_files",
    "sources",
    "methods",
]

BIOGRID_A_COLUMNS = [
    "Systematic Name Interactor A",
    "Official Symbol Interactor A",
    "Synonyms Interactor A",
    "Aliases Interactor A",
    "Interactor A",
    "Protein A",
    "Gene A",
    "Locus A",
    "AGI A",
    "TAIR A",
    "protein_a",
]

BIOGRID_B_COLUMNS = [
    "Systematic Name Interactor B",
    "Official Symbol Interactor B",
    "Synonyms Interactor B",
    "Aliases Interactor B",
    "Interactor B",
    "Protein B",
    "Gene B",
    "Locus B",
    "AGI B",
    "TAIR B",
    "protein_b",
]

TABLE_S1_A_COLUMNS = [
    "Protein A",
    "Protein 1",
    "Protein1",
    "Interactor A",
    "Interactor 1",
    "Interactor1",
    "Gene A",
    "Gene 1",
    "Gene1",
    "Locus A",
    "Locus 1",
    "AGI A",
    "AGI 1",
    "TAIR A",
    "TAIR 1",
    "ORF A",
    "ORF 1",
    "Bait",
    "protein_a",
]

TABLE_S1_B_COLUMNS = [
    "Protein B",
    "Protein 2",
    "Protein2",
    "Interactor B",
    "Interactor 2",
    "Interactor2",
    "Gene B",
    "Gene 2",
    "Gene2",
    "Locus B",
    "Locus 2",
    "AGI B",
    "AGI 2",
    "TAIR B",
    "TAIR 2",
    "ORF B",
    "ORF 2",
    "Prey",
    "protein_b",
]

METHOD_COLUMNS = [
    "Experimental System",
    "Method",
    "method",
    "Assay",
    "Technique",
    "Detection Method",
    "Experimental Method",
    "Interaction Detection Method",
]

SOURCE_COLUMNS = [
    "Publication Source",
    "Source",
    "Reference",
    "Citation",
    "Publication",
    "Pubmed ID",
    "PubMed ID",
    "PMID",
    "Author",
]

LABEL_COLUMNS = [
    "label",
    "Label",
    "Evidence class",
    "Evidence Class",
    "Interaction type",
    "Interaction Type",
    "Type",
]

IN_VIVO_METHOD_PATTERNS = [
    "affinity capture ms",
    "affinity capture western",
    "co fractionation",
    "co localization",
    "proximity label ms",
    "co immunoprecipitation",
    "co ip",
    "coip",
    "co purifications",
    "co purification",
    "immunoaffinity purification",
]

IN_VITRO_METHOD_PATTERNS = [
    "two hybrid",
    "2 hybrid",
    "y2h",
    "pca",
    "protein complementation assay",
    "reconstituted complex",
    "biochemical activity",
    "fret",
    "in vitro",
]


def normalize_column_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def normalize_method_text(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).lower()).strip()


def clean_text(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    cleaned = " ".join(str(value).replace("\n", " ").split())
    if not cleaned or cleaned.lower() == "nan":
        return None
    return cleaned


def parse_column_list(value: str | None) -> list[str] | None:
    if value is None:
        return None
    return [item.strip() for item in value.split(",") if item.strip()]


def find_column(df: pd.DataFrame, name: str) -> str | None:
    if name in df.columns:
        return name

    normalized_lookup = {
        normalize_column_name(column): column for column in df.columns
    }
    return normalized_lookup.get(normalize_column_name(name))


def resolve_columns(
    df: pd.DataFrame,
    candidates: Iterable[str],
    override: str | None,
    side: str | None = None,
) -> list[str]:
    requested = parse_column_list(override)
    source = requested if requested is not None else candidates

    columns: list[str] = []
    missing: list[str] = []
    for name in source:
        column = find_column(df, name)
        if column is None:
            missing.append(name)
            continue
        if column not in columns:
            columns.append(column)

    if requested is not None and missing:
        raise ValueError(f"Missing requested columns: {missing}")

    if not columns and side is not None:
        columns.extend(discover_endpoint_columns(df, side))

    return columns


def resolve_optional_column(
    df: pd.DataFrame,
    candidates: Iterable[str],
    override: str | None,
) -> str | None:
    requested = parse_column_list(override)
    source = requested if requested is not None else candidates
    for name in source:
        column = find_column(df, name)
        if column is not None:
            return column
    if requested is not None:
        raise ValueError(f"None of the requested columns were found: {requested}")
    return None


def discover_endpoint_columns(df: pd.DataFrame, side: str) -> list[str]:
    side = side.upper()
    side_markers = ("a", "1", "one") if side == "A" else ("b", "2", "two")
    side_words = (
        f"interactor{side.lower()}",
        f"protein{side.lower()}",
        f"gene{side.lower()}",
        f"locus{side.lower()}",
        f"agi{side.lower()}",
        f"tair{side.lower()}",
    )
    keywords = ("interactor", "protein", "gene", "locus", "agi", "tair", "orf")

    columns: list[str] = []
    for column in df.columns:
        normalized = normalize_column_name(column)
        if any(word in normalized for word in side_words):
            columns.append(column)
            continue
        if not any(keyword in normalized for keyword in keywords):
            continue
        if normalized.endswith(side_markers):
            columns.append(column)

    return columns


def extract_tair_loci(value: object) -> list[str]:
    text = clean_text(value)
    if text is None:
        return []
    matches = TAIR_LOCUS_RE.findall(text.upper())
    return sorted(set(matches))


def row_loci(row: pd.Series, columns: list[str]) -> list[str]:
    for column in columns:
        loci = extract_tair_loci(row.get(column))
        if loci:
            return loci
    return []


def classify_method(method: str | None, label_hint: str | None) -> str | None:
    method_text = normalize_method_text(method or "")
    hint_text = normalize_method_text(label_hint or "")
    combined = " ".join(part for part in [method_text, hint_text] if part)

    if any(pattern in combined for pattern in IN_VIVO_METHOD_PATTERNS):
        return "in vivo"
    if any(pattern in combined for pattern in IN_VITRO_METHOD_PATTERNS):
        return "in vitro"

    if "in vivo" in combined:
        return "in vivo"
    if "in vitro" in combined:
        return "in vitro"

    return None


def parse_sheet_name(value: str | None) -> int | str:
    if value is None:
        return 0
    try:
        return int(value)
    except ValueError:
        return value


def read_source_table(path: Path, sheet_name: str | None, header_row: int) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)

    suffix = path.suffix.lower()
    try:
        if suffix in {".xlsx", ".xls"}:
            return pd.read_excel(
                path,
                sheet_name=parse_sheet_name(sheet_name),
                header=header_row,
                dtype=str,
            )

        sep = "\t" if suffix in {".tsv", ".tab", ".txt"} or ".tab" in path.name.lower() else ","
        return pd.read_csv(path, sep=sep, header=header_row, dtype=str, low_memory=False)
    except ImportError as exc:
        raise RuntimeError(
            "Reading Excel inputs requires the optional pandas Excel engine "
            "(usually openpyxl). Install it in the environment and rerun."
        ) from exc


def iter_evidence_records(
    df: pd.DataFrame,
    source_file: Path,
    protein_a_columns: list[str],
    protein_b_columns: list[str],
    method_column: str | None,
    source_column: str | None,
    label_column: str | None,
    fallback_source: str,
    unknown_method_policy: str,
) -> tuple[list[dict[str, str]], pd.DataFrame]:
    records: list[dict[str, str]] = []
    unknown_methods: list[dict[str, str]] = []

    for row_index, row in df.iterrows():
        left_loci = row_loci(row, protein_a_columns)
        right_loci = row_loci(row, protein_b_columns)
        if not left_loci or not right_loci:
            continue

        method = clean_text(row.get(method_column)) if method_column else None
        source = clean_text(row.get(source_column)) if source_column else None
        label_hint = clean_text(row.get(label_column)) if label_column else None
        label = classify_method(method, label_hint)

        if label is None:
            unknown_methods.append(
                {
                    "source_file": source_file.name,
                    "row_number": str(row_index + 2),
                    "method": method or "",
                    "label_hint": label_hint or "",
                }
            )
            continue

        source_value = source or fallback_source
        method_value = method or label_hint or label

        for left in left_loci:
            for right in right_loci:
                protein_a, protein_b = sorted([left, right])
                records.append(
                    {
                        "protein_a": protein_a,
                        "protein_b": protein_b,
                        "label": label,
                        "source_file": source_file.name,
                        "source": source_value,
                        "method": method_value,
                    }
                )

    if unknown_methods and unknown_method_policy == "error":
        methods = pd.DataFrame(unknown_methods)
        examples = methods.drop_duplicates(["method", "label_hint"]).head(10)
        raise ValueError(
            "Encountered evidence rows with methods that could not be classified. "
            f"Examples:\n{examples.to_string(index=False)}"
        )

    return records, pd.DataFrame(unknown_methods)


def join_unique(values: pd.Series) -> str:
    unique = sorted({str(value) for value in values if clean_text(value) is not None})
    return "; ".join(unique)


def aggregate_records(records: list[dict[str, str]]) -> pd.DataFrame:
    if not records:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    df = pd.DataFrame(records)
    aggregated = (
        df.groupby(["protein_a", "protein_b", "label"], as_index=False)
        .agg(
            evidence_count=("source_file", "size"),
            source_files=("source_file", join_unique),
            sources=("source", join_unique),
            methods=("method", join_unique),
        )
        .sort_values(["protein_a", "protein_b", "label"])
    )

    return aggregated[OUTPUT_COLUMNS]


def curate_one_source(
    path: Path,
    sheet_name: str | None,
    header_row: int,
    a_candidates: list[str],
    b_candidates: list[str],
    a_override: str | None,
    b_override: str | None,
    method_override: str | None,
    source_override: str | None,
    label_override: str | None,
    fallback_source: str,
    unknown_method_policy: str,
) -> tuple[list[dict[str, str]], pd.DataFrame]:
    df = read_source_table(path, sheet_name, header_row)
    df = df.dropna(how="all")

    protein_a_columns = resolve_columns(df, a_candidates, a_override, side="A")
    protein_b_columns = resolve_columns(df, b_candidates, b_override, side="B")
    if not protein_a_columns or not protein_b_columns:
        raise ValueError(
            f"Could not identify protein endpoint columns in {path}. "
            "Use --*-protein-a-columns and --*-protein-b-columns to set them."
        )

    method_column = resolve_optional_column(df, METHOD_COLUMNS, method_override)
    source_column = resolve_optional_column(df, SOURCE_COLUMNS, source_override)
    label_column = resolve_optional_column(df, LABEL_COLUMNS, label_override)

    if method_column is None and label_column is None:
        raise ValueError(
            f"Could not identify a method or label column in {path}. "
            "Use --*-method-column or --*-label-column to set it."
        )

    print(f"{path.name}: endpoint A columns = {protein_a_columns}")
    print(f"{path.name}: endpoint B columns = {protein_b_columns}")
    print(f"{path.name}: method column = {method_column}")
    print(f"{path.name}: source column = {source_column}")
    print(f"{path.name}: label column = {label_column}")

    return iter_evidence_records(
        df=df,
        source_file=path,
        protein_a_columns=protein_a_columns,
        protein_b_columns=protein_b_columns,
        method_column=method_column,
        source_column=source_column,
        label_column=label_column,
        fallback_source=fallback_source,
        unknown_method_policy=unknown_method_policy,
    )


def write_output(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def curate_evidence(args: argparse.Namespace) -> None:
    all_records: list[dict[str, str]] = []
    unknown_reports: list[pd.DataFrame] = []

    biogrid_records, biogrid_unknown = curate_one_source(
        path=args.biogrid,
        sheet_name=args.biogrid_sheet,
        header_row=args.biogrid_header_row,
        a_candidates=BIOGRID_A_COLUMNS,
        b_candidates=BIOGRID_B_COLUMNS,
        a_override=args.biogrid_protein_a_columns,
        b_override=args.biogrid_protein_b_columns,
        method_override=args.biogrid_method_column,
        source_override=args.biogrid_source_column,
        label_override=args.biogrid_label_column,
        fallback_source="BioGRID",
        unknown_method_policy=args.unknown_method_policy,
    )
    all_records.extend(biogrid_records)
    unknown_reports.append(biogrid_unknown)

    table_s1_records, table_s1_unknown = curate_one_source(
        path=args.table_s1,
        sheet_name=args.table_s1_sheet,
        header_row=args.table_s1_header_row,
        a_candidates=TABLE_S1_A_COLUMNS,
        b_candidates=TABLE_S1_B_COLUMNS,
        a_override=args.table_s1_protein_a_columns,
        b_override=args.table_s1_protein_b_columns,
        method_override=args.table_s1_method_column,
        source_override=args.table_s1_source_column,
        label_override=args.table_s1_label_column,
        fallback_source="Arabidopsis consortium",
        unknown_method_policy=args.unknown_method_policy,
    )
    all_records.extend(table_s1_records)
    unknown_reports.append(table_s1_unknown)

    curated = aggregate_records(all_records)
    in_vivo = curated[curated["label"] == "in vivo"].copy()
    in_vitro = curated[curated["label"] == "in vitro"].copy()

    write_output(in_vivo, args.in_vivo_output)
    write_output(in_vitro, args.in_vitro_output)

    unknown = pd.concat(unknown_reports, ignore_index=True)
    if args.unknown_methods_output is not None and not unknown.empty:
        write_output(unknown, args.unknown_methods_output)

    print(f"Raw evidence records retained: {len(all_records):,}")
    print(f"In-vivo curated edges:         {len(in_vivo):,}")
    print(f"In-vitro curated edges:        {len(in_vitro):,}")
    print(f"Rows skipped for unknown methods: {len(unknown):,}")
    if args.unknown_methods_output is not None and not unknown.empty:
        print(f"Unknown-method report: {args.unknown_methods_output}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--biogrid",
        type=Path,
        default=Path("data-files/BIOGRID-ARABIDOPSIS-THALIANA-5.0.258.tab3.xlsx"),
        help="Raw BioGRID Arabidopsis table, usually TAB3 exported as xlsx/tsv.",
    )
    parser.add_argument(
        "--table-s1",
        type=Path,
        default=Path("data-files/Table S1.xlsx"),
        help="Raw Table S1 workbook/table from the Arabidopsis interactome source.",
    )
    parser.add_argument(
        "--in-vivo-output",
        type=Path,
        default=Path("data-files/ppi_in_vivo_edges.csv"),
        help="Output curated in-vivo evidence table.",
    )
    parser.add_argument(
        "--in-vitro-output",
        type=Path,
        default=Path("data-files/ppi_in_vitro_edges.csv"),
        help="Output curated in-vitro evidence table.",
    )
    parser.add_argument(
        "--unknown-methods-output",
        type=Path,
        default=None,
        help="Optional CSV report of rows skipped because their method was unknown.",
    )
    parser.add_argument(
        "--unknown-method-policy",
        choices=["skip", "error"],
        default="skip",
        help="How to handle rows whose method cannot be classified.",
    )

    parser.add_argument("--biogrid-sheet", default="0", help="Excel sheet index/name.")
    parser.add_argument("--table-s1-sheet", default="0", help="Excel sheet index/name.")
    parser.add_argument("--biogrid-header-row", type=int, default=0)
    parser.add_argument("--table-s1-header-row", type=int, default=0)

    parser.add_argument("--biogrid-protein-a-columns")
    parser.add_argument("--biogrid-protein-b-columns")
    parser.add_argument("--biogrid-method-column")
    parser.add_argument("--biogrid-source-column")
    parser.add_argument("--biogrid-label-column")

    parser.add_argument("--table-s1-protein-a-columns")
    parser.add_argument("--table-s1-protein-b-columns")
    parser.add_argument("--table-s1-method-column")
    parser.add_argument("--table-s1-source-column")
    parser.add_argument("--table-s1-label-column")

    return parser.parse_args()


def main() -> None:
    curate_evidence(parse_args())


if __name__ == "__main__":
    main()
