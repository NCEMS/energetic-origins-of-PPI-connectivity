#!/usr/bin/env python3

"""
Add protein-entanglement annotations to a pickled node table.

Entanglement result files are expected to be pipe-delimited CSV files with
names such as:

    AF-Q9ZWT3-F1-model_v6.csv

The UniProt accession is extracted from the filename. A file containing only
the header row is interpreted as not entangled. A file containing one or more
nonempty data rows is interpreted as entangled.

The resulting Boolean column is:

    is_entangled
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List

import pandas as pd


ALPHAFOLD_FILENAME_PATTERN = re.compile(
    r"^AF-(?P<uniprot>.+?)-F\d+-model_v\d+\.csv$",
    flags=re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Determine whether each AlphaFold structure has entanglement "
            "records and merge the result into a pickled node table."
        )
    )

    parser.add_argument(
        "--nodes",
        required=True,
        type=Path,
        help="Input pickled pandas node table.",
    )
    parser.add_argument(
        "--ent_data_dir",
        required=True,
        type=Path,
        help="Directory containing AlphaFold entanglement CSV files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output path for the annotated pickled node table.",
    )
    parser.add_argument(
        "--uniprot_col",
        default="UniProtKB-AC",
        help=(
            "UniProt accession column in the node table. "
            "Default: UniProtKB-AC"
        ),
    )
    parser.add_argument(
        "--audit_output",
        type=Path,
        default=None,
        help=(
            "Optional CSV file containing one row per processed entanglement "
            "file."
        ),
    )
    parser.add_argument(
        "--unmatched_output",
        type=Path,
        default=None,
        help=(
            "Optional CSV file listing entanglement files whose UniProt "
            "accessions are absent from the node table."
        ),
    )

    return parser.parse_args()


def normalize_uniprot(value: object) -> str | None:
    """
    Normalize a UniProt accession.

    Terminal isoform suffixes are removed so that, for example, P12345-1
    matches P12345.
    """

    if pd.isna(value):
        return None

    accession = str(value).strip().upper()

    if not accession:
        return None

    accession = re.sub(r"-\d+$", "", accession)

    return accession or None


def extract_uniprot_from_filename(path: Path) -> str:
    """Extract the UniProt accession from an AlphaFold result filename."""

    match = ALPHAFOLD_FILENAME_PATTERN.match(path.name)

    if match is None:
        raise ValueError(
            "Could not extract a UniProt accession from entanglement filename "
            f"{path.name!r}. Expected a name such as "
            "'AF-Q9ZWT3-F1-model_v6.csv'."
        )

    accession = normalize_uniprot(match.group("uniprot"))

    if accession is None:
        raise ValueError(
            f"Filename did not contain a usable UniProt accession: {path}"
        )

    return accession


def file_has_data_rows(path: Path) -> tuple[bool, int]:
    """
    Determine whether a pipe-delimited file contains nonempty data rows.

    Returns
    -------
    tuple
        ``(is_entangled, data_row_count)``.
    """

    data_row_count = 0

    with path.open("r", encoding="utf-8-sig") as handle:
        header_seen = False

        for raw_line in handle:
            line = raw_line.strip()

            # Ignore completely blank lines.
            if not line:
                continue

            if not header_seen:
                header_seen = True
                continue

            data_row_count += 1

    if not header_seen:
        raise ValueError(f"Entanglement file is completely empty: {path}")

    return data_row_count > 0, data_row_count


def build_entanglement_table(
    ent_data_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build a one-row-per-UniProt entanglement lookup and file-level audit.

    Multiple files for the same accession are allowed only when they agree.
    """

    if not ent_data_dir.is_dir():
        raise NotADirectoryError(
            f"Entanglement data directory does not exist: {ent_data_dir}"
        )

    csv_paths = sorted(ent_data_dir.glob("*.csv"))

    if not csv_paths:
        raise FileNotFoundError(
            f"No CSV files were found in entanglement directory: "
            f"{ent_data_dir}"
        )

    audit_rows: List[Dict[str, object]] = []

    for path in csv_paths:
        uniprot = extract_uniprot_from_filename(path)
        is_entangled, data_row_count = file_has_data_rows(path)

        audit_rows.append(
            {
                "filename": path.name,
                "UniProtKB-AC-normalized": uniprot,
                "data_row_count": data_row_count,
                "is_entangled": is_entangled,
            }
        )

    audit_df = pd.DataFrame(audit_rows)

    disagreement_counts = (
        audit_df.groupby("UniProtKB-AC-normalized")["is_entangled"]
        .nunique()
    )
    conflicting_accessions = disagreement_counts[
        disagreement_counts > 1
    ].index.tolist()

    if conflicting_accessions:
        raise ValueError(
            "Multiple entanglement files for the same UniProt accession "
            "disagree about entanglement status. Conflicting accessions: "
            f"{conflicting_accessions[:20]}"
        )

    # If duplicate files agree, collapse them to one accession-level result.
    entanglement_df = (
        audit_df.groupby(
            "UniProtKB-AC-normalized",
            as_index=False,
            sort=False,
        )
        .agg(
            is_entangled=("is_entangled", "first"),
            entanglement_file_count=("filename", "size"),
            entanglement_record_count=("data_row_count", "sum"),
        )
    )

    return entanglement_df, audit_df


def main() -> None:
    """Run the entanglement integration."""

    args = parse_args()

    if not args.nodes.is_file():
        raise FileNotFoundError(f"Node table does not exist: {args.nodes}")

    print(f"Reading node table: {args.nodes}")
    nodes_df = pd.read_pickle(args.nodes)

    if args.uniprot_col not in nodes_df.columns:
        raise KeyError(
            f"Node table does not contain required column "
            f"{args.uniprot_col!r}."
        )

    if "is_entangled" in nodes_df.columns:
        raise ValueError(
            "Node table already contains an 'is_entangled' column. "
            "Refusing to overwrite it."
        )

    entanglement_df, audit_df = build_entanglement_table(
        args.ent_data_dir
    )

    nodes_df = nodes_df.copy()
    nodes_df["_entanglement_uniprot_key"] = nodes_df[
        args.uniprot_col
    ].map(normalize_uniprot)

    duplicate_mapping = (
        nodes_df.dropna(subset=["_entanglement_uniprot_key"])
        .groupby("_entanglement_uniprot_key")["node"]
        .nunique()
    )

    ambiguous_accessions = duplicate_mapping[
        duplicate_mapping > 1
    ].index.tolist()

    if ambiguous_accessions:
        print(
            "WARNING: "
            f"{len(ambiguous_accessions):,} UniProt accession(s) occur for "
            "multiple nodes. The same entanglement status will be associated "
            "with each corresponding node because UniProtKB-AC is the "
            "pipeline-selected structure accession."
        )

    original_row_count = len(nodes_df)
    original_index = nodes_df.index.copy()

    annotated_df = nodes_df.merge(
        entanglement_df[
            ["UniProtKB-AC-normalized", "is_entangled"]
        ],
        how="left",
        left_on="_entanglement_uniprot_key",
        right_on="UniProtKB-AC-normalized",
        validate="many_to_one",
        sort=False,
    )

    if len(annotated_df) != original_row_count:
        raise RuntimeError(
            "Entanglement merge unexpectedly changed the node-row count "
            f"from {original_row_count:,} to {len(annotated_df):,}."
        )

    annotated_df.index = original_index

    matched_mask = annotated_df["is_entangled"].notna()

    matched_count = int(matched_mask.sum())
    entangled_count = int(
        annotated_df.loc[matched_mask, "is_entangled"].sum()
    )
    not_entangled_count = matched_count - entangled_count
    missing_count = int((~matched_mask).sum())

    # Use pandas' nullable Boolean type. Nodes without a corresponding
    # entanglement file remain <NA>, rather than being incorrectly labeled
    # False.
    annotated_df["is_entangled"] = annotated_df[
        "is_entangled"
    ].astype("boolean")

    annotated_df = annotated_df.drop(
        columns=[
            "_entanglement_uniprot_key",
            "UniProtKB-AC-normalized",
        ]
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    annotated_df.to_pickle(args.output)

    if args.audit_output is not None:
        args.audit_output.parent.mkdir(parents=True, exist_ok=True)
        audit_df.to_csv(args.audit_output, index=False)

    node_accessions = set(
        nodes_df["_entanglement_uniprot_key"].dropna()
    )

    unmatched_df = audit_df[
        ~audit_df["UniProtKB-AC-normalized"].isin(node_accessions)
    ].copy()

    if args.unmatched_output is not None:
        args.unmatched_output.parent.mkdir(parents=True, exist_ok=True)
        unmatched_df.to_csv(args.unmatched_output, index=False)

    print("")
    print("Entanglement integration summary")
    print("--------------------------------")
    print(f"Node rows:                    {len(annotated_df):,}")
    print(f"Entanglement CSV files:       {len(audit_df):,}")
    print(
        "Unique file accessions:      "
        f"{audit_df['UniProtKB-AC-normalized'].nunique():,}"
    )
    print(f"Nodes matched to files:       {matched_count:,}")
    print(f"Entangled nodes:              {entangled_count:,}")
    print(f"Not-entangled nodes:          {not_entangled_count:,}")
    print(f"Nodes without result file:    {missing_count:,}")
    print(f"Files absent from node table: {len(unmatched_df):,}")
    print(f"Output written to:            {args.output}")


if __name__ == "__main__":
    main()
