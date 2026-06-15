#!/usr/bin/env python3

"""Extract A. thaliana rows from a BioGRID ALL TAB3 file.

BioGRID ALL files are large. This helper streams the file and keeps rows where
both interactors have organism ID 3702, producing a small Arabidopsis-only TAB3
file that can be passed to run-discrepancy-test.py.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def find_column(header: list[str], name: str) -> int:
    try:
        return header.index(name)
    except ValueError as exc:
        raise ValueError(f"Missing required BioGRID column: {name}") from exc


def filter_biogrid(input_path: Path, output_path: Path, tax_id: str) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    scanned = 0
    kept = 0

    with input_path.open("r", encoding="utf-8", newline="") as in_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        header = next(reader)
        organism_a = find_column(header, "Organism ID Interactor A")
        organism_b = find_column(header, "Organism ID Interactor B")

        with output_path.open("w", encoding="utf-8", newline="") as out_handle:
            writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for row in reader:
                scanned += 1
                if len(row) <= max(organism_a, organism_b):
                    continue
                if row[organism_a] == tax_id and row[organism_b] == tax_id:
                    writer.writerow(row)
                    kept += 1

    print(f"Scanned rows: {scanned:,}")
    print(f"Kept rows for taxon {tax_id}: {kept:,}")
    print(f"Output: {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="BioGRID ALL TAB3 input file.")
    parser.add_argument("output", type=Path, help="Arabidopsis-only TAB3 output file.")
    parser.add_argument("--tax-id", default="3702", help="NCBI taxonomy ID to keep.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    filter_biogrid(args.input, args.output, args.tax_id)


if __name__ == "__main__":
    main()
