#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


AGI_PATTERN = re.compile(r"^AT[1-5CM]G\d{5}$", re.IGNORECASE)
UNIPROT_ISOFORM_PATTERN = re.compile(r"^([A-Z0-9]+)-\d+$", re.IGNORECASE)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Map mixed Arabidopsis protein identifiers in an edge list to AGI locus IDs "
            "using a UniProt idmapping.dat file. Self-interactions are retained."
        )
    )

    parser.add_argument(
        "--edges",
        required=True,
        help="Input edge CSV containing protein_a/protein_b or source/target columns.",
    )

    parser.add_argument(
        "--idmapping",
        required=True,
        help="UniProt idmapping.dat file for Arabidopsis.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output mapped edge CSV with protein_a and protein_b columns.",
    )

    parser.add_argument(
        "--report",
        required=True,
        help="Output CSV report describing unmapped and ambiguous identifiers.",
    )

    parser.add_argument(
        "--drop-unmapped",
        action="store_true",
        help="Drop edges where either endpoint cannot be uniquely mapped to an AGI.",
    )

    return parser.parse_args()


def normalize_identifier(identifier):
    """
    Normalize an identifier for lookup.

    Handles:
      - AGI IDs: AT1G01010, ATCG00500, ATMG01190
      - AGI isoforms: AT1G01010.1 -> AT1G01010
      - KEGG IDs: ath:ArthCp058 -> ArthCp058
      - UniProt isoforms: Q8VZJ1-1 -> Q8VZJ1
    """
    if pd.isna(identifier):
        return None

    ident = str(identifier).strip()

    if ident == "":
        return None

    # Remove KEGG species prefix, e.g. ath:ArthCp058
    if ":" in ident:
        ident = ident.split(":")[-1]

    # Normalize AGI transcript/protein isoform suffix, e.g. AT1G01010.1
    if "." in ident:
        possible_base = ident.split(".")[0]
        if AGI_PATTERN.match(possible_base):
            ident = possible_base

    # Normalize UniProt isoform suffix, e.g. Q8VZJ1-1 -> Q8VZJ1
    m = UNIPROT_ISOFORM_PATTERN.match(ident)
    if m:
        ident = m.group(1)

    # AGI IDs should be uppercase
    if AGI_PATTERN.match(ident):
        ident = ident.upper()

    return ident


def extract_agis(value):
    """
    Extract AGI locus IDs from a field in the idmapping file.
    """
    if value is None:
        return set()

    text = str(value)
    candidates = re.split(r"[;,\s]+", text)

    agis = set()
    for candidate in candidates:
        candidate = normalize_identifier(candidate)
        if candidate is not None and AGI_PATTERN.match(candidate):
            agis.add(candidate.upper())

    return agis


def build_identifier_to_agi_map(idmapping_path):
    """
    Build a mapping from mixed identifiers to AGI IDs.

    Strategy:
      1. Group all idmapping.dat entries by UniProt accession.
      2. For each accession, collect AGI IDs from Gene_OrderedLocusName-like rows.
      3. Collect aliases including UniProt accessions and KEGG IDs.
      4. If an accession maps to exactly one AGI, map all aliases to that AGI.
      5. If an alias maps to multiple AGIs, classify it as ambiguous.
    """
    accession_to_agis = defaultdict(set)
    accession_to_aliases = defaultdict(set)

    idmapping_path = Path(idmapping_path)

    with idmapping_path.open("r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            parts = line.split("\t")

            if len(parts) != 3:
                continue

            accession, db, value = parts

            accession_norm = normalize_identifier(accession)
            value_norm = normalize_identifier(value)

            if accession_norm:
                accession_to_aliases[accession].add(accession_norm)

            if value_norm:
                accession_to_aliases[accession].add(value_norm)

            # Any AGI-looking value in any mapping field is useful.
            agis = extract_agis(value)
            accession_to_agis[accession].update(agis)

            # Keep KEGG identifiers such as ath:ArthCp058 and ArthCp058 as aliases.
            if db.upper() == "KEGG":
                if value_norm:
                    accession_to_aliases[accession].add(value_norm)
                if str(value).startswith("ath:"):
                    accession_to_aliases[accession].add(
                        str(value).replace("ath:", "", 1)
                    )

    alias_to_agis = defaultdict(set)

    for accession, agis in accession_to_agis.items():
        if not agis:
            continue

        aliases = accession_to_aliases[accession]

        for alias in aliases:
            alias_norm = normalize_identifier(alias)
            if alias_norm:
                alias_to_agis[alias_norm].update(agis)

    identifier_to_agi = {}
    ambiguous = {}

    for alias, agis in alias_to_agis.items():
        if len(agis) == 1:
            identifier_to_agi[alias] = next(iter(agis))
        else:
            ambiguous[alias] = sorted(agis)

    return identifier_to_agi, ambiguous


def map_identifier(identifier, identifier_to_agi, ambiguous):
    """
    Map a single identifier to an AGI ID.

    Returns:
        mapped_id, status
    """
    ident_norm = normalize_identifier(identifier)

    if ident_norm is None:
        return pd.NA, "missing"

    # Already an AGI.
    if AGI_PATTERN.match(ident_norm):
        return ident_norm.upper(), "already_agi"

    if ident_norm in ambiguous:
        return pd.NA, "ambiguous"

    if ident_norm in identifier_to_agi:
        return identifier_to_agi[ident_norm], "mapped"

    return pd.NA, "unmapped"


def standardize_edge_columns(df):
    """
    Accept either protein_a/protein_b or source/target input columns.
    Return a DataFrame with protein_a/protein_b endpoint columns.
    """
    cols = set(df.columns)

    if {"protein_a", "protein_b"}.issubset(cols):
        return df.copy()

    if {"source", "target"}.issubset(cols):
        df = df.rename(columns={"source": "protein_a", "target": "protein_b"})
        return df.copy()

    raise ValueError(
        "Input edge file must contain either protein_a/protein_b or source/target columns."
    )


def build_mapping_report(edges):
    """
    Build endpoint-level report before dropping unmapped edges.
    """
    report_rows = []

    for endpoint_col, status_col, original_col in [
        ("protein_a", "protein_a_mapping_status", "original_protein_a"),
        ("protein_b", "protein_b_mapping_status", "original_protein_b"),
    ]:
        tmp = (
            edges[[original_col, endpoint_col, status_col]]
            .drop_duplicates()
            .rename(
                columns={
                    original_col: "original_identifier",
                    endpoint_col: "mapped_identifier",
                    status_col: "mapping_status",
                }
            )
        )

        tmp["endpoint"] = endpoint_col
        report_rows.append(tmp)

    report = pd.concat(report_rows, ignore_index=True)

    report = report[
        ["endpoint", "original_identifier", "mapped_identifier", "mapping_status"]
    ].drop_duplicates()

    return report


def count_self_edges(edges):
    """
    Count mapped self-interactions where protein_a == protein_b.

    Missing endpoints are ignored.
    """
    complete_edges = edges.dropna(subset=["protein_a", "protein_b"])
    return int((complete_edges["protein_a"] == complete_edges["protein_b"]).sum())


def canonicalize_mapped_undirected_edges(edges):
    """
    Canonicalize mapped undirected edges after identifier mapping.

    This makes A-B and B-A equivalent while preserving A-A self-interactions.

    Rows with missing endpoints are left unchanged. In the usual pipeline,
    --drop-unmapped is used before this function, so missing endpoints should
    already be gone.
    """
    edges = edges.copy()

    complete_mask = edges["protein_a"].notna() & edges["protein_b"].notna()

    complete_edges = edges.loc[complete_mask].copy()
    incomplete_edges = edges.loc[~complete_mask].copy()

    if not complete_edges.empty:
        edge_min = complete_edges[["protein_a", "protein_b"]].min(axis=1)
        edge_max = complete_edges[["protein_a", "protein_b"]].max(axis=1)

        complete_edges["protein_a"] = edge_min
        complete_edges["protein_b"] = edge_max

    if incomplete_edges.empty:
        return complete_edges

    return pd.concat([complete_edges, incomplete_edges], ignore_index=True)


def main():
    args = parse_args()

    identifier_to_agi, ambiguous = build_identifier_to_agi_map(args.idmapping)

    input_edges = pd.read_csv(args.edges)
    edges = standardize_edge_columns(input_edges)

    mapped_a = edges["protein_a"].apply(
        lambda x: map_identifier(x, identifier_to_agi, ambiguous)
    )
    mapped_b = edges["protein_b"].apply(
        lambda x: map_identifier(x, identifier_to_agi, ambiguous)
    )

    edges["original_protein_a"] = edges["protein_a"]
    edges["original_protein_b"] = edges["protein_b"]

    edges["protein_a"] = [x[0] for x in mapped_a]
    edges["protein_b"] = [x[0] for x in mapped_b]

    edges["protein_a_mapping_status"] = [x[1] for x in mapped_a]
    edges["protein_b_mapping_status"] = [x[1] for x in mapped_b]

    # Build report before dropping anything so unmapped/ambiguous identifiers
    # remain visible in the report.
    report = build_mapping_report(edges)

    if args.drop_unmapped:
        before_drop = len(edges)

        edges = edges.dropna(subset=["protein_a", "protein_b"])

        edges = edges[
            ~edges["protein_a_mapping_status"].isin(["unmapped", "ambiguous", "missing"])
            & ~edges["protein_b_mapping_status"].isin(["unmapped", "ambiguous", "missing"])
        ].copy()

        dropped_mapping_issues = before_drop - len(edges)
    else:
        dropped_mapping_issues = 0

    # Self-interactions are intentionally retained.
    self_edges_retained_before_dedup = count_self_edges(edges)

    # Canonicalize after mapping so A-B and B-A collapse to the same edge.
    # This preserves A-A self-interactions.
    edges = canonicalize_mapped_undirected_edges(edges)

    # Output clean edge list. Preserve extra metadata columns except temporary
    # mapping columns.
    drop_cols = [
        "original_protein_a",
        "original_protein_b",
        "protein_a_mapping_status",
        "protein_b_mapping_status",
    ]

    output_cols = [c for c in edges.columns if c not in drop_cols]
    edges_for_output = edges[output_cols].copy()

    before_dedup = len(edges_for_output)

    # Remove duplicate mapped topological edges.
    # This removes duplicate A-B rows and duplicate A-A rows, but does not
    # remove self-interactions as a class.
    edges_out = edges_for_output.drop_duplicates(
        subset=["protein_a", "protein_b"]
    ).copy()

    duplicate_edges_removed = before_dedup - len(edges_out)
    self_edges_retained_after_dedup = count_self_edges(edges_out)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)

    edges_out.to_csv(args.output, index=False)
    report.to_csv(args.report, index=False)

    print(f"Input edges:                                  {len(input_edges):,}")
    print(f"Output edges:                                 {len(edges_out):,}")
    print(f"Edges dropped due to mapping issues:          {dropped_mapping_issues:,}")
    print(f"Duplicate mapped edges removed:               {duplicate_edges_removed:,}")
    print(f"Self-edges retained before deduplication:     {self_edges_retained_before_dedup:,}")
    print(f"Self-edges retained after deduplication:      {self_edges_retained_after_dedup:,}")

    print()
    print("Endpoint mapping status counts:")
    print(report["mapping_status"].value_counts(dropna=False).to_string())

    print()
    print(f"Wrote mapped edges to:                        {args.output}")
    print(f"Wrote mapping report to:                      {args.report}")


if __name__ == "__main__":
    main()
