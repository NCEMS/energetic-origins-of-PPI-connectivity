#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner


AGI_PATTERN = re.compile(r"^AT[1-5CM]G\d{5}$", re.IGNORECASE)
AGI_ISOFORM_PATTERN = re.compile(r"^(AT[1-5CM]G\d{5})(?:\.(\d+))?$", re.IGNORECASE)
UNIPROT_ISOFORM_PATTERN = re.compile(r"^([A-Z0-9]+)-\d+$", re.IGNORECASE)


def normalize_agi(value):
    """
    Normalize AGI locus or AGI isoform identifiers.

    Examples:
        AT1G01010   -> AT1G01010
        AT1G01010.1 -> AT1G01010
    """
    if pd.isna(value):
        return None

    value = str(value).strip()

    if value == "":
        return None

    if "." in value:
        base = value.split(".")[0]
    else:
        base = value

    if AGI_PATTERN.match(base):
        return base.upper()

    return None


def normalize_uniprot_accession(value):
    """
    Normalize UniProt accession identifiers.

    Examples:
        Q8VZJ1   -> Q8VZJ1
        Q8VZJ1-1 -> Q8VZJ1
    """
    if pd.isna(value):
        return None

    value = str(value).strip()

    if value == "":
        return None

    match = UNIPROT_ISOFORM_PATTERN.match(value)
    if match:
        value = match.group(1)

    return value


def get_node_column(nodes_df):
    """
    Identify the node identifier column.

    Earlier parsing scripts may emit either:
        node
    or:
        name
    """
    if "node" in nodes_df.columns:
        return "node"

    if "name" in nodes_df.columns:
        return "name"

    raise ValueError("Nodes file must contain either a 'node' or 'name' column.")


def parse_tair_isoform_record_id(record_id: str):
    """
    Parse a TAIR protein FASTA record ID into base AGI locus and isoform number.

    Examples:
        AT1G01010.1 -> ("AT1G01010", 1)
        AT1G01010.2 -> ("AT1G01010", 2)
        AT1G01010   -> ("AT1G01010", None)

    Returns:
        locus_id, isoform_number
    """
    record_id = str(record_id).strip().split()[0].upper()

    match = AGI_ISOFORM_PATTERN.match(record_id)

    if match is None:
        return None, None

    locus_id = match.group(1).upper()
    isoform_number = int(match.group(2)) if match.group(2) is not None else None

    return locus_id, isoform_number


def get_base_locus_id(record_id: str) -> str:
    """
    Convert a TAIR protein isoform ID to a base locus ID.

    Example:
        AT1G01010.1 -> AT1G01010
    """
    return record_id.split(".")[0].upper()


def add_sequences(nodes_df: pd.DataFrame, fasta_file: str) -> pd.DataFrame:
    """
    Add TAIR12 protein sequence information to the node table.

    Node names are expected to be base AGI locus IDs without isoform suffixes,
    for example AT1G01010.

    Sequence selection rule:
      1. Use AT1G01010.1 if present.
      2. If .1 is absent, fall back to the lowest available numbered isoform:
         .2, then .3, etc.
      3. If only an unsuffixed record exists, use it as a final fallback.

    The selected FASTA record ID is stored in:
        tair12_isoform_id
    """
    node_col = get_node_column(nodes_df)

    nodes_df[node_col] = nodes_df[node_col].astype(str).str.upper()

    # Read TAIR12 sequences and group them by base locus ID.
    seqs_by_locus = defaultdict(list)
    malformed_record_ids = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        locus_id, isoform_number = parse_tair_isoform_record_id(record.id)

        if locus_id is None:
            malformed_record_ids.append(record.id)
            continue

        seqs_by_locus[locus_id].append(
            {
                "record_id": record.id.upper(),
                "isoform_number": isoform_number,
                "sequence": str(record.seq).rstrip("*"),
            }
        )

    if malformed_record_ids:
        print("Warning: some FASTA record IDs could not be parsed as TAIR AGI isoforms.")
        print(f"Number of malformed FASTA record IDs: {len(malformed_record_ids):,}")
        print("Examples:")
        for record_id in malformed_record_ids[:10]:
            print(f"  {record_id}")

    # Select the lowest available numbered isoform for each locus.
    seqs = {}
    selected_isoform_ids = {}
    selected_isoform_numbers = {}
    duplicate_numbered_isoforms = {}

    for locus_id, records in seqs_by_locus.items():
        numbered_records = [
            record for record in records
            if record["isoform_number"] is not None
        ]

        # Check for duplicate records for the same locus + numbered isoform.
        records_by_isoform = defaultdict(list)

        for record in numbered_records:
            records_by_isoform[record["isoform_number"]].append(record["record_id"])

        duplicate_isoforms = {
            isoform_number: record_ids
            for isoform_number, record_ids in records_by_isoform.items()
            if len(record_ids) > 1
        }

        if duplicate_isoforms:
            duplicate_numbered_isoforms[locus_id] = duplicate_isoforms
            continue

        if numbered_records:
            selected = sorted(
                numbered_records,
                key=lambda record: (
                    record["isoform_number"],
                    record["record_id"],
                ),
            )[0]
        else:
            # Final fallback for unusual FASTA records that have a valid AGI
            # but no isoform suffix.
            selected = sorted(
                records,
                key=lambda record: record["record_id"],
            )[0]

        seqs[locus_id] = selected["sequence"]
        selected_isoform_ids[locus_id] = selected["record_id"]
        selected_isoform_numbers[locus_id] = selected["isoform_number"]

    if duplicate_numbered_isoforms:
        examples = list(duplicate_numbered_isoforms.items())[:10]
        example_text = []

        for locus_id, dup_info in examples:
            isoform_text = "; ".join(
                f".{isoform_number}: {', '.join(record_ids)}"
                for isoform_number, record_ids in dup_info.items()
            )
            example_text.append(f"{locus_id}: {isoform_text}")

        raise ValueError(
            "Duplicate TAIR12 records were found for the same locus and isoform. "
            "FASTA record IDs should be unique.\n"
            "Examples:\n"
            + "\n".join(example_text)
        )

    nodes_df["has_verified_sequence"] = nodes_df[node_col].isin(seqs.keys())

    nodes_df["sequence"] = nodes_df[node_col].apply(
        lambda x: seqs[x] if x in seqs else np.nan
    )

    nodes_df["tair12_isoform_id"] = nodes_df[node_col].apply(
        lambda x: selected_isoform_ids[x] if x in selected_isoform_ids else np.nan
    )

    nodes_df["tair12_sequence_length"] = nodes_df["sequence"].apply(
        lambda x: len(x) if isinstance(x, str) else np.nan
    )

    nodes_with_isoform_1 = nodes_df[node_col].apply(
        lambda x: selected_isoform_numbers.get(x) == 1
    )

    nodes_with_fallback_isoform = nodes_df[node_col].apply(
        lambda x: (
            x in selected_isoform_numbers
            and selected_isoform_numbers.get(x) != 1
        )
    )

    if nodes_with_fallback_isoform.any():
        print("Warning: some nodes lacked TAIR12 .1 records and used fallback isoforms.")
        print(f"Number of affected nodes: {nodes_with_fallback_isoform.sum():,}")
        print("Examples:")

        fallback_nodes = nodes_df.loc[nodes_with_fallback_isoform, node_col].head(10)

        for node in fallback_nodes:
            available_ids = sorted(
                record["record_id"]
                for record in seqs_by_locus.get(node, [])
            )
            selected_id = selected_isoform_ids.get(node)
            print(
                f"{node}: expected {node}.1; "
                f"selected {selected_id}; "
                f"available {', '.join(available_ids)}"
            )

    print(f"Input nodes: {len(nodes_df):,}")
    print(f"Nodes with selected TAIR12 sequence: {nodes_df['has_verified_sequence'].sum():,}")
    print(f"Nodes using TAIR12 .1 sequence: {nodes_with_isoform_1.sum():,}")
    print(f"Nodes using fallback TAIR12 isoform: {nodes_with_fallback_isoform.sum():,}")
    print(f"Nodes without any TAIR12 sequence: {(~nodes_df['has_verified_sequence']).sum():,}")

    return nodes_df


def extract_agis(value):
    """
    Extract AGI locus IDs from an idmapping field.
    """
    if pd.isna(value):
        return set()

    tokens = re.split(r"[;,\s]+", str(value).strip())

    agis = set()
    for token in tokens:
        agi = normalize_agi(token)
        if agi is not None:
            agis.add(agi)

    return agis


def build_agi_to_uniprot_candidates(map_file: str) -> dict:
    """
    Build AGI -> set(UniProt accessions) from UniProt idmapping.dat.

    This focuses on Gene_OrderedLocusName rows, which are the AGI mappings.
    """
    column_names = ["UniProtKB-AC", "ID_type", "ID"]

    cross_df = pd.read_csv(
        map_file,
        names=column_names,
        sep="\t",
        dtype=str,
    )

    cross_df = cross_df.dropna(subset=["UniProtKB-AC", "ID_type", "ID"]).copy()

    # Use Gene_OrderedLocusName rows as the primary AGI mapping source.
    cross_df = cross_df[cross_df["ID_type"] == "Gene_OrderedLocusName"].copy()

    agi_to_uniprot = defaultdict(set)

    for _, row in cross_df.iterrows():
        accession = normalize_uniprot_accession(row["UniProtKB-AC"])
        agis = extract_agis(row["ID"])

        if accession is None:
            continue

        for agi in agis:
            agi_to_uniprot[agi].add(accession)

    n_multi_mapping_ids = sum(
        len(accessions) > 1 for accessions in agi_to_uniprot.values()
    )

    print(f"AGIs with at least one UniProtKB-AC candidate: {len(agi_to_uniprot):,}")
    print(f"AGIs with multiple UniProtKB-AC candidates: {n_multi_mapping_ids:,}")

    return dict(agi_to_uniprot)


def parse_uniprot_accession_from_header(header: str):
    """
    Parse UniProt FASTA headers.

    Typical headers:
        sp|P93033|...
        tr|A0A178U7Y4|...

    Returns:
        accession, reviewed
    """
    first_token = header.split()[0]

    if "|" in first_token:
        parts = first_token.split("|")
        if len(parts) >= 2:
            database_prefix = parts[0].lower()
            accession = normalize_uniprot_accession(parts[1])
            reviewed = database_prefix == "sp"
            return accession, reviewed

    accession = normalize_uniprot_accession(first_token)
    reviewed = False

    return accession, reviewed


def parse_uniprot_fasta(uniprot_fasta: str) -> dict:
    """
    Parse UniProt proteome FASTA into accession -> record metadata.

    If duplicate normalized accessions are encountered, the longest sequence is
    retained. This is conservative for the usual accession-level proteome FASTA.
    """
    accession_to_record = {}

    for record in SeqIO.parse(uniprot_fasta, "fasta"):
        accession, reviewed = parse_uniprot_accession_from_header(record.description)

        if accession is None:
            continue

        sequence = str(record.seq).rstrip("*")

        new_record = {
            "accession": accession,
            "sequence": sequence,
            "reviewed": reviewed,
            "description": record.description,
            "length": len(sequence),
        }

        if accession not in accession_to_record:
            accession_to_record[accession] = new_record
        else:
            old_record = accession_to_record[accession]

            # Deterministic tie handling if the same accession appears more than once.
            if (
                len(sequence) > old_record["length"]
                or (
                    len(sequence) == old_record["length"]
                    and reviewed
                    and not old_record["reviewed"]
                )
            ):
                accession_to_record[accession] = new_record

    print(f"UniProt FASTA accessions parsed: {len(accession_to_record):,}")

    return accession_to_record


def global_sequence_identity(seq_a: str, seq_b: str) -> float:
    """
    Estimate global sequence identity for two protein sequences.

    Exact matches are handled before alignment. For non-exact matches, this uses
    Bio.Align.PairwiseAligner with match-only scoring and normalizes by the
    longer sequence length.

    This is sufficient for ranking candidate UniProt accessions for a single
    AGI without adding an external MUSCLE dependency.
    """
    if not isinstance(seq_a, str) or not isinstance(seq_b, str):
        return np.nan

    if len(seq_a) == 0 and len(seq_b) == 0:
        return 1.0

    if len(seq_a) == 0 or len(seq_b) == 0:
        return 0.0

    if seq_a == seq_b:
        return 1.0

    if len(seq_a) == len(seq_b):
        matches = sum(a == b for a, b in zip(seq_a, seq_b))
        return matches / len(seq_a)

    aligner = PairwiseAligner()
    aligner.mode = "global"

    # With these scores, aligner.score is the number of matched residues in the
    # best global alignment. Normalizing by max length penalizes gaps.
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = 0.0
    aligner.extend_gap_score = 0.0

    matches = aligner.score(seq_a, seq_b)
    return float(matches) / max(len(seq_a), len(seq_b))


def choose_best_uniprot_for_node(
    agi: str,
    tair_sequence,
    agi_to_uniprot_candidates: dict,
    uniprot_records: dict,
) -> dict:
    """
    Choose the best UniProt accession for one AGI node.

    Rules:
      1. Get UniProt candidates from idmapping.dat.
      2. Retrieve candidate sequences from UniProt FASTA.
      3. Prefer exact TAIR12 .1 sequence matches.
      4. If one or more exact matches exist, choose among them using reviewed
         Swiss-Prot status and accession string for deterministic tie-breaking.
      5. If no exact match exists, retain the best nonexact candidate only for
         diagnostics and leave UniProtKB-AC unpopulated.
    """
    candidates = sorted(agi_to_uniprot_candidates.get(agi, set()))

    result = {
        "UniProtKB-AC": np.nan,
        "uniprot_best_candidate_accession": np.nan,
        "UniProtKB-AC_candidates": ";".join(candidates) if candidates else np.nan,
        "uniprot_candidate_count": len(candidates),
        "uniprot_candidate_count_with_sequence": 0,
        "uniprot_exact_sequence_match": pd.NA,
        "uniprot_best_sequence_identity": np.nan,
        "uniprot_match_status": "no_uniprot_candidates",
        "uniprot_best_sequence_length": np.nan,
        "uniprot_best_reviewed": pd.NA,
    }

    if not candidates:
        return result

    if not isinstance(tair_sequence, str) or tair_sequence == "":
        result["uniprot_match_status"] = "no_tair12_sequence"
        return result

    candidate_records = []

    for accession in candidates:
        record = uniprot_records.get(accession)

        if record is None:
            continue

        candidate_records.append(record)

    result["uniprot_candidate_count_with_sequence"] = len(candidate_records)

    if not candidate_records:
        result["uniprot_match_status"] = "no_candidate_sequence_in_uniprot_fasta"
        return result

    scored_candidates = []

    for record in candidate_records:
        uniprot_sequence = record["sequence"]
        exact_match = tair_sequence == uniprot_sequence
        identity = global_sequence_identity(tair_sequence, uniprot_sequence)
        length_delta = abs(len(tair_sequence) - len(uniprot_sequence))

        scored_candidates.append(
            {
                "accession": record["accession"],
                "sequence": uniprot_sequence,
                "reviewed": record["reviewed"],
                "length": record["length"],
                "exact_match": exact_match,
                "identity": identity,
                "length_delta": length_delta,
            }
        )

    scored_candidates = sorted(
        scored_candidates,
        key=lambda x: (
            not x["exact_match"],       # exact matches first
            -float(x["identity"]),      # highest identity first
            not bool(x["reviewed"]),    # reviewed Swiss-Prot first
            x["length_delta"],          # closest length first
            x["accession"],             # deterministic final tie-break
        ),
    )

    best = scored_candidates[0]

    # Always retain information about the best-ranked candidate for auditing,
    # even when it is not an exact sequence match.
    result["uniprot_best_candidate_accession"] = best["accession"]
    result["uniprot_exact_sequence_match"] = bool(best["exact_match"])
    result["uniprot_best_sequence_identity"] = round(float(best["identity"]), 6)
    result["uniprot_best_sequence_length"] = best["length"]
    result["uniprot_best_reviewed"] = bool(best["reviewed"])

    # Only assign the canonical UniProt accession when its protein sequence
    # exactly matches the selected TAIR12 sequence for this AGI.
    if best["exact_match"]:
        result["UniProtKB-AC"] = best["accession"]
        result["uniprot_match_status"] = "exact_sequence_match"
    else:
        result["UniProtKB-AC"] = np.nan
        result["uniprot_match_status"] = "no_exact_sequence_match"

    return result


def add_sequence_vetted_uniprot_mappings(
    nodes_df: pd.DataFrame,
    map_file: str,
    uniprot_fasta: str,
) -> pd.DataFrame:
    """
    Add one exact-sequence-matched UniProt accession per AGI node.

    UniProtKB-AC is populated only when a candidate UniProt sequence is exactly
    identical to the selected TAIR12 sequence.

    The best-ranked candidate accession, including nonexact candidates, is
    retained separately in:
        uniprot_best_candidate_accession

    All mapped candidates are written to:
        UniProtKB-AC_candidates
    """
    node_col = get_node_column(nodes_df)

    n_before = len(nodes_df)

    agi_to_uniprot_candidates = build_agi_to_uniprot_candidates(map_file)
    uniprot_records = parse_uniprot_fasta(uniprot_fasta)

    mapping_rows = []

    for _, row in nodes_df.iterrows():
        agi = normalize_agi(row[node_col])
        tair_sequence = row["sequence"]

        if agi is None:
            choice = {
                "UniProtKB-AC": np.nan,
                "uniprot_best_candidate_accession": np.nan,
                "UniProtKB-AC_candidates": np.nan,
                "uniprot_candidate_count": 0,
                "uniprot_candidate_count_with_sequence": 0,
                "uniprot_exact_sequence_match": pd.NA,
                "uniprot_best_sequence_identity": np.nan,
                "uniprot_match_status": "invalid_agi_node_id",
                "uniprot_best_sequence_length": np.nan,
                "uniprot_best_reviewed": pd.NA,
            }
        else:
            choice = choose_best_uniprot_for_node(
                agi=agi,
                tair_sequence=tair_sequence,
                agi_to_uniprot_candidates=agi_to_uniprot_candidates,
                uniprot_records=uniprot_records,
            )

        choice[node_col] = row[node_col]
        mapping_rows.append(choice)

    mapping_df = pd.DataFrame(mapping_rows)

    nodes_df = nodes_df.merge(
        mapping_df,
        on=node_col,
        how="left",
        validate="one_to_one",
    )

    n_after = len(nodes_df)

    if n_after != n_before:
        raise ValueError(
            f"Row count changed during UniProt mapping merge: "
            f"{n_before:,} -> {n_after:,}. This should not happen."
        )

    print()
    print("UniProt sequence-mapping status counts:")
    print(nodes_df["uniprot_match_status"].value_counts(dropna=False).to_string())

    print()
    print("UniProt exact sequence-match counts:")
    print(nodes_df["uniprot_exact_sequence_match"].value_counts(dropna=False).to_string())

    return nodes_df


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Add TAIR12 protein sequences to PPI network nodes and assign "
            "a UniProt accession only when its sequence exactly matches the "
            "selected TAIR12 sequence."
        )
    )

    parser.add_argument(
        "--nodes",
        required=True,
        help="Path to the nodes CSV file.",
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to TAIR12 protein FASTA file.",
    )

    parser.add_argument(
        "--uniprot_fasta",
        required=True,
        help="Path to UniProt proteome FASTA file, e.g. UP000006548_3702.fasta.",
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Path to output directory.",
    )

    parser.add_argument(
        "--output_prefix",
        required=True,
        help="Prefix to be applied to output file.",
    )

    parser.add_argument(
        "--output_suffix",
        required=True,
        help="Suffix to be applied to output file.",
    )

    parser.add_argument(
        "--id_mappings",
        required=True,
        help="UniProt idmapping.dat file.",
    )

    parser.add_argument(
        "--organism_tag",
        required=True,
        help="Tag to label the organism for this run.",
    )

    args = parser.parse_args()

    nodes_df = pd.read_csv(args.nodes)

    nodes_df = add_sequences(nodes_df, args.fasta)

    nodes_df = add_sequence_vetted_uniprot_mappings(
        nodes_df=nodes_df,
        map_file=args.id_mappings,
        uniprot_fasta=args.uniprot_fasta,
    )

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    output_file = (
        f"{args.output_dir}/"
        f"{args.output_prefix}-{args.organism_tag}-seqs-{args.output_suffix}.csv"
    )

    nodes_df.to_csv(output_file, index=False)

    print()
    print(f"Wrote sequence-annotated nodes to: {output_file}")


if __name__ == "__main__":
    main()
