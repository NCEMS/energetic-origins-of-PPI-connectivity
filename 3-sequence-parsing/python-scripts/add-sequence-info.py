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


def get_base_locus_id(record_id: str) -> str:
    """
    Convert a TAIR protein isoform ID to a base locus ID.

    Example:
        AT1G01010.1 -> AT1G01010
    """
    return record_id.split(".")[0].upper()


def add_sequences(nodes_df: pd.DataFrame, fasta_file: str) -> pd.DataFrame:
    """
    Add TAIR12 .1 protein sequence information to the node table.

    Node names are expected to be base AGI locus IDs without isoform suffixes,
    for example AT1G01010. The interactome is treated as representing the .1
    isoform, so this function assigns sequence from AT1G01010.1 when available.
    """
    node_col = get_node_column(nodes_df)

    nodes_df[node_col] = nodes_df[node_col].astype(str).str.upper()

    # Read TAIR12 sequences and group them by base locus ID.
    seqs_by_locus = defaultdict(list)

    for record in SeqIO.parse(fasta_file, "fasta"):
        locus_id = get_base_locus_id(record.id)
        seqs_by_locus[locus_id].append(record)

    # Select the .1 isoform for each base locus ID.
    seqs = {}
    loci_without_isoform_1 = {}
    duplicate_isoform_1 = {}

    for locus_id, records in seqs_by_locus.items():
        expected_isoform_1_id = f"{locus_id}.1"

        isoform_1_records = [
            record for record in records
            if record.id.upper() == expected_isoform_1_id
        ]

        if len(isoform_1_records) == 1:
            seqs[locus_id] = str(isoform_1_records[0].seq).rstrip("*")

        elif len(isoform_1_records) > 1:
            duplicate_isoform_1[locus_id] = [record.id for record in isoform_1_records]

        else:
            loci_without_isoform_1[locus_id] = [record.id for record in records]

    if duplicate_isoform_1:
        examples = list(duplicate_isoform_1.items())[:10]
        example_text = "\n".join(
            f"{locus}: {', '.join(ids)}" for locus, ids in examples
        )
        raise ValueError(
            "Multiple TAIR12 .1 records were found for at least one locus. "
            "FASTA record IDs should be unique.\n"
            f"Examples:\n{example_text}"
        )

    # Report nodes with FASTA records but no matching .1 isoform.
    nodes_missing_isoform_1 = {
        node: loci_without_isoform_1[node]
        for node in nodes_df[node_col]
        if node in loci_without_isoform_1
    }

    if nodes_missing_isoform_1:
        print("Warning: some nodes had TAIR12 sequence records but no matching .1 isoform.")
        print(f"Number of affected nodes: {len(nodes_missing_isoform_1):,}")
        print("Examples:")

        for node, isoform_ids in list(nodes_missing_isoform_1.items())[:10]:
            print(f"{node}: expected {node}.1; found {', '.join(isoform_ids)}")

    nodes_df["has_verified_sequence"] = nodes_df[node_col].isin(seqs.keys())

    nodes_df["sequence"] = nodes_df[node_col].apply(
        lambda x: seqs[x] if x in seqs else np.nan
    )

    nodes_df["tair12_isoform_id"] = nodes_df[node_col].apply(
        lambda x: f"{x}.1" if x in seqs else np.nan
    )

    nodes_df["tair12_sequence_length"] = nodes_df["sequence"].apply(
        lambda x: len(x) if isinstance(x, str) else np.nan
    )

    print(f"Input nodes: {len(nodes_df):,}")
    print(f"Nodes with TAIR12 .1 sequence: {nodes_df['has_verified_sequence'].sum():,}")
    print(f"Nodes without TAIR12 .1 sequence: {(~nodes_df['has_verified_sequence']).sum():,}")

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
      4. If no exact match, choose the candidate with highest global sequence identity.
      5. Break ties by reviewed Swiss-Prot status, smaller length difference,
         then accession string for determinism.
    """
    candidates = sorted(agi_to_uniprot_candidates.get(agi, set()))

    result = {
        "UniProtKB-AC": np.nan,
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

    result["UniProtKB-AC"] = best["accession"]
    result["uniprot_exact_sequence_match"] = bool(best["exact_match"])
    result["uniprot_best_sequence_identity"] = round(float(best["identity"]), 6)
    result["uniprot_best_sequence_length"] = best["length"]
    result["uniprot_best_reviewed"] = bool(best["reviewed"])

    if best["exact_match"]:
        result["uniprot_match_status"] = "exact_sequence_match"
    else:
        result["uniprot_match_status"] = "best_nonexact_sequence_match"

    return result


def add_sequence_vetted_uniprot_mappings(
    nodes_df: pd.DataFrame,
    map_file: str,
    uniprot_fasta: str,
) -> pd.DataFrame:
    """
    Add one sequence-vetted UniProt accession per AGI node.

    The selected accession is written to:
        UniProtKB-AC

    All candidates are written to:
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
            "Add TAIR12 .1 sequences to PPI network nodes and select one "
            "sequence-vetted UniProt accession per AGI node."
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
