#!/usr/bin/env python3

import os
import argparse
import typing
from typing import Optional

import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa


AA_MAP = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


class CleavageSelect(Select):
    """
    Bio.PDB residue selector used to remove N-terminal residues through cut_pos.

    Residues are retained when residue.id[1] > cut_pos.
    """

    def __init__(self, cut_pos: int):
        self.cut_pos = cut_pos

    def accept_residue(self, residue):
        return residue.id[1] > self.cut_pos


def is_missing(value) -> bool:
    """
    Safely test for missing scalar values.
    """
    if value is None:
        return True

    try:
        missing = pd.isna(value)
    except TypeError:
        return False

    if isinstance(missing, (bool, np.bool_)):
        return bool(missing)

    return False


def clean_sequence(seq) -> Optional[str]:
    """
    Return a clean uppercase protein sequence string or None.
    """
    if not isinstance(seq, str):
        return None

    seq = seq.strip().upper().rstrip("*")

    if seq == "":
        return None

    return seq


def parse_uniprot_accessions(value) -> list[str]:
    """
    Parse a UniProtKB-AC field.

    The updated sequence-parsing step should provide a single selected accession
    in UniProtKB-AC. This parser remains tolerant of older semicolon-delimited
    files for backward compatibility.

    Examples:
        Q9SP32 -> ["Q9SP32"]
        F4HQG6;Q56Y50;Q9SP32 -> ["F4HQG6", "Q56Y50", "Q9SP32"]
        np.nan -> []
    """
    if is_missing(value):
        return []

    accessions = []

    for accession in str(value).split(";"):
        accession = accession.strip()

        if not accession:
            continue

        if accession.lower() in {"nan", "none", "na", "null", "<na>"}:
            continue

        accessions.append(accession)

    return accessions


def choose_selected_uniprot_accession(value) -> Optional[str]:
    """
    Choose the selected UniProt accession for structure lookup.

    In the updated pipeline, UniProtKB-AC should contain one sequence-vetted
    accession. If an older semicolon-delimited value is encountered, this
    function uses the first accession and diagnostics will report that case.
    """
    accessions = parse_uniprot_accessions(value)

    if not accessions:
        return None

    return accessions[0]


def truncate_fasta(seq: str, cut: Optional[int]) -> Optional[str]:
    """
    Create a truncated sequence based on a SignalP predicted cleavage site.

    Args:
        seq: sequence to be trimmed
        cut: number of residues to remove from the N-terminus, or None/NaN

    Returns:
        Truncated sequence, original sequence, or None if seq/cut are invalid.
    """
    seq = clean_sequence(seq)

    if seq is None:
        return None

    if is_missing(cut):
        return seq

    try:
        cut = int(cut)
    except (TypeError, ValueError):
        return None

    if cut < 0:
        return None

    return seq[cut:]


def truncate_structure(pdb_path: str, cut_pos: int, output_path: str) -> None:
    """
    Write a cleaved PDB file by removing residues through cut_pos.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, CleavageSelect(cut_pos))


def truncate_row_structure(row) -> Optional[str]:
    """
    Create a cleaved PDB file if a structure exists and a SignalP cleavage site
    is available.
    """
    if row.get("structure_exists") != 1:
        return None

    if is_missing(row.get("cleavage_site_start")):
        return None

    input_path = row.get("structure_path")

    if not isinstance(input_path, str) or not os.path.isfile(input_path):
        return None

    try:
        cut_pos = int(row["cleavage_site_start"])
    except (TypeError, ValueError):
        return None

    output_path = input_path.replace(".pdb", "-cleaved.pdb")

    try:
        truncate_structure(input_path, cut_pos, output_path)
        return output_path
    except Exception as e:
        print(f"Error cleaving structure {input_path}: {e}")
        return None


def extract_sequence_from_pdb(pdb_path: str) -> Optional[str]:
    """
    Extract the full protein sequence from a PDB by walking residues in order.

    This is used as a fallback when the AF2 FASTA file is unavailable.
    """
    if not isinstance(pdb_path, (str, bytes, os.PathLike)):
        return None

    if not os.path.isfile(pdb_path):
        return None

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)

    model = structure[0]

    seq_chars = []
    seen = set()

    for chain in model:
        for residue in chain:
            hetflag, resseq, icode = residue.id
            key = (chain.id, resseq, icode)

            if key in seen:
                continue

            seen.add(key)

            if not is_aa(residue, standard=False):
                continue

            resname = residue.get_resname().upper()
            aa = AA_MAP.get(resname, "X")
            seq_chars.append(aa)

    if not seq_chars:
        return None

    return "".join(seq_chars)


def read_fasta_sequence(fasta_path: str) -> Optional[str]:
    """
    Read the first sequence stored in a FASTA file.
    """
    if is_missing(fasta_path):
        return None

    if not isinstance(fasta_path, str):
        return None

    if not os.path.isfile(fasta_path):
        return None

    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        return clean_sequence(str(record.seq))
    except Exception as e:
        print(f"Error reading FASTA {fasta_path}: {e}")
        return None


def locate_structure(nodes_df: pd.DataFrame, structure_dir: str) -> pd.DataFrame:
    """
    Add AF2-EBI structure path information using the selected UniProtKB-AC.

    The updated pipeline treats UniProtKB-AC as the selected, sequence-vetted
    accession chosen during sequence parsing.
    """
    nodes_df = nodes_df.copy()

    if "UniProtKB-AC" not in nodes_df.columns:
        raise ValueError(
            "Input nodes file does not contain UniProtKB-AC. "
            "Run the updated sequence-parsing step before prepare-structures.py."
        )

    nodes_df["n_UniProtKB_AC"] = nodes_df["UniProtKB-AC"].apply(
        lambda x: len(parse_uniprot_accessions(x))
    )

    nodes_df["selected_UniProtKB_AC_for_structure"] = nodes_df["UniProtKB-AC"].apply(
        choose_selected_uniprot_accession
    )

    nodes_df["structure_path"] = nodes_df["selected_UniProtKB_AC_for_structure"].apply(
        lambda accession: (
            os.path.join(structure_dir, f"AF-{accession}-F1-model_v6.pdb")
            if accession is not None
            else None
        )
    )

    nodes_df["structure_exists"] = nodes_df["structure_path"].apply(
        lambda path: 1 if isinstance(path, str) and os.path.isfile(path) else 0
    )

    n_missing_uniprot = nodes_df["n_UniProtKB_AC"].eq(0).sum()
    n_multiple_uniprot = nodes_df["n_UniProtKB_AC"].gt(1).sum()
    n_existing_structures = nodes_df["structure_exists"].sum()

    print("Structure accession diagnostics:")
    print(f"Nodes with no selected UniProtKB-AC:        {n_missing_uniprot:,}")
    print(f"Nodes with >1 UniProtKB-AC value:           {n_multiple_uniprot:,}")
    print(f"Nodes with existing selected AF2-EBI PDB:   {n_existing_structures:,}")

    if n_multiple_uniprot > 0:
        print()
        print(
            "Warning: some rows still contain multiple UniProtKB-AC values. "
            "This script will use the first accession only."
        )

        example_rows = nodes_df.loc[
            nodes_df["n_UniProtKB_AC"].gt(1),
            ["node", "UniProtKB-AC", "selected_UniProtKB_AC_for_structure"],
        ].head(10)

        for _, row in example_rows.iterrows():
            print(
                f"{row['node']}: UniProtKB-AC={row['UniProtKB-AC']}; "
                f"selected={row['selected_UniProtKB_AC_for_structure']}"
            )

    return nodes_df


def locate_structure_sequence(nodes_df: pd.DataFrame, structure_dir: str) -> pd.DataFrame:
    """
    Add AF2-derived FASTA path and structure sequence information.

    Sequence source priority:
      1. AF-{accession}-F1-model_v6.fasta if present
      2. extracted PDB sequence if FASTA is absent
    """
    nodes_df = nodes_df.copy()

    nodes_df["structure_fasta_path"] = nodes_df[
        "selected_UniProtKB_AC_for_structure"
    ].apply(
        lambda accession: (
            os.path.join(structure_dir, f"AF-{accession}-F1-model_v6.fasta")
            if accession is not None
            else None
        )
    )

    structure_sequences = []
    sequence_sources = []

    for _, row in nodes_df.iterrows():
        fasta_path = row.get("structure_fasta_path")
        pdb_path = row.get("structure_path")

        fasta_seq = read_fasta_sequence(fasta_path)

        if fasta_seq is not None:
            structure_sequences.append(fasta_seq)
            sequence_sources.append("fasta")
            continue

        pdb_seq = extract_sequence_from_pdb(pdb_path)

        if pdb_seq is not None:
            structure_sequences.append(pdb_seq)
            sequence_sources.append("pdb")
            continue

        structure_sequences.append(None)
        sequence_sources.append(None)

    nodes_df["structure_sequence"] = structure_sequences
    nodes_df["structure_sequence_source"] = sequence_sources

    return nodes_df


def expected_final_sequence(row) -> Optional[str]:
    """
    Determine the sequence the final structure should represent.

    Prefer signalP_trimmed_sequence if present. Otherwise, fall back to sequence
    and cleavage_site_start.
    """
    signalp_seq = clean_sequence(row.get("signalP_trimmed_sequence"))

    if signalp_seq is not None:
        return signalp_seq

    full_seq = clean_sequence(row.get("sequence"))

    if full_seq is None:
        return None

    return truncate_fasta(full_seq, row.get("cleavage_site_start"))


def expected_full_sequence(row) -> Optional[str]:
    """
    Determine the full pre-cleavage sequence for a node.
    """
    return clean_sequence(row.get("sequence"))


def first_mismatch_index(seq_a: str, seq_b: str) -> Optional[int]:
    """
    Return the first zero-based mismatch index between two strings.
    """
    if not isinstance(seq_a, str) or not isinstance(seq_b, str):
        return None

    for i, (a, b) in enumerate(zip(seq_a, seq_b)):
        if a != b:
            return i

    if len(seq_a) != len(seq_b):
        return min(len(seq_a), len(seq_b))

    return None


def compare_structure_to_expected(row) -> pd.Series:
    """
    Compare the candidate structure sequence to the expected final sequence.

    For cleaved proteins, this compares cleaved_structure_sequence against the
    expected final sequence. For uncleaved proteins, it compares structure_sequence
    directly against the expected final sequence.
    """
    expected_seq = expected_final_sequence(row)

    if is_missing(row.get("cleavage_site_start")):
        struct_seq = clean_sequence(row.get("structure_sequence"))
    else:
        struct_seq = clean_sequence(row.get("cleaved_structure_sequence"))

    if struct_seq is None and expected_seq is None:
        return pd.Series([False, "missing_structure_and_expected_sequence", np.nan, np.nan, np.nan])

    if struct_seq is None:
        return pd.Series([False, "missing_structure_sequence", np.nan, len(expected_seq), np.nan])

    if expected_seq is None:
        return pd.Series([False, "missing_expected_sequence", len(struct_seq), np.nan, np.nan])

    match = struct_seq == expected_seq

    if match:
        return pd.Series([True, "match", len(struct_seq), len(expected_seq), np.nan])

    mismatch_idx = first_mismatch_index(struct_seq, expected_seq)

    return pd.Series(
        [
            False,
            "mismatch",
            len(struct_seq),
            len(expected_seq),
            mismatch_idx,
        ]
    )


def extract_plddt_ca_only(pdb_filename):
    """
    Extract pLDDT values from CA atoms in the B-factor field of an AF2 PDB file.
    """
    if is_missing(pdb_filename):
        return np.nan

    if not isinstance(pdb_filename, (str, bytes, os.PathLike)):
        return np.nan

    if not os.path.isfile(pdb_filename):
        return np.nan

    plddt_values = []

    with open(pdb_filename, "r") as handle:
        for line in handle:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    b_factor = float(line[60:66].strip())
                    plddt_values.append(b_factor)
                except ValueError:
                    continue

    if len(plddt_values) == 0:
        return np.nan

    return float(np.mean(plddt_values))


def compute_mean_plddt(row):
    """
    Compute mean pLDDT from the final selected structure.
    """
    path = row.get("final_structure_path")

    if is_missing(path):
        return np.nan

    if not isinstance(path, (str, bytes, os.PathLike)):
        return np.nan

    if os.path.isfile(path):
        return extract_plddt_ca_only(path)

    return np.nan


def locate_rescue_structure(
    node: str,
    af2_dir: str,
) -> typing.Tuple[Optional[str], Optional[str]]:
    """
    Locate an additional/custom AF2 rescue structure and extract sequence.

    Expected layout:
        {af2_dir}/{node}/ranked_0.pdb
    """
    if not af2_dir:
        return None, None

    if not isinstance(node, str):
        return None, None

    pdb_path = os.path.join(af2_dir, node, "ranked_0.pdb")

    if not os.path.isfile(pdb_path):
        return None, None

    try:
        af2_seq = extract_sequence_from_pdb(pdb_path)

        if af2_seq is None:
            print(f"No sequence could be extracted from rescue AF2 structure for {node}")
            return None, None

        return pdb_path, af2_seq

    except Exception as e:
        print(f"Error reading rescue AF2 structure for {node}: {e}")
        return None, None


def cleave_rescue_structure_if_needed(
    pdb_path: str,
    node: str,
    cut_pos,
) -> Optional[str]:
    """
    Cleave a rescue AF2 PDB if it represents the full sequence but the expected
    final sequence is SignalP-trimmed.
    """
    if is_missing(cut_pos):
        return pdb_path

    try:
        cut_pos = int(cut_pos)
    except (TypeError, ValueError):
        return None

    output_path = pdb_path.replace(".pdb", "-cleaved.pdb")

    try:
        truncate_structure(pdb_path, cut_pos, output_path)
        return output_path
    except Exception as e:
        print(f"[{node}] Error cleaving rescue AF2 structure {pdb_path}: {e}")
        return None


def select_final_structure(row, af2_dir):
    """
    Decide which structure to use.

    Priority:
      1. Selected UniProt AF2-EBI structure if its checked sequence matches the
         expected final node sequence.
      2. Rescue/custom AF2 structure if available and sequence-compatible.

    For selected UniProt AF2-EBI structures:
      - uncleaved proteins use structure_path
      - cleaved proteins use cleaved_structure_path
    """
    node = row.get("node")

    expected_final = expected_final_sequence(row)
    expected_full = expected_full_sequence(row)

    # Rule 1: selected UniProt AF2-EBI structure is acceptable.
    if row.get("structure_exists") == 1 and row.get("sequence_matches_structure") is True:
        if is_missing(row.get("cleavage_site_start")):
            return pd.Series(
                [
                    row.get("structure_path"),
                    row.get("structure_sequence"),
                    "EBI-selected-UniProt",
                ]
            )

        cleaved_path = row.get("cleaved_structure_path")
        cleaved_seq = clean_sequence(row.get("cleaved_structure_sequence"))

        if isinstance(cleaved_path, str) and os.path.isfile(cleaved_path):
            return pd.Series(
                [
                    cleaved_path,
                    cleaved_seq,
                    "EBI-selected-UniProt-cleaved",
                ]
            )

    # Rule 2: try rescue/custom AF2 structure.
    if af2_dir and isinstance(node, str):
        af2_path, af2_seq = locate_rescue_structure(node, af2_dir)

        af2_seq = clean_sequence(af2_seq)

        if af2_path is None or af2_seq is None:
            return pd.Series([None, None, "None"])

        if expected_final is not None and af2_seq == expected_final:
            source = (
                "AF2-rescue-trimmed"
                if not is_missing(row.get("cleavage_site_start"))
                else "AF2-rescue"
            )

            return pd.Series([af2_path, af2_seq, source])

        if (
            not is_missing(row.get("cleavage_site_start"))
            and expected_full is not None
            and af2_seq == expected_full
        ):
            cleaved_rescue_path = cleave_rescue_structure_if_needed(
                af2_path,
                node,
                row.get("cleavage_site_start"),
            )

            cleaved_rescue_seq = truncate_fasta(
                af2_seq,
                row.get("cleavage_site_start"),
            )

            if (
                isinstance(cleaved_rescue_path, str)
                and os.path.isfile(cleaved_rescue_path)
                and expected_final is not None
                and cleaved_rescue_seq == expected_final
            ):
                return pd.Series(
                    [
                        cleaved_rescue_path,
                        cleaved_rescue_seq,
                        "AF2-rescue-cleaved",
                    ]
                )

        print(
            f"[{node}] Rescue AF2 sequence did not match expected sequence. "
            f"AF2 length={len(af2_seq) if af2_seq is not None else 'NA'}; "
            f"expected_final length={len(expected_final) if expected_final is not None else 'NA'}; "
            f"expected_full length={len(expected_full) if expected_full is not None else 'NA'}"
        )

    return pd.Series([None, None, "None"])


def print_structure_diagnostics(nodes_df: pd.DataFrame) -> None:
    """
    Print concise diagnostics for structure preparation.
    """
    print()
    print("Structure preparation diagnostics:")
    print(f"Total nodes:                              {len(nodes_df):,}")
    print(f"Nodes with selected UniProtKB-AC:          {nodes_df['selected_UniProtKB_AC_for_structure'].notna().sum():,}")
    print(f"Nodes with existing AF2-EBI structures:    {nodes_df['structure_exists'].sum():,}")
    print(
        "Nodes with structure sequence strings:     "
        f"{nodes_df['structure_sequence'].apply(lambda x: isinstance(x, str)).sum():,}"
    )

    if "signalP_trimmed_sequence" in nodes_df.columns:
        print(
            "Nodes with SignalP-trimmed sequences:      "
            f"{nodes_df['signalP_trimmed_sequence'].apply(lambda x: isinstance(x, str)).sum():,}"
        )

    if "cleavage_site_start" in nodes_df.columns:
        print(
            "Nodes with SignalP cleavage sites:         "
            f"{nodes_df['cleavage_site_start'].notna().sum():,}"
        )

    print()
    print("Structure sequence-check status counts:")
    print(nodes_df["structure_sequence_check_status"].value_counts(dropna=False).to_string())

    print()
    print("Final structure source counts:")
    print(nodes_df["final_structure_source"].value_counts(dropna=False).to_string())


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare selected UniProt AF2 structures for downstream stability "
            "calculations with explicit sequence checking."
        )
    )

    parser.add_argument("--nodes", required=True)

    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where output files and cleaved PDBs will be saved.",
    )

    parser.add_argument("--output_prefix", required=True)

    parser.add_argument(
        "--input_dir",
        default="processed-data",
        help="Directory containing AF2-EBI structures and optional FASTA files.",
    )

    parser.add_argument("--organism_tag", required=True)

    parser.add_argument(
        "--AF2_dir",
        help="Directory containing additional/custom AF2 predictions.",
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    nodes_df = pd.read_pickle(args.nodes)

    required_columns = ["node", "UniProtKB-AC", "sequence"]

    for col in required_columns:
        if col not in nodes_df.columns:
            raise ValueError(
                f"Input nodes file is missing required column '{col}'. "
                "This script expects output from the updated sequence-parsing step."
            )

    if "signalP_trimmed_sequence" not in nodes_df.columns:
        print(
            "Warning: signalP_trimmed_sequence column is absent. "
            "The script will use sequence plus cleavage_site_start where possible."
        )

    if "cleavage_site_start" not in nodes_df.columns:
        print(
            "Warning: cleavage_site_start column is absent. "
            "All structures will be treated as uncleaved."
        )
        nodes_df["cleavage_site_start"] = np.nan

    # Locate selected UniProt AF2-EBI structures.
    nodes_df = locate_structure(nodes_df, args.input_dir)

    # Read structure sequence from AF2 FASTA when available, otherwise from PDB.
    nodes_df = locate_structure_sequence(nodes_df, args.input_dir)

    # Cleave selected UniProt AF2-EBI PDBs when SignalP predicts a cleavage site.
    nodes_df["cleaved_structure_path"] = nodes_df.apply(
        truncate_row_structure,
        axis=1,
    )

    nodes_df["cleaved_structure_sequence"] = nodes_df.apply(
        lambda row: truncate_fasta(
            row.get("structure_sequence"),
            row.get("cleavage_site_start"),
        ),
        axis=1,
    )

    # Compare selected structure sequence to the expected final sequence.
    nodes_df[
        [
            "sequence_matches_structure",
            "structure_sequence_check_status",
            "structure_sequence_length",
            "expected_structure_sequence_length",
            "structure_first_mismatch_index",
        ]
    ] = nodes_df.apply(compare_structure_to_expected, axis=1)

    # Select the final structure path/sequence/source.
    nodes_df[["final_structure_path", "final_sequence", "final_structure_source"]] = (
        nodes_df.apply(lambda row: select_final_structure(row, args.AF2_dir), axis=1)
    )

    # Compute mean pLDDT for the final selected structure.
    nodes_df["mean_plddt"] = nodes_df.apply(compute_mean_plddt, axis=1)

    print_structure_diagnostics(nodes_df)

    output_pkl = f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-temp.pkl"
    output_csv = f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-temp.csv"

    nodes_df.to_pickle(output_pkl)

    nodes_df = nodes_df.replace({None: np.nan})

    nodes_df.to_csv(
        output_csv,
        index=False,
    )

    print()
    print(f"Wrote structure-prepared pickle to: {output_pkl}")
    print(f"Wrote structure-prepared CSV to:    {output_csv}")


if __name__ == "__main__":
    main()
