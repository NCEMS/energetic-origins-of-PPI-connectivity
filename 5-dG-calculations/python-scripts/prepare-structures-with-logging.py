#!/usr/bin/env python3
"""Structure preparation with rescue-first cleavage selection and decision logging.

Existing CLI flags and -temp.pkl/-temp.csv names are preserved. Dependencies:
pandas, numpy, Biopython; pint/pint-pandas if needed by the input pickle.
Relative paths use the process working directory, as in the original pipeline.
Cleavage retains original residue numbers > cleavage_site_start. Existing
cleaved PDBs are reused; fallback acceptance checks their actual sequence.
Rescue structures are never modified. No class/sequence-verification eligibility
filter is added here; downstream steps retain responsibility for those filters.
"""
import os, sys
import logging
from pathlib import Path
from collections import Counter
import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
from typing import Optional
import argparse
# Register extension types when installed (some input pickles use pint).
try:
    import pint
    import pint_pandas
except ImportError:
    pass
import numpy as np
LOG = logging.getLogger("prepare-structures")

AA_MAP = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

class CleavageSelect(Select):
    def __init__(self, cut_pos: int):
        self.cut_pos = cut_pos
    def accept_residue(self, residue):
        return residue.id[1] > self.cut_pos

def truncate_structure(pdb_path: str, cut_pos: int, output_path: str) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, CleavageSelect(cut_pos))

def truncate_row_structure(row) -> Optional[str]:
    if row["structure_exists"] == 0 or pd.isna(row["cleavage_site_start"]):
        return None
    node = row["node"]
    input_path = row["structure_path"]
    output_path = str(Path(input_path).with_name(Path(input_path).stem + "-cleaved.pdb"))
    try:
        cut_value = float(row["cleavage_site_start"])
        if not np.isfinite(cut_value) or not cut_value.is_integer() or cut_value < 0:
            raise ValueError(f"Invalid cleavage position: {cut_value}")
        cut_pos = int(cut_value)
        if os.path.isfile(output_path):
            LOG.info("[%s] CLEAVED_EBI_REUSED path=%s; sequence will be checked if selected",
                     node, output_path)
        else:
            truncate_structure(input_path, cut_pos, output_path)
            LOG.info("[%s] CLEAVED_EBI_CREATED cut=%s path=%s", node, cut_pos, output_path)
        return output_path
    except Exception as exc:
        LOG.error("[%s] CLEAVED_EBI_PREPARATION_FAILED input=%s error=%s", node, input_path, exc)
        return None


def locate_structure(nodes_df: pd.DataFrame, structure_dir: str) -> pd.DataFrame:
    """
    Adds information to nodes_df regarding if a structure exists for a node and, if so, its path
    Args:
        nodes_df (pd.DataFrame): DataFrame to which structure information will be added
        structure_dir (str): path to the directory containing AlphaFold2 structures
    Returns:
        Updated nodes_df with structure_path and structure_exists columns inserted
    """
    # locate and add structures to dataframe
    nodes_df["structure_path"] = nodes_df["UniProtKB-AC"].apply(
        lambda id: f"{structure_dir}/AF-{id}-F1-model_v4.pdb" if pd.notna(id) else None
    )
    # create the structure_exists column by checking if the file actually exists
    nodes_df["structure_exists"] = nodes_df["structure_path"].apply(
        lambda path: 1 if path is not None and os.path.exists(path) else 0
    )
    return nodes_df

def locate_structure_fasta(nodes_df: pd.DataFrame, fasta_dir: str) -> pd.DataFrame:
    """
    Adds information to nodes_df regarding the fasta sequence extracted from the AF2 structure
    Args:
        nodes_df (pd.DataFrame): DataFrame to which structure information will be added
        fasta_dir (str): path to the directory containing AlphaFold2 structure-based sequencers
    Returns:
        Updated nodes_df with AF2_fasta_path column inserted
    """
    # locate the fasta file and add to dataframe
    nodes_df["structure_fasta_path"] = nodes_df["UniProtKB-AC"].apply(
        lambda id: f"{fasta_dir}/AF-{id}-F1-model_v4.fasta" if pd.notna(id) else None
    )
    # add the sequence from this fasta file if it exists
    nodes_df["structure_sequence"] = nodes_df["structure_fasta_path"].apply(
        read_fasta_sequence
    )
    return nodes_df

def extract_sequence_from_af2_pdb(pdb_path: str) -> Optional[str]:
    """
    Extract the full protein sequence from an AlphaFold2 PDB by walking residues
    in order, independent of geometric chain breaks.
    Returns:
        str or None
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("af2", pdb_path)
    # AF2 outputs a single model
    model = structure[0]
    seq_chars = []
    seen = set()  # avoid duplicates if any weird altlocs/duplicates exist
    for chain in model:
        for residue in chain:
            # residue.id is a tuple: (hetflag, resseq, icode)
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
    Read in the sequence stored in a fasta file and return it
    Args:
        fasta_path (str): path to the fasta file
    Returns: either the sequence as a str or None
    """
    if fasta_path is None:
        return None
    if not os.path.isfile(fasta_path):
        LOG.info("EBI_FASTA_MISSING path=%s", fasta_path)
        return None
    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        return str(record.seq)
    except Exception as e:
        LOG.warning("EBI_FASTA_READ_FAILED path=%s error=%s", fasta_path, e)
        return None

def truncate_fasta(seq: str, cut: Optional[int]) -> Optional[str]:
    """
    Creates a truncated fasta sequence based on the original AF2 sequence and the SignalP predicted cleavage site
    Args:
        seq (str): sequence to be trimmed
        cut (Optional[int]): the integer value of the first value to keep, or None
    Returns:
        The truncated sequence
    """
    if seq is None or pd.isna(seq):
        return None
    if pd.isna(cut) or cut is None:
        return seq
    return seq[int(cut) :]

def compare_sequence(row) -> bool:
    # Diagnostic comparison retained for compatibility. Final selection separately
    # checks the sequence extracted from the cleaved PDB or rescue PDB.
    struct_seq = (row.get("structure_sequence") if pd.isna(row.get("cleavage_site_start"))
                  else row.get("cleaved_structure_sequence"))
    expected = row.get("signalP_trimmed_sequence")
    if not isinstance(struct_seq, str) or not isinstance(expected, str):
        LOG.info("[%s] EBI_SEQUENCE_COMPARISON_UNAVAILABLE", row["node"])
        return False
    match = struct_seq.strip().upper() == expected.strip().upper()
    if not match:
        LOG.info("[%s] EBI_SEQUENCE_MISMATCH structure_length=%s expected_length=%s",
                 row["node"], len(struct_seq.strip()), len(expected.strip()))
    return match


def extract_plddt_ca_only(pdb_filename):
    """
    Extracts pLDDT values from CA atoms in B-factor field of an AlphaFold2 PDB file.
    Parameters:
    - pdb_filename (str): Path to the PDB file.
    Returns:
    - avg_plddt (float): Mean pLDDT value across all residues.
    """
    plddt_values = []
    with open(pdb_filename, "r") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    b_factor = float(line[60:66].strip())
                    plddt_values.append(b_factor)
                except ValueError:
                    continue  # skip lines with malformed B-factors
    plddt_array = np.array(plddt_values)
    avg_plddt = np.mean(plddt_array) if len(plddt_array) > 0 else float("nan")
    return avg_plddt

def compute_mean_plddt(row):
    path = row.get("final_structure_path")
    if isinstance(path, str) and os.path.isfile(path):
        return extract_plddt_ca_only(path)
    return np.nan

def select_final_structure(row, af2_dir):
    """Preserve matching rescue predictions; fall back only when rescue is absent.

    Noncleaved: matching EBI -> rescue.
    Cleaved: matching rescue -> cleaved EBI only if rescue file is absent.
    Existing but invalid/mismatched rescue files are explicitly flagged; no fallback.
    The AF2-cleaved label means the rescue matches the mature sequence. This
    function does NOT truncate the rescue prediction.
    """
    node = row["node"]
    cleavage = pd.notna(row.get("cleavage_site_start"))
    expected = row.get("signalP_trimmed_sequence")

    def selected(path, sequence, source, reason):
        LOG.info("[%s] SELECTED source=%s reason=%s length=%s path=%s",
                 node, source, reason, len(sequence), path)
        return pd.Series([path, sequence, source, reason])

    def rejected(reason, detail):
        LOG.warning("[%s] NOT_SELECTED reason=%s detail=%s", node, reason, detail)
        return pd.Series([None, None, "None", reason])

    LOG.info("[%s] SELECTION_START cleavage=%s DeepTMHMM_class=%s",
             node, cleavage, row.get("DeepTMHMM_class"))
    if not isinstance(expected, str) or not expected.strip():
        return rejected("invalid_expected_sequence", "signalP_trimmed_sequence is missing/empty")
    expected_clean = expected.strip().upper()

    # Preserve the original priority for noncleaved, matching EBI structures.
    if (not cleavage and row.get("structure_exists") == 1
            and row.get("sequence_matches_structure", False)):
        path, seq = row.get("structure_path"), row.get("structure_sequence")
        if isinstance(path, str) and os.path.isfile(path):
            return selected(path, seq, "EBI", "matching_uncleaved_ebi")

    rescue_path = None
    if af2_dir and isinstance(node, str):
        if Path(node).name != node or node in {".", ".."}:
            return rejected("invalid_node_path", "Node cannot be used as a rescue-directory name")
        rescue_path = os.path.join(af2_dir, node, "ranked_0.pdb")

    # lexists also detects a broken symlink, which must be reviewed, not ignored.
    if rescue_path is not None and os.path.lexists(rescue_path):
        LOG.info("[%s] RESCUE_FOUND path=%s", node, rescue_path)
        if not os.path.isfile(rescue_path):
            return rejected("rescue_unreadable", f"Not a readable regular file: {rescue_path}; no fallback")
        try:
            rescue_seq = extract_sequence_from_af2_pdb(rescue_path)
        except Exception as exc:
            return rejected("rescue_unreadable", f"{rescue_path}: {exc}; no fallback")
        if not isinstance(rescue_seq, str) or not rescue_seq.strip():
            return rejected("rescue_empty_sequence", f"{rescue_path}; no fallback")
        if rescue_seq.strip().upper() != expected_clean:
            return rejected("rescue_sequence_mismatch",
                            f"expected_length={len(expected_clean)} actual_length={len(rescue_seq)} "
                            f"path={rescue_path}; no fallback")
        return selected(rescue_path, rescue_seq, "AF2-cleaved" if cleavage else "AF2",
                        "matching_rescue_preferred")

    LOG.info("[%s] RESCUE_ABSENT path=%s", node, rescue_path or "AF2_dir not configured")
    if not cleavage:
        return rejected("no_matching_structure", "No matching uncleaved EBI structure and no rescue")

    # New fallback: validate the actual cleaved PDB, not only its companion FASTA.
    path = row.get("cleaved_structure_path")
    LOG.info("[%s] TRY_CLEAVED_EBI_FALLBACK path=%s", node, path)
    if not isinstance(path, str) or not os.path.isfile(path):
        return rejected("cleaved_ebi_missing", "No rescue and no usable cleaved EBI file")
    try:
        sequence = extract_sequence_from_af2_pdb(path)
    except Exception as exc:
        return rejected("cleaved_ebi_unreadable", f"{path}: {exc}")
    if not isinstance(sequence, str) or not sequence.strip():
        return rejected("cleaved_ebi_empty_sequence", path)
    if sequence.strip().upper() != expected_clean:
        return rejected("cleaved_ebi_sequence_mismatch",
                        f"expected_length={len(expected_clean)} actual_length={len(sequence)} path={path}")
    return selected(path, sequence, "EBI-cleaved", "rescue_absent_used_cleaved_ebi")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare structures, preferring matching rescue predictions over cleaved EBI fallback.")
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--input_dir", default="processed-data")
    parser.add_argument("--organism_tag", required=True)
    parser.add_argument("--AF2_dir")
    parser.add_argument("--compare_nodes", help="Optional previous annotated pickle for selection comparison")
    parser.add_argument("--log_file", help="Optional additional log file (existing shell redirection also works)")
    args = parser.parse_args()
    handlers = [logging.StreamHandler(sys.stdout)]
    if args.log_file:
        Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.log_file, mode="w"))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", handlers=handlers)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{args.output_prefix}-{args.organism_tag}"
    LOG.info("START nodes=%s working_directory=%s", args.nodes, os.getcwd())
    LOG.info("POLICY noncleaved=matching EBI then rescue; cleaved=matching rescue then EBI-cleaved if rescue absent")
    LOG.info("POLICY existing invalid rescue is flagged with no fallback; existing cleaved PDBs are reused")
    LOG.info("PATHS input_dir=%s AF2_dir=%s output_dir=%s", args.input_dir, args.AF2_dir, output_dir)
    nodes_df = pd.read_pickle(args.nodes)
    required = {"node", "UniProtKB-AC", "cleavage_site_start", "signalP_trimmed_sequence"}
    if not nodes_df.columns.is_unique or not required.issubset(nodes_df.columns):
        raise ValueError(f"Input needs unique columns including {sorted(required)}")
    if nodes_df["node"].isna().any() or nodes_df["node"].duplicated().any():
        raise ValueError("Input node identifiers must be non-null and unique")
    baseline = pd.read_pickle(args.compare_nodes) if args.compare_nodes else None
    if baseline is not None:
        fields = {"node", "final_structure_path", "final_sequence", "final_structure_source"}
        if not baseline.columns.is_unique or not fields.issubset(baseline.columns):
            raise ValueError(f"Comparison file needs unique columns including {sorted(fields)}")
        if baseline["node"].isna().any() or baseline["node"].duplicated().any():
            raise ValueError("Comparison file node identifiers must be non-null and unique")
    original_index = nodes_df.index.copy()
    nodes_df = locate_structure(nodes_df, args.input_dir)
    nodes_df = locate_structure_fasta(nodes_df, args.input_dir)
    LOG.info("INPUT rows=%s original_EBI_files=%s cleavage_sites=%s", len(nodes_df),
             int(nodes_df["structure_exists"].sum()), int(nodes_df["cleavage_site_start"].notna().sum()))
    nodes_df["cleaved_structure_path"] = nodes_df.apply(truncate_row_structure, axis=1)
    nodes_df["cleaved_structure_sequence"] = nodes_df.apply(
        lambda row: truncate_fasta(row["structure_sequence"], row["cleavage_site_start"]), axis=1)
    nodes_df["sequence_matches_structure"] = nodes_df.apply(compare_sequence, axis=1)
    columns = ["final_structure_path", "final_sequence", "final_structure_source", "structure_selection_reason"]
    decisions = [select_final_structure(row, args.AF2_dir).tolist() for _, row in nodes_df.iterrows()]
    for i, column in enumerate(columns):
        nodes_df[column] = [decision[i] for decision in decisions]
    nodes_df["mean_plddt"] = nodes_df.apply(compute_mean_plddt, axis=1)
    if not nodes_df.index.equals(original_index):
        raise AssertionError("Input row index changed")
    sources = Counter(nodes_df["final_structure_source"])
    reasons = Counter(nodes_df["structure_selection_reason"])
    for source, count in sorted(sources.items()):
        LOG.info("SUMMARY source=%s count=%s", source, count)
    for reason, count in sorted(reasons.items()):
        LOG.info("SUMMARY reason=%s count=%s", reason, count)
    report_cols = [c for c in ["node", "DeepTMHMM_class", "cleavage_site_start", "structure_exists",
                  "structure_path", "cleaved_structure_path", "sequence_matches_structure",
                  *columns, "mean_plddt"] if c in nodes_df.columns]
    report = nodes_df[report_cols].copy()
    if baseline is not None:
        old = baseline.set_index("node")
        def value(v):
            return None if pd.isna(v) else v
        changes = []
        for _, row in nodes_df.iterrows():
            node = row["node"]
            if node not in old.index:
                status = "not_in_baseline"
            else:
                previous = old.loc[node]
                before = tuple(value(previous[c]) for c in columns[:3])
                after = tuple(value(row[c]) for c in columns[:3])
                if before == after:
                    status = "unchanged"
                elif before[0] is None and after[0] is not None:
                    status = "recovered"
                elif before[0] is not None and after[0] is None:
                    status = "lost_selection"
                else:
                    status = "changed_selection"
                if status != "unchanged":
                    LOG.warning("[%s] COMPARISON status=%s old_source=%s old_path=%s new_source=%s new_path=%s",
                                node, status, before[2], before[0], after[2], after[0])
            changes.append(status)
        report["comparison_status"] = changes
        for c in columns[:3]:
            report["previous_" + c] = report["node"].map(old[c])
        for status, count in sorted(Counter(changes).items()):
            LOG.info("COMPARISON_SUMMARY status=%s count=%s", status, count)
        LOG.info("COMPARISON_SUMMARY baseline_nodes_absent_from_input=%s",
                 len(set(old.index) - set(nodes_df["node"])))
    for suffix, data in [("structure-selection.csv", report),
                         ("structure-unselected.csv", report.loc[nodes_df["final_structure_path"].isna()])]:
        path = output_dir / f"{prefix}-{suffix}"
        data.to_csv(path, index=False)
        LOG.info("WROTE %s rows=%s", path, len(data))
    pkl_path = output_dir / f"{prefix}-temp.pkl"
    csv_path = output_dir / f"{prefix}-temp.csv"
    nodes_df.to_pickle(pkl_path)
    nodes_df.to_csv(csv_path, index=False)
    LOG.info("WROTE %s", pkl_path)
    LOG.info("WROTE %s", csv_path)
    LOG.info("COMPLETE selected=%s unselected=%s",
             int(nodes_df["final_structure_path"].notna().sum()),
             int(nodes_df["final_structure_path"].isna().sum()))

if __name__ == "__main__":
    try:
        main()
    except Exception:
        logging.basicConfig(level=logging.INFO)
        LOG.exception("Structure preparation failed")
        sys.exit(1)
