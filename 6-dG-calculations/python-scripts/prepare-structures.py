import os, sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
import typing
from typing import Optional
import argparse
import pint
import pint_pandas
import numpy as np

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


def truncate_fasta(seq: str, cut: Optional[int]) -> Optional[str]:
    """
    Creates a truncated FASTA sequence based on the original AF2 sequence and
    the SignalP predicted cleavage site.

    Args:
        seq (str): sequence to be trimmed
        cut (Optional[int]): number of residues to remove from the N-terminus,
            or None/NaN if no trimming should be performed

    Returns:
        The truncated sequence, the original sequence, or None if seq is missing
    """

    if not isinstance(seq, str):
        return None

    seq = seq.strip().rstrip("*")

    if pd.isna(cut) or cut is None:
        return seq

    try:
        cut = int(cut)
    except (TypeError, ValueError):
        return None

    return seq[cut:]


def truncate_structure(pdb_path: str, cut_pos: int, output_path: str) -> None:
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

    if pd.isna(row.get("cleavage_site_start")):
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
        print(f"Error processing {input_path}: {e}")
        return None


def parse_uniprot_accessions(value) -> list[str]:
    """
    Parse a UniProtKB-AC field that may contain one or more semicolon-delimited
    accessions. Missing values return an empty list.

    Examples:
        Q9SP32 -> ["Q9SP32"]
        F4HQG6;Q56Y50;Q9SP32 -> ["F4HQG6", "Q56Y50", "Q9SP32"]
        np.nan -> []
        "nan" -> []
    """

    if value is None or pd.isna(value):
        return []

    accessions = []

    for accession in str(value).split(";"):
        accession = accession.strip()

        if not accession:
            continue

        if accession.lower() in {"nan", "none", "na", "null"}:
            continue

        accessions.append(accession)

    return accessions


def choose_first_uniprot_accession(value) -> Optional[str]:
    """
    Choose the first UniProt accession from a potentially semicolon-delimited
    field. Return None if no valid accession is available.
    """

    accessions = parse_uniprot_accessions(value)
    return accessions[0] if accessions else None


def locate_structure(nodes_df: pd.DataFrame, structure_dir: str) -> pd.DataFrame:
    """
    Adds AF2-EBI structure path information to nodes_df.

    Handles UniProtKB-AC values that are:
        - missing
        - single accessions, e.g. Q9SP32
        - semicolon-delimited accessions, e.g. F4HQG6;Q56Y50;Q9SP32

    If multiple accessions have existing structures, this is logged. The first
    accession in the UniProtKB-AC field is chosen deterministically.
    """

    nodes_df = nodes_df.copy()

    def existing_structure_accessions(uniprot_field):
        existing = []

        for accession in parse_uniprot_accessions(uniprot_field):
            path = os.path.join(
                structure_dir,
                f"AF-{accession}-F1-model_v6.pdb",
            )

            if os.path.isfile(path):
                existing.append(accession)

        return existing

    nodes_df["n_UniProtKB_AC"] = nodes_df["UniProtKB-AC"].apply(
        lambda x: len(parse_uniprot_accessions(x))
    )

    nodes_df["existing_AF2_EBI_structure_accessions"] = nodes_df["UniProtKB-AC"].apply(
        existing_structure_accessions
    )

    nodes_df["n_existing_AF2_EBI_structures"] = nodes_df[
        "existing_AF2_EBI_structure_accessions"
    ].apply(len)

    n_missing_uniprot = nodes_df["n_UniProtKB_AC"].eq(0).sum()
    n_multiple_uniprot = nodes_df["n_UniProtKB_AC"].gt(1).sum()
    n_multiple_structures = nodes_df["n_existing_AF2_EBI_structures"].gt(1).sum()

    print(f"Nodes with no UniProtKB-AC mapping: {n_missing_uniprot:,}")
    print(f"Nodes with multiple UniProtKB-AC mappings: {n_multiple_uniprot:,}")
    print(f"Nodes with multiple existing AF2-EBI structures: {n_multiple_structures:,}")

    if n_multiple_structures > 0:
        print("Examples of nodes with multiple existing AF2-EBI structures:")

        example_rows = nodes_df.loc[
            nodes_df["n_existing_AF2_EBI_structures"].gt(1),
            ["node", "UniProtKB-AC", "existing_AF2_EBI_structure_accessions"],
        ].head(10)

        for _, row in example_rows.iterrows():
            print(
                f"{row['node']}: "
                f"UniProtKB-AC={row['UniProtKB-AC']}; "
                f"existing={','.join(row['existing_AF2_EBI_structure_accessions'])}"
            )

    nodes_df["selected_UniProtKB_AC_for_structure"] = nodes_df["UniProtKB-AC"].apply(
        choose_first_uniprot_accession
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

    return nodes_df


def locate_structure_fasta(nodes_df: pd.DataFrame, fasta_dir: str) -> pd.DataFrame:
    """
    Adds AF2-derived FASTA path and sequence information.

    Uses selected_UniProtKB_AC_for_structure if present. Missing UniProt IDs
    produce None rather than AF-nan paths.
    """

    nodes_df = nodes_df.copy()

    if "selected_UniProtKB_AC_for_structure" not in nodes_df.columns:
        nodes_df["selected_UniProtKB_AC_for_structure"] = nodes_df["UniProtKB-AC"].apply(
            choose_first_uniprot_accession
        )

    nodes_df["structure_fasta_path"] = nodes_df[
        "selected_UniProtKB_AC_for_structure"
    ].apply(
        lambda accession: (
            os.path.join(fasta_dir, f"AF-{accession}-F1-model_v6.fasta")
            if accession is not None
            else None
        )
    )

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
    Read the sequence stored in a FASTA file and return it.
    """

    if fasta_path is None or pd.isna(fasta_path):
        return None

    if not isinstance(fasta_path, str):
        return None

    if not os.path.isfile(fasta_path):
        print(f"Missing file: {fasta_path}")
        return None

    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        return str(record.seq).rstrip("*")

    except Exception as e:
        print(f"Error reading {fasta_path}: {e}")
        return None


def compare_sequence(row) -> bool:

    if pd.isna(row["cleavage_site_start"]):
        struct_seq = row.get("structure_sequence")
    else:
        struct_seq = row.get("cleaved_structure_sequence")

    signalp_seq = row.get("signalP_trimmed_sequence")

    # basic checks
    if not isinstance(struct_seq, str) or not isinstance(signalp_seq, str):
        print("Problem with", row["node"])
        print(f"Either {struct_seq} or {signalp_seq} or both is not a str object")
        print(f"type of struct_seq", type(struct_seq))
        print(f"type of signalp_seq", type(signalp_seq), "\n")
        return False

    # strip and uppercase
    struct_seq_clean = struct_seq.strip().upper()
    signalp_seq_clean = signalp_seq.strip().upper()

    if struct_seq_clean != signalp_seq_clean:
        print(
            f"Mismatch for",
            row["node"],
            f":\n  STRUCT:   {struct_seq_clean}\n  SIGNALP: {signalp_seq_clean}\n",
        )

    return struct_seq_clean == signalp_seq_clean


def extract_plddt_ca_only(pdb_filename):
    """
    Extract pLDDT values from CA atoms in the B-factor field of an AlphaFold2 PDB file.
    """

    if pdb_filename is None or pd.isna(pdb_filename):
        return np.nan

    if not isinstance(pdb_filename, (str, bytes, os.PathLike)):
        return np.nan

    if not os.path.isfile(pdb_filename):
        return np.nan

    plddt_values = []

    with open(pdb_filename, "r") as f:
        for line in f:
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
    path = row.get("final_structure_path")

    if path is None or pd.isna(path):
        return np.nan

    if not isinstance(path, (str, bytes, os.PathLike)):
        return np.nan

    if os.path.isfile(path):
        return extract_plddt_ca_only(path)

    return np.nan


def select_final_structure(row, af2_dir):
    """
    Decide which structure to use: EBI, AF2 replacement, or cleaved AF2.
    For cleaved sequences, confirms AF2 structure matches the expected trimmed sequence.
    """
    node = row["node"]
    sgd_seq = row.get("signalP_trimmed_sequence")

    # Rule 1: EBI structure is acceptable
    if (
        row.get("structure_exists") == 1
        and pd.isna(row.get("cleavage_site_start"))
        and row.get("sequence_matches_structure", False)
    ):
        return pd.Series(
            [row.get("structure_path"), row.get("structure_sequence"), "EBI"]
        )

    # Rules 2 & 3: Try AF2 structure
    if af2_dir and isinstance(node, str):
        af2_path, af2_seq = locate_rescue_structure(node, af2_dir)

        if af2_path is None or af2_seq is None:
            return pd.Series([None, None, "None"])

        # make sure both sequences are valid strings
        if not (isinstance(af2_seq, str) and isinstance(sgd_seq, str)):
            print(f"[{node}] One or both sequences are not valid strings.")
            return pd.Series([None, None, "None"])

        af2_clean = af2_seq.strip().upper()
        sgd_clean = sgd_seq.strip().upper()
        print (node)
        print (f"af2_clean: {af2_clean}")
        print (f"sgd_clean: {sgd_clean}")

        if af2_clean == sgd_clean:
            source = (
                "AF2-cleaved" if not pd.isna(row.get("cleavage_site_start")) else "AF2"
            )
            return pd.Series([af2_path, af2_seq, source])
        else:
            if not pd.isna(row.get("cleavage_site_start")):
                print(
                    f"[{node}] Cleaved AF2 structure does not match expected trimmed sequence."
                )
            else:
                print(f"[{node}] AF2 structure does not match expected full sequence.")
            return pd.Series([None, None, "None"])

    return pd.Series([None, None, "None"])


def locate_rescue_structure(
    node: str, af2_dir: str
) -> typing.Tuple[Optional[str], Optional[str]]:
    """
    Locate the AlphaFold2 rescue structure (ranked_0.pdb) and extract sequence from it.

    Args:
        node (str): The ORF name (e.g. YGR188C) from nodes_df["node"]
        af2_dir (str): Root directory for AF2 predictions (contains subfolders for each ORF)

    Returns:
        Tuple[str or None, str or None]: path to ranked_0.pdb and extracted sequence, or (None, None)

    """
    pdb_path = os.path.join(af2_dir, node, "ranked_0.pdb")

    if not os.path.isfile(pdb_path):
        return None, None

    try:
        af2_seq = extract_sequence_from_af2_pdb(pdb_path)
        if af2_seq is None:
            print(f"No sequence could be extracted from AF2 structure for {node}")
            return None, None

        return pdb_path, af2_seq

    except Exception as e:
        print(f"Error reading AF2 structure for {node}: {e}")
        return None, None


def main():

    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", required=True)
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where truncated PDBs will be saved",
    )
    parser.add_argument("--output_prefix")
    parser.add_argument(
        "--input_dir",
        default="processed-data",
        help="Directory containing predicted structures for the proteome under consideration",
    )
    parser.add_argument("--organism_tag")
    parser.add_argument(
        "--AF2_dir", help="Directory containing additional AF2 predictions"
    )
    args = parser.parse_args()

    # read in the nodes_df from the previous step
    nodes_df = pd.read_pickle(args.nodes)

    # check to see if a structure exists for each node and add to column "structure_exists" with path in "structure_path"
    nodes_df = locate_structure(nodes_df, args.input_dir)

    # add AF2 fasta information as well; adds "structure_fasta_path" and "structure_sequence" to the pd.DataFrame
    nodes_df = locate_structure_fasta(nodes_df, args.input_dir)

    print("Structure preparation diagnostics:")
    print(f"Total nodes: {len(nodes_df):,}")
    print(f"Nodes with UniProtKB-AC: {nodes_df['UniProtKB-AC'].notna().sum():,}")
    print(
        "Nodes with selected UniProt accession: "
        f"{nodes_df['selected_UniProtKB_AC_for_structure'].notna().sum():,}"
    )
    print(f"Nodes with existing structures: {nodes_df['structure_exists'].sum():,}")
    print(
        "Nodes with structure_sequence strings: "
        f"{nodes_df['structure_sequence'].apply(lambda x: isinstance(x, str)).sum():,}"
    )
    print(
        "Nodes with SignalP cleavage sites: "
        f"{nodes_df['cleavage_site_start'].notna().sum():,}"
    )

    # for proteins with a cleavage site predicted by SignalP, create a truncated structure and update structure_path value
    nodes_df["cleaved_structure_path"] = nodes_df.apply(truncate_row_structure, axis=1)

    # for proteins with a cleavage site predicted by SignalP, create a truncated FASTA based on the sequence from the AF2 structure
    nodes_df["cleaved_structure_sequence"] = nodes_df.apply(
        lambda row: truncate_fasta(
            row.get("structure_sequence"),
            row.get("cleavage_site_start"),
        ),
        axis=1,
    )
    # compare sequences between structures and trimmed sequences and add this information to a Boolean column named "sequence_matches_structure"
    nodes_df["sequence_matches_structure"] = nodes_df.apply(compare_sequence, axis=1)

    # determine the final structure to use - includes accounting for new AF2 structures
    nodes_df[["final_structure_path", "final_sequence", "final_structure_source"]] = (
        nodes_df.apply(lambda row: select_final_structure(row, args.AF2_dir), axis=1)
    )

    # compute the mean of per-residue pLDDT for each AF2 structure
    nodes_df["mean_plddt"] = nodes_df.apply(compute_mean_plddt, axis=1)

    # save a temporary output file that has structure information; this is the input to the cagiada-stability.py calculations in the next rule
    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-temp.pkl"
    )

    # repace None for np.nan for simplicity
    nodes_df = nodes_df.replace({None: np.nan})

    # save a csv file as well (for testing purposes)
    nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-temp.csv",
        index=False,
    )


if __name__ == "__main__":

    main()
