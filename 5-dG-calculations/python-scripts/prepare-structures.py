import os, sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser, PDBIO, Select
import typing
from typing import Optional
import argparse
import pint
import pint_pandas
import numpy as np


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

    input_path = row["structure_path"]
    cut_pos = int(row["cleavage_site_start"])

    output_path = input_path.replace(".pdb", "-cleaved.pdb")

    try:
        truncate_structure(input_path, cut_pos, output_path)
        return output_path
    except Exception as e:
        print(f"Error processing {input_path}: {e}")
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
        print(f"Missing file: {fasta_path}")
        return None

    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        return str(record.seq)

    except Exception as e:
        print(f"Error reading {fasta_path}: {e}")
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

    if pd.isna(cut) or cut is None:
        return seq

    return seq[int(cut) :]


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
    Extracts pLDDT values from CA atoms in B-factor field of an AlphaFold2 PDB file.

    Parameters:
    - pdb_filename (str): Path to the PDB file.

    Returns:
    - avg_plddt (float): Mean pLDDT value across all residues.
    """
    plddt_values = []

    with open(pdb_filename, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    b_factor = float(line[60:66].strip())
                    plddt_values.append(b_factor)
                except ValueError:
                    continue  # skip lines with malformed B-factors

    plddt_array = np.array(plddt_values)
    avg_plddt = np.mean(plddt_array) if len(plddt_array) > 0 else float('nan')
    return avg_plddt


def compute_mean_plddt(row):
    if row["structure_exists"] == 1:
        if pd.isna(row["cleavage_site_start"]):
            return extract_plddt_ca_only(row["structure_path"])
        else:
            return extract_plddt_ca_only(row["cleaved_structure_path"])
    else:
        return np.nan

"""
def compute_mean_plddt(row):

    if row["structure_exists"] == 1:
        return extract_plddt_ca_only(row["structure_path"])
    else:
        return np.nan
"""
"""
def prepare_alphafold_fasta_dirs(nodes_df: pd.DataFrame, output_dir: str, organism_tag: str, output_prefix: str) -> None:
    #
    #Creates AlphaFold2 prediction directories and writes a list of target nodes to a file.

    #Outputs a text file: {output_prefix}-{organism_tag}-alphafold-targets.txt
    #Each line: one node ID
    #

    af_base_dir = os.path.join(output_dir, "alphafold")
    os.makedirs(af_base_dir, exist_ok=True)

    prepared_nodes = set()

    for _, row in nodes_df.iterrows():
        node_id = str(row["node"])
        seq = row.get("signalP_trimmed_sequence")
        cut = row.get("cleavage_site_start")
        struct_ok = row.get("structure_exists") == 1
        matches = row.get("sequence_matches_structure")

        subdir = os.path.join(af_base_dir, node_id)
        output_pdb = os.path.join(subdir, "result_model_1_pred_0.pdb")

        if os.path.exists(output_pdb):
            continue  # Skip if prediction already exists

        # CASE 1: Cleaved prediction needed
        if struct_ok and pd.notna(cut) and isinstance(seq, str):
            prepared_nodes.add(node_id)
            os.makedirs(subdir, exist_ok=True)
            fasta_path = os.path.join(subdir, f"{node_id}.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">{node_id}\n{seq.strip().upper()}\n")

        # CASE 2: Full-length prediction due to mismatch
        elif matches is False and isinstance(seq, str) and node_id not in prepared_nodes:
            prepared_nodes.add(node_id)
            os.makedirs(subdir, exist_ok=True)
            fasta_path = os.path.join(subdir, f"{node_id}.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">{node_id}\n{seq.strip().upper()}\n")

    # Write target node IDs to file
    targets_file = os.path.join(output_dir, f"{output_prefix}-{organism_tag}-alphafold-targets.txt")
    with open(targets_file, "w") as f:
        for node_id in sorted(prepared_nodes):
            f.write(f"{node_id}\n")
"""

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
    args = parser.parse_args()

    # read in the nodes_df from the previous step
    nodes_df = pd.read_pickle(args.nodes)

    # check to see if a structure exists for each node and add to column "structure_exists" with path in "structure_path"
    nodes_df = locate_structure(nodes_df, args.input_dir)

    # add AF2 fasta information as well; adds "structure_fasta_path" and "structure_sequence" to the pd.DataFrame
    nodes_df = locate_structure_fasta(nodes_df, args.input_dir)

    # for proteins with a cleavage site predicted by SignalP, create a truncated structure and update structure_path value
    nodes_df["cleaved_structure_path"] = nodes_df.apply(truncate_row_structure, axis=1)

    # for proteins with a cleavage site predicted by SignalP, create a truncated FASTA based on the sequence from the AF2 structure
    nodes_df["cleaved_structure_sequence"] = nodes_df.apply(
        lambda row: truncate_fasta(
            row["structure_sequence"], row["cleavage_site_start"]
        ),
        axis=1,
    )

    # compare sequences between structures and trimmed sequences and add this information to a Boolean column named "sequence_matches_structure"
    nodes_df["sequence_matches_structure"] = nodes_df.apply(compare_sequence, axis=1)

    # compute the mean of per-residue pLDDT for each AF2 structure
    nodes_df["mean_plddt"] = nodes_df.apply(compute_mean_plddt, axis=1)

    # call function to prepare AF2 input fasta and directory (not currently used; AF2 run in separate pipeline)
    #prepare_alphafold_fasta_dirs(nodes_df, args.output_dir, args.organism_tag, args.output_prefix)

    # save a temporary output file that has structure information; this is the input to the cagiada-stability.py calculations in the next rule
    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-temp.pkl"
    )


if __name__ == "__main__":

    main()
