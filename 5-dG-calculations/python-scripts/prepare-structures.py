import os, sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser, PDBIO, Select
import typing
from typing import Optional

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
    if not row["structure_exists"] or pd.isna(row["cleavage_start_site"]):
        return None

    input_path = row["structure_path"]
    cut_pos = int(row["cleavage_start_site"])

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
        lambda id: f"{structure_dir}/AF-{id}-F1-model_v4.pdb"
    )

    # create the structure_exists column by checking if the file actually exists
    nodes_df["structure_exists"] = nodes_df["structure_path"].apply(
        lambda path: 1 if os.path.exists(path) else None
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
    nodes_df["structure_fasta_path"] = nodes_df"UniProtKB-AC"].apply(
        lambda id: f"{structure_dir}/AF-{id}-F1-model_v4.fasta"
    )

    # add the sequence from this fasta file if it exists
    nodes_df["structure_sequence"] = nodes_df["structure_fasta_path"].apply(read_fasta_sequence)

    return nodes_df


def read_fasta_sequence(fasta_path: str) -> Optional(str):
    """
    Read in the sequence stored in a fasta file and return it

    Args:
        fasta_path (str): path to the fasta file

    Returns: either the sequence as a str or None
    """

    if not os.path.isfile(fasta_path):
        print(f"Missing file: {fasta_path}")
        return None

    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        return str(record.seq)

    except Exception as e:
        print(f"Error reading {fasta_path}: {e}")
        return None


def truncate_fasta(seq: str, cut: Optional[int]) -> Optional(str):
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

    return seq[int(cut):]


def compare_sequence(row) -> bool:

    struct_seq = row.get("cleaved_structure_sequence")
    signalp_seq = row.get("signalP_trimmed_sequence")

    # only compare if both sequences are present and not null
    if pd.isna(struct_seq) or pd.isna(signalp_seq) or struct_seq is None or signalp_seq is None:
        return False

    return str(struct_seq) == str(signalp_seq)


def main():

    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--nodes", required=True, help="Output from network-analysis.py"
    )
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where truncated PDBs will be saved",
    )
    parser.add_argument(
        "--input_dir",
        default="processed-data",
        help="Directory containing predicted structures for the proteome under consideration",
    )
    parser.add_argument(
        "--organism_tag")
    args = parser.parse_args()

    # read in the nodes_df from the previous step
    nodes_df = pd.read_pickle(args.nodes)

    # check to see if a structure exists for each node and add to column "structure_exists"
    nodes_df = locate_structure(nodes_df, args.input_dir)

    # add AF2 fasta information as well
    nodes_df = locate_structure_fasta(nodes_df, args.input_dir)

    # for proteins with a cleavage site predicted by SignalP, create a truncated structure and update structure_path value
    nodes_df["cleaved_structure_path"] = nodes_df.apply(truncate_row_structure, axis=1)

    # for proteins with a cleavage site predicted by SignalP, create a truncated FASTA based on the sequence from the AF2 structure
    nodes_df["cleaved_structure_sequence"] = nodes_df.apply(
        lambda row: truncate_sequence(row["structure_sequence"], row["cleavage_site_start"]),
        axis=1
    )

    # compare sequences between structures and trimmed sequences and add this information to a Boolean column named "sequence_matches_structure"
    nodes_df["sequence_matches_structure"] = nodes_df.apply(compare_sequences, axis=1)

    # save a temporary output file that has structure information; this is the input to the cagiada-stability.py calculations in the next rule
    nodes_df.to_pickle("{args.output_dir}/{args.output_prefix}-{args.organism_tag}-temp.pkl")

if __name__ == "__main__":

    main()
