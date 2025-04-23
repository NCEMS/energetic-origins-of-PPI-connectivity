import os, sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import typing
from typing import Optional

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


def truncate_fasta(nodes_df:pd.DataFrame, fasta_dir: str) -> pd.DataFrame:
    """
    Creates a truncated fasta sequence based on the original AF2 sequence and the SignalP
    Args:
        nodes_df (pd.DataFrame): DataFrame to which structure information will be added
        fasta_dir (str): path to the directory containing AlphaFold2 structure-based sequencers

    Returns:
        Updated nodes_df with AF2_fasta_path column inserted
    """


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
    

    # for proteins with a cleavage site predicted by SignalP, create a truncated FASTA based on the sequence
    # present in the truncated AF2 structure

    # compare sequences between structures and trimmed sequences and add this information to a Boolean column named "sequence_matches_structure"

    # save a temporary output file that has structure information; this is the input to the cagiada-stability.py calculations in the next rule
    nodes_df.to_pickle("{args.output_dir}/{args.output_prefix}-{args.organism_tag}-temp.pkl")

if __name__ == "__main__":

    main()
