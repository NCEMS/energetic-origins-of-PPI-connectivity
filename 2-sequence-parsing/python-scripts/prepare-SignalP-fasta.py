import pandas as pd
import typing
from typing import List
import argparse

def write_SignalP_fasta(nodes_df: pd.DataFrame, DeepTMHMM_classes_to_use: List[str], output_file: str) -> None:

    """
    Creates a fasta file that contains the sequences for which SignalP cleavage site predictions need to be run

    Args:
        nodes_df (pd.DataFrame): the input nodes_df; must contain sequence information
        DeepTMHMM_classes_to_use (List[str]): list of the DeepTMHMM classes on which SignalP predictions will be run
        output_file (str): path to the output file to be written by the function
    """

    

def main():

    parser = argparse.ArgumentParser(description="Create input fasta file for SignalP6.0; only proteins with classification SP will be included.")
    parser.add_argument(
        "--nodes",
        help="Path to the input nodes CSV file",
    )
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Path to output directory",
    )
    parser.add_argument(
        "--output_prefix",
        default="0",
        help="Prefix to be applied to output file",
    )
    parser.add_argument(
        "--organism_tag", default="s288c", help="Tag to label the organism for this run"
    )
    args = parser.parse_args()

    DeepTMHMM_classes_to_use = ["SP"]

    write_SignalP_fasta(nodes_df, DeepTMHMM_classes_to_use)

if __name__ == "__main__":
    main()
