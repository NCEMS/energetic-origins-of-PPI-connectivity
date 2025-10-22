import pandas as pd
import typing
from typing import List
import argparse


def write_SignalP_fasta(
    nodes_df: pd.DataFrame, DeepTMHMM_classes_to_use: List[str], output_file: str
) -> None:
    """
    Creates a fasta file that contains the sequences for which SignalP cleavage site predictions need to be run

    Args:
        nodes_df (pd.DataFrame): the input nodes_df; must contain sequence information
        DeepTMHMM_classes_to_use (List[str]): list of the DeepTMHMM classes on which SignalP predictions will be run
        output_file (str): path to the output file to be written by the function

    Returns:
        Nothing, but writes a file to the path output_file
    """

    with open(output_file, "w") as f:
        for _, row in nodes_df.iterrows():
            if row["DeepTMHMM_class"] in DeepTMHMM_classes_to_use:
                f.write(f">{row['node']}\n{row['sequence']}\n")

    return


def main():

    parser = argparse.ArgumentParser(
        description="Create input fasta file for SignalP6.0"
    )
    parser.add_argument(
        "--nodes",
        help="Path to the input nodes CSV file",
    )
    parser.add_argument(
        "--output_dir",
        help="Path to output directory",
    )
    parser.add_argument(
        "--output_prefix",
        help="Prefix to be applied to output file",
    )
    parser.add_argument("--output_fasta", help="Path to the fasta file to be written")
    parser.add_argument("--organism_tag", help="Tag to label the organism for this run")
    args = parser.parse_args()

    # create a list of the DeepTMHMM output tags indicating proteins on which SignalP will be run
    DeepTMHMM_classes_to_use = ["SP", "TM", "GLOB", "BETA", "SP+TM"]

    # read in the nodes information
    nodes_df = pd.read_csv(args.nodes)

    # create the fasta file needed as input for SignalP
    write_SignalP_fasta(nodes_df, DeepTMHMM_classes_to_use, args.output_fasta)


if __name__ == "__main__":
    main()
