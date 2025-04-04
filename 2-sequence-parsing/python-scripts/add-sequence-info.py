import pandas as pd
from Bio import SeqIO
import typing
import argparse


def add_sequences(nodes_df: pd.DataFrame, fasta_file: str) -> pd.DataFrame:
    """
    Adds protein sequence information to input DataFrame by matching gene names

    Args:
        nodes_df (pd.DataFrame): nodes information, loaded based on output from step 1
        fasta_file (str): path to file containing protein sequence information in fasta format

    Returns:
        Updated nodes_df containing sequence information where available and "None" where not available
    """

    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    nodes_df["has_verified_sequence"] = nodes_df["node"].isin(seqs.keys())

    nodes_df["sequence"] = nodes_df["node"].apply(
        lambda x: str(seqs[x].seq).rstrip("*") if x in seqs else None
    )

    return nodes_df


def main():

    parser = argparse.ArgumentParser(description="Add sequences to PPI network.")
    parser.add_argument(
        "--nodes",
        default="../1-network-centrality/processed-data/0_nodes-centrality.csv",
        help="Path to the nodes CSV file",
    )
    parser.add_argument(
        "--fasta",
        default="data-files/orf_trans.fasta",
        help="Path to open reading frame FASTA file",
    )
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Path to output directory",
    )
    parser.add_argument(
        "--output_prefix",
        default="0_",
        help="Prefix to be applied to output file",
    )
    parser.add_argument(
        "--organism_tag", default="s288c", help="Tag to label the organism for this run"
    )
    args = parser.parse_args()

    # load the nodes csv file
    nodes_df = pd.read_csv(args.nodes)

    # insert sequence information into nodes_df
    nodes_df = add_sequences(nodes_df, args.fasta)

    # replace missing values/nan with "None"
    nodes_df = nodes_df.replace("", "None")
    nodes_df = nodes_df.fillna("None")

    # write the updated nodes_df to file
    nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs.csv", index=False
    )


if __name__ == "__main__":

    main()
