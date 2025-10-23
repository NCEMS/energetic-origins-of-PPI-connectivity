import pandas as pd
from Bio import SeqIO
import typing
import argparse
import numpy as np


def add_sequences(nodes_df: pd.DataFrame, fasta_file: str) -> pd.DataFrame:
    """
    Adds protein sequence information to input DataFrame by matching gene names

    Args:
        nodes_df (pd.DataFrame): nodes information, loaded based on output from step 1
        fasta_file (str): path to file containing protein sequence information in fasta format

    Returns:
        Updated nodes_df containing sequence information where available and np.nan where not available
    """

    # read sequences
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # add column stating which nodes have sequence info
    nodes_df["has_verified_sequence"] = nodes_df["node"].isin(seqs.keys())

    # add the sequence information
    nodes_df["sequence"] = nodes_df["node"].apply(
        lambda x: str(seqs[x].seq).rstrip("*") if x in seqs else np.nan
    )

    return nodes_df


def add_mappings(nodes_df: pd.DataFrame, map_file: str) -> pd.DataFrame:
    """
    Add UniProtKB-AC identifiers to each node as possible

    Args:
        nodes_df (pd.DataFrame): current nodes_df to be updated
        map_file (str): path to the *_idmapping.dat file to be used (e.g., YEAST_559292_idmapping.dat)

    Returns:
        Updated nodes_df (pd.DataFrame) containing mapped names in column 'UniProtKB-AC'
    """

    # read in the mapping file as a pd.DataFrame
    column_names = ["UniProtKB-AC", "ID_type", "ID"]
    cross_df = pd.read_csv(map_file, names=column_names, sep="\t")
    cross_df = cross_df[cross_df["ID_type"] == "Gene_OrderedLocusName"]

    # add information to nodes_df with a left merge
    nodes_df = nodes_df.merge(
        cross_df[["UniProtKB-AC", "ID"]], left_on="node", right_on="ID", how="left"
    )

    # drop the unneeded "ID" column left over from the merge
    nodes_df.drop(columns=["ID"], inplace=True)

    return nodes_df


def main():

    parser = argparse.ArgumentParser(description="Add sequences to PPI network.")
    parser.add_argument(
        "--nodes",
        help="Path to the nodes CSV file",
    )
    parser.add_argument(
        "--fasta",
        help="Path to open reading frame FASTA file",
    )
    parser.add_argument(
        "--output_dir",
        help="Path to output directory",
    )
    parser.add_argument(
        "--output_prefix",
        help="Prefix to be applied to output file",
    )
    parser.add_argument(
        "--output_suffix",
        help="Suffix to be applied to output file",
    )
    parser.add_argument(
        "--id_mappings",
        help="File containing UniProt ID mappings to current identifier",
    )
    parser.add_argument("--organism_tag", help="Tag to label the organism for this run")
    args = parser.parse_args()

    # load the nodes csv file
    nodes_df = pd.read_csv(args.nodes)

    # insert sequence information into nodes_df
    nodes_df = add_sequences(nodes_df, args.fasta)

    # add UniProt IDs
    nodes_df = add_mappings(nodes_df, args.id_mappings)

    # write the updated nodes_df to file
    nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-seqs-{args.output_suffix}.csv",
        index=False,
    )


if __name__ == "__main__":

    main()
