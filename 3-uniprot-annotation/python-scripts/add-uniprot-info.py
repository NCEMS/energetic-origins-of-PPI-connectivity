import pandas as pd
import typing
import argparse


def add_UniProt_info(nodes_df: pd.DataFrame, uniprot_data: str) -> pd.DataFrame:
    """
    Function that reads in a pre-processed annotation file from UniProt and adds selected information it to nodes_df

    Args:
        nodes_df (pd.DataFrame): nodes DataFrame so far; must contain mapped UniProt IDs for this to work (added by process_nodes)
        uniprot_data (str): path to the pre-processed UniProt information for the organism currently under study

    Returns:
        pd.DataFrame
    """

    uniprot_df = pd.read_csv(uniprot_data)

    to_add = [
        "PrimaryAccession",
        "ProteinName",
        "GO_terms",
        "GO_terms_human_readable",
        "localization_keywords",
        "parsed_functions",
        "parsed_PTMs",
    ]

    return nodes_df.merge(
        uniprot_df[to_add],
        left_on="UniProtKB-AC",
        right_on="PrimaryAccession",
        how="left",
    ).drop(columns=["PrimaryAccession"])


def main():

    parser = argparse.ArgumentParser(
        description="Process UniProt database to extract entries matching organism name."
    )
    parser.add_argument(
        "--nodes",
        default="../2-sequence-parsing/processed-data/0_nodes-centrality-seqs-DeepTMHMM.csv",
        type=str,
        help="Path to input nodes csv file",
    )
    parser.add_argument(
        "--uniprot",
        default="../0-download-inputs/data-files/uniprot_sprot.dat",
        type=str,
        help="Path to input UniProt.dat file",
    )
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        type=str,
        help="Path to the directory to which data will be written",
    )
    parser.add_argument(
        "--output_prefix",
        default="0_",
        type=str,
        help="Prefix to be appended to the output file",
    )
    parser.add_argument(
        "--organism_tag", default="s288c", help="Tag to label the organism for this run"
    )
    args = parser.parse_args()

    # read in the nodes_df file to which annotations will be added
    nodes_df = pd.read_csv(args.nodes)

    # update nodes_df with annotations from UniProt
    nodes_df = add_UniProt_info(nodes_df, args.uniprot)

    # save the updated nodes_df to file
    nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
