import pandas as pd
import typing
import argparse
import re


def add_UniProt_info(nodes_df: pd.DataFrame, uniprot_data: str) -> pd.DataFrame:
    """
    Function that reads in a pre-processed annotation file from UniProt and adds selected information from it to nodes_df

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


def add_PTM_exchange_info(nodes_df: pd.DataFrame, ptm_file: str) -> pd.DataFrame:
    """
    Adds Gold, Silver, and Bronze PTM annotations from PTMeXchange to nodes_df.

    Args:
        nodes_df (pd.DataFrame): Existing nodes DataFrame with UniProtKB-AC column.
        ptm_file (str): Path to the PTM Exchange CSV file.

    Returns:
        pd.DataFrame: nodes_df with added columns: ptm_gold, ptm_silver, ptm_bronze
    """
    ptm_df = pd.read_csv(ptm_file)

    # drop rows with missing PTM data
    ptm_df = ptm_df.dropna(subset=["additional_PTMs", "UniProtKB-AC"])

    # split each PTM entry into individual site entries with category
    expanded_rows = []
    for _, row in ptm_df.iterrows():
        accession = row["UniProtKB-AC"]
        ptms = row["additional_PTMs"].split(";")
        for ptm in ptms:
            match = re.match(
                r"([STYACDEFGHIKLMNPQRVW]{1}\d+)\((Gold|Silver|Bronze)\)", ptm.strip()
            )
            if match:
                site, category = match.groups()
                expanded_rows.append((accession, site, category))

    expanded_df = pd.DataFrame(
        expanded_rows, columns=["UniProtKB-AC", "site", "category"]
    )

    # aggregate sites by UniProt ID and category
    categorized = (
        expanded_df.groupby(["UniProtKB-AC", "category"])["site"]
        .apply(lambda x: ";".join(sorted(set(x))))
        .unstack(fill_value="")
    )

    # rename columns
    categorized = categorized.rename(
        columns={"Gold": "ptm_gold", "Silver": "ptm_silver", "Bronze": "ptm_bronze"}
    ).reset_index()

    merged = nodes_df.merge(categorized, on="UniProtKB-AC", how="left")

    return merged


def main():

    parser = argparse.ArgumentParser(
        description="Add UniProt database information to PPI network"
    )
    parser.add_argument(
        "--nodes",
        type=str,
        help="Path to input nodes csv file",
    )
    parser.add_argument(
        "--uniprot",
        type=str,
        help="Path to input UniProt.dat file",
    )
    parser.add_argument("--ptmexchange")
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path to the directory to which data will be written",
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        help="Prefix to be appended to the output file",
    )
    parser.add_argument(
        "--output_suffix",
        type=str,
        help="Suffix to be appended to the output file",
    )
    parser.add_argument("--organism_tag", help="Tag to label the organism for this run")
    args = parser.parse_args()

    # read in the nodes_df file to which annotations will be added
    nodes_df = pd.read_csv(args.nodes)

    # update nodes_df with annotations from UniProt
    nodes_df = add_UniProt_info(nodes_df, args.uniprot)

    # add additional PTM information from PTMeXchange
    #nodes_df = add_PTM_exchange_info(nodes_df, args.ptmexchange)

    # save the updated nodes_df to file
    nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
