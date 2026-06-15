import pandas as pd
from Bio import SeqIO
import typing
import argparse
import numpy as np
from collections import defaultdict

def get_base_locus_id(record_id: str) -> str:
    """
    Convert a TAIR protein isoform ID to a base locus ID.

    Example:
        AT1G01010.1 -> AT1G01010
    """
    return record_id.split(".")[0]


def add_sequences(nodes_df: pd.DataFrame, fasta_file: str) -> pd.DataFrame:
    """
    Adds protein sequence information to input DataFrame by matching base gene/locus names.

    Node names are expected to be base locus IDs with isoform suffixes removed,
    for example ATG00001. If multiple isoforms are present in the FASTA file,
    this function assigns the sequence from isoform 1, for example ATG00001.1.

    Args:
        nodes_df (pd.DataFrame): nodes information, loaded based on output from step 1
        fasta_file (str): path to file containing protein sequence information in fasta format

    Returns:
        Updated nodes_df containing sequence information where isoform .1 is available
        and np.nan where not available.
    """

    # Read sequences and group them by base locus ID
    seqs_by_locus = defaultdict(list)

    for record in SeqIO.parse(fasta_file, "fasta"):
        locus_id = get_base_locus_id(record.id)
        seqs_by_locus[locus_id].append(record)

    # Select the .1 isoform for each base locus ID
    seqs = {}
    loci_without_isoform_1 = {}

    for locus_id, records in seqs_by_locus.items():
        expected_isoform_1_id = f"{locus_id}.1"

        isoform_1_records = [
            record for record in records
            if record.id == expected_isoform_1_id
        ]

        if len(isoform_1_records) == 1:
            seqs[locus_id] = isoform_1_records[0]

        elif len(isoform_1_records) > 1:
            raise ValueError(
                f"Multiple records found with ID {expected_isoform_1_id}. "
                "FASTA IDs should be unique."
            )

        else:
            loci_without_isoform_1[locus_id] = [record.id for record in records]

    # Report nodes with FASTA records but no matching .1 isoform
    nodes_missing_isoform_1 = {
        node: loci_without_isoform_1[node]
        for node in nodes_df["node"]
        if node in loci_without_isoform_1
    }

    if nodes_missing_isoform_1:
        print("Warning: some nodes had sequence records but no matching .1 isoform.")
        print(f"Number of affected nodes: {len(nodes_missing_isoform_1):,}")
        print("Examples:")

        for node, isoform_ids in list(nodes_missing_isoform_1.items())[:10]:
            print(f"{node}: expected {node}.1; found {', '.join(isoform_ids)}")

    # Add column stating which nodes have isoform .1 sequence info
    nodes_df["has_verified_sequence"] = nodes_df["node"].isin(seqs.keys())

    # Add sequence information while keeping the node name unchanged
    nodes_df["sequence"] = nodes_df["node"].apply(
        lambda x: str(seqs[x].seq).rstrip("*") if x in seqs else np.nan
    )

    return nodes_df


def add_mappings(nodes_df: pd.DataFrame, map_file: str) -> pd.DataFrame:
    """
    Add UniProtKB-AC identifiers to each node without increasing the number of rows.

    The UniProt idmapping file can contain multiple UniProtKB accessions for the
    same ordered locus name. To preserve one row per network node, this function
    collapses all accessions for a locus into a semicolon-delimited string before
    merging.

    Args:
        nodes_df (pd.DataFrame): current nodes_df to be updated
        map_file (str): path to the *_idmapping.dat file

    Returns:
        Updated nodes_df containing mapped names in column 'UniProtKB-AC'
    """

    n_before = len(nodes_df)

    column_names = ["UniProtKB-AC", "ID_type", "ID"]
    cross_df = pd.read_csv(map_file, names=column_names, sep="\t")

    cross_df = cross_df[cross_df["ID_type"] == "Gene_OrderedLocusName"].copy()
    cross_df["ID"] = cross_df["ID"].str.upper()

    # Collapse to one row per locus ID to avoid one-to-many merge expansion
    cross_df_collapsed = (
        cross_df
        .dropna(subset=["ID", "UniProtKB-AC"])
        .drop_duplicates(subset=["ID", "UniProtKB-AC"])
        .sort_values(["ID", "UniProtKB-AC"])
        .groupby("ID", as_index=False)["UniProtKB-AC"]
        .agg(lambda accessions: ";".join(accessions))
    )

    # Optional diagnostic: how many locus IDs had multiple UniProt mappings?
    n_multi_mapping_ids = (
        cross_df
        .drop_duplicates(subset=["ID", "UniProtKB-AC"])
        .groupby("ID")["UniProtKB-AC"]
        .nunique()
        .gt(1)
        .sum()
    )

    print(f"Input nodes: {n_before:,}")
    print(f"Locus IDs with multiple UniProtKB-AC mappings: {n_multi_mapping_ids:,}")

    nodes_df = nodes_df.merge(
        cross_df_collapsed,
        left_on="node",
        right_on="ID",
        how="left",
        validate="one_to_one",
    )

    nodes_df.drop(columns=["ID"], inplace=True)

    n_after = len(nodes_df)

    if n_after != n_before:
        raise ValueError(
            f"Row count changed during UniProt mapping merge: "
            f"{n_before:,} -> {n_after:,}. This should not happen."
        )

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
