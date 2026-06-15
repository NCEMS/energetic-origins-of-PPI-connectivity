import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


def normalize_gene_id(value):
    """
    Normalize Arabidopsis gene/locus identifiers.

    Examples:
        AT3G22530 -> AT3G22530
        At3g22530 -> AT3G22530
        at3g22530 -> AT3G22530
        NaN       -> None
    """
    if pd.isna(value):
        return None

    value = str(value).strip()

    if not value:
        return None

    return value.upper()


def get_chaperone_metadata(prot, interaction_map, chaperone_info):
    """
    Return metadata for chaperones that interact with a given protein.
    """
    chaps = interaction_map.get(prot, set())

    return [
        {"Gene locus": c, **chaperone_info[c]}
        for c in sorted(chaps)
        if c in chaperone_info
    ]


def main():
    parser = argparse.ArgumentParser(
        description="Add chaperone interaction annotations to network nodes."
    )

    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_suffix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="athaliana")
    parser.add_argument("--input_nodes", required=True)
    parser.add_argument("--chap_data", required=True)
    parser.add_argument("--edges", required=True)

    args = parser.parse_args()

    # Read input files
    chap_df = pd.read_csv(args.chap_data, sep="\t")
    edges_df = pd.read_csv(args.edges, usecols=["source", "target"])
    nodes_df = pd.read_pickle(args.input_nodes)

    # Required columns in the new Arabidopsis chaperone table
    required_chap_cols = {
        "Family",
        "Gene locus",
        "Protein names",
        "Gene name",
        "Presumed subcellular location",
    }

    missing_chap_cols = required_chap_cols - set(chap_df.columns)

    if missing_chap_cols:
        raise ValueError(
            "Chaperone data file is missing required column(s): "
            f"{', '.join(sorted(missing_chap_cols))}. "
            f"Observed columns: {list(chap_df.columns)}"
        )

    required_edge_cols = {"source", "target"}
    missing_edge_cols = required_edge_cols - set(edges_df.columns)

    if missing_edge_cols:
        raise ValueError(
            "Edges file is missing required column(s): "
            f"{', '.join(sorted(missing_edge_cols))}"
        )

    if "node" not in nodes_df.columns:
        raise ValueError("Input nodes DataFrame must contain a 'node' column.")

    n_nodes_before = len(nodes_df)

    # Normalize identifiers
    nodes_df = nodes_df.copy()
    nodes_df["node"] = nodes_df["node"].apply(normalize_gene_id)

    edges_df = edges_df.copy()
    edges_df["source"] = edges_df["source"].apply(normalize_gene_id)
    edges_df["target"] = edges_df["target"].apply(normalize_gene_id)

    chap_df = chap_df.copy()
    chap_df["Gene locus"] = chap_df["Gene locus"].apply(normalize_gene_id)

    # Drop rows without a valid chaperone locus
    chap_df = chap_df.dropna(subset=["Gene locus"])

    # If duplicate chaperone rows exist for the same Gene locus, collapse them
    # to one row per locus. This prevents ambiguous metadata lookup.
    metadata_cols = [
        "Family",
        "Protein names",
        "Gene name",
        "Reference",
        "Entrez GeneID (NCBI)",
        "Presumed subcellular location",
    ]

    metadata_cols = [col for col in metadata_cols if col in chap_df.columns]

    def collapse_unique_values(series):
        values = [
            str(x).strip()
            for x in series.dropna()
            if str(x).strip()
        ]

        unique_values = sorted(set(values))

        if len(unique_values) == 0:
            return None

        return ";".join(unique_values)

    chap_df_collapsed = (
        chap_df[["Gene locus"] + metadata_cols]
        .groupby("Gene locus", as_index=False)
        .agg(collapse_unique_values)
    )

    # Make lookup for chaperone metadata
    chaperone_info = chap_df_collapsed.set_index("Gene locus")[
        metadata_cols
    ].to_dict(orient="index")

    chaperone_loci = set(chaperone_info.keys())

    print(f"Input nodes: {n_nodes_before:,}")
    print(f"Unique chaperone loci in chaperone table: {len(chaperone_loci):,}")

    # Build map of protein -> set of interacting chaperones
    interaction_map = defaultdict(set)

    for _, row in edges_df.iterrows():
        src = row["source"]
        tgt = row["target"]

        if src is None or tgt is None:
            continue

        if tgt in chaperone_loci:
            interaction_map[src].add(tgt)

        if src in chaperone_loci:
            interaction_map[tgt].add(src)

    # Add chaperone interaction information to nodes_df
    nodes_df["interacting_chaperones"] = nodes_df["node"].apply(
        lambda prot: sorted(interaction_map.get(prot, []))
    )

    nodes_df["n_interacting_chaperones"] = nodes_df["interacting_chaperones"].apply(len)

    nodes_df["interacting_chaperone_families"] = nodes_df["interacting_chaperones"].apply(
        lambda chaps: sorted(
            {
                chaperone_info[c].get("Family")
                for c in chaps
                if c in chaperone_info and chaperone_info[c].get("Family") is not None
            }
        )
    )

    nodes_df["n_interacting_chaperone_families"] = nodes_df[
        "interacting_chaperone_families"
    ].apply(len)

    nodes_df["interacting_chaperone_info"] = nodes_df["node"].apply(
        lambda prot: get_chaperone_metadata(prot, interaction_map, chaperone_info)
    )

    n_nodes_after = len(nodes_df)

    if n_nodes_after != n_nodes_before:
        raise ValueError(
            f"Row count changed during chaperone annotation: "
            f"{n_nodes_before:,} -> {n_nodes_after:,}. "
            "This should not happen."
        )

    n_nodes_with_chaperones = nodes_df["n_interacting_chaperones"].gt(0).sum()

    print(f"Nodes with at least one interacting chaperone: {n_nodes_with_chaperones:,}")

    print("Top chaperone family counts among annotated node interactions:")
    family_counts = (
        nodes_df["interacting_chaperone_families"]
        .explode()
        .dropna()
        .value_counts()
    )

    if len(family_counts) > 0:
        print(family_counts.to_string())
    else:
        print("No chaperone family interactions found.")

    # Save output
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    output_path = (
        f"{args.output_dir}/"
        f"{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )

    nodes_df.to_pickle(output_path)

    print(f"Wrote output: {output_path}")


if __name__ == "__main__":
    main()
