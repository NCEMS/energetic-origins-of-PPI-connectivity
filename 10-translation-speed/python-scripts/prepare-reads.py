import pandas as pd
import numpy as np
from typing import List
import argparse

def normalize_counts(row):
    """Normalize raw_counts string to CPM and return list of floats."""
    counts = [float(x) for x in row["raw_counts"].split(",")]
    total = sum(counts)
    if total > 0:
        norm_counts = [x / total * 1e6 for x in counts]
    else:
        norm_counts = [0.0] * len(counts)
    return norm_counts


def pooled_ribo_profile(study_dfs: List[pd.DataFrame], return_codon_counts=True):
    """
    Pools normalized ribosome profiles across multiple studies.

    Parameters:
    - study_dfs: list of DataFrames, each with a 'gene' column and 'normalized_counts' column
    - return_codon_counts: whether to reduce to canonical frame (default: True)

    Returns:
    - pooled_df: DataFrame with gene and pooled codon-level counts
    """
    # Get union of all genes
    all_genes = set()
    for df in study_dfs:
        all_genes.update(df["gene"].unique())

    pooled_data = []

    for gene in sorted(all_genes):
        # Collect all normalized profiles for this gene across studies
        gene_profiles = []
        for df in study_dfs:
            row = df[df["gene"] == gene]
            if not row.empty:
                gene_profiles.append(row.iloc[0]["normalized_counts"])

        if not gene_profiles:
            continue

        # Pad to same length
        max_len = max(len(c) for c in gene_profiles)
        padded = [np.pad(c, (0, max_len - len(c)), constant_values=0.0) for c in gene_profiles]

        # Sum across studies
        pooled = np.sum(padded, axis=0)

        if return_codon_counts:
            pooled = pooled[::3]

        pooled_data.append({"gene": gene, "pooled_counts": pooled})

    return pd.DataFrame(pooled_data)


def main():

    parser = argparse.ArgumentParser(description="Pool ribosome profiles.")
    parser.add_argument("--ribo_seq_data_files", nargs="+", required=True, help="Input TSV files with ribo seq data")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--output_prefix", required=True, help="Prefix appended to output files")
    parser.add_argument("--organism_tag", required=True, help="Organism label to be applied to output files")
    args = parser.parse_args()

    study_dfs = []
    for file_path in args.ribo_seq_data_files:
        df = pd.read_csv(file_path, sep="\t", header=None, names=["gene", "num_ncs", "raw_counts"])
        df["normalized_counts"] = df.apply(normalize_counts, axis=1)
        study_dfs.append(df)

    pooled_df = pooled_ribo_profile(study_dfs)

    # save the output file
    pooled_df.to_csv(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-pooled-ribo-seq-data.csv", index=False)


if __name__ == "__main__":
    main()
