import pandas as pd
import numpy as np
from typing import List
import argparse

def parse_raw_counts(row):
    """Convert raw_counts string to a list of floats."""
    return [float(x) for x in row["raw_counts"].split(",")]


def normalize_globally(study_dfs: List[pd.DataFrame]) -> List[pd.DataFrame]:
    """
    Normalize raw A-site counts across all genes and studies to global CPM.

    Returns updated study_dfs with a 'normalized_counts' column.
    """
    all_counts = []
    for df in study_dfs:
        df["raw_counts_list"] = df.apply(parse_raw_counts, axis=1)
        all_counts.extend([x for sublist in df["raw_counts_list"] for x in sublist])

    global_total = sum(all_counts)

    for df in study_dfs:
        df["normalized_counts"] = df["raw_counts_list"].apply(
            lambda counts: [x / global_total * 1e6 for x in counts]
        )

    return study_dfs


def pooled_ribo_profile(study_dfs: List[pd.DataFrame], return_codon_counts=True):
    """
    Pools globally normalized ribosome profiles across studies.
    """
    all_genes = set()
    for df in study_dfs:
        all_genes.update(df["gene"].unique())

    pooled_data = []

    for gene in sorted(all_genes):
        gene_profiles = []
        for df in study_dfs:
            row = df[df["gene"] == gene]
            if not row.empty:
                gene_profiles.append(row.iloc[0]["normalized_counts"])

        if not gene_profiles:
            continue

        lengths = [len(c) for c in gene_profiles]
        if len(set(lengths)) > 1:
            print(f"Warning: Length mismatch for gene {gene}: {lengths}")

        max_len = max(lengths)
        padded = [np.pad(c, (0, max_len - len(c)), constant_values=0.0) for c in gene_profiles]

        pooled = np.sum(padded, axis=0)

        if return_codon_counts:
            pooled = pooled[::3]

        pooled_data.append({"gene": gene, "pooled_counts": pooled})

    return pd.DataFrame(pooled_data)


def main():
    parser = argparse.ArgumentParser(description="Pool ribosome profiles with global CPM normalization.")
    parser.add_argument("--ribo_seq_data_files", nargs="+", required=True, help="Input TSV files with ribo seq data")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--output_prefix", required=True, help="Prefix appended to output files")
    parser.add_argument("--organism_tag", required=True, help="Organism label to be applied to output files")
    args = parser.parse_args()

    study_dfs = []
    for file_path in args.ribo_seq_data_files:
        df = pd.read_csv(file_path, sep="\t", header=None, names=["gene", "num_ncs", "raw_counts"])
        study_dfs.append(df)

    study_dfs = normalize_globally(study_dfs)

    pooled_df = pooled_ribo_profile(study_dfs)

    pooled_df["avg_dwell"] = pooled_df["pooled_counts"].apply(np.mean)
    mean_dwell = pooled_df["avg_dwell"].mean()
    centered = pooled_df["avg_dwell"] - mean_dwell
    max_dev = np.mean(np.abs(centered))
    pooled_df["translation_speed_score"] = -1.0 * (centered / max_dev)

    pooled_df.to_pickle(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-pooled-ribo-seq-data.pkl")

if __name__ == "__main__":
    main()
