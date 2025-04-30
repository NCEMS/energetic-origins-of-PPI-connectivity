import os, sys
import numpy as np
import pandas as pd
import argparse
import typing
import pint
import pint_pandas
from pathlib import Path
import re

def select_structure(row):
    if row.get("structure_exists", 0) == 1:
        if pd.notna(row["cleaved_structure_path"]):
            return row["cleaved_structure_path"]
        else:
            return row["structure_path"]
    else:
        return None

def get_log_file(row, score_dir="processed-data/scores"):
    structure_path = select_structure(row)
    if pd.notna(structure_path):
        stem = Path(structure_path).stem
        return score_dir / f"{stem}.log"
    else:
        return None

def parse_foldx_stdout(log_file_path):
    """
    Parse FoldX .foldx.log file to extract energy terms.
    """
    if not log_file_path or not Path(log_file_path).exists():
        return None  # missing file

    try:
        with open(log_file_path, "r") as f:
            content = f.read()
    except Exception as e:
        print(f"Warning: Failed to read {log_file_path}: {e}")
        return None

    scores = {}
    pattern = re.compile(r"(.+?)\s*=\s*(-?\d+\.\d+)")
    for match in pattern.finditer(content):
        key_raw, value = match.groups()
        key_clean = key_raw.strip().lower().replace(" ", "_")
        if key_clean == "total":
            key_clean = "total_energy"
        key_clean = "FoldX_" + key_clean
        scores[key_clean] = float(value)

    return scores

def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Parse FoldX output logs and add energy scores to nodes_df"
    )
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--output_dir", default="processed-data", help="Directory where results will be saved")
    parser.add_argument("--input_dir")
    parser.add_argument("--organism_tag")
    parser.add_argument("--output_prefix", default="0")
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)

    score_dir = Path(args.input_dir)

    # 1. Point to .foldx.log files instead of .fxout
    nodes_df["FoldX_score_file"] = nodes_df.apply(lambda row: get_log_file(row, score_dir), axis=1)

    # 2. Parse .foldx.log files
    nodes_df["FoldX_scores"] = nodes_df["FoldX_score_file"].apply(parse_foldx_stdout)

    # 3. Expand the dictionary into columns
    foldx_scores_df = nodes_df["FoldX_scores"].apply(pd.Series)
    nodes_df = pd.concat([nodes_df, foldx_scores_df], axis=1)

    # 4. Drop intermediate column
    nodes_df = nodes_df.drop(columns=["FoldX_scores"])

    # 5. Save the updated DataFrame
    output_path = Path(args.output_dir) / f"{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX.pkl"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    nodes_df.to_pickle(output_path)

if __name__ == "__main__":
    main()
