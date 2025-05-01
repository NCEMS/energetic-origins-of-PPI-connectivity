import os, sys
import numpy as np
import pandas as pd
import argparse
import typing
import pint
import pint_pandas
from pathlib import Path

def select_structure(row):
    if row.get("structure_exists", 0) == 1:
        if pd.notna(row["cleaved_structure_path"]):
            return row["cleaved_structure_path"]
        else:
            return row["structure_path"]
    else:
        return None

def get_score_file(row, score_dir="processed-data/scores"):
    structure_path = select_structure(row)
    if pd.notna(structure_path):
        stem = Path(structure_path).stem
        return score_dir / f"{stem}_0001.sc"
    else:
        return None

def parse_rosetta_score_file(score_file_path):

    if score_file_path is None or not score_file_path.exists():
        return None  # missing file

    with open(score_file_path, "r") as f:
        lines = f.readlines()

    header = None
    for line in lines:
        if line.startswith("SCORE:") and "total_score" in line:
            # This is the header line
            header = line.strip().split()[1:]  # Skip "SCORE:" keyword
            continue

        if line.startswith("SCORE:") and header:
            # This is the first data line
            values = line.strip().split()[1:]  # Skip "SCORE:" keyword
            score_dict = dict(zip(header, map(float_or_str, values)))
            if "total_score" not in score_dict:
                print(f"Warning: 'total_score' not found in {score_file_path}")
            return score_dict

def float_or_str(x):
    try:
        return float(x)
    except ValueError:
        return x

def main():

    # parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Extract Rosetta energy scoring of protein structures and add to nodes_df"
    )
    parser.add_argument(
        "--nodes", required=True
    )
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where results will be saved (default: 'outputs')",
    )
    parser.add_argument(
        "--input_dir"
    )
    parser.add_argument("--organism_tag")
    parser.add_argument("--output_prefix", default="0", help="Prefix for output files")
    args = parser.parse_args()

    # read in the nodes_df
    nodes_df = pd.read_pickle(args.nodes)

    # add path to the Rosetta score file
    score_dir = Path(args.input_dir)
    nodes_df["Rosetta_score_file"] = nodes_df.apply(lambda row: get_score_file(row, score_dir), axis=1)

    # add Rosetta scoring function information to the nodes_df
    nodes_df["Rosetta_scores"] = nodes_df["Rosetta_score_file"].apply(parse_rosetta_score_file)

    # expand dictionary into multiple columns
    rosetta_scores_df = nodes_df["Rosetta_scores"].apply(pd.Series)
    rosetta_scores_df = rosetta_scores_df.add_prefix("Rosetta_")
    nodes_df = pd.concat([nodes_df, rosetta_scores_df], axis=1)
    nodes_df = nodes_df.drop(columns=["Rosetta_scores"])

    nodes_df.to_pickle(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta.pkl"
    )


if __name__ == "__main__":
    main()
