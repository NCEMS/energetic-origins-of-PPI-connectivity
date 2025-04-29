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
        return score_dir / f"{stem}.fxout"
    else:
        return None

def parse_FoldX_score_file(score_file_path, file_columns):
    if not score_file_path or not Path(score_file_path).exists():
        return None  # file missing

    try:
        with open(score_file_path, "r") as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) == len(file_columns):
                    return dict(zip(file_columns, _coerce_types(fields)))
    except Exception as e:
        print(f"Warning: Failed to parse {score_file_path}: {e}")
        return None

    return None  # fallback if no line matches

def _coerce_types(fields):
    def try_float(x):
        try:
            return float(x)
        except ValueError:
            return x
    return [try_float(val) for val in fields]

def main():

    # parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Run FoldX energy scoring of protein structures"
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

    nodes_df = pd.read_pickle(args.nodes)

    score_dir = Path(args.input_dir)

    nodes_df["FoldX_score_file"] = nodes_df.apply(lambda row: get_score_file(row, score_dir), axis=1)

    foldx_columns = [
        "PDB",
        "total_energy",
        "backbone_hbond",
        "sidechain_hbond",
        "van_der_waals_clashes",
        "electrostatics",
        "solvation_polar",
        "solvation_hydrophobic",
        "van_der_waals",
        "entropy_mainchain",
        "entropy_sidechain",
        "torsional_clash",
        "backbone_clash",
        "helix_dipole",
        "water_bridges",
        "disulfide",
        "electrostatic_kon",
        "partial_covalent",
        "energy_Ionisation",
        "entropy_complex",
        "number_of_residues",
        "interface_energy"
    ]

    nodes_df["FoldX_scores"] = nodes_df["FoldX_score_file"].apply(lambda path: parse_FoldX_score_file(path, foldx_columns))

    nodes_df.to_pickle(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX.pkl"
    )

if __name__ == "__main__":

    main()
