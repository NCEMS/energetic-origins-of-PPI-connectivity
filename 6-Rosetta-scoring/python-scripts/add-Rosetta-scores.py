#!/usr/bin/env python3
import os, sys
import argparse
from pathlib import Path
from typing import Optional, List, Dict, Tuple

import numpy as np
import pandas as pd

# -----------------------------
# Helpers for structure naming
# -----------------------------

def select_structure(row) -> Optional[str]:
    """
    Original gating from your script: use cleaved_structure_path if present,
    else structure_path; only when structure_exists == 1.
    """
    if row.get("structure_exists", 0) == 1:
        if pd.notna(row.get("cleaved_structure_path")):
            return row["cleaved_structure_path"]
        else:
            return row.get("structure_path")
    return None


def candidate_stems_for_row(row) -> List[str]:
    """
    Return possible stems for this protein's outputs, to support both naming styles:
      1) <stem> from the selected structure path (historical)
      2) row["node"] (label-based runs)
    We try them in this order.
    """
    stems = []
    structure_path = select_structure(row)
    if pd.notna(structure_path):
        stems.append(Path(structure_path).stem)
    node_label = row.get("node")
    if pd.notna(node_label):
        stems.append(str(node_label))
    # Deduplicate while preserving order
    seen = set()
    uniq = []
    for s in stems:
        if s and s not in seen:
            uniq.append(s); seen.add(s)
    return uniq


# -----------------------------
# Rosetta score parsing
# -----------------------------

def float_or_str(x: str):
    try:
        return float(x)
    except ValueError:
        return x

def parse_rosetta_sc_one(score_file_path: Path) -> Optional[Dict[str, float]]:
    """
    Parse *one* Rosetta .sc file and return a dict of scores for the FIRST data line.
    Returns None if file is missing or malformed.
    """
    if score_file_path is None or not score_file_path.exists():
        return None

    with open(score_file_path, "r") as f:
        header = None
        for line in f:
            if not line.startswith("SCORE:"):
                continue
            fields = line.strip().split()
            if "total_score" in fields:
                # Header row (skip leading 'SCORE:')
                header = fields[1:]
                continue
            if header is not None:
                # First data row
                values = fields[1:]
                # Make dict; cast numerics where possible
                d = dict(zip(header, map(float_or_str, values)))
                return d
    return None


# -----------------------------
# Collect per-protein replicates
# -----------------------------

def list_replicate_files(score_dir: Path, stems: List[str], nstruct: Optional[int]) -> List[Path]:
    """
    Return the list of candidate .sc files for the first stem that yields hits.
    - If nstruct is given, expect exact files: <stem>_0001.sc .. <stem>_{nstruct:04d}.sc
    - If nstruct is None, auto-discover by globbing <stem>_*.sc and filtering for 4-digit suffix.
    We try stems in order; as soon as one stem yields any files, we use that set.
    """
    def suffix_is_4digit(p: Path) -> bool:
        # Match *_0001.sc style
        name = p.stem  # e.g., "STEM_0001"
        if "_" not in name: return False
        suf = name.split("_")[-1]
        return len(suf) == 4 and suf.isdigit()

    for stem in stems:
        files = []
        if nstruct and nstruct > 0:
            files = [score_dir / f"{stem}_{i:04d}.sc" for i in range(1, nstruct + 1)]
            files = [p for p in files if p.exists()]
        else:
            files = sorted([p for p in score_dir.glob(f"{stem}_*.sc") if suffix_is_4digit(p)])
        if files:
            return files
    return []


def best_total_score_for_row(row, score_dir: Path, nstruct: Optional[int]) -> Dict[str, object]:
    """
    For a given row, scan replicate .sc files, parse scores, and return:
      {
        'Rosetta_best_total_score': float or np.nan,
        'Rosetta_best_pose_index':  int or np.nan,
        'Rosetta_best_score_file':  str or None
      }
    If nothing found, returns NaNs/None.
    """
    stems = candidate_stems_for_row(row)
    files = list_replicate_files(score_dir, stems, nstruct)
    if not files:
        return {
            "Rosetta_best_total_score": np.nan,
            "Rosetta_best_pose_index":  np.nan,
            "Rosetta_best_score_file":  None,
        }

    best_score = None
    best_file = None
    best_pose = None

    for p in files:
        d = parse_rosetta_sc_one(p)
        if not d or "total_score" not in d:
            continue
        ts = d["total_score"]
        # Some parsers might leave non-floats; skip if not numeric
        try:
            tsf = float(ts)
        except Exception:
            continue
        if (best_score is None) or (tsf < best_score):
            best_score = tsf
            best_file = str(p)
            # Extract pose index from *_0001.sc
            try:
                pose_idx = int(Path(p).stem.split("_")[-1])
            except Exception:
                pose_idx = None
            best_pose = pose_idx

    if best_score is None:
        return {
            "Rosetta_best_total_score": np.nan,
            "Rosetta_best_pose_index":  np.nan,
            "Rosetta_best_score_file":  None,
        }

    return {
        "Rosetta_best_total_score": best_score,
        "Rosetta_best_pose_index":  best_pose if best_pose is not None else np.nan,
        "Rosetta_best_score_file":  best_file,
    }


# -----------------------------
# Main
# -----------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extract best (lowest) Rosetta total_score across replicates and add to nodes_df"
    )
    parser.add_argument("--nodes", required=True, help="Path to pickled nodes DataFrame")
    parser.add_argument("--output_dir", default="processed-data",
                        help="Directory to write the updated pickle (default: processed-data)")
    parser.add_argument("--input_dir", required=True,
                        help="Directory containing Rosetta .sc files (scores)")
    parser.add_argument("--organism_tag", required=True)
    parser.add_argument("--output_prefix", default="0", help="Prefix for output filename")
    parser.add_argument("--nstruct", type=int, default=None,
                        help="Number of replicates per protein. If omitted, auto-discover by glob.")
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)
    score_dir = Path(args.input_dir)
    score_dir.mkdir(parents=True, exist_ok=True)

    # Compute best score info per row
    best_df = nodes_df.apply(
        lambda row: best_total_score_for_row(row, score_dir, args.nstruct),
        axis=1, result_type="expand"
    )

    # Merge into nodes_df
    nodes_df = pd.concat([nodes_df, best_df], axis=1)

    out_path = (
        Path(args.output_dir)
        / f"{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta.pkl"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    nodes_df.to_pickle(out_path)
    print(f"Wrote: {out_path}")

if __name__ == "__main__":
    main()
