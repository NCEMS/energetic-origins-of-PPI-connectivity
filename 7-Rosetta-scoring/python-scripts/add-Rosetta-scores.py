#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Optional, List, Dict, Tuple

import numpy as np
import pandas as pd

# -----------------------------
# Naming helpers
# -----------------------------


def candidate_stems_for_row(row) -> List[str]:
    """
    Prefer the new labeling scheme (row['node']), but also fall back to
    final_structure_path / cleaved_structure_path / structure_path stems.
    """
    stems: List[str] = []
    node_label = row.get("node")
    if pd.notna(node_label):
        stems.append(str(node_label))

    for key in ("final_structure_path", "cleaved_structure_path", "structure_path"):
        p = row.get(key)
        if pd.notna(p) and isinstance(p, str) and p.strip():
            stems.append(Path(p).stem)

    # dedupe preserving order
    seen = set()
    out = []
    for s in stems:
        if s and s not in seen:
            out.append(s)
            seen.add(s)
    return out


# -----------------------------
# Rosetta score parsing
# -----------------------------


def float_or_str(x: str):
    try:
        return float(x)
    except Exception:
        return x


def parse_rosetta_sc_first(score_file_path: Path) -> Optional[Dict[str, float]]:
    """
    Parse one Rosetta .sc file and return dict for the FIRST data row.
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
                header = fields[1:]
                continue
            if header is not None:
                values = fields[1:]
                return dict(zip(header, map(float_or_str, values)))
    return None


def _suffix_is_4digit(p: Path) -> Optional[int]:
    """
    Return 4-digit replicate index if name matches *_dddd.sc, else None.
    """
    stem = p.stem  # e.g., STEM_0001
    if "_" not in stem:
        return None
    suf = stem.split("_")[-1]
    return int(suf) if len(suf) == 4 and suf.isdigit() else None


def list_replicate_files(
    score_dir: Path, stems: List[str], nstruct: Optional[int]
) -> List[Path]:
    """
    Return the list of .sc files for the first stem that yields hits.
    If nstruct is provided, look for exact STEM_0001..STEM_NNNN and keep existing.
    Otherwise, glob STEM_*.sc and filter *_dddd.sc.
    """
    for stem in stems:
        if nstruct and nstruct > 0:
            files = [score_dir / f"{stem}_{i:04d}.sc" for i in range(1, nstruct + 1)]
            files = [p for p in files if p.exists()]
        else:
            files = sorted(
                [
                    p
                    for p in score_dir.glob(f"{stem}_*.sc")
                    if _suffix_is_4digit(p) is not None
                ]
            )
        if files:
            return files
    return []


# -----------------------------
# Per-row aggregation
# -----------------------------


def per_replicate_scores_for_row(
    row, score_dir: Path, nstruct: Optional[int]
) -> Dict[str, object]:
    """
    Build a dict with:
      - Rosetta_total_score_0001 .. Rosetta_total_score_NNNN (NaN if missing)
      - Rosetta_best_total_score
      - Rosetta_best_pose_index
      - Rosetta_best_score_file
      - Rosetta_n_found
    """
    stems = candidate_stems_for_row(row)
    files = list_replicate_files(score_dir, stems, nstruct)

    # Collect totals by pose index
    totals: Dict[int, float] = {}
    filemap: Dict[int, str] = {}

    for p in files:
        idx = _suffix_is_4digit(p)
        if idx is None:
            continue
        d = parse_rosetta_sc_first(p)
        if not d or "total_score" not in d:
            continue
        try:
            ts = float(d["total_score"])
        except Exception:
            continue
        totals[idx] = ts
        filemap[idx] = str(p)

    out: Dict[str, object] = {}

    # Decide which indices to materialize as columns
    if nstruct and nstruct > 0:
        pose_indices = list(range(1, nstruct + 1))
    else:
        pose_indices = sorted(totals.keys())

    # One column per replicate
    for i in pose_indices:
        out[f"Rosetta_total_score_{i:04d}"] = totals.get(i, np.nan)

    # Best over available replicates
    if totals:
        best_idx = min(totals, key=lambda k: totals[k])
        out["Rosetta_best_total_score"] = totals[best_idx]
        out["Rosetta_best_pose_index"] = best_idx
        out["Rosetta_best_score_file"] = filemap.get(best_idx)
    else:
        out["Rosetta_best_total_score"] = np.nan
        out["Rosetta_best_pose_index"] = np.nan
        out["Rosetta_best_score_file"] = None

    out["Rosetta_n_found"] = int(len(totals))
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Add per-replicate Rosetta total_score columns and best score to nodes_df"
    )
    parser.add_argument(
        "--nodes", required=True, help="Path to pickled nodes DataFrame"
    )
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory to write the updated pickle (default: processed-data)",
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing Rosetta .sc files (scores)",
    )
    parser.add_argument("--organism_tag", required=True)
    parser.add_argument(
        "--output_prefix", default="0", help="Prefix for output filename"
    )
    parser.add_argument(
        "--output_suffix", required=True, help="Suffix for output filename"
    )
    parser.add_argument(
        "--nstruct",
        type=int,
        default=None,
        help="Number of replicates per protein. If omitted, auto-discover by glob.",
    )
    parser.add_argument("--nstruc", type=int, dest="nstruct")  # legacy alias

    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)
    score_dir = Path(args.input_dir)
    score_dir.mkdir(parents=True, exist_ok=True)

    # Build per-row dicts, then concat
    perrow = nodes_df.apply(
        lambda row: per_replicate_scores_for_row(row, score_dir, args.nstruct),
        axis=1,
        result_type="expand",
    )
    out_df = pd.concat([nodes_df, perrow], axis=1)

    out_path = (
        Path(args.output_dir)
        / f"{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_pickle(out_path)
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
