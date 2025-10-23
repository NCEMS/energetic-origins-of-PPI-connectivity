#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple

import numpy as np
import pandas as pd


# -----------------------------
# Stem discovery (aligns with Rosetta/FoldX naming)
# -----------------------------


def candidate_stems_for_row(row: pd.Series) -> List[str]:
    """Prefer 'node' label; fall back to structure-path stems."""
    stems: List[str] = []
    node_label = row.get("node")
    if pd.notna(node_label):
        stems.append(str(node_label))
    for key in ("final_structure_path", "cleaved_structure_path", "structure_path"):
        p = row.get(key)
        if pd.notna(p) and isinstance(p, str) and p.strip():
            stems.append(Path(p).stem)
    # dedupe preserving order
    seen, out = set(), []
    for s in stems:
        if s and s not in seen:
            out.append(s)
            seen.add(s)
    return out


# -----------------------------
# File discovery for replicates
# -----------------------------


def first_existing_for_index(score_dir: Path, stem: str, idx: int) -> Optional[Path]:
    """
    For a given stem and replicate index, return the first existing file
    in this priority: _ST.fxout -> .fxout -> _Summary.fxout -> .log
    """
    base = f"{stem}_{idx:04d}"
    candidates = [
        score_dir / f"{base}_ST.fxout",
        score_dir / f"{base}.fxout",
        score_dir / f"{base}_Summary.fxout",
        score_dir / f"{base}.log",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def list_replicate_files(
    score_dir: Path, stems: List[str], nstruct: Optional[int]
) -> List[Tuple[int, Path]]:
    """
    Try stems in order; return list of (pose_index, file_path) for the first stem that yields hits.
    If nstruct is given, check exactly indices 1..nstruct; otherwise, glob and infer indices.
    """

    def infer_idx_from_name(p: Path) -> Optional[int]:
        m = re.search(
            r"_(\d{4})(?:\.\w+|_ST\.fxout|_Summary\.fxout|\.fxout|\.log)$", p.name
        )
        return int(m.group(1)) if m else None

    for stem in stems:
        found: List[Tuple[int, Path]] = []
        if nstruct and nstruct > 0:
            for i in range(1, nstruct + 1):
                p = first_existing_for_index(score_dir, stem, i)
                if p is not None:
                    found.append((i, p))
        else:
            # discover by glob
            pats = [f"{stem}_*.fxout", f"{stem}_*.log"]
            hits = []
            for pat in pats:
                hits.extend(score_dir.glob(pat))
            # unique by inferred index
            seen = set()
            for h in sorted(hits):
                i = infer_idx_from_name(h)
                if i is not None and i not in seen:
                    seen.add(i)
                    found.append((i, h))
        if found:
            return sorted(found, key=lambda t: t[0])
    return []


# -----------------------------
# Parsing FoldX outputs
# -----------------------------


def parse_foldx_fxout(path: Path) -> Optional[float]:
    """
    Parse FoldX *.fxout to get total energy.
    Handles both:
      1) Headered (TSV) files with a 'Total energy' column.
      2) Headerless single-line files: '<pdb>\t<total_energy>\t...'
    """
    try:
        with open(path, "r") as f:
            # keep non-empty lines
            lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    except Exception:
        return None
    if not lines:
        return None

    # --- Try headered format first (any line that looks like a header) ---
    header_idx = None
    headers = []
    for i, ln in enumerate(lines):
        parts = re.split(r"\s*\t\s*|\s{2,}", ln.strip())
        if len(parts) > 1 and any("energy" in h.lower() for h in parts):
            header_idx = i
            headers = parts
            break

    if header_idx is not None:
        norm = [re.sub(r"[\s_]+", " ", h).strip().lower() for h in headers]
        try:
            col_idx = norm.index("total energy")
        except ValueError:
            col_idx = next(
                (j for j, h in enumerate(norm) if "total" in h and "energy" in h), None
            )
            if col_idx is None:
                return None
        # take first data row with a numeric value
        for ln in lines[header_idx + 1 :]:
            row = re.split(r"\s*\t\s*|\s{2,}", ln.strip())
            if col_idx < len(row):
                try:
                    return float(row[col_idx])
                except Exception:
                    continue
        return None

    # --- Headerless format: first non-empty line like "<pdb>\t<num>\t..." ---
    # e.g.: ./YKL048C_0008.pdb\t560.749\t...
    first = lines[0].strip()
    toks = re.split(r"\s*\t\s*|\s{2,}", first)
    if len(toks) >= 2:
        # Usually toks[0] is a *.pdb path; find the first numeric token after it
        for tok in toks[1:]:
            try:
                return float(tok)
            except Exception:
                continue

    # If it's multi-line but still headerless, scan all lines for the first numeric after a *.pdb token
    for ln in lines:
        toks = re.split(r"\s*\t\s*|\s{2,}", ln.strip())
        if len(toks) >= 2 and toks[0].lower().endswith(".pdb"):
            for tok in toks[1:]:
                try:
                    return float(tok)
                except Exception:
                    continue

    return None


def parse_foldx_log(path: Path) -> Optional[float]:
    """
    Parse our per-replicate .log to extract a 'Total' or 'Total energy' style key if present.
    Fallback regex for '... = <float>'.
    """
    try:
        txt = path.read_text(errors="ignore")
    except Exception:
        return None

    # Try direct 'Total energy' capture
    m = re.search(r"(?i)\btotal\s*energy\b\s*=\s*(-?\d+(?:\.\d+)?)", txt)
    if m:
        try:
            return float(m.group(1))
        except Exception:
            pass

    # Generic key=value scan; prefer a key that looks like total/energy
    candidates = []
    for key, val in re.findall(r"([A-Za-z0-9 _-]+?)\s*=\s*(-?\d+(?:\.\d+)?)", txt):
        k = key.strip().lower().replace("_", " ")
        try:
            v = float(val)
        except Exception:
            continue
        score = 0
        if "total" in k:
            score += 1
        if "energy" in k:
            score += 1
        candidates.append((score, v))
    if candidates:
        # best match (most total/energy tokens); if tie, take first
        candidates.sort(key=lambda t: (-t[0], abs(t[1])))
        return candidates[0][1]

    return None


def parse_foldx_file(path: Path) -> Optional[float]:
    """Dispatch parser by extension; prefer fxout semantics when available."""
    name = path.name.lower()
    if (
        name.endswith(".fxout")
        or name.endswith("_st.fxout")
        or name.endswith("_summary.fxout")
    ):
        v = parse_foldx_fxout(path)
        if v is not None:
            return v
    # fallback to log
    return parse_foldx_log(path)


# -----------------------------
# Per-row aggregation
# -----------------------------


def per_row_foldx_scores(
    row: pd.Series, score_dir: Path, nstruct: Optional[int]
) -> Dict[str, object]:
    """
    Produce:
      - FoldX_total_energy_0001 .. _NNNN (NaN if missing)
      - FoldX_best_total_energy
      - FoldX_best_pose_index
      - FoldX_best_file
      - FoldX_n_found
    """
    stems = candidate_stems_for_row(row)
    idx_files = list_replicate_files(score_dir, stems, nstruct)

    # If nstruct fixed, ensure columns 1..nstruct, else use discovered indices.
    if nstruct and nstruct > 0:
        indices = list(range(1, nstruct + 1))
        # map index -> file (if found)
        fmap = {i: p for (i, p) in idx_files}
    else:
        indices = [i for (i, _) in idx_files]
        fmap = {i: p for (i, p) in idx_files}

    out: Dict[str, object] = {}
    totals: Dict[int, float] = {}
    filemap: Dict[int, str] = {}

    for i in indices:
        p = fmap.get(i)
        val = np.nan
        if p is not None:
            t = parse_foldx_file(p)
            if t is not None:
                val = float(t)
                totals[i] = val
                filemap[i] = str(p)
        out[f"FoldX_total_energy_{i:04d}"] = val

    if totals:
        best_idx = min(totals, key=lambda k: totals[k])
        out["FoldX_best_total_energy"] = totals[best_idx]
        out["FoldX_best_pose_index"] = best_idx
        out["FoldX_best_file"] = filemap.get(best_idx)
    else:
        out["FoldX_best_total_energy"] = np.nan
        out["FoldX_best_pose_index"] = np.nan
        out["FoldX_best_file"] = None

    out["FoldX_n_found"] = int(len(totals))
    return out


# -----------------------------
# Main
# -----------------------------


def main():
    ap = argparse.ArgumentParser(
        description="Add per-replicate FoldX total energies and best to nodes_df"
    )
    ap.add_argument("--nodes", required=True)
    ap.add_argument(
        "--input_dir", required=True, help="Directory with FoldX outputs (.fxout/.log)"
    )
    ap.add_argument("--output_dir", default="processed-data")
    ap.add_argument("--organism_tag", required=True)
    ap.add_argument("--output_prefix", default="0")
    ap.add_argument("--output_suffix", required=True)
    ap.add_argument(
        "--nstruct", type=int, default=None, help="Expected # replicates per protein"
    )
    ap.add_argument("--nstruc", type=int, dest="nstruct")  # legacy alias
    args = ap.parse_args()

    nodes_df = pd.read_pickle(args.nodes)
    score_dir = Path(args.input_dir)

    # build per-row dicts, concat to nodes_df
    perrow = nodes_df.apply(
        lambda row: per_row_foldx_scores(row, score_dir, args.nstruct),
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
