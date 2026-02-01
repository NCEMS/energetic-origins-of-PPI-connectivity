#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd


def populate_structures(
    nodes_pkl: Path,
    out_dir: Path,
    *,
    node_col: str = "node",
    path_col: str = "final_structure_path",
    mode: str = "copy",          # "copy" | "symlink" | "hardlink"
    overwrite: bool = False,
) -> dict:
    """
    Populate out_dir with PDBs renamed to {node}.pdb using paths in nodes_pkl[path_col].

    Returns a small summary dict.
    """
    nodes_pkl = nodes_pkl.resolve()
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_pickle(nodes_pkl)

    if node_col not in df.columns:
        raise KeyError(f"Missing column '{node_col}' in {nodes_pkl}. Columns: {list(df.columns)}")
    if path_col not in df.columns:
        raise KeyError(f"Missing column '{path_col}' in {nodes_pkl}. Columns: {list(df.columns)}")

    # Normalize node IDs to strings
    df[node_col] = df[node_col].astype(str)

    # Track duplicates (same node appears multiple times)
    dup_nodes = df.loc[df[node_col].duplicated(), node_col].unique().tolist()
    if dup_nodes:
        # This is usually a data issue; you can choose to raise instead.
        print(f"WARNING: {len(dup_nodes)} duplicate node IDs detected. "
              f"Later rows may overwrite earlier ones if overwrite=True.")

    n_total = len(df)
    n_missing_path = 0
    n_missing_file = 0
    n_written = 0
    n_skipped_existing = 0

    for _, row in df.iterrows():
        node_id = str(row[node_col]).strip()
        raw_path = row[path_col]

        if not node_id:
            continue

        if raw_path is None or (isinstance(raw_path, float) and pd.isna(raw_path)) or str(raw_path).strip() == "":
            n_missing_path += 1
            continue

        raw = str(raw_path).strip()

        # Special-case: paths stored relative to repo root starting with "0-download-inputs"
        # Make them work from the script/notebook working directory by prefixing "../../"
        if raw.startswith("0-download-inputs"):
            raw = "../../" + raw

        src = Path(raw).expanduser()

        # If still relative, interpret relative to the nodes.pkl location (usually ./data)
        if not src.is_absolute():
            src = (nodes_pkl.parent / src).resolve()
        else:
            src = src.resolve()

        if not src.exists():
            n_missing_file += 1
            continue

        dst = out_dir / f"{node_id}.pdb"

        if dst.exists():
            if overwrite:
                dst.unlink()
            else:
                n_skipped_existing += 1
                continue

        if mode == "copy":
            shutil.copy2(src, dst)
        elif mode == "symlink":
            dst.symlink_to(src)
        elif mode == "hardlink":
            # Hardlink only works when src/dst are on same filesystem
            dst.hardlink_to(src)
        else:
            raise ValueError("mode must be one of: copy, symlink, hardlink")

        n_written += 1

    return {
        "nodes_pkl": str(nodes_pkl),
        "out_dir": str(out_dir),
        "mode": mode,
        "overwrite": overwrite,
        "n_total_rows": n_total,
        "n_written": n_written,
        "n_skipped_existing": n_skipped_existing,
        "n_missing_path": n_missing_path,
        "n_missing_file": n_missing_file,
        "n_duplicate_nodes": len(dup_nodes),
    }


def main() -> None:
    p = argparse.ArgumentParser(
        description="Populate a structures directory from nodes.pkl final_structure_path, renaming to {node}.pdb"
    )
    p.add_argument("--nodes-pkl", default="data/nodes.pkl", help="Path to nodes.pkl")
    p.add_argument("--out-dir", default="structures", help="Output structures directory")
    p.add_argument("--node-col", default="node", help="Column with node IDs")
    p.add_argument("--path-col", default="final_structure_path", help="Column with source PDB paths")
    p.add_argument("--mode", choices=["copy", "symlink", "hardlink"], default="copy",
                   help="How to populate: copy (recommended), symlink, or hardlink")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing {node}.pdb files")

    args = p.parse_args()

    summary = populate_structures(
        Path(args.nodes_pkl),
        Path(args.out_dir),
        node_col=args.node_col,
        path_col=args.path_col,
        mode=args.mode,
        overwrite=args.overwrite,
    )

    print("\n=== Populate structures summary ===")
    for k, v in summary.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
