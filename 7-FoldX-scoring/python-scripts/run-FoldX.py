#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
import shutil
import subprocess
import multiprocessing as mp
from pathlib import Path
from typing import Optional, Tuple, List, Dict

import pandas as pd

# -----------------------------
# Selection & naming utilities
# -----------------------------

def should_include(row: pd.Series) -> bool:
    return (
        row.get("has_verified_sequence") is True
        and row.get("DeepTMHMM_class") in ["GLOB", "SP"]
        and pd.notna(row.get("final_structure_path"))
        and str(row.get("final_structure_source")) != "None"
    )

def candidate_stems_for_row(row: pd.Series) -> List[str]:
    """
    Prefer the new label (row['node']) used by your Rosetta step, but
    also fall back to stems from structure paths for backward compat.
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
            out.append(s); seen.add(s)
    return out

# -----------------------------
# Filesystem helpers
# -----------------------------

def ensure_symlink(src: Path, dst: Path, copy_instead: bool = False, overwrite: bool = False) -> Path:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if overwrite:
            try:
                dst.unlink()
            except IsADirectoryError:
                raise
        else:
            if dst.is_symlink():
                try:
                    if os.path.realpath(dst) == os.path.realpath(src):
                        return dst
                    dst.unlink()
                except OSError:
                    pass
            else:
                return dst
    if copy_instead:
        shutil.copy2(src, dst)
        return dst
    rel = os.path.relpath(src, start=dst.parent)
    dst.symlink_to(rel)
    return dst

def foldx_outputs_exist(output_dir: Path, base: str) -> bool:
    candidates = [
        f"{base}_ST.fxout",
        f"{base}.fxout",
        f"{base}_Summary.fxout",
        f"{base}.fxout_ST.fxout",
    ]
    return any((output_dir / c).exists() for c in candidates)

# -----------------------------
# Rosetta replicate discovery
# -----------------------------

def list_relaxed_pdbs_for_row(row: pd.Series, rosetta_dir: Path, nstruct: Optional[int]) -> List[Tuple[Path, str]]:
    """
    Return list of (pdb_path, base) where base is STEM_#### without extension.
    Try stems in order; stop at the first that yields any files.
    If nstruct is provided, only those indices (1..nstruct) are considered.
    """
    stems = candidate_stems_for_row(row)
    for stem in stems:
        found: List[Tuple[Path, str]] = []
        if nstruct and nstruct > 0:
            for i in range(1, nstruct + 1):
                base = f"{stem}_{i:04d}"
                p = rosetta_dir / f"{base}.pdb"
                if p.exists():
                    found.append((p, base))
        else:
            for p in sorted(rosetta_dir.glob(f"{stem}_*.pdb")):
                name = p.stem
                if "_" not in name:
                    continue
                suf = name.split("_")[-1]
                if len(suf) == 4 and suf.isdigit():
                    found.append((p, name))
        if found:
            return found
    return []

# -----------------------------
# FoldX worker
# -----------------------------

def score_with_foldx(args: Tuple[str, str, str, str, bool, bool]) -> Tuple[str, bool]:
    """
    Worker to run FoldX Stability on a single PDB replicate.

    args:
        pdb_path_str: path to input PDB
        base:         STEM used for outputs (e.g., LABEL_0007)
        foldx_exec:   path/binary
        output_dir:   where outputs/logs go
        copy_instead: copy PDB instead of symlink
        overwrite:    re-run even if outputs exist
    """
    pdb_path_str, base, foldx_exec, output_dir, copy_instead, overwrite = args
    src_pdb = Path(pdb_path_str)
    output_dir = Path(output_dir)

    if not src_pdb.exists():
        sys.stderr.write(f"[WARN] Missing PDB: {src_pdb}\n")
        return str(src_pdb), False

    link_pdb = output_dir / f"{base}.pdb"
    ensure_symlink(src_pdb, link_pdb, copy_instead=copy_instead, overwrite=overwrite)

    if not overwrite and foldx_outputs_exist(output_dir, base):
        print(f"[SKIP] {base}: outputs present.")
        return str(src_pdb), True

    cmd = [
        foldx_exec,
        "--command=Stability",
        f"--pdb={link_pdb.name}",
        f"--output-file={base}",
    ]

    print(f"[FoldX] {base}")
    try:
        result = subprocess.run(cmd, check=False, text=True, capture_output=True, cwd=output_dir)
        log_path = output_dir / f"{base}.log"
        with open(log_path, "w") as fh:
            fh.write("=== CMD ===\n" + " ".join(cmd) + "\n\n")
            fh.write("=== STDOUT ===\n" + (result.stdout or "") + "\n\n")
            fh.write("=== STDERR ===\n" + (result.stderr or "") + "\n")

        if result.returncode != 0:
            sys.stderr.write(f"[ERROR] FoldX failed for {base}; see {log_path}\n")
            return str(src_pdb), False

    except Exception as e:
        sys.stderr.write(f"[ERROR] Exception running FoldX for {base}: {e}\n")
        return str(src_pdb), False

    ok = foldx_outputs_exist(output_dir, base)
    if not ok:
        sys.stderr.write(f"[ERROR] No recognized FoldX outputs for {base}\n")
    return str(src_pdb), ok

# -----------------------------
# Main
# -----------------------------

def main():
    p = argparse.ArgumentParser(description="Run FoldX Stability on Rosetta-relaxed PDB replicates (STEM_####.pdb).")
    p.add_argument("--nodes", required=True, help="Path to nodes DataFrame (.pkl/.pickle or .csv)")
    p.add_argument("--rosetta_dir", required=True, help="Directory containing Rosetta <stem>_####.pdb files")
    p.add_argument("--output_dir", required=True, help="Directory where FoldX outputs will be written")
    p.add_argument("--foldx_executable", default="foldx", help="FoldX binary (default: foldx in PATH)")
    p.add_argument("--nprocessors", type=int, default=0, help="Parallel workers (0 = cpu_count)")
    p.add_argument("--copy-instead-of-symlink", action="store_true",
                   help="Copy PDBs into output_dir instead of creating symlinks")
    p.add_argument("--overwrite", action="store_true",
                   help="Re-run even if outputs already exist for a target")
    p.add_argument("--dry-run", action="store_true", help="Plan only; do not run FoldX")
    p.add_argument("--nstruct", type=int, default=None, help="Number of Rosetta replicates per protein")

    args = p.parse_args()

    # Load DataFrame
    if str(args.nodes).endswith((".pkl", ".pickle")):
        nodes = pd.read_pickle(args.nodes)
    else:
        nodes = pd.read_csv(args.nodes)

    rosetta_dir = Path(args.rosetta_dir).resolve()
    output_dir  = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if shutil.which(args.foldx_executable) is None and not Path(args.foldx_executable).exists():
        sys.stderr.write(f"[ERROR] FoldX executable not found: {args.foldx_executable}\n")
        sys.exit(1)

    # Build (row, replicate) jobs
    jobs: List[Tuple[str, str, str, str, bool, bool]] = []
    selected_rows = 0
    rows_without_inputs = 0

    for _, row in nodes.iterrows():
        if not should_include(row):
            continue
        selected_rows += 1

        pairs = list_relaxed_pdbs_for_row(row, rosetta_dir, args.nstruct)
        if not pairs:
            rows_without_inputs += 1
            continue

        for pdb_path, base in pairs:
            jobs.append((
                str(pdb_path),
                base,  # STEM_#### (also used for output-file)
                args.foldx_executable,
                str(output_dir),
                args.copy_instead_of_symlink,
                args.overwrite
            ))

    if selected_rows == 0:
        print("[INFO] No rows selected by gating; nothing to do.")
        return

    if rows_without_inputs:
        print(f"[WARN] {rows_without_inputs} selected rows had no matching Rosetta replicates in {rosetta_dir}")

    # If not overwriting, drop jobs whose outputs already exist
    if not args.overwrite:
        pre = []
        for (pdb_path, base, fx, outdir, copyflag, overwrite) in jobs:
            if not foldx_outputs_exist(Path(outdir), base):
                pre.append((pdb_path, base, fx, outdir, copyflag, overwrite))
        jobs = pre

    if not jobs:
        print("[INFO] Nothing to run (all outputs present).")
        (output_dir / ".all_foldx_done").touch()
        return

    if args.dry_run:
        print("[DRY-RUN] Planned FoldX jobs:")
        for _, base, _, outdir, _, _ in jobs:
            print(f" - {base}  (cwd={outdir})")
        return

    # Parallel execution (replicate-level)
    n = args.nprocessors or mp.cpu_count()
    with mp.Pool(processes=n) as pool:
        results = pool.map(score_with_foldx, jobs)

    ok = [p for (p, s) in results if s]
    bad = [p for (p, s) in results if not s]

    print(f"\n[SUMMARY] FoldX complete. OK replicates: {len(ok)}  Failed replicates: {len(bad)}")
    if bad:
        print("Failed PDBs:")
        for p in bad:
            print(" -", p)
        with open(output_dir / "foldx_failed.txt", "w") as fh:
            for p in bad:
                fh.write(f"{p}\n")
    else:
        (output_dir / ".all_foldx_done").touch()

if __name__ == "__main__":
    main()
