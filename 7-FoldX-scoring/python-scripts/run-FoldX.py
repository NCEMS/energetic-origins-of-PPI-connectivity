#!/usr/bin/env python3
"""
Run FoldX Score on Rosetta-relaxed PDBs that follow the <label>_0001.pdb naming scheme.

- Loads a nodes DataFrame (pickle or CSV)
- Selects targets with the same gating used elsewhere:
    has_verified_sequence == True
    DeepTMHMM_class in ["GLOB", "SP"]
    final_structure_path is not NA
    final_structure_source != "None"
- For each selected target, expects a Rosetta-relaxed PDB named: <label>_0001.pdb
  where label defaults to `row["node"]`, else `row["UniProtKB-AC"]`, else "UNKNOWN".
- Creates a symlink in --output_dir named <label>_0001.pdb pointing to the Rosetta PDB
- Runs FoldX Score in cwd=--output_dir so all outputs appear next to that symlink
- Skips targets if <label>_0001.fxout or <label>_0001_Summary.fxout already exists
"""

from __future__ import annotations

import argparse
import os
import sys
import shutil
import subprocess
import multiprocessing as mp
from pathlib import Path
from typing import Optional, Tuple, List

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


def select_label(row: pd.Series) -> str:
    return str(row.get("node") or row.get("UniProtKB-AC") or "UNKNOWN")


def find_rosetta_pdb_for_label(label: str, rosetta_dir: Path) -> Path:
    """
    Expect Rosetta-relaxed file named <label>_0001.pdb in rosetta_dir.
    """
    return rosetta_dir / f"{label}_0001.pdb"


# -----------------------------
# Filesystem helpers
# -----------------------------

def ensure_symlink(src: Path, dst: Path, copy_instead: bool = False) -> Path:
    """
    Ensure `dst` points to `src`:
      - If copy_instead=False (default): create/update a symlink at dst → src
      - If copy_instead=True: copy src to dst (overwrites an existing symlink to elsewhere)
    If dst is an existing *regular file*, leave it as-is.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.is_symlink():
        # If correct already, do nothing; else replace the symlink.
        if os.path.realpath(dst) == os.path.realpath(src):
            return dst
        dst.unlink()

    if dst.exists() and not dst.is_symlink():
        # Regular file present; keep it.
        return dst

    if copy_instead:
        shutil.copy2(src, dst)
        return dst

    # Create a relative symlink for nicer paths
    rel = os.path.relpath(src, start=dst.parent)
    dst.symlink_to(rel)
    return dst


def foldx_outputs_exist(output_dir: Path, base: str) -> bool:
    """
    Treat presence of either file as 'already scored'.
    Different FoldX builds may produce <base>.fxout OR <base>_Summary.fxout.
    """
    return ((output_dir / f"{base}.fxout").exists() or
            (output_dir / f"{base}_Summary.fxout").exists())


# -----------------------------
# FoldX worker
# -----------------------------

def score_with_foldx(args: Tuple[str, str, str, str, bool, bool]) -> Tuple[str, bool]:
    """
    Worker function to run FoldX Score on a single PDB.

    args:
        pdb_path_str: source Rosetta PDB (absolute or relative)
        label:        label used as stem (<label>_0001)
        foldx_exec:   path to FoldX executable
        output_dir:   directory where symlink and outputs live
        copy_instead: copy PDB instead of symlink (if fs doesn't support symlinks)
        overwrite:    re-run even if outputs already exist
    """
    pdb_path_str, label, foldx_exec, output_dir, copy_instead, overwrite = args
    src_pdb = Path(pdb_path_str)
    output_dir = Path(output_dir)

    base = f"{label}_0001"
    link_pdb = output_dir / f"{base}.pdb"  # target file in output_dir (symlink or copy)

    # Ensure we have an input in output_dir
    if not src_pdb.exists():
        sys.stderr.write(f"[WARN] Source PDB missing: {src_pdb}\n")
        return str(src_pdb), False

    ensure_symlink(src_pdb, link_pdb, copy_instead=copy_instead)

    # Skip if already scored and not overwriting
    if not overwrite and foldx_outputs_exist(output_dir, base):
        print(f"[SKIP] {base}: FoldX outputs present.")
        return str(src_pdb), True

    # Build FoldX command; run with cwd=output_dir so outputs appear next to the input PDB
    cmd = [
        foldx_exec,
        "--command=Score",
        f"--pdb={link_pdb.name}",  # just basename; cwd=output_dir
        f"--output-file={base}",
        # Optional extras if you use them:
        # "--water=CRYSTAL",
        # "--ionStrength=0.05",
        # "--pH=7",
        # "--vdwDesign=2",
        # "--rotabase=/path/to/rotabase.txt",
    ]

    print(f"[FoldX] Scoring: {base}")
    try:
        # capture_output=True keeps stdout/stderr; switch to False if you want real-time logs
        subprocess.run(cmd, check=True, text=True, capture_output=True, cwd=output_dir)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"[ERROR] FoldX failed for {src_pdb}:\n{e.stderr}\n")
        return str(src_pdb), False

    ok = foldx_outputs_exist(output_dir, base)
    if not ok:
        sys.stderr.write(f"[ERROR] FoldX produced no expected outputs for {base}\n")
    return str(src_pdb), ok


# -----------------------------
# Main
# -----------------------------

def main():
    p = argparse.ArgumentParser(description="Run FoldX Score on Rosetta-relaxed PDBs (<label>_0001.pdb).")
    p.add_argument("--nodes", required=True, help="Path to nodes DataFrame (.pkl/.pickle or .csv)")
    p.add_argument("--rosetta_dir", required=True, help="Directory containing <label>_0001.pdb from Rosetta")
    p.add_argument("--output_dir", required=True, help="Directory where FoldX outputs will be written")
    p.add_argument("--foldx_executable", default="foldx", help="FoldX binary (default: foldx in PATH)")
    p.add_argument("--nprocessors", type=int, default=0, help="Parallel workers (0 = cpu_count)")
    p.add_argument("--copy-instead-of-symlink", action="store_true",
                   help="Copy PDBs into output_dir instead of creating symlinks")
    p.add_argument("--overwrite", action="store_true",
                   help="Re-run even if FoldX outputs exist for a target")
    p.add_argument("--dry-run", action="store_true", help="Plan only; do not run FoldX")
    args = p.parse_args()

    # Load DataFrame
    if args.nodes.endswith((".pkl", ".pickle")):
        nodes = pd.read_pickle(args.nodes)
    else:
        nodes = pd.read_csv(args.nodes)

    rosetta_dir = Path(args.rosetta_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Optional: fast-fail on missing FoldX
    if shutil.which(args.foldx_executable) is None and not Path(args.foldx_executable).exists():
        sys.stderr.write(f"[ERROR] FoldX executable not found: {args.foldx_executable}\n")
        sys.exit(1)

    # Build job list from DataFrame
    jobs: List[Tuple[str, str, str, str, bool, bool]] = []
    missing_inputs = 0

    for _, row in nodes.iterrows():
        if not should_include(row):
            continue
        label = select_label(row)
        src_pdb = find_rosetta_pdb_for_label(label, rosetta_dir)
        if not src_pdb.exists():
            missing_inputs += 1
            continue

        # We'll symlink/copy to output_dir as <label>_0001.pdb
        jobs.append((
            str(src_pdb),
            label,
            args.foldx_executable,
            str(output_dir),
            args.copy_instead_of_symlink,
            args.overwrite
        ))

    if missing_inputs:
        print(f"[WARN] {missing_inputs} selected targets had no Rosetta PDB in {rosetta_dir}")

    if not jobs:
        print("[INFO] No jobs to run (nothing selected or inputs missing).")
        return

    # Prefilter jobs unless overwriting
    if not args.overwrite:
        pre = []
        for (pdb_path, label, fx, outdir, copyflag, overwrite) in jobs:
            base = f"{label}_0001"
            if not foldx_outputs_exist(Path(outdir), base):
                pre.append((pdb_path, label, fx, outdir, copyflag, overwrite))
        jobs = pre

    if not jobs:
        print("[INFO] Nothing to do; all FoldX outputs already present.")
        (output_dir / ".all_foldx_done").touch()
        return

    if args.dry_run:
        print("[DRY-RUN] Planned FoldX Score jobs:")
        for (_, label, _, outdir, _, _) in jobs:
            print(f" - {label}_0001  (cwd={outdir})")
        return

    # Parallel execution
    n = args.nprocessors or mp.cpu_count()
    with mp.Pool(processes=n) as pool:
        results = pool.map(score_with_foldx, jobs)

    ok = [p for (p, s) in results if s]
    bad = [p for (p, s) in results if not s]

    print(f"\n[SUMMARY] FoldX Score complete. OK: {len(ok)}  Failed: {len(bad)}")
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
