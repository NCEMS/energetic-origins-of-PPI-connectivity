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
        and row.get("DeepTMHMM_class") in ["GLOB", "SP"]  # your newer gating
        and pd.notna(row.get("final_structure_path"))
        and str(row.get("final_structure_source")) != "None"
    )


def select_label(row: pd.Series) -> str:
    return str(row.get("node") or row.get("UniProtKB-AC") or "UNKNOWN")


def find_relaxed_pdb(row: pd.Series, rosetta_dir: Path) -> Optional[Path]:
    """
    Prefer the original (old-script) resolution: derive <stem> from final_structure_path
    and look for <stem>_0001.pdb in rosetta_dir. If that fails, fall back to the new
    label-based name (<label>_0001.pdb).
    """
    # 1) Old reliable method: from final_structure_path stem
    fsp = row.get("final_structure_path")
    if pd.notna(fsp):
        try:
            stem = Path(str(fsp)).stem
            candidate = rosetta_dir / f"{stem}_0001.pdb"
            if candidate.exists():
                return candidate
        except Exception:
            pass

    # 2) Fallback: new label mapping
    label = select_label(row)
    candidate = rosetta_dir / f"{label}_0001.pdb"
    return candidate if candidate.exists() else None


# -----------------------------
# Filesystem helpers
# -----------------------------

def ensure_symlink(src: Path, dst: Path, copy_instead: bool = False, overwrite: bool = False) -> Path:
    """
    Ensure `dst` points to `src` in dst.parent:
      - If copy_instead=False: create/update a symlink at dst → src
      - If copy_instead=True: copy src to dst
    If dst is an existing *regular file*:
      - overwrite=True: replace it with symlink/copy
      - overwrite=False: keep it
    """
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists() or dst.is_symlink():
        if overwrite:
            try:
                dst.unlink()
            except IsADirectoryError:
                raise
        else:
            # Keep existing file/symlink if not overwriting
            # But if it's a symlink to a different target, refresh it
            if dst.is_symlink():
                try:
                    if os.path.realpath(dst) == os.path.realpath(src):
                        return dst
                    dst.unlink()
                except OSError:
                    pass
            else:
                return dst  # keep existing regular file

    if copy_instead:
        shutil.copy2(src, dst)
        return dst

    rel = os.path.relpath(src, start=dst.parent)
    dst.symlink_to(rel)
    return dst


def foldx_outputs_exist(output_dir: Path, base: str) -> bool:
    """
    Recognize common FoldX Stability outputs across builds/configs.
    `base` is the *stem* passed via --output-file (we pass <label>_0001).

    Common patterns observed:
      <base>_ST.fxout                (classic Stability output)
      <base>.fxout                   (some builds when --output-file given)
      <base>_Summary.fxout           (summary variants)
      <base>.fxout_ST.fxout          (when user passes '.fxout' and FoldX appends '_ST.fxout')
    """
    candidates = [
        f"{base}_ST.fxout",
        f"{base}.fxout",
        f"{base}_Summary.fxout",
        f"{base}.fxout_ST.fxout",
    ]
    return any((output_dir / c).exists() for c in candidates)


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
    link_pdb = output_dir / f"{base}.pdb"  # input in output_dir

    # Ensure we have an input in output_dir
    if not src_pdb.exists():
        sys.stderr.write(f"[WARN] Source PDB missing: {src_pdb}\n")
        return str(src_pdb), False

    ensure_symlink(src_pdb, link_pdb, copy_instead=copy_instead, overwrite=overwrite)

    # Skip if already scored and not overwriting
    if not overwrite and foldx_outputs_exist(output_dir, base):
        print(f"[SKIP] {base}: FoldX outputs present.")
        return str(src_pdb), True

    # Build FoldX command; run with cwd=output_dir so outputs appear next to the input PDB
    # IMPORTANT: no trailing space in '--command=Stability', and pass a STEM to --output-file
    cmd = [
        foldx_exec,
        "--command=Stability",
        f"--pdb={link_pdb.name}",     # basename only; cwd=output_dir
        f"--output-file={base}",      # STEM; let FoldX append its suffix(es)
        # Optional:
        # "--water=CRYSTAL",
        # "--ionStrength=0.05",
        # "--pH=7",
        # "--vdwDesign=2",
        # "--rotabase=/path/to/rotabase.txt",
    ]

    print(f"[FoldX] Scoring: {base}")
    try:
        # Capture out/err and always write a log
        result = subprocess.run(cmd, check=False, text=True, capture_output=True, cwd=output_dir)
        log_path = output_dir / f"{base}.log"
        with open(log_path, "w") as fh:
            fh.write("=== CMD ===\n" + " ".join(cmd) + "\n\n")
            fh.write("=== STDOUT ===\n" + (result.stdout or "") + "\n\n")
            fh.write("=== STDERR ===\n" + (result.stderr or "") + "\n")

        if result.returncode != 0:
            sys.stderr.write(f"[ERROR] FoldX failed for {src_pdb}; see {log_path}\n")
            return str(src_pdb), False

    except Exception as e:
        sys.stderr.write(f"[ERROR] Exception running FoldX for {src_pdb}: {e}\n")
        return str(src_pdb), False

    ok = foldx_outputs_exist(output_dir, base)
    if not ok:
        sys.stderr.write(f"[ERROR] FoldX produced no recognized outputs for {base}\n")
    return str(src_pdb), ok


# -----------------------------
# Main
# -----------------------------

def main():
    p = argparse.ArgumentParser(description="Run FoldX Score on Rosetta-relaxed PDBs (<stem>_0001.pdb).")
    p.add_argument("--nodes", required=True, help="Path to nodes DataFrame (.pkl/.pickle or .csv)")
    p.add_argument("--rosetta_dir", required=True, help="Directory containing <stem}_0001.pdb from Rosetta")
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

        src_pdb = find_relaxed_pdb(row, rosetta_dir)
        if src_pdb is None:
            missing_inputs += 1
            continue

        # Determine the label used for output naming (align with our base=<label>_0001)
        # If the file came from final_structure_path stem, keep that stem for naming;
        # else use select_label(row).
        fsp = row.get("final_structure_path")
        if pd.notna(fsp) and Path(str(fsp)).stem in src_pdb.name:
            label = Path(str(fsp)).stem
        else:
            label = select_label(row)

        jobs.append((
            str(src_pdb),
            label,
            args.foldx_executable,
            str(output_dir),
            args.copy_instead_of_symlink,
            args.overwrite
        ))

    if missing_inputs:
        print(f"[WARN] {missing_inputs} selected targets had no Rosetta relaxed PDB in {rosetta_dir}")

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
        for (pdb_path, label, _, outdir, _, _) in jobs:
            print(f" - {label}_0001  (src={pdb_path}, cwd={outdir})")
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
