#!/usr/bin/env python3
import os, sys
import numpy as np
import pandas as pd
import typing
import multiprocessing as mp
from pathlib import Path
import argparse
import subprocess
import tempfile
import shutil

# ---------------------------
# Selection
# ---------------------------

def select_structure(row):
    """
    Select the structure path decided upstream and attach a unique label for outputs.
    We mirror the gating logic and use `node` as the unique base name.
    """
    if (
        row.get("has_verified_sequence") is True
        and row.get("DeepTMHMM_class") in ["GLOB", "SP"]
        and pd.notna(row.get("final_structure_path"))
        and str(row.get("final_structure_source")) != "None"
    ):
        path = row["final_structure_path"]
        label = str(row.get("node") or row.get("UniProtKB-AC") or Path(path).parent.name)
        return (path, label)
    return None

# ---------------------------
# Helpers
# ---------------------------

def expected_outputs_exist(output_dir: Path, label: str, nstruct: int) -> bool:
    """Return True if ALL expected pose files exist for this label."""
    for i in range(1, nstruct + 1):
        stem = f"{label}_{i:04d}"
        if not ( (output_dir / f"{stem}.pdb").exists() and (output_dir / f"{stem}.sc").exists() ):
            return False
    return True

# ---------------------------
# Worker
# ---------------------------

def relax_pdb(args):
    """
    Run Rosetta FastRelax for N poses serially for a single input PDB.

    Args tuple:
      pdb_path_str, label, rosetta_exec, output_dir, nstruct
    """
    pdb_path_str, label, rosetta_exec, output_dir, nstruct = args
    pdb_path   = Path(pdb_path_str)
    output_dir = Path(output_dir)

    # Quick skip: if all expected outputs exist, do nothing
    if expected_outputs_exist(output_dir, label, nstruct):
        print(f"Skipping {pdb_path.name} (already have {label}_0001..{label}_{nstruct:04d}.*)")
        return str(pdb_path), True

    # One scratch root per dataset; per-pose subdirs to keep runs clean/isolated
    scratch_root = output_dir / ".scratch"
    scratch_root.mkdir(parents=True, exist_ok=True)

    all_ok = True
    for i in range(1, nstruct + 1):
        final_sc  = output_dir / f"{label}_{i:04d}.sc"
        final_pdb = output_dir / f"{label}_{i:04d}.pdb"

        # Skip this pose if both outputs already exist
        if final_sc.exists() and final_pdb.exists():
            print(f"  Pose {i:04d}: already present; skipping")
            continue

        pose_scratch = Path(
            tempfile.mkdtemp(prefix=f"rosetta_{label}_{i:04d}_", dir=str(scratch_root))
        )

        # Build Rosetta command (nstruct=1 per serial replicate)
        relax_cmd = [
            rosetta_exec,
            "-s", str(pdb_path),
            "-relax:fast",
            "-relax:constrain_relax_to_start_coords",
            "-nstruct", "1",
            "-score:weights", "ref2015",
            "-out:path:all", str(pose_scratch),
            # Explicit scorefile path; Rosetta may still append extra info, so we glob below
            "-out:file:scorefile", str(pose_scratch / f"{label}_{i:04d}.sc"),
            "-out:pdb", "true",
        ]

        print(f"[{label}] Pose {i:04d}: running FastRelax in scratch")
        try:
            res = subprocess.run(relax_cmd, check=False, capture_output=True, text=True)
        except Exception as e:
            print(f"[{label}] Pose {i:04d}: exception while running Rosetta: {e}")
            all_ok = False
            # Clean and continue to next pose
            shutil.rmtree(pose_scratch, ignore_errors=True)
            continue

        # Always write a small log for debugging this pose
        try:
            with open(output_dir / f"{label}_{i:04d}.log", "w") as fh:
                fh.write("=== CMD ===\n" + " ".join(relax_cmd) + "\n\n")
                fh.write("=== STDOUT ===\n" + (res.stdout or "") + "\n\n")
                fh.write("=== STDERR ===\n" + (res.stderr or "") + "\n")
        except Exception:
            pass

        if res.returncode != 0:
            print(f"[{label}] Pose {i:04d}: Rosetta returned nonzero; see log.")
            all_ok = False
            shutil.rmtree(pose_scratch, ignore_errors=True)
            continue

        # Find produced files (be lenient about exact names)
        produced_pdbs = sorted(pose_scratch.glob("*.pdb"))
        produced_scs  = sorted(pose_scratch.glob("*.sc"))

        if not produced_pdbs or not produced_scs:
            print(f"[{label}] Pose {i:04d}: missing PDB or scorefile in scratch.")
            all_ok = False
            shutil.rmtree(pose_scratch, ignore_errors=True)
            continue

        # Choose the first of each (nstruct=1 → singletons expected)
        src_pdb = produced_pdbs[0]
        src_sc  = produced_scs[0]

        # Move to final labeled names (atomic on same FS)
        try:
            os.replace(src_pdb, final_pdb)
            os.replace(src_sc,  final_sc)
        except Exception as e:
            print(f"[{label}] Pose {i:04d}: failed to move outputs: {e}")
            all_ok = False
        finally:
            shutil.rmtree(pose_scratch, ignore_errors=True)

        ok = final_pdb.exists() and final_sc.exists()
        print(f"[{label}] Pose {i:04d}: {'OK' if ok else 'FAILED'}")
        all_ok = all_ok and ok

    return str(pdb_path), all_ok

# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser(description="Run Rosetta FastRelax with N poses per structure.")
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--output_dir", default="processed-data",
                        help="Directory where results will be saved (default: processed-data)")
    parser.add_argument("--nprocessors", type=int)
    parser.add_argument("--organism_tag")
    parser.add_argument("--relax_executable")
    parser.add_argument("--output_prefix", default="0", help="Prefix for output files")
    parser.add_argument("--nstruct", type=int, default=1,
                        help="Number of poses (replicates) to generate per input structure (default: 1)")
    args = parser.parse_args()

    rosetta_exec = args.relax_executable
    if not rosetta_exec or not Path(rosetta_exec).is_file():
        raise FileNotFoundError(f"Rosetta executable not found: {rosetta_exec}")

    nodes_df = pd.read_pickle(args.nodes)

    # just use a single small protein for testing purposes
    # nodes_df = nodes_df[nodes_df["node"] == "YOR167C"]

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build job list from final_structure_path
    selected = nodes_df.apply(select_structure, axis=1).tolist()
    jobs = []
    for item in selected:
        if item is None:
            continue
        path, label = item
        if isinstance(path, str) and path and os.path.isfile(path):
            jobs.append((path, label, rosetta_exec, str(output_dir), int(args.nstruct)))

    if not jobs:
        print("No structures selected (verify final_structure_path/source filters).")
        return

    # Prefilter: skip items that already have ALL expected replicates
    jobs_to_run = []
    for (path, label, rexec, outdir, nstruct) in jobs:
        if not expected_outputs_exist(Path(outdir), label, nstruct):
            jobs_to_run.append((path, label, rexec, outdir, nstruct))

    if not jobs_to_run:
        print("Nothing to do; all replicate files present.")
        (output_dir / ".all_scores_done").touch()
        return

    ncpus = args.nprocessors or mp.cpu_count()
    with mp.Pool(processes=ncpus) as pool:
        results = pool.map(relax_pdb, jobs_to_run)

    successful_paths = [r[0] for r in results if r[1] is True]
    failed_paths     = [r[0] for r in results if r[1] is False]

    # Final report
    print(f"\nSelected structures: {len(jobs)}")
    print(f"Ran (needed work):  {len(jobs_to_run)}")
    print(f"OK:                 {len(successful_paths)}")
    print(f"Failed:             {len(failed_paths)}")

    if failed_paths:
        with open(output_dir / "missing_structures.txt", "w") as f:
            for p in failed_paths:
                f.write(f"{p}\n")
        print(f"Details written to {output_dir/'missing_structures.txt'}")
    else:
        (output_dir / ".all_scores_done").touch()

if __name__ == "__main__":
    main()
