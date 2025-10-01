import os, sys
import numpy as np
import pandas as pd
import pint
import pint_pandas
import typing
import multiprocessing as mp
from pathlib import Path
import argparse
import subprocess
import tempfile
import shutil

def select_structure(row):
    """
    Select the structure path decided upstream and attach a unique label for outputs.
    We mirror the cagiada gating logic and use `node` as the unique base name.
    """
    if (
        row.get("has_verified_sequence") is True
        and row.get("DeepTMHMM_class") in ["GLOB", "SP"]
        and pd.notna(row.get("final_structure_path"))
        and str(row.get("final_structure_source")) != "None"
    ):
        path = row["final_structure_path"]
        # Use a stable, unique label per protein:
        label = str(row.get("node") or row.get("UniProtKB-AC") or Path(path).parent.name)
        return (path, label)
    return None


def score_pdb(args):
    pdb_path_str, rosetta_exec, output_dir = args
    pdb_path = Path(pdb_path_str)
    output_scorefile = output_dir / (pdb_path.stem + ".sc")

    cmd = [
        rosetta_exec,
        "-in:file:s",
        str(pdb_path),
        "-score:weights",
        "ref2015",
        "-out:file:scorefile",
        str(output_scorefile),
        "-out:pdb",
        "false",
    ]

    print(f"Scoring {pdb_path.name} ...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\nRosetta failed on {pdb_path.name}:\n{e.stderr}")
        return str(pdb_path), False
    if not output_scorefile.exists():
        print(
            f"\nNo score file produced for {pdb_path.name} (no crash, but silent failure?)"
        )
        return str(pdb_path), False
    return str(pdb_path), True


def relax_pdb(args):
    pdb_path_str, label, rosetta_exec, output_dir = args
    pdb_path = Path(pdb_path_str)

    base = label  # unique per-target stem (e.g., node or UniProtKB-AC)
    final_sc  = output_dir / f"{base}_0001.sc"
    final_pdb = output_dir / f"{base}_0001.pdb"

    # If both final outputs exist, we’re done
    if final_sc.exists() and final_pdb.exists():
        print(f"Skipping {pdb_path.name} (already have {base}_0001.*)")
        return str(pdb_path), True

    # Make a per-job scratch dir inside output_dir/.scratch
    scratch_root = output_dir / ".scratch"
    scratch_root.mkdir(parents=True, exist_ok=True)
    scratch_dir = Path(tempfile.mkdtemp(prefix=f"rosetta_{base}_", dir=str(scratch_root)))

    try:
        # Use explicit names inside scratch to be extra clear
        scratch_sc  = scratch_dir / f"{base}_0001.sc"
        scratch_pdb = scratch_dir / f"{base}_0001.pdb"

        # Build Rosetta command to write into scratch
        relax_cmd = [
            rosetta_exec,
            "-s", str(pdb_path),
            "-relax:fast",
            "-relax:constrain_relax_to_start_coords",
            "-nstruct", "1",
            "-score:weights", "ref2015",
            "-out:path:all", str(scratch_dir),
            "-out:file:scorefile", str(scratch_sc),
            "-out:pdb", "true",
            "-out:file:o", str(scratch_pdb),  # some builds may ignore, but we’re isolated anyway
        ]

        print(f"Running FastRelax on {pdb_path.name} -> {base}_0001.* (in scratch)")
        try:
            subprocess.run(relax_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"FastRelax failed on {pdb_path.name}:\n{e.stderr}")
            return str(pdb_path), False

        # Locate outputs (handle Rosetta ignoring -out:file:o)
        # Expect exactly one .pdb and one .sc in scratch (nstruct=1)
        produced_pdbs = list(scratch_dir.glob("*.pdb"))
        produced_scs  = list(scratch_dir.glob("*.sc"))

        if not produced_pdbs:
            print(f"Missing PDB in scratch for {pdb_path.name}")
            return str(pdb_path), False
        if not produced_scs:
            print(f"Missing scorefile in scratch for {pdb_path.name}")
            return str(pdb_path), False

        # Pick the (only) files – or the first if multiple, which shouldn’t happen with nstruct=1
        src_pdb = produced_pdbs[0]
        src_sc  = produced_scs[0]

        # Atomically move into final location with the label-based stem
        # (os.replace is atomic on the same filesystem)
        os.replace(src_pdb, final_pdb)
        os.replace(src_sc,  final_sc)

        print(f"FastRelax completed for {pdb_path.name} -> {base}_0001.*")
        return str(pdb_path), (final_pdb.exists() and final_sc.exists())

    finally:
        # Clean up scratch dir
        try:
            shutil.rmtree(scratch_dir, ignore_errors=True)
        except Exception:
            pass


def main():

    # parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Run Rosetta energy scoring of protein structures"
    )
    parser.add_argument("--nodes", required=True)
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where results will be saved (default: 'outputs')",
    )
    parser.add_argument("--nprocessors", type=int)
    parser.add_argument("--organism_tag")
    parser.add_argument("--relax_executable")
    parser.add_argument("--output_prefix", default="0", help="Prefix for output files")
    args = parser.parse_args()

    rosetta_exec = args.relax_executable
    if not rosetta_exec or not Path(rosetta_exec).is_file():
        raise FileNotFoundError(f"Rosetta executable not found: {rosetta_exec}")

    # read in the nodes_df from file
    nodes_df = pd.read_pickle(args.nodes)

    # setup output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # create list of jobs (path + unique label), from final_structure_path only
    selected = nodes_df.apply(select_structure, axis=1).tolist()
    jobs = []
    for item in selected:
        if item is None:
            continue
        path, label = item
        if isinstance(path, str) and path and os.path.isfile(path):
            jobs.append((path, label, rosetta_exec, output_dir))

    if not jobs:
        print("No structures selected (verify final_structure_path/source filters).")
        return

    # default CPU count if not provided
    ncpus = args.nprocessors or mp.cpu_count()

    # (optional) prefilter: skip already-scored items to avoid launching trivial workers
    jobs_to_run = []
    for (path, label, rexec, outdir) in jobs:
        sf = outdir / f"{label}_0001.sc"
        pf = outdir / f"{label}_0001.pdb"
        if not (sf.exists() and pf.exists()):
            jobs_to_run.append((path, label, rexec, outdir))

    if not jobs_to_run:
        print("Nothing to do; all score files present.")
        (Path(args.output_dir) / ".all_scores_done").touch()
        return

    # run in parallel
    with mp.Pool(processes=ncpus) as pool:
        results = pool.map(relax_pdb, jobs_to_run)

    successful_paths = [r[0] for r in results if r[1] is True]
    failed_paths = [r[0] for r in results if r[1] is False]

    # check all expected score files are present
    score_files = list(output_dir.glob("*.sc"))
    scored_ids = {f.stem for f in score_files}

    expected_ids = {f"{label}_0001" for (_, label, _, _) in jobs}
    missing_ids = expected_ids - scored_ids

    # Map missing ids back to their original PDB paths for a helpful report
    label_by_path = {path: label for (path, label, _, _) in jobs}
    missing_paths = [path for (path, label, _, _) in jobs if f"{label}_0001" in missing_ids]

    print(f"\nNumber of structures selected: {len(jobs)}")
    print(f"Number of score files present: {len(score_files)}")
    print(f"Number of missing scorefiles : {len(missing_paths)}")

    if missing_paths:
        print(f"\nMissing score files for {len(missing_paths)} structures:")
        for m in missing_paths:
            print(f" - {m}  (label={label_by_path[m]})")
        with open(output_dir / "missing_structures.txt", "w") as f:
            for m in missing_paths:
                f.write(f"{m}\n")
    else:
        print("All structures were scored successfully.")
        (Path(args.output_dir) / ".all_scores_done").touch()

if __name__ == "__main__":
    main()
