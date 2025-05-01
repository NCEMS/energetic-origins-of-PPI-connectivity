import os, sys
import numpy as np
import pandas as pd
import multiprocessing as mp
from pathlib import Path
import argparse
import subprocess

def select_relaxed_structure(row, relaxed_dir):
    if (
        row.get("structure_exists", 0) == 1 and
        row.get("DeepTMHMM_class", 0) not in ["TM", "SP+TM", "BETA"] and
        row.get("sequence_matches_structure", 0) == True and
        row.get("has_verified_sequence", 0) == True
    ):
        structure_path = row["cleaved_structure_path"] if pd.notna(row["cleaved_structure_path"]) else row["structure_path"]
        stem = Path(structure_path).stem
        relaxed_path = Path(relaxed_dir) / f"{stem}_0001.pdb"
        print ("TEST:", relaxedf_path)
        if relaxed_path.exists():
            return relaxed_path
    return None

def score_pdb(args):
    pdb_path_str, FoldX_exec, output_dir = args
    pdb_path = Path(pdb_path_str).resolve()
    pdb_filename = pdb_path.name

    # Symlink into output_dir if not already there
    destination_pdb = output_dir / pdb_filename
    if not destination_pdb.exists():
        try:
            destination_pdb.symlink_to(pdb_path)
        except FileExistsError:
            pass
        except Exception as e:
            print(f"Warning: Could not symlink {pdb_path} → {destination_pdb}: {e}")
            return str(pdb_path), False

    # FoldX scoring command
    cmd = (
        f"{FoldX_exec} --command=Stability "
        f"--pdb={pdb_filename} "
        f"--output-file={pdb_filename.replace('.pdb', '.fxout')}"
    )

    result = subprocess.run(
        cmd, shell=True, cwd=output_dir,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    # Save stdout/stderr
    log_path = output_dir / f"{pdb_path.stem}.log"
    with open(log_path, "w") as log_file:
        log_file.write("=== STDOUT ===\n" + result.stdout)
        log_file.write("\n\n=== STDERR ===\n" + result.stderr)

    if result.returncode != 0:
        print(f"FoldX failed on {pdb_filename}. See log: {log_path}")
        return str(pdb_path), False

    # Check for expected output file
    expected_output = output_dir / pdb_filename.replace(".pdb", "_ST.fxout")
    if not expected_output.exists():
        print(f"Missing .fxout file for {pdb_filename} (no crash, but no output)")
        return str(pdb_path), False

    print(f"FoldX completed for {pdb_filename}")
    return str(pdb_path), True

def main():
    parser = argparse.ArgumentParser(description="Run FoldX stability scoring")
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--nprocessors", type=int)
    parser.add_argument("--organism_tag")
    parser.add_argument("--output_prefix", default="0")
    parser.add_argument("--struc_dir")
    parser.add_argument("--executable")
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)
    FoldX_exec = args.executable
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Select and deduplicate structure paths
    #structure_paths = list({p for p in nodes_df.apply(select_structure, axis=1).dropna()})
    structure_paths = list({p for p in nodes_df.apply(lambda row: select_relaxed_structure(row, args.struc_dir), axis=1).dropna()})

    print(f"Total unique structures to score: {len(structure_paths)}")

    job_args = [(path, FoldX_exec, output_dir) for path in structure_paths]

    with mp.Pool(processes=args.nprocessors) as pool:
        results = pool.map(score_pdb, job_args)

    success_paths = [r[0] for r in results if r[1]]
    failed_paths = [r[0] for r in results if not r[1]]

    print(f"\nTotal scoring attempts: {len(structure_paths)}")
    print(f"Successful outputs:     {len(success_paths)}")
    print(f"Failures:               {len(failed_paths)}")

    if failed_paths:
        print("\nMissing score files for the following structures:")
        for p in failed_paths:
            print(f" - {p}")
        with open(output_dir / "missing_structures.txt", "w") as f:
            for p in failed_paths:
                f.write(f"{p}\n")
    else:
        print("All structures were scored successfully.")
        (output_dir / ".all_scores_done").touch()

if __name__ == "__main__":
    main()
