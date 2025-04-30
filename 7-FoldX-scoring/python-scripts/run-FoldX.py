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

def select_structure(row):
    if row.get("structure_exists", 0) == 1:
        if pd.notna(row["cleaved_structure_path"]):
            return row["cleaved_structure_path"]
        else:
            return row["structure_path"]
    else:
        return None

def score_pdb(args):

    pdb_path_str, FoldX_exec, output_dir = args
    pdb_path = Path(pdb_path_str).resolve()
    pdb_filename = pdb_path.name

    # Symlink the pdb into the output directory if it doesn't already exist
    destination_pdb = output_dir / pdb_filename
    if not destination_pdb.exists():
        try:
            destination_pdb.symlink_to(pdb_path)
        except FileExistsError:
            pass  # In case it was created in parallel by another process
        except Exception as e:
            print(f"Warning: Failed to create symlink for {pdb_path} -> {destination_pdb}: {e}")
            raise

    # run FoldX in output_dir with the pdb file name
    cmd = (
        f"{FoldX_exec} --command=Stability "
        f"--pdb={pdb_filename} "
        f"--output-file={pdb_filename.replace('.pdb', '.fxout')}"
    )
    #print(cmd)
    #subprocess.run(cmd, shell=True, cwd=output_dir, check=True)
    # Run FoldX, capturing stdout and stderr
    result = subprocess.run(
        cmd,
        shell=True,
        cwd=output_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # Save stdout and stderr to a .log file
    log_path = output_dir / f"{pdb_path.stem}.log"
    with open(log_path, "w") as log_file:
        log_file.write("=== STDOUT ===\n")
        log_file.write(result.stdout)
        log_file.write("\n\n=== STDERR ===\n")
        log_file.write(result.stderr)

    # Check if FoldX exited cleanly
    if result.returncode != 0:
        print(f"FoldX error on {pdb_filename}. See log: {log_path}")
        return {"pdb": pdb_filename, "error": True, "log_file": str(log_path)}

    print(f"FoldX completed for {pdb_filename}. Log saved: {log_path}")
    return {"pdb": pdb_filename, "error": False, "log_file": str(log_path)}

def main():

    # parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Run FoldX energy scoring of protein structures"
    )
    parser.add_argument(
        "--nodes", required=True
    )
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where results will be saved (default: 'outputs')",
    )
    parser.add_argument(
        "--nprocessors", type=int
    )
    parser.add_argument("--organism_tag")
    parser.add_argument("--output_prefix", default="0", help="Prefix for output files")
    parser.add_argument("--executable")
    args = parser.parse_args()

    # read in the nodes_df from file
    nodes_df = pd.read_pickle(args.nodes)

    # for testing purposes only
    #nodes_df = nodes_df.head(2)

    # Rosetta executable to use for scoring structures
    FoldX_exec = args.executable

    # number of processors over which jobs will be distributed
    ncpus = args.nprocessors

    # setup output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # create file with paths to structures to be analyzed
    structure_paths = nodes_df.apply(select_structure, axis=1).dropna().tolist()

    print (structure_paths)

    # setup arguments
    job_args = [(path, FoldX_exec, output_dir) for path in structure_paths]

    # run scoring in parallel
    with mp.Pool(processes=ncpus) as pool:
        results = pool.map(score_pdb, job_args)

    # check all expected score files are present
    score_files = list(output_dir.glob("*.fxout"))

    # print some things for debugging
    print (f"\nNumber of structures on which scoring was attempted:", len(structure_paths))
    print (f"Number of score files produced for these poses     :", len(score_files), "\n")

    if len(score_files) == len(structure_paths):
        (Path(args.output_dir) / ".all_scores_done").touch()

if __name__ == "__main__":
    main()
