import os, sys
import numpy as np
import pandas as pd
import pint
import pint_pandas
import typing
import multiprocessing as mp
from pathlib import Path
import argparse

def select_structure(row):
    if (
        row.get("structure_exists", 0) == 1 and
        row.get("DeepTMHMM_class", 0) not in ["TM", "SP+TM", "BETA"] and
        row.get("sequence_matches_structure", 0) == True and
        row.get("has_verified_sequence", 0) == True
       ):

        if pd.notna(row["cleaved_structure_path"]):
            return row["cleaved_structure_path"]
        else:
            return row["structure_path"]
    else:
        return None

def score_pdb(args):

    pdb_path_str, rosetta_exec, output_dir = args

    pdb_path = Path(pdb_path_str)
    output_scorefile = output_dir / (pdb_path.stem + ".sc")  # use filename without extension
    cmd = (
        f"{rosetta_exec} "
        f"-in:file:s {pdb_path} "
        f"-score:weights ref2015 "
        f"-out:file:scorefile {output_scorefile} "
        f"-out:pdb false"
    )
    print(f"Scoring {pdb_path.name} ...")
    os.system(cmd)
    return pdb_path.name

def main():

    # parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Run Rosetta energy scoring of protein structures"
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
    parser.add_argument("--executable")
    parser.add_argument("--output_prefix", default="0", help="Prefix for output files")
    args = parser.parse_args()

    # read in the nodes_df from file
    nodes_df = pd.read_pickle(args.nodes)

    # Rosetta executable to use for scoring structures
    rosetta_exec = args.executable

    # number of processors over which jobs will be distributed
    ncpus = args.nprocessors

    # setup output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # create file with paths to structures to be analyzed
    structure_paths = nodes_df.apply(select_structure, axis=1).dropna().tolist()

    # setup arguments
    job_args = [(path, rosetta_exec, output_dir) for path in structure_paths]

    # run scoring in parallel
    with mp.Pool(processes=ncpus) as pool:
        results = pool.map(score_pdb, job_args)

    # check all expected score files are present
    score_files = list(output_dir.glob("*.sc"))

    print (Path(output_dir / "scores"))

    # print some things for debugging
    print(f"\nNumber of structures on which scoring was attempted:", len(structure_paths))
    print(f"Number of score files produced for these poses     :", len(score_files), "\n")

    # Identify missing structures
    structure_ids = {Path(p).stem for p in structure_paths}
    score_ids = {f.stem for f in score_files}
    missing_ids = structure_ids - score_ids
    missing_paths = [p for p in structure_paths if Path(p).stem in missing_ids]

    if missing_paths:
        print(f"\nMissing score files for {len(missing_paths)} structures:")
        for m in missing_paths:
            print(f" - {m}")
        # save to file
        with open(output_dir / "missing_structures.txt", "w") as f:
            for m in missing_paths:
                f.write(f"{m}\n")
    else:
        print("All structures were scored successfully.")
        (Path(args.output_dir) / ".all_scores_done").touch()

if __name__ == "__main__":
    main()
