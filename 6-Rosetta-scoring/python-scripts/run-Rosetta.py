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
    if (
        row.get("structure_exists", 0) == 1
        and row.get("DeepTMHMM_class", 0) not in ["TM", "SP+TM", "BETA"]
        and row.get("sequence_matches_structure", 0) == True
        and row.get("has_verified_sequence", 0) == True
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

    pdb_path_str, rosetta_exec, output_dir = args
    pdb_path = Path(pdb_path_str).resolve()
    pdb_filename = pdb_path.name
    base = pdb_filename.replace(".pdb", "")

    relaxed_pdb = output_dir / f"{base}_0001.pdb"
    scorefile = output_dir / f"{base}_0001.sc"

    # run FastRelax
    relax_cmd = [
        rosetta_exec,
        "-s",
        str(pdb_path),
        "-relax:fast",
        "-relax:constrain_relax_to_start_coords",
        "-nstruct",
        "1",
        "-score:weights",
        "ref2015",
        "-out:file:scorefile",
        str(scorefile),
        "-out:pdb",
        "true",
        "-out:path:all",
        str(output_dir),
    ]
    print(f"Running FastRelax on {pdb_filename}")
    try:
        subprocess.run(relax_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"FastRelax failed on {pdb_filename}:\n{e.stderr}")
        return str(pdb_path), False

    if not scorefile.exists():
        print(f"Missing relaxed score file for {pdb_filename}")
        return str(pdb_path), False

    print(f"FastRelax completed for {pdb_filename}")
    return str(pdb_path), True


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

    # read in the nodes_df from file
    nodes_df = pd.read_pickle(args.nodes)

    # relatively quick test run
    # nodes_df = nodes_df[nodes_df["L"] < 120]

    # Rosetta executable to use for scoring structures
    rosetta_exec = args.relax_executable

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
        results = pool.map(relax_pdb, job_args)

        # filter success/failure
        successful_paths = [r[0] for r in results if r[1] is True]
        failed_paths = [r[0] for r in results if r[1] is False]

    # check all expected score files are present
    score_files = list(output_dir.glob("*.sc"))

    print(Path(output_dir / "scores"))

    # print some things for debugging
    print(
        f"\nNumber of structures on which scoring was attempted:", len(structure_paths)
    )
    print(
        f"Number of score files produced for these poses     :", len(score_files), "\n"
    )

    # check all expected score files are present
    score_files = list(output_dir.glob("*.sc"))
    scored_ids = {f.stem for f in score_files}
    expected_ids = {Path(p).stem + "_0001" for p in structure_paths}
    missing_ids = expected_ids - scored_ids

    missing_paths = [p for p in structure_paths if Path(p).stem in missing_ids]

    print(
        f"\nNumber of structures on which scoring was attempted: {len(structure_paths)}"
    )
    print(f"Number of score files produced for these poses     : {len(score_files)}")
    print(f"Number of missing score files                      : {len(missing_paths)}")

    if len(missing_paths) > 0:
        print(f"\nMissing score files for {len(missing_paths)} structures:")
        for m in missing_paths:
            print(f" - {m}")
        with open(output_dir / "missing_structures.txt", "w") as f:
            for m in missing_paths:
                f.write(f"{m}\n")
    else:
        print("All structures were scored successfully.")
        (Path(args.output_dir) / ".all_scores_done").touch()


if __name__ == "__main__":
    main()
