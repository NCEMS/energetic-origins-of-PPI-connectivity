#!/usr/bin/env python3

import argparse
import os
import shlex
from pathlib import Path

import pandas as pd


def select_structure(row: pd.Series):
    """
    Reproduce the structure-selection logic from run-Rosetta.py.
    """
    if (
        row.get("has_verified_sequence") is True
        and row.get("DeepTMHMM_class") in {"GLOB", "SP"}
        and pd.notna(row.get("final_structure_path"))
        and str(row.get("final_structure_source")) != "None"
    ):
        structure_path = Path(str(row["final_structure_path"])).resolve()
        label = str(
            row.get("node")
            or row.get("UniProtKB-AC")
            or structure_path.parent.name
        )
        return structure_path, label

    return None


def shell_quote(value) -> str:
    return shlex.quote(str(value))


def build_command(
    *,
    rosetta_exec: Path,
    rosetta_database: Path | None,
    structure_path: Path,
    output_dir: Path,
    label: str,
    replica: int,
) -> str:
    """
    Generate one restart-safe shell command for one protein replica.

    Rosetta runs in a task-specific temporary directory. The single generated
    PDB and score file are moved to deterministic final filenames.
    """
    replica_tag = f"{replica:04d}"

    final_pdb = output_dir / f"{label}_{replica_tag}.pdb"
    final_score = output_dir / f"{label}_{replica_tag}.sc"
    final_log = output_dir / f"{label}_{replica_tag}.log"

    database_argument = ""
    if rosetta_database is not None:
        database_argument = (
            f"-database {shell_quote(rosetta_database)} "
        )

    command = f"""
set -euo pipefail

FINAL_PDB={shell_quote(final_pdb)}
FINAL_SCORE={shell_quote(final_score)}
FINAL_LOG={shell_quote(final_log)}

if [[ -s "$FINAL_PDB" && -s "$FINAL_SCORE" ]]; then
    echo "{label} replica {replica_tag}: already complete"
    exit 0
fi

mkdir -p {shell_quote(output_dir)}

SCRATCH_BASE="${{TMPDIR:-/tmp}}"
TASK_SCRATCH=$(mktemp -d "$SCRATCH_BASE/rosetta_{label}_{replica_tag}.XXXXXX")

cleanup() {{
    rm -rf "$TASK_SCRATCH"
}}
trap cleanup EXIT

{shell_quote(rosetta_exec)} \
    {database_argument}\
    -s {shell_quote(structure_path)} \
    -relax:fast \
    -relax:constrain_relax_to_start_coords \
    -nstruct 1 \
    -score:weights ref2015 \
    -out:path:all "$TASK_SCRATCH" \
    -out:file:scorefile "$TASK_SCRATCH/score.sc" \
    -out:pdb true \
    > "$FINAL_LOG" 2>&1

mapfile -t PDB_FILES < <(find "$TASK_SCRATCH" -maxdepth 1 -type f -name '*.pdb')
mapfile -t SCORE_FILES < <(find "$TASK_SCRATCH" -maxdepth 1 -type f -name '*.sc')

if [[ "${{#PDB_FILES[@]}}" -ne 1 ]]; then
    echo "ERROR: Expected one PDB, found ${{#PDB_FILES[@]}}" >> "$FINAL_LOG"
    exit 1
fi

if [[ "${{#SCORE_FILES[@]}}" -ne 1 ]]; then
    echo "ERROR: Expected one score file, found ${{#SCORE_FILES[@]}}" >> "$FINAL_LOG"
    exit 1
fi

mv "${{PDB_FILES[0]}}" "$FINAL_PDB"
mv "${{SCORE_FILES[0]}}" "$FINAL_SCORE"

test -s "$FINAL_PDB"
test -s "$FINAL_SCORE"

echo "{label} replica {replica_tag}: complete"
""".strip()

    # Store each task as one physical line in commands.txt.
    return f"bash -lc {shlex.quote(command)}"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate one Rosetta FastRelax command per protein replica."
        )
    )

    parser.add_argument(
        "--nodes",
        required=True,
        help="Annotated nodes pickle file.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for final Rosetta PDB, score, and log files.",
    )
    parser.add_argument(
        "--commands",
        required=True,
        help="Output command-list file.",
    )
    parser.add_argument(
        "--relax-executable",
        required=True,
        help="Path to relax.static.linuxgccrelease.",
    )
    parser.add_argument(
        "--rosetta-database",
        default=None,
        help="Optional path to the Rosetta database directory.",
    )
    parser.add_argument(
        "--nstruct",
        type=int,
        default=10,
        help="Number of replicas per protein. Default: 10.",
    )
    parser.add_argument(
        "--include-completed",
        action="store_true",
        help=(
            "Include commands for replicas whose final PDB and score files "
            "already exist. By default, completed replicas are omitted."
        ),
    )

    return parser.parse_args()


def main():
    args = parse_args()

    nodes_path = Path(args.nodes).resolve()
    output_dir = Path(args.output_dir).resolve()
    commands_path = Path(args.commands).resolve()
    rosetta_exec = Path(args.relax_executable).resolve()

    rosetta_database = None
    if args.rosetta_database:
        rosetta_database = Path(args.rosetta_database).resolve()

    if not nodes_path.is_file():
        raise FileNotFoundError(f"Nodes file not found: {nodes_path}")

    if not rosetta_exec.is_file():
        raise FileNotFoundError(
            f"Rosetta executable not found: {rosetta_exec}"
        )

    if rosetta_database is not None and not rosetta_database.is_dir():
        raise NotADirectoryError(
            f"Rosetta database not found: {rosetta_database}"
        )

    if args.nstruct < 1:
        raise ValueError("--nstruct must be at least 1")

    output_dir.mkdir(parents=True, exist_ok=True)
    commands_path.parent.mkdir(parents=True, exist_ok=True)

    nodes_df = pd.read_pickle(nodes_path)

    selected = []
    missing_structures = []
    seen_labels = set()

    for _, row in nodes_df.iterrows():
        result = select_structure(row)

        if result is None:
            continue

        structure_path, label = result

        if label in seen_labels:
            raise ValueError(
                f"Duplicate Rosetta output label encountered: {label}"
            )
        seen_labels.add(label)

        if not structure_path.is_file():
            missing_structures.append((label, structure_path))
            continue

        selected.append((structure_path, label))

    commands = []
    completed = 0

    for structure_path, label in selected:
        for replica in range(1, args.nstruct + 1):
            replica_tag = f"{replica:04d}"
            final_pdb = output_dir / f"{label}_{replica_tag}.pdb"
            final_score = output_dir / f"{label}_{replica_tag}.sc"

            if (
                not args.include_completed
                and final_pdb.is_file()
                and final_pdb.stat().st_size > 0
                and final_score.is_file()
                and final_score.stat().st_size > 0
            ):
                completed += 1
                continue

            commands.append(
                build_command(
                    rosetta_exec=rosetta_exec,
                    rosetta_database=rosetta_database,
                    structure_path=structure_path,
                    output_dir=output_dir,
                    label=label,
                    replica=replica,
                )
            )

    with commands_path.open("w") as handle:
        for command in commands:
            handle.write(command)
            handle.write("\n")

    missing_path = commands_path.with_suffix(
        commands_path.suffix + ".missing-structures.tsv"
    )

    with missing_path.open("w") as handle:
        handle.write("node\tstructure_path\n")
        for label, structure_path in missing_structures:
            handle.write(f"{label}\t{structure_path}\n")

    expected_tasks = len(selected) * args.nstruct

    print(f"Input node rows:              {len(nodes_df):,}")
    print(f"Selected structures:          {len(selected):,}")
    print(f"Missing structure files:      {len(missing_structures):,}")
    print(f"Expected protein-replicas:    {expected_tasks:,}")
    print(f"Already completed replicas:   {completed:,}")
    print(f"Commands written:             {len(commands):,}")
    print(f"Command file:                 {commands_path}")
    print(f"Missing-structure report:     {missing_path}")


if __name__ == "__main__":
    main()
