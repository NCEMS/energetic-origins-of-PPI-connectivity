#!/usr/bin/env python3

"""
Run pipeline command scripts in the required step order.

The scripts are expected under:
    command-files/

By default, execution stops immediately when any command script fails.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


COMMAND_SCRIPTS = [
    "step1-parse-interactome-commands.sh",
    "step2-centrality-commands.sh",
    "step3-sequence-parsing-commands.sh",
    "step4-uniprot-info-commands.sh",
    "step5-idr-props-commands.sh",
    "step6-dG-calcs-commands.sh",
    "step7-Rosetta-scoring-commands.sh",
    "step8-half-life-commands.sh",
    "step9-protein-expression-commands.sh",
    "step10-translation-efficiency-commands.sh",
    "step11-predict-PTMs-commands.sh",
    "step12-chaperones-commands.sh",
    "step13-oligomers-commands.sh",
    "step14-domain-annotation-commands.sh",
    "step15-phenotype-commands.sh",
    "step16-melting-point-commands.sh",
    "step17-pathogen-commands.sh",
    "step18-entanglement-commands.sh",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Execute pipeline bash command files in numerical order."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
        help=(
            "Repository root containing command-files/. "
            "Default: current working directory."
        ),
    )
    parser.add_argument(
        "--command-dir",
        type=Path,
        default=None,
        help=(
            "Optional explicit command-files directory. "
            "Default: <repo-root>/command-files."
        ),
    )
    parser.add_argument(
        "--bash",
        default="/usr/bin/bash",
        help="Bash executable to use. Default: /usr/bin/bash",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the execution order without running the scripts.",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help=(
            "Continue running later scripts after a failure. "
            "The program still exits nonzero if any script fails."
        ),
    )
    parser.add_argument(
        "--start-step",
        type=int,
        default=1,
        choices=range(1, len(COMMAND_SCRIPTS) + 1),
        metavar=f"1-{len(COMMAND_SCRIPTS)}",
        help="Begin at this pipeline step. Default: 1",
    )
    parser.add_argument(
        "--end-step",
        type=int,
        default=len(COMMAND_SCRIPTS),
        choices=range(1, len(COMMAND_SCRIPTS) + 1),
        metavar=f"1-{len(COMMAND_SCRIPTS)}",
        help=f"Stop after this pipeline step. Default: {len(COMMAND_SCRIPTS)}",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.start_step > args.end_step:
        raise ValueError("--start-step cannot be greater than --end-step.")

    repo_root = args.repo_root.expanduser().resolve()
    command_dir = (
        args.command_dir.expanduser().resolve()
        if args.command_dir is not None
        else repo_root / "command-files"
    )

    if not repo_root.is_dir():
        raise NotADirectoryError(f"Repository root does not exist: {repo_root}")

    if not command_dir.is_dir():
        raise NotADirectoryError(
            f"Command directory does not exist: {command_dir}"
        )

    selected_scripts = COMMAND_SCRIPTS[
        args.start_step - 1 : args.end_step
    ]

    missing = [
        command_dir / script_name
        for script_name in selected_scripts
        if not (command_dir / script_name).is_file()
    ]

    if missing:
        missing_text = "\n".join(f"  {path}" for path in missing)
        raise FileNotFoundError(
            "One or more command scripts are missing:\n"
            f"{missing_text}"
        )

    print(f"Repository root: {repo_root}")
    print(f"Command directory: {command_dir}")
    print(
        f"Selected steps: {args.start_step} through {args.end_step}"
    )
    print("Mode:", "DRY RUN" if args.dry_run else "EXECUTE")
    print()

    failures: list[tuple[int, Path, int]] = []

    for step_number, script_name in enumerate(COMMAND_SCRIPTS, start=1):
        if step_number < args.start_step or step_number > args.end_step:
            continue

        script_path = command_dir / script_name

        print("=" * 80)
        print(f"STEP {step_number}: {script_name}")
        print("=" * 80)
        sys.stdout.flush()

        if args.dry_run:
            print(f"WOULD RUN: {args.bash} {script_path}")
            print()
            continue

        completed = subprocess.run(
            [args.bash, str(script_path)],
            cwd=repo_root,
            env=os.environ.copy(),
            check=False,
        )

        if completed.returncode != 0:
            failures.append(
                (step_number, script_path, completed.returncode)
            )
            print(
                f"\nERROR: step {step_number} failed with exit code "
                f"{completed.returncode}: {script_path}",
                file=sys.stderr,
            )

            if not args.continue_on_error:
                raise subprocess.CalledProcessError(
                    completed.returncode,
                    [args.bash, str(script_path)],
                )

        else:
            print(f"\nCompleted step {step_number} successfully.")

        print()

    if failures:
        print("Pipeline completed with failures:", file=sys.stderr)
        for step_number, script_path, return_code in failures:
            print(
                f"  Step {step_number}: {script_path.name} "
                f"(exit code {return_code})",
                file=sys.stderr,
            )
        raise SystemExit(1)

    print("All selected command scripts completed successfully.")


if __name__ == "__main__":
    main()
