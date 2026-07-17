#!/usr/bin/env python3

"""
Purge processed pipeline outputs from selected pipeline-step directories.

By default, this script performs a dry run and prints what would be removed.
Pass --execute to actually delete files.

Special handling:
    11-predict-PTMs/processed-data
        Only files matching:
            *predictions-processed.csv
            *.step11.pkl
        are removed. All other files and subdirectories are preserved.

For all other listed pipeline steps:
    Every file and subdirectory inside processed-data is removed, but the
    processed-data directory itself is preserved.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Iterable


PIPELINE_STEPS = [
    "1-parse-interactome",
    "2-network-centrality",
    "3-sequence-parsing",
    "4-uniprot-annotation",
    "5-idr-properties",
    "6-dG-calculations",
    "7-Rosetta-scoring",
    "8-protein-half-life",
    "9-protein-expression",
    "10-translation-efficiency",
    "11-predict-PTMs",
    "12-chaperones",
    "13-oligomers",
    "14-domain-annotations",
    "15-phenotype",
    "16-thermal-proteome-profiling",
    "17-pathogen-target",
    "18-entanglement",
]

PTM_STEP = "11-predict-PTMs"
PTM_PATTERNS = [
    "*predictions-processed.csv",
    "*step11.pkl",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Purge processed-data outputs from pipeline-step directories. "
            "Runs as a dry run unless --execute is provided."
        )
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
        help=(
            "Repository root containing the numbered pipeline directories. "
            "Default: current working directory."
        ),
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually delete files. Without this flag, only print a dry run.",
    )
    parser.add_argument(
        "--ignore-missing",
        action="store_true",
        help=(
            "Do not fail when a listed step or processed-data directory is "
            "missing."
        ),
    )
    return parser.parse_args()


def iter_all_children(directory: Path) -> Iterable[Path]:
    """Yield all immediate children of a processed-data directory."""
    yield from sorted(directory.iterdir(), key=lambda path: path.name)


def iter_ptm_targets(directory: Path) -> Iterable[Path]:
    """
    Yield only the PTM outputs that should be purged.

    Duplicate matches are removed while preserving deterministic order.
    """
    matches: set[Path] = set()

    for pattern in PTM_PATTERNS:
        matches.update(directory.glob(pattern))

    yield from sorted(matches, key=lambda path: path.name)


def remove_path(path: Path, execute: bool) -> None:
    """Print or remove one file, symlink, or directory."""
    action = "DELETE" if execute else "WOULD DELETE"
    print(f"{action}: {path}")

    if not execute:
        return

    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)
    else:
        raise FileNotFoundError(f"Path disappeared before deletion: {path}")


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.expanduser().resolve()

    if not repo_root.is_dir():
        raise NotADirectoryError(f"Repository root does not exist: {repo_root}")

    print(f"Repository root: {repo_root}")
    print("Mode:", "EXECUTE" if args.execute else "DRY RUN")
    print()

    total_targets = 0
    missing_directories: list[Path] = []

    for step in PIPELINE_STEPS:
        processed_dir = repo_root / step / "processed-data"

        if not processed_dir.is_dir():
            missing_directories.append(processed_dir)
            print(f"MISSING: {processed_dir}")
            continue

        print(f"[{step}]")

        if step == PTM_STEP:
            targets = list(iter_ptm_targets(processed_dir))
        else:
            targets = list(iter_all_children(processed_dir))

        if not targets:
            print("  No matching files or directories.")
            print()
            continue

        for target in targets:
            remove_path(target, execute=args.execute)
            total_targets += 1

        print()

    print(
        f"{'Deleted' if args.execute else 'Matched'} "
        f"{total_targets:,} file(s) or directorie(s)."
    )

    if missing_directories:
        print(f"Missing processed-data directories: {len(missing_directories):,}")

        if not args.ignore_missing:
            missing_text = "\n".join(f"  {path}" for path in missing_directories)
            raise FileNotFoundError(
                "One or more expected processed-data directories were missing:\n"
                f"{missing_text}\n"
                "Use --ignore-missing to allow missing directories."
            )


if __name__ == "__main__":
    main()
