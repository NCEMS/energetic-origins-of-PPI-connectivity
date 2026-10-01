#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


# =====================================================================
# Structure selection
# =====================================================================


def select_structure(row):
    """
    Apply the same selection criteria used by the FastRelax pipeline.

    Returns
    -------
    tuple or None
        A tuple containing:

        1. The original input-structure path.
        2. The label used to name the Rosetta output files.

        Returns None when the protein does not satisfy the selection
        criteria.
    """
    if (
        row.get("has_verified_sequence") is True
        and row.get("DeepTMHMM_class") in ["GLOB", "SP"]
        and pd.notna(row.get("final_structure_path"))
        and str(row.get("final_structure_source")) != "None"
    ):
        structure_path = str(
            row["final_structure_path"]
        )

        label = str(
            row.get("node")
            or row.get("UniProtKB-AC")
            or Path(structure_path).parent.name
        )

        return structure_path, label

    return None


# =====================================================================
# Replicate inspection
# =====================================================================


def inspect_replicates(
    output_dir,
    label,
    nstruct,
):
    """
    Inspect the expected Rosetta outputs for one protein.

    A replicate is considered complete only when both its PDB file and
    score file exist.
    """
    complete_replicates = []
    missing_pdb_replicates = []
    missing_score_replicates = []

    for replicate in range(
        1,
        nstruct + 1,
    ):
        stem = (
            f"{label}_{replicate:04d}"
        )

        pdb_path = (
            output_dir / f"{stem}.pdb"
        )

        score_path = (
            output_dir / f"{stem}.sc"
        )

        pdb_exists = pdb_path.is_file()
        score_exists = score_path.is_file()

        if pdb_exists and score_exists:
            complete_replicates.append(
                replicate
            )

        if not pdb_exists:
            missing_pdb_replicates.append(
                replicate
            )

        if not score_exists:
            missing_score_replicates.append(
                replicate
            )

    return {
        "complete_replicates": complete_replicates,
        "missing_pdb_replicates": missing_pdb_replicates,
        "missing_score_replicates": missing_score_replicates,
    }


# =====================================================================
# Formatting
# =====================================================================


def format_replicate_list(
    replicate_numbers,
):
    """
    Format replicate numbers as a comma-separated string.
    """
    return ",".join(
        str(number)
        for number in replicate_numbers
    )


# =====================================================================
# Main
# =====================================================================


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Identify proteins missing one or more expected Rosetta "
            "FastRelax replicate files."
        )
    )

    parser.add_argument(
        "--nodes",
        required=True,
        help=(
            "Pickle file containing the nodes DataFrame used by the "
            "FastRelax pipeline."
        ),
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help=(
            "Directory containing the Rosetta PDB and score files."
        ),
    )

    parser.add_argument(
        "--nstruct",
        type=int,
        default=10,
        help=(
            "Expected number of replicates per protein "
            "(default: 10)."
        ),
    )

    parser.add_argument(
        "--report",
        default=None,
        help=(
            "Path for the output CSV report. By default, the report "
            "is written to "
            "<output-dir>/incomplete_rosetta_structures.csv."
        ),
    )

    args = parser.parse_args()

    nodes_path = Path(
        args.nodes
    )

    output_dir = Path(
        args.output_dir
    )

    if args.report is None:
        report_path = (
            output_dir
            / "incomplete_rosetta_structures.csv"
        )
    else:
        report_path = Path(
            args.report
        )

    if not nodes_path.is_file():
        raise FileNotFoundError(
            f"Nodes pickle not found: {nodes_path}"
        )

    if not output_dir.is_dir():
        raise NotADirectoryError(
            f"Rosetta output directory not found: {output_dir}"
        )

    if args.nstruct < 1:
        raise ValueError(
            "--nstruct must be at least 1."
        )

    # -----------------------------------------------------------------
    # Read nodes table and reproduce the pipeline selection
    # -----------------------------------------------------------------

    nodes_df = pd.read_pickle(
        nodes_path
    )

    required_columns = {
        "has_verified_sequence",
        "DeepTMHMM_class",
        "final_structure_path",
        "final_structure_source",
    }

    missing_columns = (
        required_columns
        - set(nodes_df.columns)
    )

    if missing_columns:
        raise KeyError(
            "The nodes DataFrame is missing required columns: "
            + ", ".join(
                sorted(missing_columns)
            )
        )

    selected = nodes_df.apply(
        select_structure,
        axis=1,
    )

    selected = [
        item
        for item in selected
        if item is not None
    ]

    if not selected:
        print(
            "No structures passed the selection criteria."
        )
        return

    # Do not require the original input PDB to still exist. The purpose
    # of this program is to audit outputs that have already been
    # generated, potentially on a different computer.
    expected_structures = selected

    # -----------------------------------------------------------------
    # Ensure output labels are unique
    # -----------------------------------------------------------------

    label_counts = pd.Series(
        [
            label
            for _, label in expected_structures
        ],
        dtype="object",
    ).value_counts()

    duplicate_labels = label_counts[
        label_counts > 1
    ]

    if not duplicate_labels.empty:
        duplicate_text = ", ".join(
            (
                f"{label} "
                f"({count} selected rows)"
            )
            for label, count
            in duplicate_labels.items()
        )

        raise ValueError(
            "Multiple selected structures use the same output label. "
            "The output files cannot be assigned unambiguously. "
            f"Duplicated labels: {duplicate_text}"
        )

    # -----------------------------------------------------------------
    # Inspect expected outputs
    # -----------------------------------------------------------------

    incomplete_records = []
    complete_protein_count = 0
    proteins_with_no_outputs = 0

    for (
        structure_path,
        label,
    ) in expected_structures:
        inspection = inspect_replicates(
            output_dir=output_dir,
            label=label,
            nstruct=args.nstruct,
        )

        complete_replicates = inspection[
            "complete_replicates"
        ]

        missing_pdb_replicates = inspection[
            "missing_pdb_replicates"
        ]

        missing_score_replicates = inspection[
            "missing_score_replicates"
        ]

        complete_count = len(
            complete_replicates
        )

        if complete_count == args.nstruct:
            complete_protein_count += 1
            continue

        if (
            len(missing_pdb_replicates) == args.nstruct
            and len(missing_score_replicates) == args.nstruct
        ):
            proteins_with_no_outputs += 1

        incomplete_records.append(
            {
                "label": label,
                "input_structure_path": structure_path,
                "complete_replicate_count": complete_count,
                "expected_replicate_count": args.nstruct,
                "complete_replicates": format_replicate_list(
                    complete_replicates
                ),
                "missing_pdb_replicates": format_replicate_list(
                    missing_pdb_replicates
                ),
                "missing_score_replicates": format_replicate_list(
                    missing_score_replicates
                ),
            }
        )

    # -----------------------------------------------------------------
    # Write report
    # -----------------------------------------------------------------

    report_columns = [
        "label",
        "input_structure_path",
        "complete_replicate_count",
        "expected_replicate_count",
        "complete_replicates",
        "missing_pdb_replicates",
        "missing_score_replicates",
    ]

    report_df = pd.DataFrame(
        incomplete_records,
        columns=report_columns,
    )

    if not report_df.empty:
        report_df = report_df.sort_values(
            by=[
                "complete_replicate_count",
                "label",
            ],
            ascending=[
                True,
                True,
            ],
        ).reset_index(
            drop=True
        )

    report_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    report_df.to_csv(
        report_path,
        index=False,
    )

    # -----------------------------------------------------------------
    # Print summary
    # -----------------------------------------------------------------

    incomplete_count = len(
        report_df
    )

    partially_complete_count = (
        incomplete_count
        - proteins_with_no_outputs
    )

    print()
    print(
        "Rosetta replicate audit"
    )
    print(
        "======================="
    )
    print(
        f"Expected proteins:                 "
        f"{len(expected_structures):,}"
    )
    print(
        f"Expected replicates per protein:   "
        f"{args.nstruct:,}"
    )
    print(
        f"Proteins with all replicates:      "
        f"{complete_protein_count:,}"
    )
    print(
        f"Proteins with partial outputs:     "
        f"{partially_complete_count:,}"
    )
    print(
        f"Proteins with no outputs:          "
        f"{proteins_with_no_outputs:,}"
    )
    print(
        f"Total proteins missing replicates: "
        f"{incomplete_count:,}"
    )
    print(
        f"Report written to: {report_path}"
    )

    if not report_df.empty:
        display_columns = [
            "label",
            "complete_replicate_count",
            "missing_pdb_replicates",
            "missing_score_replicates",
        ]

        print()
        print(
            "Incomplete proteins"
        )
        print(
            "==================="
        )
        print(
            report_df[
                display_columns
            ].to_string(
                index=False
            )
        )


if __name__ == "__main__":
    main()
