#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# Important for reading pickles that contain pint / pint-pandas objects.
import pint
import pint_pandas


def main():
    parser = argparse.ArgumentParser(
        description="Write dummy Rosetta score columns without running Rosetta."
    )

    parser.add_argument("--nodes", required=True)
    parser.add_argument("--output_nodes", required=True)
    parser.add_argument("--done_signal", required=True)
    parser.add_argument("--nstruct", type=int, required=True)

    args = parser.parse_args()

    output_nodes = Path(args.output_nodes)
    done_signal = Path(args.done_signal)

    output_nodes.parent.mkdir(parents=True, exist_ok=True)
    done_signal.parent.mkdir(parents=True, exist_ok=True)

    nodes_df = pd.read_pickle(args.nodes)

    # Edit these names to exactly match the columns normally produced by
    # add-Rosetta-scores.py, if downstream steps expect specific names.
    dummy_columns = [
        "Rosetta_total_score",
        "Rosetta_best_total_score",
        "Rosetta_mean_total_score",
        "Rosetta_median_total_score",
        "Rosetta_std_total_score",
        "Rosetta_nstruct_completed",
        "Rosetta_scores_available",
        "Rosetta_dummy_run",
    ]

    for col in dummy_columns:
        if col == "Rosetta_nstruct_completed":
            nodes_df[col] = 0
        elif col == "Rosetta_scores_available":
            nodes_df[col] = False
        elif col == "Rosetta_dummy_run":
            nodes_df[col] = True
        else:
            nodes_df[col] = np.nan

    nodes_df.to_pickle(output_nodes)

    with open(done_signal, "w") as fh:
        fh.write(
            "Dummy Rosetta scoring completed.\n"
            "No Rosetta jobs were executed.\n"
            f"Input nodes: {args.nodes}\n"
            f"Output nodes: {output_nodes}\n"
            f"nstruct placeholder: {args.nstruct}\n"
        )


if __name__ == "__main__":
    main()
