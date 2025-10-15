import os, sys
import argparse
import pandas as pd
from functools import reduce


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--recovery")
    args = parser.parse_args()

    # load data into memory
    df = pd.read_csv(
        args.recovery,
        usecols=["Protein ID", "No. of Significant Peptides (Adj. P-value)"],
    )

    # if at least this many peptides are significantly different, consider protein different from control
    Npep_nonrefoldable = 2

    # determine if each protein has significant differences between recovery and control
    df["recovery_diff"] = (
        df["No. of Significant Peptides (Adj. P-value)"] >= Npep_nonrefoldable
    ).astype(int)
    df.rename(
        columns={"No. of Significant Peptides (Adj. P-value)": "recovery_NsigPep"},
        inplace=True,
    )

    # save the output intermediate file
    df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-processed-recovery.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
