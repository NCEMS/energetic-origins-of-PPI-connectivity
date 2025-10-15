import os, sys
import argparse
import pandas as pd
from functools import reduce


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--refolding1_1")
    parser.add_argument("--refolding5_1")
    parser.add_argument("--refolding120_1")
    parser.add_argument("--refolding1_2")
    parser.add_argument("--refolding5_2")
    parser.add_argument("--refolding1_2_2")
    parser.add_argument("--refolding5_2_2")
    parser.add_argument("--refolding120_2")
    args = parser.parse_args()

    # load data into memory
    df_1_1 = pd.read_csv(args.refolding1_1)
    df_5_1 = pd.read_csv(args.refolding5_1)
    df_120_1 = pd.read_csv(args.refolding120_1)
    df_1_2 = pd.read_csv(args.refolding1_2)
    df_5_2 = pd.read_csv(args.refolding5_2)
    df_1_2_2 = pd.read_csv(args.refolding1_2_2)
    df_5_2_2 = pd.read_csv(args.refolding5_2_2)
    df_120_2 = pd.read_csv(args.refolding120_2)

    datasets = {
        "refolding_1min_20210721_NsigPep": df_1_1,
        "refolding_5min_20210721_NsigPep": df_5_1,
        "refolding_120min_20210721_NsigPep": df_120_1,
        "refolding_1min_20220505_NsigPep": df_1_2,
        "refolding_5min_20220505_NsigPep": df_5_2,
        "refolding_1min_20220505_NsigPep_rep2": df_1_2_2,
        "refolding_5min_20220505_NsigPep_rep2": df_5_2_2,
        "refolding_120min_20220505_NsigPep": df_120_2,
    }

    all_columns = list(datasets.keys())
    one_min_columns = [
        "refolding_1min_20210721_NsigPep",
        "refolding_1min_20220505_NsigPep",
        "refolding_1min_20220505_NsigPep_rep2",
    ]
    five_min_columns = [
        "refolding_5min_20210721_NsigPep",
        "refolding_5min_20220505_NsigPep",
        "refolding_5min_20220505_NsigPep_rep2",
    ]
    two_hour_columns = [
        "refolding_120min_20210721_NsigPep",
        "refolding_120min_20220505_NsigPep",
    ]

    # if at least this many peptides are significantly different, consider protein nonrefoldable
    Npep_nonrefoldable = 2

    processed = [
        df[["Protein ID", "No. of Significant Peptides (Adj. P-value)"]].rename(
            columns={"No. of Significant Peptides (Adj. P-value)": name}
        )
        for name, df in datasets.items()
    ]

    combined_df = reduce(
        lambda left, right: pd.merge(left, right, on="Protein ID", how="outer"),
        processed,
    )

    ### Across all six datasets

    # add column representing whether any of the six datasets have Npep_nonrefoldable threshold met
    combined_df["nonrefoldable_any"] = (
        combined_df[all_columns].ge(Npep_nonrefoldable).any(axis=1).astype(int)
    )

    # add column representing whether all of the six datasets have Npep_nonrefoldable threshold met
    combined_df["nonrefoldable_all"] = (
        combined_df[all_columns].ge(Npep_nonrefoldable).all(axis=1).astype(int)
    )

    ### Across 1-min datasets

    # add column representing whether any of the two 1-min timepoints have Npep_nonrefoldable threshold met
    combined_df["nonrefoldable_1min_any"] = (
        combined_df[one_min_columns].ge(Npep_nonrefoldable).any(axis=1).astype(int)
    )

    # add column representing whether both of the two 1-min timepoints have Npep_nonrefoldable threshold met
    combined_df["nonrefoldable_1min_all"] = (
        combined_df[one_min_columns].ge(Npep_nonrefoldable).all(axis=1).astype(int)
    )

    ### Across 5-min datasets

    # add column representing whether any of the two 5-min timepoints have Npep_nonrefoldable threshold met
    combined_df["nonrefoldable_5min_any"] = (
        combined_df[five_min_columns].ge(Npep_nonrefoldable).any(axis=1).astype(int)
    )

    # add column representing whether both of the two 5-min timepoints have Npep_nonrefoldable threshold met
    combined_df["nonrefoldable_5min_all"] = (
        combined_df[five_min_columns].ge(Npep_nonrefoldable).all(axis=1).astype(int)
    )

    ### Across 120-min datasets

    # add column representing whether any of the two 120-min timepoints have Npep_nonrefoldable threshold met
    combined_df["nonrefoldable_120min_any"] = (
        combined_df[two_hour_columns].ge(Npep_nonrefoldable).any(axis=1).astype(int)
    )

    # add column representing whether both of the two 120-min timepoints have Npep_nonrefoldable threshold met
    combined_df["nonrefoldable_120min_all"] = (
        combined_df[two_hour_columns].ge(Npep_nonrefoldable).all(axis=1).astype(int)
    )

    ### Summarize cross-referenced results
    for col in [
        "nonrefoldable_any",
        "nonrefoldable_all",
        "nonrefoldable_1min_any",
        "nonrefoldable_1min_all",
        "nonrefoldable_5min_any",
        "nonrefoldable_5min_all",
        "nonrefoldable_120min_any",
        "nonrefoldable_120min_all",
    ]:
        print(col, combined_df[col].value_counts(), "\n")

    # save the output intermediate file
    combined_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-processed-refolding.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
