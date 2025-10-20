import pandas as pd
import numpy as np
import pint
import pint_pandas
import argparse
import typing
from typing import List, Dict, Optional

def add_halflife_db1(
    nodes_df: pd.DataFrame,
    halflife_df: pd.DataFrame,
    merge_col1: str,                      # column in nodes_df
    merge_col2: str,                      # ORIGINAL column name in halflife_df
    halflife_df_cols: List[str],         # columns to keep (OLD names)
    rename_map: Optional[Dict[str, str]] = None,  # {old_name: new_name}
) -> pd.DataFrame:
    """
    Christiano et al. 2014-like merge.
    Select columns by OLD names, then rename, then merge on the (renamed) right key.
    """
    # 1) select by OLD names (so your current lists work)
    halflife_df = halflife_df[halflife_df_cols]

    # 2) rename to new schema (if provided)
    if rename_map:
        halflife_df = halflife_df.rename(columns=rename_map)

    # 3) compute the RIGHT-ON key (after rename)
    right_on_key = rename_map.get(merge_col2, merge_col2) if rename_map else merge_col2

    # 4) drop duplicates on the (renamed) ID column
    if right_on_key in halflife_df.columns:
        halflife_df = halflife_df.drop_duplicates(subset=right_on_key, keep="first")

    # 5) merge
    nodes_df = nodes_df.merge(
        halflife_df, how="left", left_on=merge_col1, right_on=right_on_key
    )
    return nodes_df


def add_halflife_db2(
    nodes_df: pd.DataFrame,
    halflife_df: pd.DataFrame,
    merge_col1: str,                      # column in nodes_df
    merge_col2: str,                      # ORIGINAL column name in halflife_df
    halflife_df_cols: List[str],         # columns to keep (OLD names)
    rename_map: Optional[Dict[str, str]] = None,  # {old_name: new_name}
) -> pd.DataFrame:
    """
    Martin-Perez & Villén 2017-like merge.
    Select by OLD names, then rename, then merge on the (renamed) right key.
    """
    # 1) select by OLD names
    halflife_df = halflife_df[halflife_df_cols]

    # 2) rename (if provided)
    if rename_map:
        halflife_df = halflife_df.rename(columns=rename_map)

    # 3) compute right key after rename
    right_on_key = rename_map.get(merge_col2, merge_col2) if rename_map else merge_col2

    # 4) merge
    nodes_df = nodes_df.merge(
        halflife_df, how="left", left_on=merge_col1, right_on=right_on_key
    )
    return nodes_df

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", required=True, help="Input nodes file")
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where results will be saved",
    )
    parser.add_argument("--output_prefix", default="0", help="Prefix for output files")
    parser.add_argument(
        "--output_suffix", default="step8", help="Suffix for output files"
    )
    parser.add_argument("--organism_tag", help="Organism label for this run")
    parser.add_argument("--halflife_db1")
    parser.add_argument("--halflife_db2")
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)

    # ADD FIRST HALF-LIFE DATABASE HERE
    halflife_df1 = pd.read_csv(args.halflife_db1)
    halflife_df1_cols = [
        "ENSG",
        "Degradation rates (min-1)",
        "R2 (quality of curve fitting)",
        "t1/2 (min)",
    ]
    rename_map1 = {"ENSG":"ENSG",
        "Degradation rates (min-1)":"Christiano_degradation_rate(min-1)",
        "R2 (quality of curve fitting)":"Christiano_degradation_R2",
        "t1/2 (min)":"Christiano_halflife_min"}
    merge_col1 = "node"
    merge_col2 = "ENSG"
    nodes_df = add_halflife_db1(
        nodes_df, halflife_df1, merge_col1, merge_col2, halflife_df1_cols, rename_map1
    )

    # ADD SECOND HALF-LIFE DATABASE HERE
    halflife_df2 = pd.read_csv(args.halflife_db2)
    halflife_df2_cols = ["Protein IDs", "t_12_avg", "t_12_sd", "t_12_cv"]
    rename_map2 = {"Protein IDs":"Protein IDs", "t_12_avg":"Villen_halflife_hours", "t_12_sd":"Villen_halflife_SD_hours", "t_12_cv":"Villen_halflife_CV_hours"}
    merge_col1 = "node"
    merge_col2 = "Protein IDs"
    nodes_df = add_halflife_db2(
        nodes_df, halflife_df2, merge_col1, merge_col2, halflife_df2_cols, rename_map2
    )

    # convert Martin-Perez & Villen numbers from hours to minutes
    nodes_df["Villen_halflife_min"] = nodes_df["Villen_halflife_hours"]*60.0
    nodes_df["Villen_halflife_SD_min"] = nodes_df["Villen_halflife_SD_hours"]*60.0

    # output the results to file
    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )


if __name__ == "__main__":
    main()
