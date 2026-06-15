import os, sys
import argparse
import pandas as pd
import pint
import pint_pandas


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_suffix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--input_nodes")
    parser.add_argument("--AF2_data")
    parser.add_argument("--exp_data")
    args = parser.parse_args()

    # column names for AF2 data
    af2_columns = [
        "AF2_gene",
        "AF2_PDB",
        "AF2_chain",
        "AF2_ENT-ID",
        "AF2_gn",
        "AF2_N_term_thread",
        "AF2_gc",
        "AF2_C_term_thread",
        "AF2_i",
        "AF2_j",
        "AF2_NC",
        "AF2_NC_wbuff",
        "AF2_NC_region",
        "AF2_crossings",
        "AF2_crossings_wbuff",
        "AF2_crossings_region",
        "AF2_ent_region",
        "AF2_loopsize",
        "AF2_num_zipper_nc",
        "AF2_perc_bb_loop",
        "AF2_num_loop_contacting_res",
        "AF2_num_cross_nearest_neighbors",
        "AF2_ent_coverage",
        "AF2_min_N_prot_depth_left",
        "AF2_min_N_thread_depth_left",
        "AF2_min_N_thread_slippage_left",
        "AF2_min_C_prot_depth_right",
        "AF2_min_C_thread_depth_right",
        "AF2_min_C_thread_slippage_right",
        "AF2_prot_size",
        "AF2_ACO",
        "AF2_RCO",
        "AF2_CCBond",
    ]

    # column names for experimental data
    exp_columns = [
        "exp_gene",
        "exp_PDB",
        "exp_chain",
        "exp_ENT-ID",
        "exp_gn",
        "exp_N_term_thread",
        "exp_gc",
        "exp_C_term_thread",
        "exp_i",
        "exp_j",
        "exp_NC",
        "exp_NC_wbuff",
        "exp_NC_region",
        "exp_crossings",
        "exp_crossings_wbuff",
        "exp_crossings_region",
        "exp_ent_region",
        "exp_loopsize",
        "exp_num_zipper_nc",
        "exp_perc_bb_loop",
        "exp_num_loop_contacting_res",
        "exp_num_cross_nearest_neighbors",
        "exp_ent_coverage",
        "exp_min_N_prot_depth_left",
        "exp_min_N_thread_depth_left",
        "exp_min_N_thread_slippage_left",
        "exp_min_C_prot_depth_right",
        "exp_min_C_thread_depth_right",
        "exp_min_C_thread_slippage_right",
        "exp_prot_size",
        "exp_ACO",
        "exp_RCO",
        "exp_CCBond",
    ]

    # load the annotated nodes file from the previous pipeline step
    nodes_df = pd.read_pickle(args.input_nodes)

    # load the AF2 structure derived data and merge with nodes_df
    af2 = pd.read_csv(args.AF2_data, sep="|", names=af2_columns)

    group_sizes = af2.groupby("AF2_gene")["AF2_gene"].transform("count")

    def AF2_compute_contains_entanglement(row):
        if group_sizes.loc[row.name] > 1:
            return 1
        else:
            return 0 if pd.isna(row["AF2_ENT-ID"]) else 1

    af2["AF2_contains_entanglement"] = af2.apply(
        AF2_compute_contains_entanglement, axis=1
    )
    af2.drop_duplicates(subset="AF2_gene", inplace=True)

    # merge AF2 entanglement information into nodes_df
    nodes_df = nodes_df.merge(
        af2[["AF2_gene", "AF2_contains_entanglement"]],
        how="left",
        right_on="AF2_gene",
        left_on="UniProtKB-AC",
    )

    # load the experimental structure derived data and merge with nodes_df
    exp = pd.read_csv(args.exp_data, sep="|", names=exp_columns)

    group_sizes = exp.groupby("exp_gene")["exp_gene"].transform("count")

    def exp_compute_contains_entanglement(row):
        if group_sizes.loc[row.name] > 1:
            return 1
        else:
            return 0 if pd.isna(row["exp_ENT-ID"]) else 1

    exp["exp_contains_entanglement"] = exp.apply(
        exp_compute_contains_entanglement, axis=1
    )
    exp.drop_duplicates(subset="exp_gene", inplace=True)

    # merge into nodes_df
    nodes_df = nodes_df.merge(
        exp[["exp_gene", "exp_contains_entanglement"]],
        how="left",
        right_on="exp_gene",
        left_on="UniProtKB-AC",
    )

    # remove unneeded columns
    nodes_df.drop(columns=["exp_gene", "AF2_gene"], inplace=True)

    # save the results to file
    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"
    )


if __name__ == "__main__":
    main()
