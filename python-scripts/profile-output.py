import os, sys
import argparse
import pandas as pd


# function to summarize the results for this DataFrame
def summarize_results(df, out_path):

    # determine the total number of unique nodes
    N_nodes = df["node"].nunique()

    # determine the number of nodes with a "verified" sequence
    N_seqs = df[df["has_verified_sequence"] == True]["node"].nunique()

    # determine the number of nodes WITHOUT a mapped UniProtKB-AC ID
    N_no_ID = df[df["UniProtKB-AC"].isna()]["node"].nunique()

    # determine the number of nodes WITHOUT an available structure
    N_no_struc = df[df["structure_exists"].isna()]["node"].nunique()

    # print (df["DeepTMHMM_class"])

    # count the number of times each DeepTMHMM output class appears
    N_TM = df[df["DeepTMHMM_class"] == "TM"]["node"].nunique()
    N_GLOB = df[df["DeepTMHMM_class"] == "GLOB"]["node"].nunique()
    N_SP = df[df["DeepTMHMM_class"] == "SP"]["node"].nunique()
    N_SPTM = df[df["DeepTMHMM_class"] == "SP+TM"]["node"].nunique()
    N_BETA = df[df["DeepTMHMM_class"] == "BETA"]["node"].nunique()
    N_none = df[df["DeepTMHMM_class"].isna()]["node"].nunique()

    # count the number of proteins that are disordered
    N_dis = df[df["is_disordered"] == 1]["node"].nunique()

    # count the number of proteins that have a non-None stability prediction
    N_stab = df[df["cagiada_stability"].notna()]["node"].nunique()

    print(df["cagiada_stability"])

    # print (N_nodes, N_seqs, N_no_map, N_no_struc, N_TM, N_GLOB, N_SP, N_SPTM, N_BETA, N_none, N_dis, N_stab)

    print("debugging")
    print(df[df["UniProtKB-AC"].str.strip() == "None"])
    print("Count:", df[df["UniProtKB-AC"].str.strip() == "None"]["node"].nunique())

    # write data to file
    with open(out_path, "w") as f:
        f.write("Overall summary information:\n\n")
        f.write(f"N_nodes   : {N_nodes}\n")
        f.write(f"N_seqs    : {N_seqs}\n")
        f.write(f"N_no_ID   : {N_no_ID}\n")
        f.write(f"N_no_struc: {N_no_struc}\n")
        f.write("\nDeepTMHMM results:\n\n")
        f.write(f"N_TM      : {N_TM}\n")
        f.write(f"N_GLOB    : {N_GLOB}\n")
        f.write(f"N_SP      : {N_SP}\n")
        f.write(f"N_SPTM    : {N_SPTM}\n")
        f.write(f"N_BETA    : {N_BETA}\n")
        f.write(f"N_none    : {N_none}\n")
        f.write("\nBiophysical annotation metrics:\n\n")
        f.write(f"N_dis     : {N_dis}\n")
        f.write(f"N_stab    : {N_stab}\n")

    return


# function carried out when called from the command line
def main():

    # setup arguments from the command line
    parser = argparse.ArgumentParser(description="Process Yeast interactome network.")
    parser.add_argument("--output_prefix", default="0_", help="Prefix for output files")
    parser.add_argument(
        "--output_dir", default="processed-data", help="Output directory"
    )
    parser.add_argument(
        "--input_annotated_network",
        required=True,
        help="Path to .csv file representing a DataFrame output by cagiada-stability.py",
    )
    args = parser.parse_args()

    # read in the DataFrame
    df_nodes = pd.read_csv(args.input_annotated_network)

    # get some summary information on the DataFrame
    summarize_results(
        df_nodes, f"{args.output_dir}/{args.output_prefix}annotated_network_summary.csv"
    )


# entry point
if __name__ == "__main__":

    main()
