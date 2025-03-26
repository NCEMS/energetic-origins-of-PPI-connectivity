 """

# profile the output DataFrame
rule profile_output:
    input:
        input_annotated_network = f"{PROCESSED_DIR}/{OUTPUT_PREFIX}network_nodes_with_annotation_and_stability.csv"
    params:
        output_prefix = OUTPUT_PREFIX,
        output_dir    = PROCESSED_DIR
    output:
        f"{PROCESSED_DIR}/{OUTPUT_PREFIX}annotated_network_summary.csv"
    shell:
        #
        conda run -n network-analysis python python-scripts/profile-output.py \
        --input_annotated_network {input.input_annotated_network} \
        --output_prefix           {params.output_prefix} \
        --output_dir              {params.output_dir}
        #
"""

import os, sys
import pandas as pd

# function to summarize the results for this DataFrame
def summarize_results(df, out_path):

	# determine the total number of unique nodes
	N_nodes = df["node"].nunique()

	# determine the number of nodes with a "verified" sequence
	N_seqs  = df[df["has_verified_sequence"] == True]["node"].nunique()

	# determine the number of nodes WITHOUT a mapped UniProtKB-AC ID
	N_no_map = df[df["UniProtKB-AC"] == "None"]["node"].nunique()

	# determine the number of nodes WITHOUT an available structure
	N_no_struc = df[df["structure_exists"] == "None"]["node"].nunique()

	# count the number of times each DeepTMHMM output class appears
	N_TM   = df[df["DeepTMHMM_class"] == "TM"]["node"].nunique()
	N_GLOB = df[df["DeepTMHMM_class"] == "GLOB"]["node"].nunique()
	N_SP   = df[df["DeepTMHMM_class"] == "SP"]["node"].nunique()
	N_SPTM = df[df["DeepTMHMM_class"] == "SP+TM"]["node"].nunique()
	N_BETA = df[df["DeepTMHMM_class"] == "BETA"]["node"].nunique()
	N_none = df[df["DeepTMHMM_class"] == "None"]["node"].nunique()

	# count the number of proteins that are disordered
	N_dis = df[df["is_disordered"] == 1]["node"].nunique()

	# count the number of proteins that have a non-None stability prediction
	N_stab = df[df["cagiada_stability"] != "None"]["node"].nunique()

	

# function carried out when called from the command line
def main():

	# setup arguments from the command line
	parser   = argparse.ArgumentParser(description="Process Yeast interactome network.")
	parser.add_argument("--output_prefix", default="0_", help="Prefix for output files")
	parser.add_argument("--output_dir", default="processed-data", help="Output directory")
	parser.add_argument("--input_annotated_network", required=True, help="Path to .csv file representing a DataFrame output by cagiada-stability.py")
        args     = parser.parse_args()

	# read in the DataFrame
	df_nodes = pd.read_csv(args.input_annotated_network)

	# get some summary information on the DataFrame
	summarize_results(df_nodes)

# entry point
if __name__ == "__main__":

        main()
