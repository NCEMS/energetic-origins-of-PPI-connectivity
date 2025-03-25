import os
import sys
import time
import torch
import numpy as np
import pandas as pd
import argparse
import esm
import pytest
from datetime import datetime
from esm.inverse_folding.util import load_structure, extract_coords_from_structure, CoordBatchConverter
from esm.inverse_folding.multichain_util import extract_coords_from_complex, _concatenate_coords, load_complex_coords

# class to store parameters for each dG calculation to be carried out
class cagiada:

	def __init__(self, structure, chainID, output_name):

		self.structure = structure
		self.chainID = chainID
		self.output_name = output_name

# reads values from last two lines of a file as output by predict_dG function
def read_values_from_file(f_path):

	with open(f_path, "r") as file:
		lines = file.readlines()
		val1  = lines[-2].strip().split(",")[1] # sum of likelihoods
		val2  = lines[-1].strip().split(",")[1] # dG prediction

	return float(val1), float(val2)

# test that results match expectation
def test_P78285_results_within_tolerance(f_path):

	# test tolerance; need to allow for some difference due to floating point differences
	tol = 1e-4

	# expected values from https://colab.research.google.com/github/KULL-Centre/_2024_cagiada_stability/blob/main/stab_ESM_IF.ipynb
	# calculated on March 19, 2025 by Dan Nissley using AF2 input structure AF-P78285-F1-model_v4.pdb
	ref_likelihood_sum = 108.69234741592663 # kcal/mol
	ref_dG_predicted   = 11.934800287565977 # kcal/mol

	# extract predicted values from the new run of the test protein
	new_likelihood_sum, new_dG_predicted = read_values_from_file(f_path)

	# check assertions
	assert abs(ref_likelihood_sum - new_likelihood_sum) < tol, "Reference and calculated likelihood sums do not match"
	assert abs(ref_dG_predicted - new_dG_predicted) < tol, "Reference and calculated ΔG do not match"

# test that sequences match between structure used for dG prediction and sequence used for all other predictions (e.g., metapredict)
def test_sequences_match(AF2_fasta_path, seq2):

	with open(AF2_fasta_path, "r") as f:
		temp = f.readlines()

	seq1 = temp[1].strip()

	assert seq1 == seq2

# function to check CUDA memory being used by script
def print_gpu_memory_usage():

	allocated = torch.cuda.memory_allocated()
	reserved  = torch.cuda.memory_reserved()
	print(f"GPU Memory Allocated: {allocated/1e6:.2f} MB")
	print(f"GPU Memory Reserved:  {reserved/1e6:.2f} MB")

# function to run the model for a set of inputs
def run_model(coords, sequence, model, alphabet, chain_target='A'):
	"""
	Run the ESM inverse folding model on the input sequence and coordinates.

	Args:
		coords (tensor): Input coordinates.
		sequence (str): Target sequence.
		model (ESM model): Pre-trained inverse folding model.
		alphabet (Alphabet): Model alphabet.
		chain_target (str): Target chain identifier.

	Returns:
		tensor: Probability scores for each residue position.
	"""
	device = next(model.parameters()).device
	batch_converter = CoordBatchConverter(alphabet)
	batch = [(coords, None, sequence)]
	coords, confidence, strs, tokens, padding_mask = batch_converter(batch, device=device)

	prev_output_tokens = tokens[:, :-1].to(device)
	target = tokens[:, 1:]
	target_padding_mask = (target == alphabet.padding_idx)

	logits, _ = model.forward(coords, padding_mask, confidence, prev_output_tokens)
	logits_swapped = torch.swapaxes(logits, 1, 2)
	token_probs = torch.softmax(logits_swapped, dim=-1)

	return token_probs

# compute likelihood scores for each residue in the sequence
def score_variants(sequence, token_probs, alphabet):
	"""
	Compute wild-type likelihood scores for each residue in the sequence.

	Args:
		sequence (str): Target sequence.
		token_probs (tensor): Output probabilities from the model.
		alphabet (Alphabet): Model alphabet.

	Returns:
		list, list: List of residues and their corresponding scores.
	"""
	aa_list = []
	wt_scores = []

	alphabetAA_L_D = {
		'-': 0, '_': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7,
		'I': 8, 'K': 9, 'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15,
		'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20
	}
	alphabetAA_D_L = {v: k for k, v in alphabetAA_L_D.items()}

	for i, n in enumerate(sequence):
		aa_list.append(n + str(i + 1))
		score_pos = [masked_absolute(alphabetAA_D_L[j], i, token_probs, alphabet) for j in range(1, 21)]
		wt_scores.append(score_pos[alphabetAA_L_D[n] - 1])

	return aa_list, wt_scores

# compute masked absolute probability score
def masked_absolute(mut, idx, token_probs, alphabet):
	"""
	Compute masked absolute probability score.

	Args:
		mut (str): Mutated residue.
		idx (int): Position in sequence.
		token_probs (tensor): Output probabilities from the model.
		alphabet (Alphabet): Model alphabet.

	Returns:
		float: Probability score for mutation.
	"""
	mt_encoded = alphabet.get_idx(mut)
	return token_probs[0, idx, mt_encoded].item()

# function to carry out various steps in model pipeline
def predict_dG(cagiada_info, output_dir, model, alphabet, create_file = False):

	# cagiada_info is a class object of the type cagiada
	# load structure
	structure                            = load_structure(cagiada_info.structure, cagiada_info.chainID)
	coords_structure, sequence_structure = extract_coords_from_structure(structure)

	# compute model probabilities
	prob_tokens        = run_model(coords_structure, sequence_structure, model, alphabet, chain_target=cagiada_info.chainID)
	aa_list, wt_scores = score_variants(sequence_structure, prob_tokens, alphabet)

	# compute ΔG predictions
	a, b       = 0.10413378327743603, 0.6162549378400894  # Fitting parameters from manuscript
	dg_IF      = np.nansum(wt_scores)
	dg_kcalmol = a * dg_IF + b

	# only save a file is specifically requested by the user
	if create_file:

		# save results
		output_file = os.path.join(output_dir, cagiada_info.output_name+"-cagiada-dG.csv")
		output_df   = pd.DataFrame({'Residue': aa_list + ['dG_IF', 'dG_kcalmol'], 'score': wt_scores + [dg_IF, dg_kcalmol]})
		output_df.to_csv(output_file, sep=',', index=False)

	# return the absolute free energy estimate
	return dg_kcalmol

# MAIN
def main():

	# check if CUDA is available
	if torch.cuda.is_available():
		print("CUDA is available")
		# Print the current device index
		device_index = torch.cuda.current_device()
		print("Current GPU device index:", device_index)
		# Print the name of the GPU
		print("GPU Name:", torch.cuda.get_device_name(device_index))
	else:
		print("CUDA is not available; this script will take a few days to run. Exiting")
		sys.exit()

	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Run stability predictions using ESM inverse folding.")
	parser.add_argument("--input_node_file", required=True, help="Output from network-analysis.py")
	parser.add_argument("--output_dir", default="processed-data", help="Directory where results will be saved (default: 'outputs')")
	parser.add_argument("--output_prefix", default="0_", help="Prefix for output files")
	args = parser.parse_args()

	# load esm model
	IF_model_name   = "data-files/esm_if1_gvp4_t16_142M_UR50.pt"
	model, alphabet = esm.pretrained.load_model_and_alphabet(IF_model_name)
	model.to("cuda")
	model.eval().cuda().requires_grad_(False)

	# load network node information and prepare set of commands to be run
	nodes_df = pd.read_csv(args.input_node_file)

	# testing purposes only
	#nodes_df = nodes_df.head(50)

	# all protein structure predictions from EBI for S288C contain a single chain with name A
	chainID = "A"

	# run a test to make sure results match expectation
	# run on the default protein P78285 from Cagiada Google Colab notebook
	predict_dG(cagiada("test/P78285/AF-P78285-F1-model_v4.pdb", chainID, "AF-P78285-F1-model_v4"), args.output_dir, model, alphabet, create_file = True)
	torch.cuda.empty_cache()
	test_P78285_results_within_tolerance(os.path.join(args.output_dir, "AF-P78285-F1-model_v4-cagiada-dG.csv"))

	# add empty column to hold dG prediction
	nodes_df["cagiada_stability"] = "None"

	# run predictions in series using CUDA
	start = datetime.now()
	for i, r in nodes_df.iterrows():

		if r['has_verified_sequence'] == True and r['DeepTMHMM_class'] == "GLOB" and r['structure_exists'] == 1:

			curr_cagiada = cagiada(r['structure_path'], chainID, r["UniProtKB-AC"])

			# test that the sequences between the orf_trans.fasta file from SGD and AF2 structures from EBI match
			test_sequences_match(r['structure_path'].split('.pdb')[0]+'.fasta', r["sequence"])

			# generate the prediction and add it to the DataFrame
			abs_dG = predict_dG(curr_cagiada, args.output_dir, model, alphabet)
			nodes_df.at[i, "cagiada_stability"] = abs_dG

			# clean up GPU memory after each iteration
			torch.cuda.empty_cache()

			print ("Done with ΔG prediction for:", r["UniProtKB-AC"])

		else:
			pass

	print ("Total execution time is:", datetime.now() - start) # total time for all dG predictions

	nodes_df   = nodes_df.replace("", "None")
	nodes_df   = nodes_df.fillna("None")

	# save updated nodes_df to file with a new name
	nodes_df.to_csv(f"{args.output_dir}/{args.output_prefix}network_nodes_with_annotation_and_stability.csv", index=False, na_rep=None)

# execute main when run from command line
if __name__ == "__main__":
	main()
