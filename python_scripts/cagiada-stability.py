import os
import sys
import time
import torch
import numpy as np
import pandas as pd
import argparse
import esm
from esm.inverse_folding.util import load_structure, extract_coords_from_structure, CoordBatchConverter
from esm.inverse_folding.multichain_util import extract_coords_from_complex, _concatenate_coords, load_complex_coords

# class to store parameters for each dG calculation to be carried out
class cagiada:

	def __init__(self, structure, chainID, output_name):

		self.structure = structure
		self.chainID = chainID
		self.output_name = output_name
		

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
def predict_dG(cagiada_info, output_dir, model, alphabet):

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

	#print(f"ΔG predicted (likelihood sum): {dg_IF}")
	#print(f"ΔG predicted (kcal/mol): {dg_kcalmol}")

	# save results
	output_file = os.path.join(output_dir, cagiada_info.output_name+"-cagiada-dG.csv")
	output_df   = pd.DataFrame({'Residue': aa_list + ['dG_IF', 'dG_kcalmol'], 'score': wt_scores + [dg_IF, dg_kcalmol]})
	output_df.to_csv(output_file, sep=',', index=False)
	#print(f"Results saved to {output_file}"

# MAIN
def main():

	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Run stability predictions using ESM inverse folding.")
	#parser.add_argument("--structure", required=True, help="Path to the input PDB structure file")
	#parser.add_argument("--chain", required=True, help="PDB chain identifier (e.g., 'A')")
	#parser.add_argument("--output_name", required=True, help="Output filename for the processed PDB")
	parser.add_argument("--input_node_file", required=True, help="Output from network-analysis.py")
	parser.add_argument("--output_dir", default="processed-data", help="Directory where results will be saved (default: 'outputs')")
	args = parser.parse_args()

	# load esm model
	IF_model_name   = "data-files/esm_if1_gvp4_t16_142M_UR50.pt"
	model, alphabet = esm.pretrained.load_model_and_alphabet(IF_model_name)
	model.eval().to("cpu").requires_grad_(False)

	# load network node information and prepare set of commands to be run with multiprocessing
	nodes_df = pd.read_csv(args.input_node_file)

	nodes_df = nodes_df.head(10) # TESTING PURPOSES ONLY

	# all protein structure predictions from EBI for S288C contain a single chain with name A
	chainID = "A"

	# create cagiada class object list for each row in nodes_df
	#cagiada_jobs = nodes_df.apply(lambda row: cagiada("data-files/AF-"+row["UniProtKB-AC"]+"-F1-model_v4.pdb", chainID, row["UniProtKB-AC"]))

	test_object1 = cagiada("data-files/AF-P35725-F1-model_v4.pdb", "A", "P35725")
	test_object2 = cagiada("data-files/AF-P15700-F1-model_v4.pdb", "A", "P15700")

	predict_dG(test_object1, args.output_dir, model, alphabet)
	predict_dG(test_object2, args.output_dir, model, alphabet)

# execute main when run from command line
if __name__ == "__main__":
	main()

