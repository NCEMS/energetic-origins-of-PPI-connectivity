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


def parse_arguments():
	"""
	Parses command-line arguments using argparse.

	Returns:
		args: Parsed command-line arguments.
	"""
	parser = argparse.ArgumentParser(description="Run stability predictions using ESM inverse folding.")
	parser.add_argument("--structure", required=True, help="Path to the input PDB structure file")
	parser.add_argument("--chain", required=True, help="PDB chain identifier (e.g., 'A')")
	parser.add_argument("--output_name", required=True, help="Output filename for the processed PDB")
	parser.add_argument("--output_dir", default="processed-data", help="Directory where results will be saved (default: 'outputs')")

	return parser.parse_args()

def main():

	# parse command-line arguments
	args = parse_arguments()

	# load esm model
	IF_model_name   = "esm_if1_gvp4_t16_142M_UR50.pt"
	model, alphabet = esm.pretrained.load_model_and_alphabet(IF_model_name)
	model.eval().to("cpu").requires_grad_(False)

	# load structure
	structure                            = load_structure(args.structure, args.chain)
	coords_structure, sequence_structure = extract_coords_from_structure(structure)

	# compute model probabilities
	prob_tokens        = run_model(coords_structure, sequence_structure, model, alphabet, chain_target=args.chain)
	aa_list, wt_scores = score_variants(sequence_structure, prob_tokens, alphabet)

	# Compute ΔG predictions
	a, b       = 0.10413378327743603, 0.6162549378400894  # Fitting parameters from manuscript
	dg_IF      = np.nansum(wt_scores)
	dg_kcalmol = a * dg_IF + b

	print(f"ΔG predicted (likelihood sum): {dg_IF}")
	print(f"ΔG predicted (kcal/mol): {dg_kcalmol}")

	# save results
	output_file = os.path.join(args.output_dir, f"{args.output_name}_dG_pos_scores_and_total.csv")
	output_df   = pd.DataFrame({'Residue': aa_list + ['dG_IF', 'dG_kcalmol'], 'score': wt_scores + [dg_IF, dg_kcalmol]})
	output_df.to_csv(output_file, sep=',', index=False)
	print(f"Results saved to {output_file}")

if __name__ == "__main__":
	main()

