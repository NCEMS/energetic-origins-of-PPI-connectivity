import os
from pathlib import Path
import torch
import evaluate
import numpy as np
import pandas as pd
import argparse
from tqdm import tqdm
tqdm.pandas()
from torch.utils.data import DataLoader, Dataset
from transformers import DataCollatorForSeq2Seq
from transformers import (
    AutoTokenizer,
    GPT2LMHeadModel,
    TrainingArguments,
    Trainer,
    GPT2Config,
)
from sklearn.metrics import (
    average_precision_score,
    matthews_corrcoef,
    f1_score,
    precision_score,
    recall_score,
    balanced_accuracy_score,
)

# N.B., functions in this code were copied from https://github.com/pallucs/PTMGPT2/GPT2-Inference.ipynb


def find_subsequences(sequence: str, chars: list, left=10, right=10):

    subsequences = []
    length = len(sequence)

    # Iterate through the sequence to find the character
    for i, c in enumerate(sequence):
        if c in chars:

            # Calculate the start and end indices for the subsequence
            start = max(0, i - left)  # Ensure start is not less than 0
            end = min(
                length, i + right + 1
            )  # Ensure end does not exceed the sequence length

            # Append the subsequence to the list
            subsequences.append(
                {
                    "Seq": sequence[start:end],
                    "Pos": i + 1,
                    "text": f"<startoftext>SEQUENCE:{sequence[start:end]}\nLABEL:",
                }
            )
    return subsequences


def read_fasta(file_path):
    """
    Reads a FASTA file and returns a dictionary with sequence identifiers as keys
    and sequences as values.

    :param file_path: str, path to the FASTA file
    :return: dict, dictionary with sequence IDs as keys and sequences as values
    """
    sequences = {}
    sequence_id = None
    sequence_data = []

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_id is not None:
                    sequences[sequence_id] = "".join(sequence_data)
                sequence_id = line[1:]
                sequence_data = []
            else:
                sequence_data.append(line)

        # Add the last sequence
        if sequence_id is not None:
            sequences[sequence_id] = "".join(sequence_data)

    return sequences


def load_model(mdl_pth):
    """
    Loads a pre-trained GPT-2 model from the specified path.

    :param mdl_pth: str, path to the model directory.
    :return: GPT2LMHeadModel, the loaded GPT-2 model in evaluation mode on the CPU.
    """
    model_config = GPT2Config.from_pretrained(mdl_pth)
    model = GPT2LMHeadModel.from_pretrained(
        mdl_pth, config=model_config, ignore_mismatched_sizes=True
    )
    return model.cpu().eval()


def tokenize(sub_sequences, tokenizer):
    """
    Tokenizes the given subsequences using the specified tokenizer.

    :param sub_sequences: list of dicts, each containing a 'text' field with the subsequence to tokenize.
    :param tokenizer: AutoTokenizer, the tokenizer to use for tokenizing the subsequences.
    :return: dict, the tokenized subsequences with padding applied.
    """
    sub_sequences = [x["text"] for x in sub_sequences]
    encoded = tokenizer(sub_sequences, return_tensors="pt", padding="longest")
    return encoded


def inference(input_seq, tokenizer_pth, model_pth, chars: list):
    """
    Performs inference on the input sequence using a specified tokenizer and model, and extracts labels.

    :param input_seq: str, the input sequence to process.
    :param tokenizer_pth: str, path to the tokenizer directory.
    :param model_pth: str, path to the model directory.
    :param chars: list of str, characters to find subsequences for.
    :return: dict, a JSON-like dictionary containing the input sequence, model type, and labeled results.
    """
    tokenizer = AutoTokenizer.from_pretrained(tokenizer_pth, padding_side="left")
    model = load_model(model_pth)
    sub_sequences = find_subsequences(input_seq, chars=chars)
    if not sub_sequences:
        return {"Sequence": input_seq, "Type": model_pth, "Results": []}
    inputs_encode = tokenize(sub_sequences=sub_sequences, tokenizer=tokenizer)
    predicted = model.generate(
        inputs_encode["input_ids"],
        attention_mask=inputs_encode["attention_mask"],
        do_sample=False,
        top_k=50,
        max_new_tokens=2,
        top_p=0.15,
        temperature=0.1,
        num_return_sequences=0,
        pad_token_id=50259,
    )
    predicted_text = tokenizer.batch_decode(predicted, skip_special_tokens=True)
    predicted_labels = [x.split("LABEL:")[-1] for x in predicted_text]
    json_results = {"Sequence": input_seq, "Type": model_pth, "Results": []}
    for label, sub_seq in zip(predicted_labels, sub_sequences):
        json_results["Results"].append({sub_seq["Pos"]: label})
    return json_results


def run_inference_on_row(row, model_path, tokenizer_path, target_residues, sequence_column):

    result = inference(
        input_seq=row[sequence_column],
        tokenizer_pth=str(tokenizer_path),
        model_pth=str(model_path),
        chars=target_residues,
    )

    return {
        "node": row["node"],
        "model": model_path.name,
        "residue": target_residues,
        "raw_result": result,
    }


def main():

    # setup command-line arguments
    parser = argparse.ArgumentParser(
        description="Predict post-translational modifications"
    )
    parser.add_argument(
        "--nodes",
        help="Path to the nodes CSV file",
    )
    parser.add_argument("--output_prefix", help="Prefix for output files")
    parser.add_argument(
        "--output_dir", default="processed-data", help="Output directory"
    )
    parser.add_argument(
        "--organism_tag", default="s288c", help="Tag to label the organism for this run"
    )
    parser.add_argument("--gpt_model_path", help="Path to directory containing models")
    parser.add_argument(
        "--tokenizer_path",
        help="Path to the Tokenizer/ directory from PTMGPT2 GitHub repo",
    )
    args = parser.parse_args()

    # load the nodes_df from the previous step
    nodes_df = pd.read_pickle(args.nodes)

    # name of the column with sequence information to use
    sequence_column = "signalP_trimmed_sequence_x"

    # take only the columns required for inference
    inference_df = nodes_df[["node", sequence_column]].copy()

    # make a small subset of inference_df for testing purposes
    inference_df = inference_df.head(10)

    # a complete list of the models available from PTMGPT2; downloaded from:
    # Part 1. https://doi.org/10.5281/zenodo.11371883
    # Part 2. https://zenodo.org/records/11362322
    model_list = [
        "Acetylation (K)",
        "Amidation (V)",
        "Formylation (K)",
        "Glutarylation (K)",
        "Glutathionylation (C)",
        "Hydroxylation (K)",
        "Hydroxylation (P)",
        "Malonylation (K)",
        "Methylation (K)",
        "Methylation (R)",
        "N-linked Glycosylation (N)",
        "O-linked Glycosylation (S,T)",
        "Phosphorylation (S,T)",
        "Phosphorylation (Y)",
        "Succinylation (K)",
        "Sumoylation (K)",
        "S-nitrosylation (C)",
        "S-palmitoylation (C)",
        "Ubiquitination (K)",
    ]

    # dictionary mapping PTMGPT2 model names to residue types
    type_dict = {
        "Acetylation (K)": ["K"],
        "Amidation (V)": ["V"],
        "Formylation (K)": ["K"],
        "Glutarylation (K)": ["K"],
        "Glutathionylation (C)": ["C"],
        "Hydroxylation (K)": ["K"],
        "Hydroxylation (P)": ["P"],
        "Malonylation (K)": ["K"],
        "Methylation (K)": ["K"],
        "Methylation (R)": ["R"],
        "N-linked Glycosylation (N)": ["N"],
        "O-linked Glycosylation (S,T)": ["S", "T"],
        "Phosphorylation (S,T)": ["S", "T"],
        "Phosphorylation (Y)": ["Y"],
        "Succinylation (K)": ["K"],
        "Sumoylation (K)": ["K"],
        "S-nitrosylation (C)": ["C"],
        "S-palmitoylation (C)": ["C"],
        "Ubiquitination (K)": ["K"],
    }

    # path to the tokenizer to use
    tokenizer_path = Path(args.tokenizer_path)

    # list to hold all inference results
    all_results = []

    # run each model in series
    for model_name in model_list:

        print (f"Running predictions for {model_name}")

        model_path = Path(args.gpt_model_path) / model_name

        target_residues = type_dict[model_name]

        results = inference_df.progress_apply(
            run_inference_on_row,
            axis=1,
            args=(model_path, tokenizer_path, target_residues, sequence_column),
        )

        all_results.extend(results.tolist())

        print ("\n")

    # convert to a pd.DataFrame
    results_df = pd.DataFrame(all_results)

    # save output from all PTMGPT2 predictions to an intermediate .csv file
    results_df.to_csv(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-PTMGPT2-predictions.csv", index=False)


if __name__ == "__main__":

    main()
