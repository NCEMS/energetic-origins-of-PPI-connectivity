import os
from pathlib import Path
import pandas as pd
import argparse
import multiprocessing as mp
from tqdm import tqdm
from transformers import AutoTokenizer, GPT2LMHeadModel, GPT2Config

tqdm.pandas()

# Note well - this program assumes that you have two CUDA-enabled GPUs with IDs {0, 1}

def find_subsequences(sequence: str, chars: list, left=10, right=10):
    subsequences = []
    length = len(sequence)
    for i, c in enumerate(sequence):
        if c in chars:
            start = max(0, i - left)
            end = min(length, i + right + 1)
            subsequences.append({
                "Seq": sequence[start:end],
                "Pos": i + 1,
                "text": f"<startoftext>SEQUENCE:{sequence[start:end]}\nLABEL:"
            })
    return subsequences


def load_model(mdl_pth):
    config = GPT2Config.from_pretrained(mdl_pth)
    model = GPT2LMHeadModel.from_pretrained(mdl_pth, config=config, ignore_mismatched_sizes=True)
    return model.cuda().eval()  # Move to GPU


def tokenize(sub_sequences, tokenizer):
    sub_sequences = [x["text"] for x in sub_sequences]
    return tokenizer(sub_sequences, return_tensors="pt", padding="longest").to("cuda")


def inference(input_seq, tokenizer, model, chars: list, model_name: str):
    sub_sequences = find_subsequences(input_seq, chars)
    if not sub_sequences:
        return {"Sequence": input_seq, "Type": model_name, "Results": []}

    inputs_encode = tokenize(sub_sequences, tokenizer)
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
    return {
        "Sequence": input_seq,
        "Type": model_name,
        "Results": [{sub["Pos"]: label} for sub, label in zip(sub_sequences, predicted_labels)]
    }


def run_inference_on_row(row, model_path, tokenizer_path, target_residues, sequence_column):
    tokenizer = AutoTokenizer.from_pretrained(tokenizer_path, padding_side="left")
    model = load_model(model_path)
    return {
        "node": row["node"],
        "model": model_path.name,
        "residue": target_residues,
        "raw_result": inference(
            input_seq=row[sequence_column],
            tokenizer=tokenizer,
            model=model,
            chars=target_residues,
            model_name=model_path.name
        )
    }


def run_inference_gpu(gpu_id, df, model_list, type_dict, tokenizer_path, sequence_column, model_root_path):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    tokenizer_path = Path(tokenizer_path)

    results_all = []
    for model_name in model_list:
        model_path = Path(model_root_path) / model_name
        target_residues = type_dict[model_name]
        print(f"[GPU {gpu_id}] Starting model: {model_name}")

        results = df.progress_apply(
            run_inference_on_row,
            axis=1,
            args=(model_path, tokenizer_path, target_residues, sequence_column)
        )
        results_all.extend(results.tolist())
        print(f"[GPU {gpu_id}] Done: {model_name}")

    return results_all


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--gpt_model_path", required=True)
    parser.add_argument("--tokenizer_path", required=True)
    args = parser.parse_args()

    model_list = [
        "Acetylation (K)", "Amidation (V)", "Formylation (K)", "Glutarylation (K)", "Glutathionylation (C)",
        "Hydroxylation (K)", "Hydroxylation (P)", "Malonylation (K)", "Methylation (K)", "Methylation (R)",
        "N-linked Glycosylation (N)", "O-linked Glycosylation (S,T)", "Phosphorylation (S,T)", "Phosphorylation (Y)",
        "Succinylation (K)", "Sumoylation (K)", "S-nitrosylation (C)", "S-palmitoylation (C)", "Ubiquitination (K)"
    ]

    type_dict = {
        "Acetylation (K)": ["K"], "Amidation (V)": ["V"], "Formylation (K)": ["K"],
        "Glutarylation (K)": ["K"], "Glutathionylation (C)": ["C"], "Hydroxylation (K)": ["K"],
        "Hydroxylation (P)": ["P"], "Malonylation (K)": ["K"], "Methylation (K)": ["K"],
        "Methylation (R)": ["R"], "N-linked Glycosylation (N)": ["N"], "O-linked Glycosylation (S,T)": ["S", "T"],
        "Phosphorylation (S,T)": ["S", "T"], "Phosphorylation (Y)": ["Y"], "Succinylation (K)": ["K"],
        "Sumoylation (K)": ["K"], "S-nitrosylation (C)": ["C"], "S-palmitoylation (C)": ["C"], "Ubiquitination (K)": ["K"]
    }

    sequence_column = "signalP_trimmed_sequence_x"
    nodes_df = pd.read_pickle(args.nodes)

    # take small subset of nodes for testing purposes
    #nodes_df = nodes_df.head(10)

    df = nodes_df[["node", sequence_column]].copy()
    df = df[df[sequence_column].apply(lambda x: isinstance(x, str) and len(x.strip()) > 0)].copy()
    print (f"Dropping {len(nodes_df) - len(df)} rows with missing or invalid sequences.")

    df1 = df.iloc[:len(df)//2].reset_index(drop=True)
    df2 = df.iloc[len(df)//2:].reset_index(drop=True)

    with mp.get_context("spawn").Pool(2) as pool:
        results_split = pool.starmap(
            run_inference_gpu,
            [
                (0, df1, model_list, type_dict, args.tokenizer_path, sequence_column, args.gpt_model_path),
                (1, df2, model_list, type_dict, args.tokenizer_path, sequence_column, args.gpt_model_path),
            ]
        )

    all_results = [item for sublist in results_split for item in sublist]
    results_df = pd.DataFrame(all_results)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{args.output_prefix}-{args.organism_tag}-PTMGPT2-predictions.csv"
    results_df.to_csv(output_file, index=False)
    print (f"Saved: {output_file}")


if __name__ == "__main__":
    main()
