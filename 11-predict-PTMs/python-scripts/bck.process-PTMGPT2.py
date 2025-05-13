import pandas as pd
import ast
import argparse

def extract_positive_sites(raw_result_str):
    try:
        result_dict = ast.literal_eval(raw_result_str)
        results = result_dict.get("Results", [])

        # find all positions with 'POSITIVE'
        positive_sites = [str(pos) for item in results for pos, label in item.items() if label == "POSITIVE"]
        return ",".join(positive_sites)
    except Exception as e:
        return ""

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    args = parser.parse_args()

    raw = pd.read_csv(args.input_file)
    raw["list_of_residue_IDs"] = raw["raw_result"].apply(extract_positive_sites)
    clean = raw[["node", "model", "list_of_residue_IDs"]]
    clean.to_csv(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-PTMGPT2-predictions-processed.csv", index=False)

if __name__ == "__main__":

    main()
