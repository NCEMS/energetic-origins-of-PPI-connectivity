import pandas as pd
import ast
import argparse

def extract_positive_sites(raw_result_str):

    try:
        result_dict = ast.literal_eval(raw_result_str)
        results = result_dict.get("Results", [])
        positive_sites = [str(pos) for item in results for pos, label in item.items() if label == "POSITIVE"]
        return ",".join(positive_sites)
    except Exception:
        return ""

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    args = parser.parse_args()

    # read in the nodes_df from the previous step in the pipeline
    nodes_df = pd.read_pickle(args.nodes)

    # load raw prediction results
    raw = pd.read_csv(args.input_file)

    # extract POSITIVE residue positions
    raw["list_of_residue_IDs"] = raw["raw_result"].apply(extract_positive_sites)

    # pivot to wide format: rows = node, columns = model, values = comma-separated residue IDs
    wide = raw.pivot(index="node", columns="model", values="list_of_residue_IDs").reset_index()

    # save interstitial output file
    output_file = f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-PTMGPT2-predictions-processed.csv"
    wide.to_csv(output_file, index=False)

    # perform merge into nodes_df and save the output from this pipeline step
    nodes_df = nodes_df.merge(wide, on="node", how="left")

    # save the result
    nodes_df.to_pickle(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife-expr-speed-PTMGPT2.pkl")

if __name__ == "__main__":
    main()
