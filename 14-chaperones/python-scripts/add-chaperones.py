import os, sys
import pandas as pd
import pint
import pint_pandas
import argparse
from collections import defaultdict

def get_chaperone_metadata(prot, interaction_map, chaperone_info):
    chaps = interaction_map.get(prot, set())
    return [
        {"ORF": c, **chaperone_info[c]}
        for c in chaps if c in chaperone_info
    ]

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--input_nodes")
    parser.add_argument("--chap_data")
    parser.add_argument("--edges")
    args = parser.parse_args()

    ## read input files

    # read in chaperone information
    chap_df = pd.read_csv(args.chap_data, sep="\t")

    # read in edges; we only need the "source" and "target" information
    edges_df = pd.read_csv(args.edges, usecols=["source", "target"])

    # read in annotated nodes
    nodes_df = pd.read_pickle(args.input_nodes)

    ## get mappings

    # make lookup for metadata
    chaperone_info = chap_df.set_index("ORF")[["CCo Family", "Name description"]].to_dict(orient="index")

    # build map of protein -> set of interacting chaperones and add to nodes_df
    interaction_map = defaultdict(set)
    for _, row in edges_df.iterrows():
        src, tgt = row["source"], row["target"]
    
        if tgt in chaperone_info:
            interaction_map[src].add(tgt)
        if src in chaperone_info:
            interaction_map[tgt].add(src)

    nodes_df["interacting_chaperones"] = nodes_df["node"].apply(
        lambda prot: list(interaction_map.get(prot, []))
    )

    # add list of chaperone metadata
    #nodes_df["interacting_chaperone_info"] = nodes_df["node"].apply(get_chaperone_metadata)
    nodes_df["interacting_chaperone_info"] = nodes_df["node"].apply(
        lambda prot: get_chaperone_metadata(prot, interaction_map, chaperone_info)
    )

    # save the output to file
    nodes_df.to_pickle(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife-expr-speed-PTMGPT2-LiPMS-entanglement-chaperones.pkl")

if __name__ == "__main__":
    main()
