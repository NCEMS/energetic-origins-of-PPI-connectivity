import os, sys
import argparse
import pandas as pd
import numpy as np
import pint
import pint_pandas
from tqdm import tqdm

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--domain_data")
    parser.add_argument("--nodes")
    args = parser.parse_args()

    # read in the nodes_df information
    nodes_df = pd.read_pickle(args.nodes)
    uniprot_ids = set(nodes_df["UniProtKB-AC"])

    with open(args.domain_data, "r") as f:
        total_lines = sum(1 for _ in f)

    output_file = f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-domain-annotations.csv"

    with open(args.domain_data, "r") as infile, open(output_file, "w") as outfile:
        for line in tqdm(infile, total=total_lines, desc="Filtering InterPro", unit="lines"):
            if not line.strip():
                continue
            uniprot = line.split("\t")[0]
            if uniprot in uniprot_ids:
                outfile.write(line)

if __name__ == "__main__":
    main()
