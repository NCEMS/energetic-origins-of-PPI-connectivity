import sys
import argparse
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np


def main():

    parser = argparse.ArgumentParser(description="Write metapredict input FASTA")
    parser.add_argument("--input_nodes", help="Path to input pd.DataFrame nodes")
    parser.add_argument("--output_fasta", help="Name/path of fasta file to be written")
    parser.add_argument(
        "--column",
        help="Column name corresponding to sequences on which IDR predictions will be run",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input_nodes)

    if args.column not in df.columns:
        print(
            f"{args.column} does not appear in the below list of column names.\n{list(df.columns)}"
        )
        sys.exit()

    df = df.dropna(subset=[args.column])
    df = df[df[args.column].astype(str).str.strip() != ""]

    records = []
    for _, row in df.iterrows():
        seq_str = str(row[args.column])
        if not seq_str or seq_str.lower() == "nan":
            continue
        records.append(SeqRecord(Seq(seq_str), id=str(row["node"]), description=""))

    if records:
        SeqIO.write(records, args.output_fasta, "fasta")
        print(f"Wrote {len(records)} sequences to {args.output_fasta}")
    else:
        print("No valid sequences found. No FASTA file written.")


if __name__ == "__main__":
    main()
