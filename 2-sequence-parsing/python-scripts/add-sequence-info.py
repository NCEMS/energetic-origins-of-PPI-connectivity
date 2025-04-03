import pandas as pd
from Bio import SeqIO
import typing

def main():

    parser = argparse.ArgumentParser(description="Add sequences and additional annotations to PPI network.")
    parser.add_argument(
        "--nodes",
        default="../1-network-centrality/processed-data/0_nodes-centrality.csv",
        help="Path to the nodes CSV file",
    )
    parser.add_argument(
        "--fasta",
        default="data-files/orf_trans.fasta",
        help="Path to open reading frame FASTA file",
    )
    parser.add_argument(
        "--seq_preds",
        default="DeepTMHMM-runs/s288c-results/all-predictions-s288c.3line",
        help="Path to 3line format prediction file from DeepTMHMM",
    )
    args = parser.parse_args()

    # load the nodes csv file
    nodes_df = pd.read_csv(args.nodes)

    # load sequence information from file
    #fasta_seqs = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

    # insert sequence information into nodes_df
    nodes_df = add_sequences(nodes_df, args.fasta)

if __name__ == "__main__":

    main()
