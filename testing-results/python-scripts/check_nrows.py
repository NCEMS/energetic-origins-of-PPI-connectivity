import os, sys
import numpy as np
import pandas as pd
import pint
import pint_pandas

output_files = ["1-network-centrality/processed-data/0-CentralityCosDist-input.csv",
                "1-network-centrality/processed-data/0-s288c-nodes-centrality.csv",
                "2-sequence-parsing/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM-SignalP.csv",
                "2-sequence-parsing/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM.csv",
                "2-sequence-parsing/processed-data/0-s288c-nodes-centrality-seqs.csv",
                "3-uniprot-annotation/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt.csv",
                "4-idr-properties/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider.pkl",
                "4-idr-properties/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross.pkl",
                "4-idr-properties/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs.pkl",
                "5-dG-calculations/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada.pkl",
                "5-dG-calculations/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill.pkl",
                "5-dG-calculations/processed-data/0-s288c-temp.pkl",
                "6-Rosetta-scoring/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta.pkl",
                "7-FoldX-scoring/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX.pkl",
                "8-protein-half-life/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife.pkl",
                "9-protein-expression/processed-data/0-s288c-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife-expr.pkl",
                "flatten/processed-data/0-s288c-nodes-final-per-node.csv"]

for file in output_files:
    file_path = "../"+file
    ext = file_path.split(".")[-1]
    if ext == "csv":
        df = pd.read_csv(file_path)
    elif ext == "pkl":
        df = pd.read_pickle(file_path)
    else:
        print ("File extension not determined for", file_path)

    print (file_path, len(df))
