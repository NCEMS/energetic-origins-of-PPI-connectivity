# load modules
import sys
from Bio import SeqIO
import networkx as nx
import numpy as np
import pandas as pd
from datetime import datetime
import metapredict as meta
import matplotlib.pyplot as plt

# function to extract first name for a node
def extract_primary_name(node):

	return node.split(";")[0] # take the first part of the name before the semi-colon

# Function to compute the fraction of residues > 0.5
def fraction_disordered(predictions):

	if predictions is None:  # Handle missing values

		return None

	return (predictions > 0.5).sum() / len(predictions)

# label applied to outputs as a prefix
output_label = "0_"

# threshold for determining if a protein is disordered
# Note well, this should be swapped for a threshold calculated from DisProt
protein_disorder_cutoff = 0.50

# edges file path
#edge_file       = "../data-files/The_Yeast_Interactome_edges.csv"
edge_file       = sys.argv[1]

# yeast fasta file path
# from http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/
s288c_fasta     = "../data-files/orf_trans.fasta"

# load S288C fasta seqs as a dictionary using Biopython
s288c_seqs      = SeqIO.to_dict(SeqIO.parse(s288c_fasta, "fasta"))

# load the file DisProt_release_2024_12_with_ambiguous_evidences.tsv
disprot         = "../data-files/DisProt_release_2024_12_with_ambiguous_evidences.tsv"
disprot_df      = pd.read_csv(disprot, sep="\t")

# load edges as dataframe
edges_df  = pd.read_csv(edge_file)

edges_df.to_csv("../processed-data/"+output_label+"network_edges.csv", columns = ["source", "target"], index=False)

## EXTRACT NODES

# make a single-column DataFrame containing the gene names of all unique nodes; should be 3927 for Yeast Interactome
nodes_df     = pd.DataFrame(pd.unique(edges_df[['source', 'target']].values.ravel()), columns = ['node'])

# check to see how many nodes lack a sequence when using the "raw" Yeast Interactome names
nodes_df['has_verified_sequence'] = nodes_df['node'].isin(s288c_seqs.keys())

print ("nodes that cannot be mapped to a sequence by raw name")
print (nodes_df[nodes_df['has_verified_sequence'] == False].info())
print (nodes_df[nodes_df['has_verified_sequence'] == False], '\n')

# rename nodes with compound names - stop gap measure, should be investigated more later (just takes the first name)
nodes_df["primary_node"] = nodes_df["node"].apply(extract_primary_name)

## ADD SEQUENCE INFORMATION

# determine how many nodes have a sequence in s288c_seqs reference fasta
nodes_df['has_verified_sequence'] = nodes_df['primary_node'].isin(s288c_seqs.keys())

print ("nodes that cannot be mapped to a sequence by first name in list")
print (nodes_df[nodes_df['has_verified_sequence'] == False].info(), '\n')
print (nodes_df[nodes_df['has_verified_sequence'] == False], '\n')

# insert sequence if available, otherwise insert None
# this line also strips off the * that indicates the stop codon location in the fasta file
nodes_df["sequence"] = nodes_df["primary_node"].apply(lambda x: str(s288c_seqs[x].seq).rstrip("*") if x in s288c_seqs else None)

print ("The following nodes_df have no sequence information:")
print (nodes_df[nodes_df['has_verified_sequence'] == False].info())
print (nodes_df[nodes_df['has_verified_sequence'] == False])

## RUN METAPREDICT

# create dictionary with node names (gene identifiers) as keys and sequences as values
# this is the input to metapredict (needs to be cleaned up a bit)
map_nodes_to_seq = dict(zip(nodes_df["primary_node"], nodes_df["sequence"]))

# create a clean version that does not include None values
map_nodes_to_seq_clean = {k: v for k, v in map_nodes_to_seq.items() if v is not None}

# perform disoder predictions with metapredict
disorder_predictions = meta.predict_disorder(map_nodes_to_seq_clean)

# add metapredict output to dataframe
nodes_df["disorder_predictions"] = nodes_df["primary_node"].map(lambda x: disorder_predictions[x][1] if x in disorder_predictions else None)

# calculate the fraction of residues in each protein predicted to be disordered
nodes_df["disorder_fraction"] = nodes_df["disorder_predictions"].apply(fraction_disordered)

# create binary classification of "disordered or not"
nodes_df["is_disordered"] = nodes_df["disorder_fraction"].apply(lambda x: 1 if x >= protein_disorder_cutoff else 0)

## NETWORK CENTRALITY CALCULATIONS

# ensure the correct column names are present
if "source" not in edges_df.columns or "target" not in edges_df.columns:
    print("Column names are incorrect. Available columns:", edges_df.columns)

# create an undirected graph
interactome_graph = nx.from_pandas_edgelist(edges_df, "source", "target")

# compute centrality metrics

# a
degree_centrality = nx.degree_centrality(interactome_graph)
nodes_df["degree_centrality"] = nodes_df["node"].map(degree_centrality)

# b
betweenness_centrality = nx.betweenness_centrality(interactome_graph)
nodes_df["betweenness_centrality"] = nodes_df["node"].map(betweenness_centrality)

# c
eigenvector_centrality = nx.eigenvector_centrality(interactome_graph, max_iter=1000)
nodes_df["eigenvector_centrality"] = nodes_df["node"].map(eigenvector_centrality)

# d
closeness_centrality = nx.closeness_centrality(interactome_graph)
nodes_df["closeness_centrality"] = nodes_df["node"].map(closeness_centrality)

# e
load_centrality = nx.load_centrality(interactome_graph)
nodes_df["load_centrality"] = nodes_df["node"].map(load_centrality)

# f
pagerank = nx.pagerank(interactome_graph)
nodes_df["pagerank"] = nodes_df["node"].map(pagerank)

# g
k_shell = nx.core_number(interactome_graph)  # k-shell decomposition
nodes_df["k_shell"] = nodes_df["node"].map(k_shell)

# save the results to file; omit the array of per-residue disorder predictions and the amino acid sequence
nodes_df.drop(columns=["disorder_predictions", "sequence"]).to_csv("../processed-data/"+output_label+"network_nodes_with_annotation.csv", index=False)


