import sys
import argparse
from Bio import SeqIO
import networkx as nx
import pandas as pd
import metapredict as meta

# function to extract primary name for a node
def extract_primary_name(node):

	# Take the first part of the name before the semi-colon
	return node.split(";")[0]

# function to compute the fraction of residues > 0.5
def fraction_disordered(predictions):

	# this is the threshold mentioned in https://www.biorxiv.org/content/10.1101/2024.11.05.622168v1
	per_residue_disorder_cutoff = 0.5

	# handle cases when no metapredict prediction is obtained
	if predictions is None:  
		return None

	# otherwise return fraction of residues with score > >0.5
	return (predictions > per_residue_disorder_cutoff).sum() / len(predictions)

# function to process nodes
def process_nodes(edges_df, s288c_seqs):

	# load nodes as a DataFrame
	nodes_df = pd.DataFrame(pd.unique(edges_df[['source', 'target']].values.ravel()), columns=["node"])

	# check to see which nodes have a sequence in s288c_seqs
	#nodes_df['has_verified_sequence'] = nodes_df['node'].isin(s288c_seqs.keys())

	# define "primary_node" as the first gene name within semi-colon delimited lists
	#nodes_df["primary_node"] = nodes_df["node"].apply(extract_primary_name)

	# remove rows in which "node" is a semi-colon separated list
	nodes_df = nodes_df[~nodes_df["node"].str.contains(";", na=False)]

	# use "node" as keys to check for sequences
	nodes_df['has_verified_sequence'] = nodes_df['node'].isin(s288c_seqs.keys())

	# grab sequences from s288c_seqs and add as a column in the DataFrame
	nodes_df["sequence"] = nodes_df["node"].apply(lambda x: str(s288c_seqs[x].seq).rstrip("*") if x in s288c_seqs else None)

	# add additional useful node label information from UniProt (required by Cagiada stability analyses 
	# to find the correct AF2 structure prediction to use for a given gene name)
	column_names = ["UniProtKB-AC", "ID_type", "ID"]
	cross_df     = pd.read_csv("data-files/YEAST_559292_idmapping.dat", names=column_names, sep="\t")
	cross_df     = cross_df[cross_df["ID_type"] == "Gene_OrderedLocusName"]
	nodes_df     = nodes_df.merge(cross_df[["UniProtKB-AC", "ID"]], left_on="node", right_on="ID", how="left")

	# return the updated DataFrame
	return nodes_df

# function to run metapredict on node sequences
def predict_disorder(nodes_df):

	# create dictionary in format needed by metapredict
	map_nodes_to_seq = {k: v for k, v in zip(nodes_df["node"], nodes_df["sequence"]) if v is not None}

	# run metapredict
	disorder_predictions = meta.predict_disorder(map_nodes_to_seq)

	# add dictionary information to DataFrame
	nodes_df["disorder_predictions"] = nodes_df["node"].map(lambda x: disorder_predictions[x][1] if x in disorder_predictions else None)

	# add column with protein's fraction of disordered residues
	nodes_df["disorder_fraction"] = nodes_df["disorder_predictions"].apply(fraction_disordered)

	# add binary classification of protein as disordered/not disordered
	nodes_df["is_disordered"] = nodes_df["disorder_fraction"].apply(lambda x: 1 if x >= 0.50 else 0)

	# return the updated DataFrame
	return nodes_df

# function to compute network centrality measures
def compute_centrality(edges_df, nodes_df):

	# make networkx style graph
	interactome_graph = nx.from_pandas_edgelist(edges_df, "source", "target")

	# perform centrality calculations
	centrality_measures = {"degree_centrality"     : nx.degree_centrality(interactome_graph),
		               "betweenness_centrality": nx.betweenness_centrality(interactome_graph),
                               "eigenvector_centrality": nx.eigenvector_centrality(interactome_graph, max_iter=1000),
                               "closeness_centrality"  : nx.closeness_centrality(interactome_graph),
                               "load_centrality"       : nx.load_centrality(interactome_graph),
                               "pagerank"              : nx.pagerank(interactome_graph),
                               "k_shell"               : nx.core_number(interactome_graph)
                              }

	# add per-node information to the DataFrame
	for key, values in centrality_measures.items():
		nodes_df[key] = nodes_df["node"].map(values)

	# return the updated DataFrame
	return nodes_df

# main function
def main():

	# setup arguments from the command line
	parser     = argparse.ArgumentParser(description="Process Yeast interactome network.")
	parser.add_argument("--edges", default="data-files/The_Yeast_Interactome_edges.csv", help="Path to the edges CSV file")
	parser.add_argument("--fasta", default="data-files/orf_trans.fasta", help="Path to the yeast protein FASTA file")
	parser.add_argument("--output_prefix", default="0_", help="Prefix for output files")
	parser.add_argument("--output_dir", default="processed-data", help="Output directory")
	#parser.add_argument("--disprot", required=True, help="Path to the DisProt TSV file")
	args       = parser.parse_args()

	# load input data
	s288c_seqs = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
	edges_df   = pd.read_csv(args.edges)
	nodes_df   = process_nodes(edges_df, s288c_seqs)

	# predict disorder using metapredict
	nodes_df   = predict_disorder(nodes_df)

	# compute network centrality measures
	nodes_df   = compute_centrality(edges_df, nodes_df)

	# add some additional identifier information; UniProtID == SGD ID == "standard gene name"

	# save outputs
	edges_df.to_csv(f"{args.output_dir}/{args.output_prefix}network_edges.csv", columns=["source", "target"], index=False)
	nodes_df.drop(columns=["disorder_predictions"]).to_csv(f"{args.output_dir}/{args.output_prefix}network_nodes_with_annotation.csv", index=False)

	# print a "DONE" statement
	print(f"Processing complete. Output saved to {args.output_dir}")

# entry point
if __name__ == "__main__":

	# run the main function
	main()

