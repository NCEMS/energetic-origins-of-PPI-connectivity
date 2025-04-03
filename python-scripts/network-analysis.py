import os, sys

sys.path.append("CentralityCosDist")
from centralitycosdist import CentralityCosDist
import argparse
from Bio import SeqIO
import networkx as nx
import pandas as pd
import metapredict as meta
import numpy as np
import typing
from typing import Optional
from typing import Union
from typing import List
import pytest


def extract_primary_name(node: str) -> str:
    """
    Returns the first name from a semi-colon delimited list

    Args:
        node (str): the original name of the node, possibly a semi-colon delimited str

    Returns:
        str: the first name from the semi-colon delimited list
    """

    return node.split(";")[0]


def fraction_disordered(predictions: np.ndarray) -> Optional[float]:
    """
    Returns the fraction of residues above the disorder threshold of 0.5

    Args:
        predictions (np.ndarray): Array of prediction values, one per residue in the protein

    Returns:
        float or None: proportion of residues in the protein that are predicted to be disordered or None if no prediction was made for this protein
    """

    # this is the threshold mentioned in https://www.biorxiv.org/content/10.1101/2024.11.05.622168v1
    per_residue_disorder_cutoff = 0.5

    if predictions is None:
        return None

    return (predictions > per_residue_disorder_cutoff).sum() / len(predictions)


def count_IDRs(
    arr: Union[np.ndarray[np.float64], float, None],
    threshold: float = 0.5,
    min_length: int = 30,
) -> int:
    """
    Counts regions in an array where values are > threshold
    for at least min_length consecutive positions

    Args:
        arr (np.ndarray or float or None): The input array or a NaN placeholder.
        threshold (float): The threshold value to define a residue as disordered
        min_length (int): Minimum number of residues a region must be to count towards total

    Returns:
        int: Number of regions satisfying the threshold and length criteria
    """

    # catch instances of NaN in the array/values
    if arr is None or isinstance(arr, float) and np.isnan(arr):
        return 0

    # create a boolean array: True where value > threshold
    mask = arr > threshold
    count = 0
    current_run = 0

    for val in mask:
        if val:
            current_run += 1
        else:
            if current_run >= min_length:
                count += 1
            current_run = 0

    # check if the last run reached the threshold
    if current_run >= min_length:
        count += 1

    return count


def process_nodes(
    edges_df: pd.DataFrame,
    s288c_seqs: dict[str, str],
    ID_mappings_path: str,
    structure_dir: str,
) -> pd.DataFrame:
    """
    Creates nodes DataFrame with node names, sequences, UniProt ID mappings, and paths to structures

    Args:
        edges_df (pd.DataFrame):
        s288c_seqs (dict[str, str]): dictionary mapping yeast gene identifiers to protein sequence
        ID_mappings_path (str): file path to the file containing ID mappings from SGD yeast gene IDs to UniProt IDs
        structure_dir (str): file path to the directory containing AlphaFold2 structure predictions

    Returns:
        pd.DataFrame: the initial nodes_df with certain annotations
    """

    # load nodes as a DataFrame
    nodes_df = pd.DataFrame(
        pd.unique(edges_df[["source", "target"]].values.ravel()), columns=["node"]
    )
    """
    # check to see which nodes have a sequence in s288c_seqs
    nodes_df['has_verified_sequence'] = nodes_df['node'].isin(s288c_seqs.keys())

    # define "primary_node" as the first gene name within semi-colon delimited lists
    nodes_df["primary_node"] = nodes_df["node"].apply(extract_primary_name)

    # remove rows in which "node" is a semi-colon separated list
    nodes_df = nodes_df[~nodes_df["node"].str.contains(";", na=False)]
    """
    # use "node" as keys to check for sequences
    nodes_df["has_verified_sequence"] = nodes_df["node"].isin(s288c_seqs.keys())

    # grab sequences from s288c_seqs and add as a column in the DataFrame
    nodes_df["sequence"] = nodes_df["node"].apply(
        lambda x: str(s288c_seqs[x].seq).rstrip("*") if x in s288c_seqs else None
    )

    # add additional useful node label information from UniProt (required by Cagiada stability analyses
    # to find the correct AF2 structure prediction to use for a given gene name)
    column_names = ["UniProtKB-AC", "ID_type", "ID"]
    cross_df = pd.read_csv(ID_mappings_path, names=column_names, sep="\t")
    cross_df = cross_df[cross_df["ID_type"] == "Gene_OrderedLocusName"]
    nodes_df = nodes_df.merge(
        cross_df[["UniProtKB-AC", "ID"]], left_on="node", right_on="ID", how="left"
    )

    # locate and add structures to dataframe
    nodes_df["structure_path"] = nodes_df["UniProtKB-AC"].apply(
        lambda id: f"{structure_dir}/AF-{id}-F1-model_v4.pdb"
    )

    # create the structure_exists column by checking if the file actually exists
    nodes_df["structure_exists"] = nodes_df["structure_path"].apply(
        lambda path: 1 if os.path.exists(path) else 0
    )
    nodes_df.loc[nodes_df["structure_exists"] == 0, "structure_path"] = None

    return nodes_df


def parse_model_output(lines: List[str]) -> dict[str, dict[str, str, str]]:
    """
    Takes list of lines output by DeepTMHMM and converts it to a dictionary

    Args:
        lines (List[str]): list of lines read in from a file


    Returns:
        dict[str, dict[str, str, str]]: dictionary; keys are gene IDs, values are dictionaries containing class_label, sequence, and a mask based on the DeepTMHMM output

    """

    parsed = {}
    i = 0

    expected_labels = ["TM", "SP", "GLOB", "SP+TM", "BETA"]

    while i < len(lines):

        header = lines[i].strip()
        seq = lines[i + 1].strip()
        mask = lines[i + 2].strip()
        i += 3

        # Extract gene_id and classification (e.g., SP, TM, GLOB)
        parts = header[1:].split("|")
        gene_id = parts[0].strip()
        class_label = parts[1].strip() if len(parts) > 1 else "UNKNOWN"
        if class_label not in expected_labels:
            print(
                "The classification of this protein took on an unexpected label:",
                gene_id,
                class_label,
            )
            sys.exit()

        parsed[gene_id] = {"class": class_label, "sequence": seq, "mask": mask}

    return parsed


def apply_model_results(
    df: pd.DataFrame, model_output_dict: dict[str, dict[str, str, str]]
) -> pd.DataFrame:
    """
    Adds DeepTMHMM information to nodes_df

    Args:
        df (pd.DataFrame): input nodes_df to be annotated
        model_output_dict (dict[str, dict[str, str, str]]): dictionary; output from parse_model_output

    Returns:
        pd.DataFrame: df annotated with DeepTMHMM information
    """

    def extract_mask(gene_id: str) -> dict[str, str]:
        return model_output_dict.get(gene_id, {}).get("mask", None)

    def trim_sequence(gene_id: str) -> Optional[str]:

        data = model_output_dict.get(gene_id)
        if not data:
            return None
        # if data['class'] == "TM":
        # 	return None
        # elif data['class'] == "SP+TM":
        # 	return None
        elif data["class"] in ["TM", "SP+TM", "BETA"]:
            return None
        elif data["class"] == "SP":
            # Remove the signal peptide: keep only residues where mask != 'S'
            trimmed = "".join(
                [aa for aa, m in zip(data["sequence"], data["mask"]) if m != "S"]
            )
        else:
            # GLOB or other cases, return original sequence
            trimmed = data["sequence"]

        return trimmed.rstrip("*")

    def extract_class(gene_id: str) -> dict[str, str]:
        return model_output_dict.get(gene_id, {}).get("class", None)

    # use apply with functions to update the input df
    df["mask"] = df["node"].apply(extract_mask)
    df["trimmed_sequence"] = df["node"].apply(trim_sequence)
    df["DeepTMHMM_class"] = df["node"].apply(extract_class)

    return df


def add_DeepTMHMM(df: pd.DataFrame, path_to_3line_file: str) -> pd.DataFrame:
    """
    Carries out steps to add DeepTMHMM information to df

    Args:
        df (pd.DataFrame): input nodes_df to be annotated
        path_to_3line_file (str): file path to DeepTMHMM output

    Returns:
        pd.DataFrame
    """
    # read in 3line file as a list of lines
    with open(path_to_3line_file) as f:
        lines = f.readlines()

    # parse the list of lines
    model_output_dict = parse_model_output(lines)

    # add information to the DataFrame
    df = apply_model_results(df, model_output_dict)

    return df


def predict_disorder(
    nodes_df: pd.DataFrame, disorder_threshold: float = 0.5
) -> pd.DataFrame:
    """
    Run metapredict on node protein sequences and add output to DataFrame

    Args:
        nodes_df (pd.DataFrame): input DataFrame to be annotated with IDR information

    Returns:
        pd.DataFrame
    """

    # create dictionary in format needed by metapredict
    map_nodes_to_seq = {
        k: v
        for k, v in zip(nodes_df["node"], nodes_df["trimmed_sequence"])
        if pd.notnull(v) and v.strip() != ""
    }

    # run metapredict
    disorder_predictions = meta.predict_disorder(map_nodes_to_seq)

    # add dictionary information to DataFrame
    nodes_df["disorder_predictions"] = nodes_df["node"].map(
        lambda x: disorder_predictions[x][1] if x in disorder_predictions else None
    )

    # add column with protein's fraction of disordered residues
    nodes_df["disorder_fraction"] = nodes_df["disorder_predictions"].apply(
        fraction_disordered
    )

    # add binary classification of protein as disordered/not disordered
    nodes_df["is_disordered"] = nodes_df["disorder_fraction"].apply(
        lambda x: 1 if x > disorder_threshold else 0
    )

    # return the updated DataFrame
    return nodes_df


def compute_centrality(
    edges_df: pd.DataFrame,
    nodes_df: pd.DataFrame,
    output_dir: str,
    output_prefix: str,
    seeds: Optional[List] = None,
) -> pd.DataFrame:
    """
    Compute network centrality measures on network

    Args:
        edges_df (pd.DataFrame): network edges
        nodes_df (pd.DataFrame): network nodes
        output_dir (str): output directory, used by add_CentralityCosDist
        output_prefix (str): output file prefix, used by add_CentralityCosDist

    Returns:
        pd.DataFrame
    """

    # make networkx style graph
    interactome_graph = nx.from_pandas_edgelist(edges_df, "source", "target")

    # perform centrality calculations
    centrality_measures = {
        "degree_centrality": nx.degree_centrality(interactome_graph),
        "betweenness_centrality": nx.betweenness_centrality(interactome_graph),
        "eigenvector_centrality": nx.eigenvector_centrality(
            interactome_graph, max_iter=1000
        ),
        "closeness_centrality": nx.closeness_centrality(interactome_graph),
        "load_centrality": nx.load_centrality(interactome_graph),
        "pagerank": nx.pagerank(interactome_graph),
        "k_shell": nx.core_number(interactome_graph),
    }

    # add per-node information to the DataFrame
    for key, values in centrality_measures.items():
        nodes_df[key] = nodes_df["node"].map(values)

    # compute CentralityCosDist; requires a file in a very specific format
    test_CentralityCosDist()
    nodes_df = add_CentralityCosDist(
        nodes_df,
        output_dir,
        output_prefix,
        ["node"]+list(centrality_measures.keys()),
        seeds=seeds,
    )

    # return the updated DataFrame
    return nodes_df


def add_CentralityCosDist(
    nodes_df: pd.DataFrame,
    output_dir: str,
    output_prefix: str,
    metrics_list: List[str],
    seeds: Optional[List] = None,
) -> pd.DataFrame:
    """
    Carries out CentralityCosDist calculations for the network input at a DataFrame

    Args:
        nodes_df (pd.DataFrame): DataFrame containing precomputed centrality metrics as required by CentralityCosDist algorithm
        output_dir (str): path to output data directory
        output_prefix (str): prefix to be appended to output file
        metrics_list (List[str]): list of the metric names to be extracted from nodes_df for cosine distance calculation
        seeds (Optional[List]): either a list of nodes to treat as seeds or None; if None, all nodes will be treated as seeds

    Returns:
        pd.DataFrame

    Note well: the CentralityCosDist program will interpret your input csv file in the following way:
               * the leftmost column will be treated as a unique identifier for nodes in the network; seeds must be cross-referenceable with this node list
               * all other columns in the DataFrame will be treated as centrality metrics and used in the cosine distance calculation
    """

    centralities_df = nodes_df[metrics_list]
    centralities_df = centralities_df.rename(columns={"node": "ID"})
    centralities_df.to_csv(
        f"{output_dir}/{output_prefix}CentralityCosDist_input.csv", index=False
    )

    # get seeds information

    # make set of the seeds to be considered
    if seeds == None:
        seeds = set(centralities_df["ID"].to_list())
    else:
        seeds = set(seeds)

    # make set of nodes
    nodes = set(centralities_df["ID"].to_list())

    # get list of seeds with centrality metric information
    seeds = list(nodes.intersection(seeds))

    # setup the algrithm
    algorithm = CentralityCosDist(
        Centrality_file=f"{output_dir}/{output_prefix}CentralityCosDist_input.csv"
    )

    # run the algorithm
    algorithm.run(seed_nodes=seeds)

    # extract results summary and rename for output DataFrame
    rank_df = algorithm.rank
    rank_df.name = "CentralityCosDist_rank"
    similarity_score_df = algorithm.similarity_score
    similarity_score_df.name = "CentralityCosDist_similarity_score"

    # add results to nodes_df
    nodes_df = pd.merge(nodes_df, rank_df, how="left", left_on="node", right_on="ID")
    nodes_df = pd.merge(
        nodes_df, similarity_score_df, how="left", left_on="node", right_on="ID"
    )

    return nodes_df


def test_CentralityCosDist():
    """
    Function to test whether or not results from CentralityCosDist match expectations
    Expected results are based on https://nilesh-iiita.github.io/CentralityCosDist/notebooks.html

    Args:
        None

    Returns:
        None
    """

    seeds = [
        "AT3G03900",
        "AT3G01850",
        "AT1G63290",
        "AT1G09100",
        "AT3G51840",
        "AT1G09770",
        "AT3G05530",
        "AT5G17310",
        "ATCG00480",
        "AT5G08670",
    ]
    nodes = pd.read_csv("test/network/data-files/Network_Centrality.csv")
    metrics_list = ["node",
        "Information_centrality",
        "Degree_centrality",
        "Betweenness_centrality",
        "Eigenvector_centrality",
        "Closeness_centrality",
        "clustering_coefficient",
        "Load_centrality",
        "Page_rank",
    ]
    result = add_CentralityCosDist(
        nodes, "test/network/processed-data", "test_", metrics_list, seeds=seeds
    )

    # expected scores from the CentralityCosDist documentation
    expected_similarity_score = {
        "AT3G03900": 0.984956,
        "AT1G09100": 0.984796,
        "AT3G51840": 0.983411,
        "AT3G05530": 0.981703,
        "AT5G17310": 0.977869,
        "ATCG00480": 0.975787,
        "AT5G08670": 0.973695,
        "AT5G08680": 0.973695,
        "AT5G08690": 0.971715,
        "AT5G19680": 0.970025,
    }

    # scores calculated in this run
    calculated_similarity_score = result.set_index("node")[
        "CentralityCosDist_similarity_score"
    ].to_dict()

    # tolerance for floating point comparisons
    tol = 1e-5

    for key, expected_val in expected_similarity_score.items():
        assert key in calculated_similarity_score, f"Missing key: {key}"
        assert (
            abs(calculated_similarity_score[key] - expected_val) < tol
        ), f"Reference and calculated values for {key} do not match"


def add_UniProt_info(nodes_df: pd.DataFrame, uniprot_data: str) -> pd.DataFrame:
    """
    Function that reads in a pre-processed annotation file from UniProt and adds selected information it to nodes_df

    Args:
        nodes_df (pd.DataFrame): nodes DataFrame so far; must contain mapped UniProt IDs for this to work (added by process_nodes)
        uniprot_data (str): path to the pre-processed UniProt information for the organism currently under study

    Returns:
        pd.DataFrame
    """

    uniprot_df = pd.read_csv(uniprot_data)

    to_add = [
        "PrimaryAccession",
        "ProteinName",
        "GO_terms",
        "localization_keywords",
        "parsed_functions",
        "parsed_PTMs",
    ]

    return nodes_df.merge(
        uniprot_df[to_add],
        left_on="UniProtKB-AC",
        right_on="PrimaryAccession",
        how="left",
    ).drop(columns=["PrimaryAccession"])

def gen_fasta(nodes_df:pd.DataFrame, out_file_path) -> None:

    """
    Takes in the nodes_df, which must contain "trimmed_sequence" column, and saves a fasta file with sequences on which to run predictions

    Args:
        nodes_df (pd.DataFrame): DataFrame containing trimmed_sequence information to be written to fasta format
        out_file_path (str): path to the fasta file to be written

    Returns:
        None, but writes a file to out_file_path
    """

    id_column = "node"
    seq_column = "trimmed_sequence"

    with open(out_file_path, "w") as f:
        for _, row in nodes_df.iterrows():
            seq = row[seq_column]
            if seq != "None":
                f.write(f">{row[id_column]}\n{sequence}\n")

def main():

    parser = argparse.ArgumentParser(description="Process Yeast interactome network.")
    parser.add_argument(
        "--edges",
        default="data-files/The_Yeast_Interactome_edges.csv",
        help="Path to the edges CSV file",
    )
    parser.add_argument(
        "--fasta",
        default="data-files/orf_trans.fasta",
        help="Path to the yeast protein FASTA file",
    )
    parser.add_argument(
        "--structure_dir",
        default="data-files/",
        help="Path to the directory containing AF2 structures for structure predictions",
    )
    parser.add_argument("--output_prefix", default="0_", help="Prefix for output files")
    parser.add_argument(
        "--output_dir", default="processed-data", help="Output directory"
    )
    parser.add_argument(
        "--seq_preds",
        default="DeepTMHMM-runs/s288c-results/all-predictions-s288c.3line",
        help="Path to 3line format prediction file from DeepTMHMM",
    )
    parser.add_argument(
        "--ID_mappings",
        default="data-files/YEAST_559292_idmapping.dat",
        help="Path to the UniProt ID mappings to be used",
    )
    parser.add_argument(
        "--uniprot_data",
        default="processed-data/uniprot_sprot-s288c.csv",
        help="Path to the UniProt data to be used for annotation",
    )
    # parser.add_argument("--disprot", required=True, help="Path to the DisProt TSV file")
    args = parser.parse_args()

    # load input data
    s288c_seqs = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
    edges_df = pd.read_csv(args.edges)
    nodes_df = process_nodes(edges_df, s288c_seqs, args.ID_mappings, args.structure_dir)

    # use DeepTMHMM results to update sequences used by metapredict
    nodes_df = add_DeepTMHMM(nodes_df, args.seq_preds)

    # predict disorder using metapredict
    nodes_df = predict_disorder(nodes_df)
    nodes_df["IDR_count"] = nodes_df["disorder_predictions"].apply(count_IDRs)

    # compute network centrality measures
    nodes_df = compute_centrality(
        edges_df, nodes_df, args.output_dir, args.output_prefix
    )

    # insert information from UniProt
    nodes_df = add_UniProt_info(nodes_df, args.uniprot_data)

    # save outputs
    edges_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}network_edges.csv",
        columns=["source", "target"],
        index=False,
    )

    nodes_df = nodes_df.drop(columns=["disorder_predictions"])
    nodes_df = nodes_df.replace("", "None")
    nodes_df = nodes_df.replace([], "None")
    nodes_df = nodes_df.fillna("None")

    # output .fasta file with trimmed_sequences for IDR property predictions with SPARROW & ALBATROSS
    gen_fasta(nodes_df, f"{args.output_dir}/{args.output_prefix}trimmed_sequence.fasta")

    # create a DataFrame with a random set of 20 rows for testing purposes
    # nodes_df   = nodes_df.sample(n=100, random_state=1991)

    nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}network_nodes_with_annotation.csv",
        index=False,
        na_rep=None,
    )

    # print a "DONE" statement with info about output locations
    print(
        f"Processing complete. Output saved to {args.output_dir}/{args.output_prefix}*"
    )


# entry point
if __name__ == "__main__":

    main()
