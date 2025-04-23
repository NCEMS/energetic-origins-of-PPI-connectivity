import os, sys
from pathlib import Path
import argparse
import networkx as nx
import pandas as pd
import numpy as np
import typing
from typing import Optional
from typing import List
from typing import Dict
import pytest


def compute_centrality(
    edges_df: pd.DataFrame,
    nodes_df: pd.DataFrame,
    output_dir: str,
    output_prefix: str,
    test_dir: str,
    CosDistPath: str,
    seeds: Optional[List] = None,
) -> pd.DataFrame:
    """
    Compute network centrality measures on network

    Args:
        edges_df (pd.DataFrame): network edges
        nodes_df (pd.DataFrame): network nodes
        output_dir (str): output directory, used by add_CentralityCosDist
        output_prefix (str): output file prefix, used by add_CentralityCosDist
        test_dir (str): path to directory with data required to test CentralityCosDist performance
        CosDistPath (str): path to CentralityCosDist code

    Returns:
        pd.DataFrame
    """

    # make networkx style graph
    interactome_graph = nx.from_pandas_edgelist(edges_df, "source", "target", create_using=nx.Graph())

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
    }

    # this version of centrality_measure contains information centrality, which currently fails as the graph is not connected...
    #centrality_measures = {
    #    "degree_centrality": nx.degree_centrality(interactome_graph),
    #    "betweenness_centrality": nx.betweenness_centrality(interactome_graph),
    #    "eigenvector_centrality": nx.eigenvector_centrality(
    #        interactome_graph, max_iter=1000
    #    ),
    #    "closeness_centrality": nx.closeness_centrality(interactome_graph),
    #    "load_centrality": nx.load_centrality(interactome_graph),
    #    "pagerank": nx.pagerank(interactome_graph),
    #    "information_centrality": nx.information_centrality(interactome_graph)
    #}

    # add per-node information to the DataFrame
    for key, values in centrality_measures.items():
        nodes_df[key] = nodes_df["node"].map(values)

    # compute CentralityCosDist; requires a file in a very specific format
    test_CentralityCosDist(test_dir, CosDistPath)
    nodes_df = add_CentralityCosDist(
        nodes_df,
        output_dir,
        output_prefix,
        ["node"] + list(centrality_measures.keys()),
        CosDistPath,
        seeds=seeds,
    )

    # return the updated DataFrame
    return nodes_df


def add_CentralityCosDist(
    nodes_df: pd.DataFrame,
    output_dir: str,
    output_prefix: str,
    metrics_list: List[str],
    CosDistPath: str,
    seeds: Optional[List] = None,
) -> pd.DataFrame:
    """
    Carries out CentralityCosDist calculations for the network input at a DataFrame

    Args:
        nodes_df (pd.DataFrame): DataFrame containing precomputed centrality metrics as required by CentralityCosDist algorithm
        output_dir (str): path to output data directory
        output_prefix (str): prefix to be appended to output file
        metrics_list (List[str]): list of the metric names to be extracted from nodes_df for cosine distance calculation
        CosDistPath (str): path to CentralityCosDist code
        seeds (Optional[List]): either a list of nodes to treat as seeds or None; if None, all nodes will be treated as seeds

    Returns:
        pd.DataFrame

    Note well: the CentralityCosDist program will interpret your input csv file in the following way:
               * the leftmost column will be treated as a unique identifier for nodes in the network; seeds must be cross-referenceable with this node list
               * all other columns in the DataFrame will be treated as centrality metrics and used in the cosine distance calculation
    """
    # load CentralityCosDist functionality
    sys.path.append(CosDistPath)
    from centralitycosdist import CentralityCosDist

    centralities_df = nodes_df[metrics_list]
    centralities_df = centralities_df.rename(columns={"node": "ID"})
    centralities_df.to_csv(
        f"{output_dir}/{output_prefix}-CentralityCosDist-input.csv", index=False
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
        Centrality_file=f"{output_dir}/{output_prefix}-CentralityCosDist-input.csv"
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


def test_CentralityCosDist(test_dir: str, CosDistPath: str):
    """
    Function to test whether or not results from CentralityCosDist match expectations
    Expected results are based on https://nilesh-iiita.github.io/CentralityCosDist/notebooks.html

    Args:
        test_dir (str): Path to directory containing input data required for the test
        CosDistPath (str): Path to CentralityCosDist code

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
    nodes = pd.read_csv(f"{test_dir}/inputs/Network_Centrality.csv")
    metrics_list = [
        "node",
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
        nodes, f"{test_dir}/outputs", "test_", metrics_list, CosDistPath, seeds=seeds
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


def main():

    # setup command-line arguments
    parser = argparse.ArgumentParser(description="Compute centrality metrics")
    parser.add_argument(
        "--edges",
        help="Path to the edges CSV file",
    )
    parser.add_argument("--output_prefix", help="Prefix for output files")
    parser.add_argument(
        "--output_dir", default="processed-data", help="Output directory"
    )
    parser.add_argument(
        "--organism_tag", default="s288c", help="Tag to label the organism for this run"
    )
    parser.add_argument(
        "--CosDistPath",
        default="1-network-centrality/python-scripts/CentralityCosDist",
        help="Path to directory containing CentralityCosDist code",
    )
    parser.add_argument(
        "--test_dir",
        default="test",
        help="Path to directory with data required to run tests for this program",
    )
    args = parser.parse_args()

    # load input data
    edges_df = pd.read_csv(args.edges)

    # extract nodes and create pd.DataFrame
    nodes_df = pd.DataFrame(
        pd.unique(edges_df[["source", "target"]].values.ravel()), columns=["node"]
    )

    # compute network centrality measures
    nodes_df = compute_centrality(
        edges_df, nodes_df, args.output_dir, args.output_prefix, args.test_dir, args.CosDistPath
    )

    # save the output to file
    nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality.csv",
        index=False,
        na_rep=None,
    )


if __name__ == "__main__":

    main()
