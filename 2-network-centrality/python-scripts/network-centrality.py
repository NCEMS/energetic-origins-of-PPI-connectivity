#!/usr/bin/env python3

import sys
from pathlib import Path
import argparse
from typing import Optional, List

import networkx as nx
import pandas as pd
import numpy as np


def validate_edge_columns(edges_df: pd.DataFrame) -> None:
    """
    Confirm that the edge file contains the expected source/target columns.
    """
    required_cols = {"source", "target"}
    missing_cols = required_cols - set(edges_df.columns)

    if missing_cols:
        raise ValueError(
            f"Edge file is missing required columns: {sorted(missing_cols)}. "
            "Expected columns: source, target."
        )


def load_authoritative_nodes(nodes_path: str) -> pd.DataFrame:
    """
    Load the input node list as the authoritative node universe.

    Expected input columns:
        name
        _wkshell

    Output columns:
        node
        _wkshell
    """
    nodes_df = pd.read_csv(nodes_path, usecols=["name", "_wkshell"])

    nodes_df = nodes_df.rename(columns={"name": "node"})

    nodes_df["node"] = nodes_df["node"].astype(str)
    nodes_df["_wkshell"] = nodes_df["_wkshell"].fillna(0.0)

    if nodes_df["node"].isna().any():
        raise ValueError("Input node file contains missing node IDs.")

    duplicated_mask = nodes_df["node"].duplicated()

    if duplicated_mask.any():
        duplicated_examples = (
            nodes_df.loc[duplicated_mask, "node"]
            .head(10)
            .tolist()
        )

        raise ValueError(
            "Input node file contains duplicated node IDs. "
            f"Examples: {duplicated_examples}"
        )

    return nodes_df


def report_node_edge_consistency(
    edges_df: pd.DataFrame,
    nodes_df: pd.DataFrame,
) -> None:
    """
    Report consistency between the edge-list node universe and the input node file.

    The input node file is treated as authoritative. Edge-list nodes absent from
    the node file are an error. Node-file nodes absent from the edge list are
    retained as isolated nodes.
    """
    edge_nodes = set(
        pd.unique(edges_df[["source", "target"]].values.ravel())
    )
    edge_nodes = {node for node in edge_nodes if pd.notna(node)}

    node_file_nodes = set(nodes_df["node"])

    missing_from_node_file = edge_nodes - node_file_nodes
    nodes_without_edges = node_file_nodes - edge_nodes

    if missing_from_node_file:
        examples = sorted(list(missing_from_node_file))[:10]

        raise ValueError(
            f"{len(missing_from_node_file):,} edge-list nodes are absent from "
            f"the input node file. Examples: {examples}"
        )

    print(f"Nodes in input node file:              {len(node_file_nodes):,}")
    print(f"Nodes appearing in edge file:          {len(edge_nodes):,}")
    print(f"Nodes in node file with no edge rows:  {len(nodes_without_edges):,}")


def build_graph_from_edges_and_nodes(
    edges_df: pd.DataFrame,
    nodes_df: pd.DataFrame,
) -> nx.Graph:
    """
    Build an undirected NetworkX graph from the edge list and then explicitly
    add all nodes from the authoritative node file.

    This preserves isolated nodes that are present in the node file but absent
    from the edge list.
    """
    interactome_graph = nx.from_pandas_edgelist(
        edges_df,
        "source",
        "target",
        create_using=nx.Graph(),
    )

    interactome_graph.add_nodes_from(nodes_df["node"])

    return interactome_graph


def get_largest_connected_subgraph(graph: nx.Graph) -> nx.Graph:
    """
    Return the largest connected component as a copied subgraph.

    Information centrality is computed only on the largest connected component,
    then nodes outside that component are assigned 0.0.
    """
    if graph.number_of_nodes() == 0:
        return graph.copy()

    if nx.is_connected(graph):
        return graph

    largest_cc = max(nx.connected_components(graph), key=len)

    return graph.subgraph(largest_cc).copy()


def compute_information_centrality_with_fallback(
    graph: nx.Graph,
    connected_graph: nx.Graph,
) -> dict:
    """
    Compute information centrality on the largest connected component and assign
    0.0 to all nodes outside that component.

    Handles very small/degenerate connected components defensively.
    """
    full_info_centrality = {node: 0.0 for node in graph.nodes}

    if connected_graph.number_of_nodes() <= 1:
        return full_info_centrality

    info_centrality = nx.information_centrality(connected_graph)
    full_info_centrality.update(info_centrality)

    return full_info_centrality


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
    Compute network centrality measures.

    The node list is authoritative. Nodes present in nodes_df but absent from
    edges_df are retained as isolated nodes in the graph.
    """
    interactome_graph = build_graph_from_edges_and_nodes(edges_df, nodes_df)

    total_nodes = interactome_graph.number_of_nodes()
    total_edges = interactome_graph.number_of_edges()
    isolated_nodes = nx.number_of_isolates(interactome_graph)

    interactome_graph_connected = get_largest_connected_subgraph(interactome_graph)
    subgraph_nodes = interactome_graph_connected.number_of_nodes()

    print(
        f"The full graph has {total_nodes:,} nodes and {total_edges:,} edges; "
        f"the largest connected subgraph has {subgraph_nodes:,} nodes"
    )
    print(f"Isolated nodes retained in graph:       {isolated_nodes:,}")

    full_info_centrality = compute_information_centrality_with_fallback(
        interactome_graph,
        interactome_graph_connected,
    )

    centrality_measures = {
        "degree_centrality": nx.degree_centrality(interactome_graph),
        "betweenness_centrality": nx.betweenness_centrality(interactome_graph),
        "eigenvector_centrality": nx.eigenvector_centrality(
            interactome_graph,
            max_iter=1000,
        ),
        "closeness_centrality": nx.closeness_centrality(interactome_graph),
        "load_centrality": nx.load_centrality(interactome_graph),
        "pagerank": nx.pagerank(interactome_graph),
        "information_centrality": full_info_centrality,
    }

    for key, values in centrality_measures.items():
        nodes_df[key] = nodes_df["node"].map(values)

    centrality_cols = list(centrality_measures.keys())
    missing_centrality = nodes_df[centrality_cols].isna().sum().sum()

    if missing_centrality != 0:
        raise ValueError(
            f"Centrality mapping produced {missing_centrality:,} missing values. "
            "This indicates a mismatch between nodes_df and the NetworkX graph."
        )

    test_CentralityCosDist(test_dir, CosDistPath)

    n_before_cosdist = len(nodes_df)

    nodes_df = add_CentralityCosDist(
        nodes_df,
        output_dir,
        output_prefix,
        ["node"] + centrality_cols + ["_wkshell"],
        CosDistPath,
        seeds=seeds,
    )

    n_after_cosdist = len(nodes_df)

    if n_after_cosdist != n_before_cosdist:
        raise ValueError(
            "Row count changed while adding CentralityCosDist results: "
            f"{n_before_cosdist:,} -> {n_after_cosdist:,}"
        )

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
    Run CentralityCosDist on precomputed centrality metrics.

    The CentralityCosDist input file is written with the node identifier in the
    leftmost column named ID.
    """
    sys.path.append(CosDistPath)
    from centralitycosdist import CentralityCosDist

    centralities_df = nodes_df[metrics_list].copy()
    centralities_df = centralities_df.rename(columns={"node": "ID"})

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    cosdist_input = f"{output_dir}/{output_prefix}-CentralityCosDist-input.csv"

    centralities_df.to_csv(cosdist_input, index=False)

    if seeds is None:
        seeds = set(centralities_df["ID"].to_list())
    else:
        seeds = set(seeds)

    nodes = set(centralities_df["ID"].to_list())
    seeds = list(nodes.intersection(seeds))

    algorithm = CentralityCosDist(
        Centrality_file=cosdist_input,
    )

    algorithm.run(seed_nodes=seeds)

    rank_df = algorithm.rank
    rank_df.name = "CentralityCosDist_rank"

    similarity_score_df = algorithm.similarity_score
    similarity_score_df.name = "CentralityCosDist_similarity_score"

    nodes_df = pd.merge(
        nodes_df,
        rank_df,
        how="left",
        left_on="node",
        right_on="ID",
    )

    nodes_df = pd.merge(
        nodes_df,
        similarity_score_df,
        how="left",
        left_on="node",
        right_on="ID",
    )

    return nodes_df


def test_CentralityCosDist(test_dir: str, CosDistPath: str):
    """
    Test whether CentralityCosDist results match expected reference values.

    Expected results are based on:
    https://nilesh-iiita.github.io/CentralityCosDist/notebooks.html
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
        nodes,
        f"{test_dir}/outputs",
        "test",
        metrics_list,
        CosDistPath,
        seeds=seeds,
    )

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

    calculated_similarity_score = result.set_index("node")[
        "CentralityCosDist_similarity_score"
    ].to_dict()

    tol = 1e-5

    for key, expected_val in expected_similarity_score.items():
        assert key in calculated_similarity_score, f"Missing key: {key}"
        assert (
            abs(calculated_similarity_score[key] - expected_val) < tol
        ), f"Reference and calculated values for {key} do not match"


def main():
    parser = argparse.ArgumentParser(description="Compute centrality metrics")

    parser.add_argument(
        "--edges",
        required=True,
        help="Path to the edges CSV file.",
    )

    parser.add_argument(
        "--nodes",
        required=True,
        help="Path to the nodes CSV file.",
    )

    parser.add_argument(
        "--output_prefix",
        required=True,
        help="Prefix for output files.",
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory.",
    )

    parser.add_argument(
        "--output_suffix",
        required=True,
        help="Output suffix for final annotated nodes from this step.",
    )

    parser.add_argument(
        "--organism_tag",
        required=True,
        help="Tag to label the organism for this run.",
    )

    parser.add_argument(
        "--CosDistPath",
        required=True,
        help="Path to directory containing CentralityCosDist code.",
    )

    parser.add_argument(
        "--test_dir",
        required=True,
        help="Path to directory with data required to run tests for this program.",
    )

    args = parser.parse_args()

    edges_df = pd.read_csv(args.edges)
    validate_edge_columns(edges_df)

    nodes_df = load_authoritative_nodes(args.nodes)

    report_node_edge_consistency(edges_df, nodes_df)

    nodes_df = compute_centrality(
        edges_df=edges_df,
        nodes_df=nodes_df,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        test_dir=args.test_dir,
        CosDistPath=args.CosDistPath,
    )

    nodes_df = nodes_df.drop(
        columns=[
            "CentralityCosDist_rank",
            "CentralityCosDist_similarity_score",
        ]
    )

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    output_file = (
        f"{args.output_dir}/"
        f"{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.csv"
    )

    nodes_df.to_csv(
        output_file,
        index=False,
        na_rep=np.nan,
    )

    print()
    print(f"Wrote centrality-annotated nodes to: {output_file}")
    print(f"Output rows:                         {len(nodes_df):,}")


if __name__ == "__main__":
    main()
