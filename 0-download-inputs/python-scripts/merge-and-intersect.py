#!/usr/bin/env python3

import argparse
import pandas as pd


def read_edges(path: str, sep: str) -> pd.DataFrame:
    """
    Read an edge file with columns protein_a and protein_b.
    """
    df = pd.read_csv(path, sep=sep, usecols=["protein_a", "protein_b"])

    df = df.rename(
        columns={
            "protein_a": "source",
            "protein_b": "target",
        }
    )

    return df


def canonicalize_undirected_edges(edges_df: pd.DataFrame) -> pd.DataFrame:
    """
    Canonicalize undirected edges so that A-B and B-A are treated as identical.
    """
    edge_min = edges_df[["source", "target"]].min(axis=1)
    edge_max = edges_df[["source", "target"]].max(axis=1)

    canonical_edges = pd.DataFrame(
        {
            "source": edge_min,
            "target": edge_max,
        }
    )

    canonical_edges = canonical_edges.dropna(subset=["source", "target"])
    canonical_edges = canonical_edges.drop_duplicates()

    return canonical_edges


def get_nodes(edges_df: pd.DataFrame) -> set:
    """
    Get unique nodes from source and target columns.
    """
    nodes = pd.unique(edges_df[["source", "target"]].values.ravel())
    nodes = {node for node in nodes if pd.notna(node)}

    return nodes


def build_node_list(nodes: set) -> pd.DataFrame:
    """
    Build node list with dummy _wkshell column.

    Output columns:
        name
        _wkshell
    """
    nodes_df = pd.DataFrame({"name": sorted(nodes)})
    nodes_df["_wkshell"] = 0

    return nodes_df


def filter_edges_to_node_set(edges_df: pd.DataFrame, node_set: set) -> pd.DataFrame:
    """
    Keep only edges for which both source and target are present in node_set.
    """
    filtered_edges = edges_df[
        edges_df["source"].isin(node_set)
        & edges_df["target"].isin(node_set)
    ].copy()

    filtered_edges = filtered_edges.drop_duplicates()

    return filtered_edges


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compute union network outputs, individual node lists, and "
            "node-intersection-filtered edge lists from two protein edge files."
        )
    )

    parser.add_argument("--file1", required=True, help="First edge file")
    parser.add_argument("--file2", required=True, help="Second edge file")

    parser.add_argument(
        "--sep",
        default=",",
        help="Input file delimiter. Default: comma. Use '\\t' for TSV files.",
    )

    parser.add_argument(
        "--union_edges",
        required=True,
        help="Output CSV file for union edge list with columns source and target.",
    )

    parser.add_argument(
        "--union_nodes",
        required=True,
        help="Output CSV file for union node list with columns name and _wkshell.",
    )

    parser.add_argument(
        "--file1_nodes",
        required=True,
        help="Output CSV file for file1 node list with columns name and _wkshell.",
    )

    parser.add_argument(
        "--file2_nodes",
        required=True,
        help="Output CSV file for file2 node list with columns name and _wkshell.",
    )

    parser.add_argument(
        "--intersection_nodes",
        required=True,
        help=(
            "Output CSV file for node intersection list. "
            "These are nodes present in both input edge files."
        ),
    )

    parser.add_argument(
        "--file1_intersection_edges",
        required=True,
        help=(
            "Output CSV file for file1 edges where both endpoints are in the "
            "intersection node list."
        ),
    )

    parser.add_argument(
        "--file2_intersection_edges",
        required=True,
        help=(
            "Output CSV file for file2 edges where both endpoints are in the "
            "intersection node list."
        ),
    )

    args = parser.parse_args()

    sep = "\t" if args.sep == "\\t" else args.sep

    # Read and standardize edge files
    edges1 = read_edges(args.file1, sep)
    edges2 = read_edges(args.file2, sep)

    # Canonicalize as undirected PPI edges
    edges1 = canonicalize_undirected_edges(edges1)
    edges2 = canonicalize_undirected_edges(edges2)

    # Compute union edge network
    union_edges = pd.concat([edges1, edges2], ignore_index=True)
    union_edges = union_edges.drop_duplicates()

    # Compute node sets
    nodes1 = get_nodes(edges1)
    nodes2 = get_nodes(edges2)

    union_nodes_set = nodes1.union(nodes2)
    intersection_nodes_set = nodes1.intersection(nodes2)

    # Build node-list DataFrames
    file1_nodes = build_node_list(nodes1)
    file2_nodes = build_node_list(nodes2)
    union_nodes = build_node_list(union_nodes_set)
    intersection_nodes = build_node_list(intersection_nodes_set)

    # Filter each original edge set to the shared node universe
    file1_intersection_edges = filter_edges_to_node_set(
        edges1,
        intersection_nodes_set,
    )

    file2_intersection_edges = filter_edges_to_node_set(
        edges2,
        intersection_nodes_set,
    )

    # Write outputs
    union_edges.to_csv(args.union_edges, index=False)
    union_nodes.to_csv(args.union_nodes, index=False)

    file1_nodes.to_csv(args.file1_nodes, index=False)
    file2_nodes.to_csv(args.file2_nodes, index=False)

    intersection_nodes.to_csv(args.intersection_nodes, index=False)
    file1_intersection_edges.to_csv(args.file1_intersection_edges, index=False)
    file2_intersection_edges.to_csv(args.file2_intersection_edges, index=False)

    # Report
    print(f"Unique edges in file1:                    {len(edges1):,}")
    print(f"Unique edges in file2:                    {len(edges2):,}")
    print(f"Union edges:                              {len(union_edges):,}")

    print(f"Nodes in file1:                           {len(file1_nodes):,}")
    print(f"Nodes in file2:                           {len(file2_nodes):,}")
    print(f"Union nodes:                              {len(union_nodes):,}")
    print(f"Intersection nodes:                       {len(intersection_nodes):,}")

    print(f"File1 edges within intersection nodes:    {len(file1_intersection_edges):,}")
    print(f"File2 edges within intersection nodes:    {len(file2_intersection_edges):,}")

    print()
    print(f"Wrote union edge list to:                 {args.union_edges}")
    print(f"Wrote union node list to:                 {args.union_nodes}")
    print(f"Wrote file1 node list to:                 {args.file1_nodes}")
    print(f"Wrote file2 node list to:                 {args.file2_nodes}")
    print(f"Wrote intersection node list to:          {args.intersection_nodes}")
    print(f"Wrote file1 intersection edge list to:    {args.file1_intersection_edges}")
    print(f"Wrote file2 intersection edge list to:    {args.file2_intersection_edges}")


if __name__ == "__main__":
    main()
