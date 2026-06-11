#!/usr/bin/env python3

"""Build A. thaliana in-vivo, in-vitro, union, and intersection PPI networks.

This script expects two processed evidence tables:

* ppi_in_vivo_edges.csv
* ppi_in_vitro_edges.csv

Each table must include protein_a and protein_b columns. Any additional columns
are preserved for the individual in-vivo/in-vitro edge outputs. The individual
network files preserve evidence-row granularity, including self-loop rows,
matching the existing A. thaliana branch data shape. Union and intersection edge
outputs are canonicalized source/target tables for downstream centrality
calculations.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


DEFAULT_OUTPUTS = {
    "in_vivo_edges": "in-vivo-edges.csv",
    "in_vivo_nodes": "in-vivo-nodes.csv",
    "in_vitro_edges": "in-vitro-edges.csv",
    "in_vitro_nodes": "in-vitro-nodes.csv",
    "union_edges": "merged-edges.csv",
    "union_nodes": "merged-nodes.csv",
    "intersection_nodes": "intersection-nodes.csv",
    "intersection_edges_in_vivo": "intersection-edges-in_vivo.csv",
    "intersection_edges_in_vitro": "intersection-edges-in_vitro.csv",
}


def clean_protein_id(value: object) -> str | None:
    if pd.isna(value):
        return None
    cleaned = str(value).strip().upper()
    return cleaned or None


def read_evidence_edges(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"protein_a", "protein_b"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")

    df = df.copy()
    source = df["protein_a"].map(clean_protein_id)
    target = df["protein_b"].map(clean_protein_id)
    df = df.drop(columns=["protein_a", "protein_b"])
    df.insert(0, "target", target)
    df.insert(0, "source", source)
    df = df.dropna(subset=["source", "target"])

    if "label" not in df.columns:
        df["label"] = label
    else:
        df["label"] = df["label"].fillna(label)

    return df


def build_nodes(edges: pd.DataFrame) -> pd.DataFrame:
    nodes = pd.unique(edges[["source", "target"]].values.ravel())
    nodes = sorted(node for node in nodes if pd.notna(node))
    nodes_df = pd.DataFrame({"name": nodes})
    nodes_df["_wkshell"] = 0.0
    return nodes_df


def bare_edges(edges: pd.DataFrame) -> pd.DataFrame:
    edge_min = edges[["source", "target"]].min(axis=1)
    edge_max = edges[["source", "target"]].max(axis=1)
    canonical = pd.DataFrame({"source": edge_min, "target": edge_max})
    return canonical.drop_duplicates().sort_values(["source", "target"])


def filter_edges_to_nodes(edges: pd.DataFrame, nodes: set[str]) -> pd.DataFrame:
    filtered = edges[edges["source"].isin(nodes) & edges["target"].isin(nodes)]
    return bare_edges(filtered)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def build_networks(args: argparse.Namespace) -> None:
    output_dir = args.output_dir
    in_vivo_edges = read_evidence_edges(args.in_vivo_input, "in vivo")
    in_vitro_edges = read_evidence_edges(args.in_vitro_input, "in vitro")

    in_vivo_nodes = build_nodes(in_vivo_edges)
    in_vitro_nodes = build_nodes(in_vitro_edges)

    in_vivo_node_set = set(in_vivo_nodes["name"])
    in_vitro_node_set = set(in_vitro_nodes["name"])
    union_node_set = in_vivo_node_set | in_vitro_node_set
    intersection_node_set = in_vivo_node_set & in_vitro_node_set

    union_edges = (
        pd.concat([bare_edges(in_vivo_edges), bare_edges(in_vitro_edges)], ignore_index=True)
        .drop_duplicates()
        .sort_values(["source", "target"])
    )
    union_nodes = build_nodes(union_edges)
    intersection_nodes = pd.DataFrame({"name": sorted(intersection_node_set)})
    intersection_nodes["_wkshell"] = 0.0

    intersection_edges_in_vivo = filter_edges_to_nodes(in_vivo_edges, intersection_node_set)
    intersection_edges_in_vitro = filter_edges_to_nodes(in_vitro_edges, intersection_node_set)

    outputs = {key: output_dir / value for key, value in DEFAULT_OUTPUTS.items()}
    write_csv(in_vivo_edges, outputs["in_vivo_edges"])
    write_csv(in_vivo_nodes, outputs["in_vivo_nodes"])
    write_csv(in_vitro_edges, outputs["in_vitro_edges"])
    write_csv(in_vitro_nodes, outputs["in_vitro_nodes"])
    write_csv(union_edges, outputs["union_edges"])
    write_csv(union_nodes, outputs["union_nodes"])
    write_csv(intersection_nodes, outputs["intersection_nodes"])
    write_csv(intersection_edges_in_vivo, outputs["intersection_edges_in_vivo"])
    write_csv(intersection_edges_in_vitro, outputs["intersection_edges_in_vitro"])

    print(f"In-vivo edges:                       {len(in_vivo_edges):,}")
    print(f"In-vitro edges:                      {len(in_vitro_edges):,}")
    print(f"Union edges:                         {len(union_edges):,}")
    print(f"In-vivo nodes:                       {len(in_vivo_nodes):,}")
    print(f"In-vitro nodes:                      {len(in_vitro_nodes):,}")
    print(f"Union nodes:                         {len(union_node_set):,}")
    print(f"Intersection nodes:                  {len(intersection_node_set):,}")
    print(f"In-vivo edges on intersection nodes: {len(intersection_edges_in_vivo):,}")
    print(f"In-vitro edges on intersection nodes:{len(intersection_edges_in_vitro):,}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--in-vivo-input",
        type=Path,
        default=Path("data-files/ppi_in_vivo_edges.csv"),
        help="Processed in-vivo evidence edge table with protein_a/protein_b columns.",
    )
    parser.add_argument(
        "--in-vitro-input",
        type=Path,
        default=Path("data-files/ppi_in_vitro_edges.csv"),
        help="Processed in-vitro evidence edge table with protein_a/protein_b columns.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data-files"),
        help="Directory where generated network CSVs should be written.",
    )
    return parser.parse_args()


def main() -> None:
    build_networks(parse_args())


if __name__ == "__main__":
    main()
