# ppi_io.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
import networkx as nx

import os

@dataclass(frozen=True)
class PPILoadResult:
    nodes_df: pd.DataFrame
    edges_df: pd.DataFrame
    graph: nx.Graph
    u_col: str
    v_col: str
    n_overlap_annot_graph: int


def _require_column(df: pd.DataFrame, col: str, df_name: str) -> None:
    if col not in df.columns:
        raise KeyError(
            f"Expected column '{col}' in {df_name}, but it was not found. "
            f"Available columns: {list(df.columns)}"
        )


def load_ppi(
    nodes_path: str | Path,
    edges_path: str | Path,
    *,
    node_col: str = "node",
    u_col: str = "source",
    v_col: str = "target",
    graph_type: str = "undirected",
    restrict_to_annotated_nodes: bool = False,
) -> PPILoadResult:
    """
    Load node annotation DataFrame (pickle) and edges (CSV), then build a NetworkX graph.

    Parameters
    ----------
    nodes_path : path to .pkl created by pandas (one row per node)
    edges_path : path to edge list CSV
    node_col   : column in nodes_df that contains node IDs
    u_col/v_col: columns in edges_df for source/target node IDs
    graph_type : 'undirected' (nx.Graph) or 'directed' (nx.DiGraph)
    restrict_to_annotated_nodes : if True, drop edges whose endpoints are not in nodes_df

    Returns
    -------
    PPILoadResult
    """
    nodes_path = Path(nodes_path)
    edges_path = Path(edges_path)

    if not nodes_path.exists():
        raise FileNotFoundError(f"Node annotation file not found: {nodes_path}")
    if not edges_path.exists():
        raise FileNotFoundError(f"Edge list file not found: {edges_path}")

    # --- Load node annotations ---
    nodes_df = pd.read_pickle(nodes_path)
    _require_column(nodes_df, node_col, "nodes_df")

    # normalize node IDs
    nodes_df = nodes_df.copy()
    nodes_df[node_col] = nodes_df[node_col].astype(str)

    # --- Load edge list ---
    edges_df = pd.read_csv(edges_path)
    _require_column(edges_df, u_col, "edges_df")
    _require_column(edges_df, v_col, "edges_df")

    edges_df = edges_df.copy()
    edges_df[u_col] = edges_df[u_col].astype(str)
    edges_df[v_col] = edges_df[v_col].astype(str)

    if restrict_to_annotated_nodes:
        annot_nodes = set(nodes_df[node_col])
        edges_df = edges_df[edges_df[u_col].isin(annot_nodes) & edges_df[v_col].isin(annot_nodes)].copy()

    # --- Build graph ---
    if graph_type not in {"undirected", "directed"}:
        raise ValueError("graph_type must be one of: {'undirected', 'directed'}")

    G = nx.from_pandas_edgelist(
        edges_df,
        source=u_col,
        target=v_col,
        create_using=(nx.DiGraph() if graph_type == "directed" else nx.Graph()),
    )

    annot_nodes = set(nodes_df[node_col])
    graph_nodes = set(G.nodes())
    overlap = annot_nodes & graph_nodes

    return PPILoadResult(
        nodes_df=nodes_df,
        edges_df=edges_df,
        graph=G,
        u_col=u_col,
        v_col=v_col,
        n_overlap_annot_graph=len(overlap),
    )
