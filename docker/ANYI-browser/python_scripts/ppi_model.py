# python_scripts/ppi_model.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Set, Tuple, Sequence

import pandas as pd
import networkx as nx


@dataclass(frozen=True)
class PPIModel:
    """
    Lightweight container for the core objects the UI will use.
    """
    nodes_df: pd.DataFrame          # indexed by node id (string)
    edges_df: pd.DataFrame          # normalized (string endpoints)
    graph: nx.Graph | nx.DiGraph
    node_col: str
    u_col: str
    v_col: str
    alignment: Dict[str, int]


def _require_columns(df: pd.DataFrame, cols: Iterable[str], df_name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(
            f"Missing required column(s) {missing} in {df_name}. "
            f"Available columns: {list(df.columns)}"
        )


def index_nodes_df(nodes_df: pd.DataFrame, *, node_col: str = "node") -> pd.DataFrame:
    """
    Ensure node IDs are strings and unique; return a copy indexed by node_col.
    """
    _require_columns(nodes_df, [node_col], "nodes_df")
    df = nodes_df.copy()
    df[node_col] = df[node_col].astype(str)

    # Defensive: disallow NA / empty IDs
    if df[node_col].isna().any():
        raise ValueError(f"{node_col} contains NA values; all node IDs must be defined.")
    if (df[node_col].str.len() == 0).any():
        raise ValueError(f"{node_col} contains empty strings; all node IDs must be non-empty.")

    # Enforce uniqueness
    if not df[node_col].is_unique:
        dupes = df.loc[df[node_col].duplicated(), node_col].head(10).tolist()
        raise ValueError(
            f"{node_col} is not unique in nodes_df. Example duplicate IDs: {dupes}"
        )

    return df.set_index(node_col, drop=False)


def normalize_edges_df(edges_df: pd.DataFrame, *, u_col: str = "source", v_col: str = "target") -> pd.DataFrame:
    """
    Ensure edge endpoints are strings; return a copy with normalized endpoints.
    """
    _require_columns(edges_df, [u_col, v_col], "edges_df")
    e = edges_df.copy()
    e[u_col] = e[u_col].astype(str)
    e[v_col] = e[v_col].astype(str)
    return e


def build_graph(
    edges_df: pd.DataFrame,
    *,
    u_col: str = "source",
    v_col: str = "target",
    graph_type: str = "undirected",
) -> nx.Graph | nx.DiGraph:
    """
    Build a NetworkX graph from an edge list DataFrame.
    """
    if graph_type not in {"undirected", "directed"}:
        raise ValueError("graph_type must be one of: {'undirected', 'directed'}")

    create_using = nx.DiGraph() if graph_type == "directed" else nx.Graph()
    G = nx.from_pandas_edgelist(edges_df, source=u_col, target=v_col, create_using=create_using)
    return G


def summarize_alignment(
    nodes_df_indexed: pd.DataFrame,
    G: nx.Graph | nx.DiGraph,
    *,
    node_col: str = "node",
) -> Dict[str, int]:
    """
    Summarize overlap between annotated nodes and graph nodes, plus simple graph diagnostics.
    """
    annot_nodes: Set[str] = set(nodes_df_indexed[node_col].astype(str))
    graph_nodes: Set[str] = set(map(str, G.nodes()))

    overlap = annot_nodes & graph_nodes
    annot_only = annot_nodes - graph_nodes
    graph_only = graph_nodes - annot_nodes

    # Graph diagnostics
    n_self_loops = nx.number_of_selfloops(G)
    n_isolates = len(list(nx.isolates(G))) if not G.is_directed() else len([n for n in G.nodes() if G.degree(n) == 0])

    return {
        "n_nodes_annot": len(annot_nodes),
        "n_nodes_graph": len(graph_nodes),
        "n_nodes_overlap": len(overlap),
        "n_nodes_annot_only": len(annot_only),
        "n_nodes_graph_only": len(graph_only),
        "n_edges_graph": int(G.number_of_edges()),
        "n_self_loops": int(n_self_loops),
        "n_isolates": int(n_isolates),
    }


def add_percentile_columns(
    nodes_df: pd.DataFrame,
    cols: Sequence[str],
    *,
    suffix: str = "_percentile",
    scale_0_100: bool = True,
) -> pd.DataFrame:
    """
    Add percentile-rank columns for numeric columns in `cols`.

    - Uses pandas rank(pct=True) computed over non-missing values.
    - Missing values remain missing in the percentile column (important for information_centrality).
    - Output is 0–100 if scale_0_100=True, else 0–1.
    """
    df = nodes_df.copy()

    for c in cols:
        if c not in df.columns:
            continue  # skip silently; you can make this strict if you prefer

        out_col = f"{c}{suffix}"
        if out_col in df.columns:
            # Avoid overwriting if caller already provided precomputed percentiles
            continue

        s = pd.to_numeric(df[c], errors="coerce")  # non-numeric -> NaN
        pct = s.rank(pct=True)  # NaNs stay NaN; ranks computed on non-NaN only

        if scale_0_100:
            pct = pct * 100.0

        # Use pandas NA-friendly dtype
        df[out_col] = pct

    return df


def prepare_model(
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
    *,
    node_col: str = "node",
    u_col: str = "source",
    v_col: str = "target",
    graph_type: str = "undirected",
    required_node_columns: Tuple[str, ...] = ("degree_centrality", "DeepTMHMM_class"),
) -> PPIModel:
    """
    End-to-end preparation for the interactive notebook.

    - indexes nodes_df by node ID
    - normalizes edge endpoints
    - builds a graph
    - validates presence of selected annotation columns (optional)
    - computes alignment summary
    - adds percentile columns for centrality metrics
    """
    nodes_ix = index_nodes_df(nodes_df, node_col=node_col)
    edges_norm = normalize_edges_df(edges_df, u_col=u_col, v_col=v_col)
    G = build_graph(edges_norm, u_col=u_col, v_col=v_col, graph_type=graph_type)

    # Validate that the annotation columns you want to expose exist
    missing = [c for c in required_node_columns if c not in nodes_ix.columns]
    if missing:
        raise KeyError(
            f"Missing required annotation column(s) in nodes_df: {missing}. "
            f"Available columns: {list(nodes_ix.columns)}"
        )

    # Add percentile columns (computed once on load)
    percentile_cols = (
        "degree_centrality",
        "betweenness_centrality",
        "eigenvector_centrality",
        "closeness_centrality",
        "load_centrality",
        "pagerank",
        "information_centrality",
        "median_molecules_per_cell",
        "Villen_halflife_min",
        "meltome-melting-point",
    )
    nodes_ix = add_percentile_columns(nodes_ix, percentile_cols, suffix="_percentile", scale_0_100=True)

    alignment = summarize_alignment(nodes_ix, G, node_col=node_col)

    return PPIModel(
        nodes_df=nodes_ix,
        edges_df=edges_norm,
        graph=G,
        node_col=node_col,
        u_col=u_col,
        v_col=v_col,
        alignment=alignment,
    )


def get_node_display_row(
    model: PPIModel,
    node_id: str,
    *,
    fields: Tuple[str, ...] = ("node", "degree_centrality", "DeepTMHMM_class"),
) -> pd.Series:
    """
    Convenience helper for the UI: return a subset of fields for a node.
    Raises KeyError if node_id is not present in the annotations.
    """
    node_id = str(node_id)
    if node_id not in model.nodes_df.index:
        raise KeyError(
            f"Node '{node_id}' not found in node annotations. "
            f"(Annotations contain {len(model.nodes_df)} nodes.)"
        )

    row = model.nodes_df.loc[node_id]
    missing = [f for f in fields if f not in row.index]
    if missing:
        raise KeyError(f"Requested fields not found in nodes_df for node '{node_id}': {missing}")
    return row[list(fields)]
