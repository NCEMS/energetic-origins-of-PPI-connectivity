# python_scripts/ppi_ui.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import pandas as pd
import ipywidgets as W
from IPython.display import display
import networkx as nx

from .ppi_model import PPIModel, get_node_display_row
from .ppi_structure import get_pdb_path_for_node, resolve_structures_dir, make_nglview_with_plddt_bins, plddt_legend_html

# ipycytoscape is optional at import time; if missing, the UI will still work
try:
    from ipycytoscape import CytoscapeWidget
    _HAVE_CYTO = True
except Exception:
    CytoscapeWidget = None  # type: ignore
    _HAVE_CYTO = False

# import guard for NGLView
try:
    import nglview as nv
    _HAVE_NGL = True
except Exception:
    nv = None  # type: ignore
    _HAVE_NGL = False


@dataclass(frozen=True)
class PPIUI:
    root: W.Widget
    node_selector: W.Widget
    status: W.HTML
    summary: W.HTML
    neighbors: W.Output
    cytoscape: Optional[object]  # CytoscapeWidget if available


def _escape_html(x: object) -> str:
    s = "" if x is None else str(x)
    return (
        s.replace("&", "&amp;")
         .replace("<", "&lt;")
         .replace(">", "&gt;")
         .replace('"', "&quot;")
         .replace("'", "&#39;")
    )


def _format_summary_table(row: pd.Series, title: str = "Selected node") -> str:
    items = [(k, row.get(k, "")) for k in row.index]
    rows_html = "\n".join(
        f"<tr>"
        f"<td style='padding:4px 10px; font-weight:600; border:1px solid #ddd;'>{_escape_html(k)}</td>"
        f"<td style='padding:4px 10px; border:1px solid #ddd;'>{_escape_html(v)}</td>"
        f"</tr>"
        for k, v in items
    )
    return f"""
    <div style="margin: 6px 0 10px 0;">
      <div style="font-weight:700; margin-bottom:6px;">{_escape_html(title)}</div>
      <table style="border-collapse:collapse; border:1px solid #ddd;">
        {rows_html}
      </table>
    </div>
    """


def _neighbors_df(G: nx.Graph | nx.DiGraph, node_id: str, max_n: int = 50) -> pd.DataFrame:
    node_id = str(node_id)
    if node_id not in G:
        return pd.DataFrame({"neighbor": []})

    # For directed graphs, "neighbors" = successors; for undirected it's the usual.
    nbrs = list(G.neighbors(node_id))
    nbrs = [str(n) for n in nbrs][:max_n]
    return pd.DataFrame({"neighbor": nbrs})


def _ego_subgraph(
    G: nx.Graph | nx.DiGraph,
    node_id: str,
    hops: int = 1,
    max_nodes: int = 200,
) -> nx.Graph | nx.DiGraph:
    """
    Build a small ego subgraph around node_id.

    - hops=1: node + first-degree neighbors
    - hops=2: include neighbors-of-neighbors (can grow quickly)
    """
    node_id = str(node_id)
    if node_id not in G:
        return G.__class__()  # empty graph of same type

    hops = int(hops)
    if hops < 1:
        hops = 1
    if hops > 3:
        hops = 3  # defensive cap

    # Start with the focal node
    visited = {node_id}
    frontier = {node_id}

    for _ in range(hops):
        nxt = set()
        for n in frontier:
            # directed: include both successors and predecessors for a more intuitive ego view
            if G.is_directed():
                nxt.update(map(str, G.successors(n)))
                nxt.update(map(str, G.predecessors(n)))
            else:
                nxt.update(map(str, G.neighbors(n)))
        nxt -= visited
        visited |= nxt
        frontier = nxt

        if len(visited) >= max_nodes:
            break

    # Create induced subgraph
    return G.subgraph(list(visited)).copy()


def _cyto_clear(cyto) -> None:
    """
    Clear cytoscape widget graph in a version-tolerant way.
    """
    # Common API: cyto.graph.clear()
    try:
        cyto.graph.clear()
        return
    except Exception:
        pass

    # Fallback: try to clear nodes/edges collections
    try:
        cyto.graph.nodes.clear()
    except Exception:
        pass
    try:
        cyto.graph.edges.clear()
    except Exception:
        pass


def _cyto_load_networkx(cyto, H: nx.Graph | nx.DiGraph) -> None:
    """
    Load a NetworkX graph into cytoscape, with minimal style.
    """
    _cyto_clear(cyto)

    # Ensure nodes have a label attribute for display
    for n in H.nodes():
        if "label" not in H.nodes[n]:
            H.nodes[n]["label"] = str(n)

    # Ensure edges have ids (helpful for some cytoscape backends)
    for i, (u, v, k) in enumerate(H.edges(keys=True) if H.is_multigraph() else [(u, v, None) for u, v in H.edges()]):
        # Only set if not already present
        if H.is_multigraph():
            if "id" not in H.edges[u, v, k]:
                H.edges[u, v, k]["id"] = f"{u}__{v}__{k}"
        else:
            if "id" not in H.edges[u, v]:
                H.edges[u, v]["id"] = f"{u}__{v}"

    # Add graph
    cyto.graph.add_graph_from_networkx(H)

    # Style: label nodes; keep other defaults minimal
    cyto.set_style([
        {"selector": "node", "style": {"label": "data(label)"}},
        {"selector": "edge", "style": {"curve-style": "bezier"}},
    ])


def build_ui(
    model: PPIModel,
    *,
    default_node: Optional[str] = None,
    display_fields: Tuple[str, ...] = ("node", "degree_centrality", "DeepTMHMM_class"),
    max_neighbors_default: int = 30,
    include_cytoscape: bool = True,
) -> PPIUI:
    node_ids = list(map(str, model.nodes_df.index))

    if default_node is None:
        default_node = node_ids[0] if node_ids else ""
    default_node = str(default_node)

    # Searchable selector (Combobox). Falls back to Dropdown if Combobox unavailable.
    try:
        node_selector = W.Combobox(
            options=node_ids,
            value=default_node if default_node in node_ids else "",
            placeholder="Type a node ID...",
            description="Node:",
            ensure_option=True,
            layout=W.Layout(width="520px"),
        )
    except Exception:
        node_selector = W.Dropdown(
            options=node_ids,
            value=default_node if default_node in node_ids else (node_ids[0] if node_ids else ""),
            description="Node:",
            layout=W.Layout(width="520px"),
        )

    max_neighbors = W.IntSlider(
        value=max_neighbors_default,
        min=0,
        max=200,
        step=10,
        description="Neighbors:",
        continuous_update=False,
        layout=W.Layout(width="520px"),
    )

    show_neighbors = W.Checkbox(
        value=True,
        description="Show neighbors table",
        indent=False,
    )

    status = W.HTML(value="")
    summary = W.HTML(value="")
    neighbors_out = W.Output()

    # Cytoscape controls + widget (optional)
    cyto = None
    cyto_panel = None

    show_network = W.Checkbox(value=True, description="Show network view", indent=False)
    hops = W.IntSlider(
        value=1, min=1, max=2, step=1,
        description="Hops:",
        continuous_update=False,
        layout=W.Layout(width="260px"),
    )
    layout_dropdown = W.Dropdown(
        options=["cose", "circle", "grid", "breadthfirst"],
        value="cose",
        description="Layout:",
        layout=W.Layout(width="260px"),
    )

    if include_cytoscape:
        if _HAVE_CYTO:
            cyto = CytoscapeWidget()
            cyto.set_layout(name=layout_dropdown.value)
            cyto.layout.width = "100%"
            cyto.layout.height = "520px"
            cyto_panel = W.VBox([
                W.HBox([hops, layout_dropdown, show_network]),
                cyto,
            ])
        else:
            cyto_panel = W.HTML(
                "<div style='color:#b00020; font-weight:600;'>"
                "ipycytoscape is not available in this environment, so the network view is disabled."
                "</div>"
            )

    # --- Structure viewer (optional) ---
    show_structure = W.Checkbox(value=True, description="Show structure", indent=False)

    structure_status = W.HTML(value="")
    structure_box = W.Box(layout=W.Layout(width="100%"))
    structure_legend = W.HTML(value=plddt_legend_html())


    # resolve once; can be overridden via env var PPI_PDB_DIR
    STRUCT_DIR = resolve_structures_dir()

    # Keep a handle to the current NGL widget so we can replace it cleanly
    ngl_widget = None

    def _render_structure(node_id: str) -> None:
        nonlocal ngl_widget

        structure_status.value = ""

        # Clear if disabled
        if not show_structure.value:
            structure_box.children = ()
            return

        if not _HAVE_NGL:
            structure_box.children = ()
            structure_status.value = (
                "<div style='color:#b00020; font-weight:600;'>"
                "nglview is not available in this environment."
                "</div>"
            )
            return

        pdb_path = get_pdb_path_for_node(
            model,
            node_id=node_id,
            structures_dir=STRUCT_DIR,
            pdb_path_col=None,  # later you can set to e.g. "pdb_path"
        )

        if pdb_path is None:
            structure_box.children = ()
            structure_status.value = (
                "<div style='color:#555;'>No PDB found for this node.</div>"
            )
            return

        # Create a fresh widget each time (simple and reliable)
        view = make_nglview_with_plddt_bins(pdb_path, chain_id="A")

        view.center()

        ngl_widget = view
        structure_box.children = (view,)
        structure_status.value = f"<div style='color:#555;'>Loaded: {_escape_html(pdb_path.name)}</div>"


    def _render(node_id: str) -> None:
        status.value = ""

        # summary
        try:
            row = get_node_display_row(model, node_id=node_id, fields=display_fields)
            summary.value = _format_summary_table(row, title="Node annotations")
        except Exception as e:
            summary.value = ""
            status.value = f"<div style='color:#b00020; font-weight:600;'>Error: {_escape_html(e)}</div>"
            with neighbors_out:
                neighbors_out.clear_output()
            # clear cyto if present
            if cyto is not None:
                _cyto_clear(cyto)
            return

        # neighbors table
        with neighbors_out:
            neighbors_out.clear_output()
            if show_neighbors.value:
                if node_id not in model.graph:
                    print(f"Node '{node_id}' is not present in the graph.")
                else:
                    df_n = _neighbors_df(model.graph, node_id=node_id, max_n=max_neighbors.value)
                    print(f"Neighbors shown: {len(df_n)} (cap={max_neighbors.value})")
                    display(df_n)

        # cytoscape network view
        if cyto is not None:
            if not show_network.value:
                _cyto_clear(cyto)
                return

            if node_id not in model.graph:
                _cyto_clear(cyto)
                status.value = (
                    status.value
                    + f"<div style='color:#b00020; font-weight:600;'>"
                      f"Node '{_escape_html(node_id)}' is not present in the graph."
                      f"</div>"
                )
                return

            H = _ego_subgraph(model.graph, node_id=node_id, hops=hops.value, max_nodes=200)

            # attach a couple annotation fields to nodes (useful now and later for styling)
            # (safe even if node not in annotations)
            for n in H.nodes():
                n_str = str(n)
                if n_str in model.nodes_df.index:
                    H.nodes[n]["degree_centrality"] = model.nodes_df.loc[n_str].get("degree_centrality", None)
                    H.nodes[n]["DeepTMHMM_class"] = model.nodes_df.loc[n_str].get("DeepTMHMM_class", None)

            _cyto_load_networkx(cyto, H)
            cyto.set_layout(name=layout_dropdown.value)

        # structure viewer
        _render_structure(node_id)

    def _on_node_change(change):
        if change.get("name") == "value":
            _render(change["new"])

    node_selector.observe(_on_node_change, names="value")
    max_neighbors.observe(lambda c: _render(node_selector.value), names="value")
    show_neighbors.observe(lambda c: _render(node_selector.value), names="value")
    show_structure.observe(lambda c: _render(node_selector.value), names="value")

    if cyto is not None:
        hops.observe(lambda c: _render(node_selector.value), names="value")
        layout_dropdown.observe(lambda c: _render(node_selector.value), names="value")
        show_network.observe(lambda c: _render(node_selector.value), names="value")

    controls = W.VBox(
        [
            node_selector,
            W.HBox([max_neighbors, show_neighbors]),
            status,
        ],
        layout=W.Layout(margin="0 0 10px 0"),
    )

    children = [controls, summary, neighbors_out]
    children.append(
        W.VBox(
            [
                W.HBox([show_structure]),
                structure_legend,
                structure_status,
                structure_box,
            ],
            layout=W.Layout(margin="10px 0 0 0"),
        )
    )
    if cyto_panel is not None:
        children.append(cyto_panel)

    root = W.VBox(children)
    _render(node_selector.value)

    return PPIUI(
        root=root,
        node_selector=node_selector,
        status=status,
        summary=summary,
        neighbors=neighbors_out,
        cytoscape=cyto,
    )
