# python_scripts/ppi_ui.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import pandas as pd
import numpy as np
import ipywidgets as W
from IPython.display import display
import networkx as nx

from .ppi_model import PPIModel, get_node_display_row
from .ppi_structure import (
    get_pdb_path_for_node,
    resolve_structures_dir,
    make_nglview_with_plddt_bins,
    plddt_legend_html,
)

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


def _hex(c: str) -> str:
    c = c.strip()
    return c if c.startswith("#") else f"#{c}"


_COLOR_BINS = [
    (0, 10,  _hex("5289C7")),
    (10, 20, _hex("7BAFDE")),
    (20, 30, _hex("4EB265")),
    (30, 40, _hex("90C987")),
    (40, 50, _hex("CAE0AB")),
    (50, 60, _hex("F7F056")),
    (60, 70, _hex("F6C141")),
    (70, 80, _hex("F1932D")),
    (80, 90, _hex("E8601C")),
    (90, 100, _hex("DC050C")),
]
_COLOR_NA = _hex("DDDDDD")
_COLOR_DEFAULT = _hex("BDBDBD")  # default node grey (non-focus)


def _color_from_percentile(p: object) -> str:
    if p is None:
        return _COLOR_NA
    try:
        if pd.isna(p):
            return _COLOR_NA
        x = float(p)
    except Exception:
        return _COLOR_NA

    # clamp into [0, 100]
    x = max(0.0, min(100.0, x))

    for i, (lo, hi, col) in enumerate(_COLOR_BINS):
        # half-open bins: [lo, hi) and last bin includes 100
        if (lo <= x < hi) or (i == len(_COLOR_BINS) - 1 and lo <= x <= hi):
            return col

    return _COLOR_NA


def _format_value(
    v: object,
    *,
    sigfig: int = 3,
    decimals: int = 2,
    sci_small: float = 1e-2,
    sci_large: float = 1e6,
) -> str:
    """
    Format values for display.
    - None/NaN -> empty string
    - Numeric -> fixed-point decimals by default, scientific only if very small/large
    - Other -> str(v)
    """
    if v is None:
        return ""

    try:
        if hasattr(pd, "isna") and pd.isna(v):
            return ""
    except Exception:
        pass

    if isinstance(v, (int, float, np.number)) and not isinstance(v, bool):
        x = float(v)
        ax = abs(x)

        if (ax != 0.0) and (ax < sci_small or ax >= sci_large):
            return f"{x:.{sigfig}e}"

        s = f"{x:.{decimals}f}"
        s = s.rstrip("0").rstrip(".")
        return s

    return str(v)


def _format_summary_table(row: pd.Series, title: str = "Selected node", field_labels: Optional[dict[str, str]] = None) -> str:
    items = [(k, row.get(k, "")) for k in row.index]
    rows_html = "\n".join(
        f"<tr>"
        f"<td style='padding:4px 10px; font-weight:600; border:1px solid #ddd;'>{_escape_html(_display_label(k, field_labels))}</td>"
        f"<td style='padding:4px 10px; border:1px solid #ddd;'>{_escape_html(_format_value(v, sigfig=3))}</td>"
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


def _display_label(field: str, field_labels: Optional[dict[str, str]] = None) -> str:
    if field_labels and field in field_labels:
        return field_labels[field]
    return field


def _format_kv_table(
    rows: list[tuple[str, object]],
    *,
    key_col_px: int = 320,
) -> str:
    trs = "\n".join(
        f"<tr>"
        f"<td style='padding:4px 10px; font-weight:600; border:1px solid #ddd; vertical-align:top; "
        f"white-space:nowrap; overflow:hidden; text-overflow:ellipsis;'>{_escape_html(k)}</td>"
        f"<td style='padding:4px 10px; border:1px solid #ddd; vertical-align:top; "
        f"overflow-wrap:anywhere; word-break:break-word;'>{_escape_html(_format_value(v, sigfig=3))}</td>"
        f"</tr>"
        for k, v in rows
    )

    return (
        "<table style='border-collapse:collapse; border:1px solid #ddd; width:100%; table-layout:fixed;'>"
        f"<colgroup><col style='width:{int(key_col_px)}px;'><col></colgroup>"
        f"{trs}"
        "</table>"
    )


def _format_annotation_sections(
    node_row: pd.Series,
    sections: tuple[tuple[str, tuple[str, ...]], ...],
    field_labels: Optional[dict[str, str]] = None,
    *,
    title: str = "NODE ANNOTATIONS",
    missing_text: str = "Not available",
) -> str:
    blocks = [f"<div style='font-weight:700; margin:6px 0 10px 0;'>{_escape_html(title)}</div>"]

    for sec_title, fields in sections:
        kv = []
        for f in fields:
            if f in node_row.index:
                val = node_row.get(f, "")
                is_missing = (
                    val is None
                    or (isinstance(val, float) and pd.isna(val))
                    or (hasattr(pd, "isna") and pd.isna(val))
                    or str(val).strip() == ""
                )
                if is_missing:
                    val = missing_text
                kv.append((_display_label(f, field_labels), val))
            else:
                kv.append((_display_label(f, field_labels), missing_text))

        blocks.append(
            f"<div style='margin:10px 0 6px 0; font-weight:700;'>{_escape_html(sec_title)}</div>"
        )
        blocks.append(_format_kv_table(kv))

    return "<div>" + "\n".join(blocks) + "</div>"


def _neighbor_ids(G: nx.Graph | nx.DiGraph, node_id: str) -> list[str]:
    node_id = str(node_id)
    if node_id not in G:
        return []

    if G.is_directed():
        nbrs = set(map(str, G.successors(node_id))) | set(map(str, G.predecessors(node_id)))
        nbrs.discard(node_id)
        return sorted(nbrs)
    else:
        return sorted(map(str, G.neighbors(node_id)))


def _neighbors_grid_df(neighbors: list[str], n_cols: int = 4) -> pd.DataFrame:
    n_cols = int(n_cols)
    if n_cols < 1:
        n_cols = 1
    pad = (-len(neighbors)) % n_cols
    items = neighbors + [""] * pad
    rows = [items[i:i + n_cols] for i in range(0, len(items), n_cols)]
    cols = [""] * n_cols
    return pd.DataFrame(rows, columns=cols)


def _ego_subgraph(
    G: nx.Graph | nx.DiGraph,
    node_id: str,
    hops: int = 1,
    max_nodes: int = 200,
) -> nx.Graph | nx.DiGraph:
    node_id = str(node_id)
    if node_id not in G:
        return G.__class__()  # empty

    hops = int(hops)
    if hops < 1:
        hops = 1
    if hops > 3:
        hops = 3

    if hops == 1:
        visited = {node_id}
        if G.is_directed():
            nbrs = list(set(map(str, G.successors(node_id))) | set(map(str, G.predecessors(node_id))))
        else:
            nbrs = list(map(str, G.neighbors(node_id)))

        # cap neighbors to avoid huge widget payloads
        if len(nbrs) > max_nodes - 1:
            nbrs = nbrs[: max_nodes - 1]

        visited |= set(nbrs)
        return G.subgraph(list(visited)).copy()

    visited = {node_id}
    frontier = {node_id}

    for _ in range(hops):
        nxt = set()
        for n in frontier:
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

    return G.subgraph(list(visited)).copy()


def _cyto_clear(cyto) -> None:
    try:
        cyto.graph.clear()
        return
    except Exception:
        pass

    try:
        cyto.graph.nodes.clear()
    except Exception:
        pass
    try:
        cyto.graph.edges.clear()
    except Exception:
        pass


def _cyto_load_networkx(cyto, H: nx.Graph | nx.DiGraph) -> None:
    _cyto_clear(cyto)

    for n in H.nodes():
        if "label" not in H.nodes[n]:
            H.nodes[n]["label"] = str(n)

    for (u, v, k) in (H.edges(keys=True) if H.is_multigraph() else [(u, v, None) for u, v in H.edges()]):
        if H.is_multigraph():
            if "id" not in H.edges[u, v, k]:
                H.edges[u, v, k]["id"] = f"{u}__{v}__{k}"
        else:
            if "id" not in H.edges[u, v]:
                H.edges[u, v]["id"] = f"{u}__{v}"

    cyto.graph.add_graph_from_networkx(H)

    cyto.set_style([
        # default nodes
        {
            "selector": "node",
            "style": {
                "label": "data(label)",
                "background-color": "data(color)",
                "border-width": 1,
                "border-color": "#4F4F4F",
                "font-size": 16,
                "width":36,
                "height":36,
            },
        },

        # focus node when NOT coloring by metric (default behavior = blue)
        {
            "selector": 'node[is_focus = "true"][color_mode = "false"]',
            "style": {
                "background-color": "#2F80ED",
                "border-width": 3,
                "border-color": "#1B4F9C",
                "font-size": 18,
                "font-weight": "bold",
                "width":48,
                "height":48,
            },
        },

        # focus node when coloring by metric (use percentile color, but keep bold)
        {
            "selector": 'node[is_focus = "true"][color_mode = "true"]',
            "style": {
                "background-color": "data(color)",
                "border-width": 3,
                "border-color": "#4F4F4F",
                "font-size": 18,
                "font-weight": "bold",
                "width":48,
                "height":48,
            },
        },

        # edges
        {
            "selector": "edge",
            "style": {
                "curve-style": "bezier",
                "line-color": "#9E9E9E",
                "width": 1,
            },
        },
    ])


def build_ui(
    model: PPIModel,
    *,
    default_node: Optional[str] = None,
    annotation_sections: Optional[tuple[tuple[str, tuple[str, ...]], ...]] = None,
    field_labels: Optional[dict[str, str]] = None,
    display_fields: Tuple[str, ...] = ("node", "degree_centrality", "DeepTMHMM_class"),
    max_neighbors_default: int = 30,
    include_cytoscape: bool = True,
) -> PPIUI:
    node_ids = list(map(str, model.nodes_df.index))

    if default_node is None:
        default_node = node_ids[0] if node_ids else ""
    default_node = str(default_node)

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

    neighbor_cols = W.Dropdown(
        options=[("3 columns", 3), ("4 columns", 4)],
        value=4,
        description="Neighbors table:",
        layout=W.Layout(width="300px"),
    )

    show_neighbors = W.Checkbox(value=True, description="Show neighbors table", indent=False)

    status = W.HTML(value="")
    summary = W.HTML(value="")
    neighbors_out = W.Output()
    neighbors_header = W.HTML("<div style='font-weight:700; margin:6px 0 6px 0;'>Neighbors</div>")

    # --- Structure viewer ---
    show_structure = W.Checkbox(value=True, description="Show structure", indent=False)
    structure_status = W.HTML(value="")
    structure_box = W.Box(layout=W.Layout(width="100%"))
    structure_legend = W.HTML(value=plddt_legend_html())
    STRUCT_DIR = resolve_structures_dir()
    ngl_widget = None

    # --- Cytoscape controls + widget (optional) ---
    network_legend = W.HTML(value="")
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
        value="grid",
        description="Layout:",
        layout=W.Layout(width="260px"),
    )

    hover_title = W.HTML("<div style='font-weight:700; margin:6px 0;'>Hover annotations</div>")
    hover_status = W.HTML("<div style='color:#777;'>Hover over a node in the network.</div>")
    hover_summary = W.HTML(value="")

    # ---- Color-by controls (mutually exclusive checkboxes) ----
    color_by_title = W.HTML("<div style='font-weight:700; margin:6px 0 6px 0;'>Color network by percentile of:</div>")

    base_cols = (
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
    percentile_map = {c: f"{c}_percentile" for c in base_cols}

    def _swatch(label: str, color: str) -> str:
        return (
            "<div style='display:flex; align-items:center; gap:8px; margin:2px 0;'>"
            f"<span style='display:inline-block; width:14px; height:14px; background:{color}; "
            "border:1px solid #666;'></span>"
            f"<span>{_escape_html(label)}</span>"
            "</div>"
        )

    def _legend_default_html() -> str:
        return (
            "<div style='margin:6px 0 0 0;'>"
            "<div style='font-weight:700; margin:0 0 4px 0;'>Legend</div>"
            f"{_swatch('Selected node', '#2F80ED')}"
            f"{_swatch('Connected nodes', _COLOR_DEFAULT)}"
            "</div>"
        )

    def _legend_percentile_html() -> str:
        # include NA as well
        items = [
            ("0–10%", _hex("5289C7")),
            ("10–20%", _hex("7BAFDE")),
            ("20–30%", _hex("4EB265")),
            ("30–40%", _hex("90C987")),
            ("40–50%", _hex("CAE0AB")),
            ("50–60%", _hex("F7F056")),
            ("60–70%", _hex("F6C141")),
            ("70–80%", _hex("F1932D")),
            ("80–90%", _hex("E8601C")),
            ("90–100%", _hex("DC050C")),
            ("Not available", _COLOR_NA),
        ]

        # two-column grid
        rows = "".join(
            "<div style='display:flex; align-items:center; gap:8px; margin:2px 0;'>"
            f"<span style='display:inline-block; width:14px; height:14px; background:{c}; border:1px solid #666;'></span>"
            f"<span>{_escape_html(lbl)}</span>"
            "</div>"
            for lbl, c in items
        )

        return (
            "<div style='margin:6px 0 0 0;'>"
            "<div style='font-weight:700; margin:0 0 4px 0;'>Legend</div>"
            "<div style='display:grid; grid-template-columns: 1fr 1fr; column-gap:18px;'>"
            f"{rows}"
            "</div>"
            "</div>"
        )


    PRETTY_COLOR_LABELS = {
        "degree_centrality": "Degree centrality",
        "betweenness_centrality": "Betweenness centrality",
        "eigenvector_centrality": "Eigenvector centrality",
        "closeness_centrality": "Closeness centrality",
        "load_centrality": "Load centrality",
        "pagerank": "PageRank",
        "information_centrality": "Information centrality",
        "median_molecules_per_cell": "Median molecules per cell",
        "villen_halflife_min": "Half-life (min)",
        "meltome-melting-point": "Melting point (°C)",
    }

    def _metric_label(base: str) -> str:
        # explicit overrides first
        if base in PRETTY_COLOR_LABELS:
            return PRETTY_COLOR_LABELS[base]
        # then config field_labels if present
        if field_labels and base in field_labels:
            return field_labels[base]
        # fallback: title-case a cleaned version
        return base.replace("_", " ").replace("-", " ").title()


    cb_none = W.Checkbox(value=True, description="None (default)", indent=False)
    color_boxes: dict[str, W.Checkbox] = {}
    for base in base_cols:
        color_boxes[base] = W.Checkbox(value=False, description=_metric_label(base), indent=False)

    def _set_exclusive(active: Optional[str]) -> None:
        if active is None:
            cb_none.value = True
            for cb in color_boxes.values():
                cb.value = False
            return
        cb_none.value = False
        for k, cb in color_boxes.items():
            cb.value = (k == active)

    def _active_color_key() -> Optional[str]:
        if cb_none.value:
            return None
        for k, cb in color_boxes.items():
            if cb.value:
                return k
        return None

    def _on_none_change(change):
        if change.get("name") == "value" and change["new"] is True:
            _set_exclusive(None)
            _render(node_selector.value)

    cb_none.observe(_on_none_change, names="value")

    def _make_box_handler(key: str):
        def _handler(change):
            if change.get("name") != "value":
                return
            if change["new"] is True:
                _set_exclusive(key)
                _render(node_selector.value)
            else:
                # if user unchecks the active one, revert to None
                if _active_color_key() is None:
                    _set_exclusive(None)
                    _render(node_selector.value)
        return _handler

    for k, cb in color_boxes.items():
        cb.observe(_make_box_handler(k), names="value")

    color_items = [cb_none] + list(color_boxes.values())  # 11 total
    grid = W.GridBox(
        children=color_items,
        layout=W.Layout(
            grid_template_columns="1fr 1fr",
            grid_gap="4px 18px",
        ),
    )

    color_by_panel = W.VBox(
        [color_by_title, grid],
        layout=W.Layout(width="100%", margin="6px 0 0 0"),
    )

    if include_cytoscape:
        if _HAVE_CYTO:
            cyto = CytoscapeWidget()
            cyto.set_layout(name=layout_dropdown.value)
            cyto.layout.width = "100%"
            cyto.layout.height = "520px"
            cyto_panel = W.VBox([
                W.HBox([hops, layout_dropdown, show_network]),
                color_by_panel,
                network_legend,
                cyto,
            ])
        else:
            cyto_panel = W.HTML(
                "<div style='color:#b00020; font-weight:600;'>"
                "ipycytoscape is not available in this environment, so the network view is disabled."
                "</div>"
            )

    # only register events if cytoscape exists
    if cyto is not None:
        def _on_cyto_mouseover(node_json):
            data = (node_json or {}).get("data", {})
            nid = data.get("id") or data.get("label")
            if nid is not None:
                _render_hover(str(nid))

        cyto.on("node", "mouseover", _on_cyto_mouseover)

    def _render_hover(node_id: str) -> None:
        try:
            if annotation_sections is None:
                row = get_node_display_row(model, node_id=node_id, fields=display_fields)
                hover_summary.value = _format_summary_table(row, title="Hovered node", field_labels=field_labels)
            else:
                all_fields = tuple(dict.fromkeys(f for _, fs in annotation_sections for f in fs))
                row2 = get_node_display_row(model, node_id=node_id, fields=all_fields)
                hover_summary.value = _format_annotation_sections(row2, annotation_sections, title="Hovered node", field_labels=field_labels)
            hover_status.value = f"<div style='color:#555;'>Hovered: {_escape_html(node_id)}</div>"
        except Exception as e:
            hover_summary.value = ""
            hover_status.value = f"<div style='color:#b00020; font-weight:600;'>Hover error: {_escape_html(e)}</div>"

    def _render_structure(node_id: str) -> None:
        nonlocal ngl_widget

        structure_status.value = ""

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
            pdb_path_col=None,
        )

        if pdb_path is None:
            structure_box.children = ()
            structure_status.value = "<div style='color:#555;'>No PDB found for this node.</div>"
            return

        view = make_nglview_with_plddt_bins(pdb_path, chain_id="A")
        view.layout.width = "100%"
        view.layout.height = "700px"
        view.center()

        ngl_widget = view
        structure_box.children = (view,)
        structure_status.value = f"<div style='color:#555;'>Loaded: {_escape_html(pdb_path.name)}</div>"

    def _render(node_id: str) -> None:
        status.value = ""

        try:
            if annotation_sections is None:
                row = get_node_display_row(model, node_id=node_id, fields=display_fields)
                summary.value = _format_summary_table(row, title="Node annotations", field_labels=field_labels)
            else:
                all_fields = tuple(dict.fromkeys(f for _, fs in annotation_sections for f in fs))
                row2 = get_node_display_row(model, node_id=node_id, fields=all_fields)
                summary.value = _format_annotation_sections(row2, annotation_sections, title="Node annotations", field_labels=field_labels)
        except Exception as e:
            summary.value = ""
            status.value = f"<div style='color:#b00020; font-weight:600;'>Error: {_escape_html(e)}</div>"
            neighbors_header.layout.display = "none"
            with neighbors_out:
                neighbors_out.clear_output()
            if cyto is not None:
                _cyto_clear(cyto)
            return

        neighbors_header.layout.display = "" if show_neighbors.value else "none"
        with neighbors_out:
            neighbors_out.clear_output()
            if show_neighbors.value:
                if node_id not in model.graph:
                    print(f"Node '{node_id}' is not present in the graph.")
                else:
                    nbrs = _neighbor_ids(model.graph, node_id)
                    print(f"Neighbors: {len(nbrs)}")
                    display(_neighbors_grid_df(nbrs, n_cols=neighbor_cols.value))

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

            active_key = _active_color_key()
            pct_col = percentile_map.get(active_key, None) if active_key else None
            have_pct = bool(pct_col) and (pct_col in model.nodes_df.columns)
            network_legend.value = _legend_percentile_html() if have_pct else _legend_default_html()
            mode = "true" if have_pct else "false"

            for n in H.nodes():
                n_str = str(n)

                H.nodes[n]["color_mode"] = mode
                H.nodes[n]["color"] = _COLOR_DEFAULT if mode == "false" else _COLOR_NA

                if n_str in model.nodes_df.index and have_pct:
                    p = model.nodes_df.loc[n_str].get(pct_col, None)
                    H.nodes[n]["color"] = _color_from_percentile(p)

            focus_id = str(node_id)
            for n in H.nodes():
                H.nodes[n]["is_focus"] = "true" if str(n) == focus_id else "false"

            _cyto_load_networkx(cyto, H)
            cyto.set_layout(name=layout_dropdown.value)

        _render_structure(node_id)

    def _on_node_change(change):
        if change.get("name") == "value":
            _render(change["new"])

    node_selector.observe(_on_node_change, names="value")
    neighbor_cols.observe(lambda c: _render(node_selector.value), names="value")
    show_neighbors.observe(lambda c: _render(node_selector.value), names="value")
    show_structure.observe(lambda c: _render(node_selector.value), names="value")

    if cyto is not None:
        hops.observe(lambda c: _render(node_selector.value), names="value")
        layout_dropdown.observe(lambda c: _render(node_selector.value), names="value")
        show_network.observe(lambda c: _render(node_selector.value), names="value")
        # note: color checkboxes call _render themselves

    controls = W.VBox(
        [
            node_selector,
            W.HBox([neighbor_cols, show_neighbors]),
            status,
        ],
        layout=W.Layout(margin="0 0 10px 0"),
    )

    structure_panel = W.VBox(
        [
            W.HBox([show_structure]),
            structure_legend,
            structure_status,
            structure_box,
        ],
        layout=W.Layout(margin="10px 0 0 0"),
    )

    children = [controls, structure_panel, summary, neighbors_header, neighbors_out]

    if cyto_panel is not None:
        children.append(
            W.VBox(
                [
                    cyto_panel,
                    W.Box(
                        [W.VBox([hover_title, hover_status, hover_summary])],
                        layout=W.Layout(width="100%", margin="10px 0 0 0"),
                    ),
                ]
            )
        )

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
