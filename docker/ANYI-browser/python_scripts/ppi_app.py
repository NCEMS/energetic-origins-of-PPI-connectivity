# python_scripts/ppi_app.py
from __future__ import annotations

from pathlib import Path

import python_scripts.ppi_paths as ppaths
import python_scripts.ppi_io as pio
import python_scripts.ppi_model as pmodel
import python_scripts.ppi_ui as pui


def _resolve_data_path(value: str, *, data_dir: Path, repo_root: Path) -> Path:
    """
    Resolve a user-configured path that may be:
      - absolute
      - relative to repo root (e.g., "data/nodes.pkl")
      - relative to data_dir (e.g., "nodes.pkl")
    Precedence for relative paths:
      1) repo_root / value if it exists
      2) data_dir  / value if it exists
      3) default to data_dir / value (for nice errors)
    """
    p = Path(value).expanduser()

    if p.is_absolute():
        return p.resolve()

    cand_repo = (repo_root / p).resolve()
    if cand_repo.exists():
        return cand_repo

    cand_data = (data_dir / p).resolve()
    if cand_data.exists():
        return cand_data

    return cand_data  # default fallback


def launch_app(
    *,
    include_cytoscape: bool = True,
    ui_config_path: str = "config/ui_config.json",
):
    import json

    config_path = Path(ui_config_path).expanduser()
    if not config_path.is_absolute():
        config_path = (Path.cwd() / config_path).resolve()

    cfg = json.loads(config_path.read_text())

    field_labels = cfg.get("field_labels", {}) if isinstance(cfg, dict) else {}
    annotation_sections = tuple(
        (sec.get("title", ""), tuple(sec.get("fields", [])))
        for sec in (cfg.get("annotation_sections", []) if isinstance(cfg, dict) else [])
    )

    # Optional: which metrics to percentile-transform AND offer for network coloring
    color_cfg = cfg.get("network_coloring", {}) if isinstance(cfg, dict) else {}
    metrics = color_cfg.get("metrics", None)

    color_metrics = None
    if isinstance(metrics, list) and len(metrics) > 0:
        color_metrics = tuple(map(str, metrics))

    # Treat repo root as the parent of the config/ directory
    # (…/ANYI-browser/config/ui_config.json -> repo_root = …/ANYI-browser)
    repo_root = config_path.parent.parent

    # Resolve data directory (uses $PPI_DATA_DIR, else ./data)
    data_dir = Path(ppaths.resolve_data_dir()).resolve()

    # Resolve nodes/edges files from config (defaults preserved)
    data_cfg = cfg.get("data", {}) if isinstance(cfg, dict) else {}
    nodes_file = data_cfg.get("nodes_file", "nodes.pkl")
    edges_file = data_cfg.get("edges_file", "edges.csv")

    nodes_path = _resolve_data_path(str(nodes_file), data_dir=data_dir, repo_root=repo_root)
    edges_path = _resolve_data_path(str(edges_file), data_dir=data_dir, repo_root=repo_root)

    io = pio.load_ppi(nodes_path, edges_path)
    model = pmodel.prepare_model(io.nodes_df, io.edges_df, percentile_columns=color_metrics)

    ui = pui.build_ui(
        model,
        include_cytoscape=include_cytoscape,
        field_labels=field_labels,
        annotation_sections=annotation_sections,
        color_metrics=color_metrics,
    )

    return ui
