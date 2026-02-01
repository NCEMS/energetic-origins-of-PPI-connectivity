from __future__ import annotations
from pathlib import Path
import python_scripts.ppi_paths as ppaths
import python_scripts.ppi_io as pio
import python_scripts.ppi_model as pmodel
import python_scripts.ppi_ui as pui

def launch_app(
    *,
    include_cytoscape: bool = True,
    ui_config_path: str = "config/ui_config.json",
):
    import json
    cfg = json.loads(Path(ui_config_path).read_text())
    field_labels = cfg.get("field_labels", {})
    annotation_sections = tuple(
        (sec["title"], tuple(sec.get("fields", [])))
        for sec in cfg.get("annotation_sections", [])
    )

    data_dir = ppaths.resolve_data_dir()
    io = pio.load_ppi(data_dir / "nodes.pkl", data_dir / "edges.csv")
    model = pmodel.prepare_model(io.nodes_df, io.edges_df)

    ui = pui.build_ui(
        model,
        include_cytoscape=include_cytoscape,
        field_labels=field_labels,
        annotation_sections=annotation_sections,
    )
    return ui
