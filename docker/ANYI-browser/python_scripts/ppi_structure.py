# python_scripts/ppi_structure.py
from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np
from Bio.PDB import PDBParser

import nglview as nv
import nglview.color as ngc

from .ppi_model import PPIModel


# bins are (min_inclusive, max_exclusive) on pLDDT/B-factor
DEFAULT_PLDDT_BINS: List[Tuple[float, float, str]] = [
    (90.0, 101.0, "#0D57D3"),  # very high
    (70.0, 90.0,  "#6ACBF1"),  # confident
    (50.0, 70.0,  "#FDD93B"),  # low
    (0.0,  50.0,  "#FD7D4D"),  # very low
]

def plddt_legend_html(
    bins: Sequence[Tuple[float, float, str]] = DEFAULT_PLDDT_BINS,
    *,
    title: str = "AlphaFold2 pLDDT",
) -> str:
    """
    Return an HTML legend for the pLDDT bin scheme.
    """
    def row(color: str, label: str) -> str:
        return (
            "<div style='display:flex; align-items:center; gap:10px; margin:3px 0;'>"
            f"<div style='width:14px; height:14px; background:{color}; border:1px solid #333;'></div>"
            f"<div style='font-family: ui-sans-serif, system-ui; font-size: 13px;'>{label}</div>"
            "</div>"
        )

    # Expect bins in descending order (90+, 70–90, 50–70, <50), but works either way
    rows = []
    for lo, hi, color in bins:
        if hi >= 101:
            label = f"{lo:.0f}–100 (very high)"
        elif lo <= 0:
            label = f"< {hi:.0f} (very low)"
        else:
            label = f"{lo:.0f}–{hi:.0f} (confidence tier)"
        rows.append(row(color, label))

    return (
        "<div style='border:1px solid #ddd; padding:10px; border-radius:8px; display:inline-block;'>"
        f"<div style='font-weight:700; margin-bottom:6px; font-family: ui-sans-serif, system-ui;'>{title}</div>"
        + "".join(rows) +
        "</div>"
    )


def resolve_structures_dir(env_var: str = "PPI_PDB_DIR", default: str = "structures") -> Path:
    """
    Resolve a directory containing PDB files.

    Priority:
      1) $PPI_PDB_DIR if set
      2) ./structures relative to current working directory
    """
    env_value = os.environ.get(env_var)
    if env_value:
        return Path(env_value).expanduser().resolve()
    return (Path.cwd() / default).resolve()


def default_pdb_path(structures_dir: Path, node_id: str) -> Path:
    """Default convention: <structures_dir>/<node_id>.pdb"""
    return structures_dir / f"{node_id}.pdb"


def get_pdb_path_for_node(
    model: PPIModel,
    node_id: str,
    *,
    structures_dir: Optional[Path] = None,
    pdb_path_col: Optional[str] = None,
) -> Optional[Path]:
    """
    Return a PDB path for node_id.

    Resolution order:
      1) If pdb_path_col is provided and exists, use that column (absolute or relative).
      2) Otherwise use default convention: <structures_dir>/<node_id>.pdb

    Returns None if no file exists.
    """
    node_id = str(node_id)
    structures_dir = resolve_structures_dir() if structures_dir is None else Path(structures_dir).resolve()

    # 1) explicit per-node path in annotations
    if pdb_path_col and pdb_path_col in model.nodes_df.columns and node_id in model.nodes_df.index:
        raw = model.nodes_df.loc[node_id].get(pdb_path_col, None)
        if raw and str(raw).strip():
            p = Path(str(raw)).expanduser()
            p = (structures_dir / p).resolve() if not p.is_absolute() else p.resolve()
            if p.exists():
                return p

    # 2) default convention
    p = default_pdb_path(structures_dir, node_id).resolve()
    return p if p.exists() else None


def _compress_int_ranges(nums: Sequence[int]) -> str:
    """
    Convert [1,2,3,7,8,10] -> "1-3 or 7-8 or 10"
    """
    if not nums:
        return ""
    nums = sorted(set(nums))

    ranges = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
        else:
            ranges.append((start, prev))
            start = prev = n
    ranges.append((start, prev))

    parts = [f"{a}-{b}" if a != b else f"{a}" for a, b in ranges]
    return " or ".join(parts)

def _compress_int_ranges_with_chain(nums: Sequence[int], chain_id: str) -> str:
    """
    Convert [1,2,3,7,8,10] into "1-3:A or 7-8:A or 10:A"
    This avoids boolean precedence issues and is robust inside selection schemes.
    """
    if not nums:
        return ""
    nums = sorted(set(nums))

    ranges = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
        else:
            ranges.append((start, prev))
            start = prev = n
    ranges.append((start, prev))

    parts = []
    for a, b in ranges:
        if a == b:
            parts.append(f"{a}:{chain_id}")
        else:
            parts.append(f"{a}-{b}:{chain_id}")
    return " or ".join(parts)

def _residue_plddt_map(pdb_path: Path, chain_id: str = "A") -> Tuple[str, List[Tuple[int, float]]]:
    """
    Return (resolved_chain_id, list of (resseq, mean_bfactor)) for residues in chain_id.
    Uses mean B-factor across atoms in each residue; for AF2, B-factor ~= pLDDT.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", str(pdb_path))
    model0 = next(structure.get_models())

    chains = list(model0.get_chains())
    chain_ids = [c.id for c in chains]

    if chain_id not in chain_ids:
        if len(chains) == 1:
            chain_id = chains[0].id
            chain = chains[0]
        else:
            raise ValueError(f"Chain '{chain_id}' not found in {pdb_path}. Chains: {chain_ids}")
    else:
        chain = model0[chain_id]

    out: List[Tuple[int, float]] = []
    for res in chain.get_residues():
        hetflag, resseq, icode = res.id
        if hetflag.strip():
            continue

        bfs = [atom.get_bfactor() for atom in res.get_atoms()]
        if bfs:
            out.append((resseq, float(np.mean(bfs))))

    return chain_id, out


def register_plddt_selection_scheme(
    pdb_path: Path,
    *,
    chain_id: str = "A",
    bins: Sequence[Tuple[float, float, str]] = DEFAULT_PLDDT_BINS,
    scheme_name: Optional[str] = None,
) -> str:
    resolved_chain, res_plddt = _residue_plddt_map(pdb_path, chain_id=chain_id)

    data_list: List[List[str]] = []
    for lo, hi, color in bins:
        resnums = [r for (r, p) in res_plddt if (p >= lo and p < hi)]
        sel_ranges = _compress_int_ranges(resnums)
        if sel_ranges:
            # IMPORTANT: omit chain qualifier (single-chain AF2; most robust in selection schemes)
            data_list.append([color, sel_ranges])

    data_list.append(["#C0C0C0", "*"])

    if scheme_name is None:
        scheme_name = f"plddt_bins__{Path(pdb_path).stem}"

    scheme_id = ngc.ColormakerRegistry.add_selection_scheme(scheme_name, data_list)

    # Robust: if the wrapper returns an id, use it; otherwise fall back to the label
    return scheme_id if scheme_id is not None else scheme_name


def make_nglview_with_plddt_bins(pdb_path: Path, *, chain_id: str = "A") -> nv.NGLWidget:
    """
    Build a widget and apply binned pLDDT coloring using multiple disjoint representations.
    This avoids ColormakerRegistry, which may not work in some Jupyter frontends.
    """
    resolved_chain, res_plddt = _residue_plddt_map(pdb_path, chain_id=chain_id)

    view = nv.show_file(str(pdb_path))
    view.clear_representations()

    # Optional fallback base (grey), will be visually overwritten where bin reps exist
    view.add_cartoon(color="#C0C0C0")

    for lo, hi, color in DEFAULT_PLDDT_BINS:
        resnums = [r for (r, p) in res_plddt if (p >= lo and p < hi)]
        sel_ranges = _compress_int_ranges(resnums)
        if not sel_ranges:
            continue

        # For single-chain AF2, residue numbers alone are usually sufficient and robust.
        # If you prefer to include chain explicitly, use: f"({sel_ranges}) and :{resolved_chain}"
        selection = sel_ranges

        view.add_cartoon(selection=selection, color=color)

    view.center()
    return view
