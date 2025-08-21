#!/usr/bin/env python3
import argparse
import re
import shlex
from pathlib import Path

import pandas as pd


def load_nodes(path: str) -> pd.DataFrame:
    """Load the nodes DataFrame into a variable named `nodes`."""
    if path.endswith((".pkl", ".pickle")):
        return pd.read_pickle(path)
    return pd.read_csv(path)


def build_acc_to_orf(nodes: pd.DataFrame) -> dict:
    """Map UniProt accession -> ORF (node)."""
    df = nodes[["UniProtKB-AC", "node"]].dropna()
    # ensure string types
    df["UniProtKB-AC"] = df["UniProtKB-AC"].astype(str)
    df["node"] = df["node"].astype(str)
    return dict(df.values)


def find_accession_in_name(name: str, acc_to_orf: dict) -> str | None:
    """
    Try to extract the UniProt accession from a filename stem.

    1) Prefer the EBI AF2 pattern: AF-<ACC>-F1-model...
    2) Fallback: scan known accessions and match as a token-like substring.
    """
    m = re.search(r"AF-([A-Za-z0-9]+)-F1-model", name)
    if m:
        return m.group(1)

    # fallback: search each accession as a token within the name
    for acc in acc_to_orf.keys():
        pat = re.compile(rf"(?<![A-Za-z0-9]){re.escape(acc)}(?![A-Za-z0-9])")
        if pat.search(name):
            return acc
    return None


def main():
    ap = argparse.ArgumentParser(description="Generate mv commands to rename AF2 outputs to ORF-based names.")
    ap.add_argument("--nodes", required=True, help="Path to nodes DataFrame (.pkl or .csv)")
    ap.add_argument("--print_lists", action="store_true", help="Also print the discovered PDB and SC lists")
    args = ap.parse_args()

    # 1) Load DataFrame into variable named `nodes`
    nodes = load_nodes(args.nodes)

    # Build mapping: UniProt accession -> ORF (node)
    acc_to_orf = build_acc_to_orf(nodes)

    cwd = Path.cwd()

    # 2) Make lists of PDB and score files (matching the *_0001.{pdb,sc} pattern)
    pdb_files = sorted(cwd.glob("*_0001.pdb"))
    sc_files = sorted(cwd.glob("*_0001.sc"))

    if args.print_lists:
        print("# PDB files:")
        for p in pdb_files:
            print(p.name)
        print("# SC files:")
        for s in sc_files:
            print(s.name)
        print()

    # Helper to extract the trailing _NNNN (e.g., _0001) to preserve it
    def extract_suffix(stem: str) -> str:
        m = re.search(r"_\d{4}$", stem)
        return m.group(0) if m else ""

    # 3) For each file, print an mv command using the ORF name as the new stem
    def emit_mv(src: Path):
        stem = src.stem
        ext = src.suffix  # includes dot
        acc = find_accession_in_name(stem, acc_to_orf)
        if not acc:
            print(f"# skipping (no accession found): {src.name}")
            return
        orf = acc_to_orf.get(acc)
        if not orf:
            print(f"# skipping (accession not in nodes): {src.name} (acc={acc})")
            return

        suffix = extract_suffix(stem) or "_0001"  # default to _0001 if missing
        dst = f"{orf}{suffix}{ext}"
        # Quote paths to be safe
        print(f"mv {shlex.quote(src.name)} {shlex.quote(dst)}")

    # Generate commands for PDBs first, then SCs (order doesn’t matter)
    for f in pdb_files:
        emit_mv(f)
    for f in sc_files:
        emit_mv(f)


if __name__ == "__main__":
    main()
