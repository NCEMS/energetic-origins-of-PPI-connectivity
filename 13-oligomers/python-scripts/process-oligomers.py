#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# Required because the input node pickle may contain pint/pint-pandas objects.
import pint
import pint_pandas


NODE_COL = "node"


HOMOMER_COLUMNS_BY_SIZE = {
    2: "homodimer",
    3: "homotrimer",
    4: "homotetramer",
}

HETEROMER_COLUMNS_BY_SIZE = {
    2: "heterodimer",
    3: "heterotrimer",
    4: "heterotetramer",
    5: "heteropentamer",
    6: "heterohexamer",
    7: "heteroheptamer",
    8: "heterooctamer",
    9: "heterononamer",
    10: "heterodecamer",
}


def normalize_node_id(value) -> str | None:
    """
    Normalize Arabidopsis AGI/locus IDs to the node-level convention.

    Examples:
        At1g01010   -> AT1G01010
        AT1G01010.1 -> AT1G01010
    """
    if value is None or pd.isna(value):
        return None

    value = str(value).strip().upper()

    if not value or value.lower() in {"nan", "none", "na", "null"}:
        return None

    # Remove isoform suffix, if present.
    value = re.sub(r"\.\d+$", "", value)

    return value


def parse_mapped_node_ids(value) -> list[str]:
    """
    Parse mapped AGI IDs from the UniProt mapping file.

    Usually there is one AGI per row, but this tolerates semicolon/comma-separated
    values if they appear.
    """
    if value is None or pd.isna(value):
        return []

    out = []
    seen = set()

    for part in re.split(r"[;,]", str(value)):
        node = normalize_node_id(part)

        if node and node not in seen:
            out.append(node)
            seen.add(node)

    return out


def build_uniprot_to_node_map(map_file: str) -> dict[str, set[str]]:
    """
    Build UniProt accession -> set of AGI node IDs from the ID mapping file.

    Expected mapping format:
        ID    ID-type    mapped-ID

    We retain rows where:
        ID-type == Gene_OrderedLocusName
    """
    id_map = pd.read_csv(
        map_file,
        sep="\t",
        header=None,
        names=["ID", "ID-type", "mapped-ID"],
        dtype=str,
    )

    id_map = id_map[id_map["ID-type"] == "Gene_OrderedLocusName"].copy()

    if id_map.empty:
        raise ValueError(
            "No Gene_OrderedLocusName rows found in mapping file. "
            f"Check map_file: {map_file}"
        )

    uniprot_to_nodes = defaultdict(set)

    for _, row in id_map.iterrows():
        uniprot_id = str(row["ID"]).strip().upper()

        if not uniprot_id:
            continue

        mapped_nodes = parse_mapped_node_ids(row["mapped-ID"])

        for node in mapped_nodes:
            uniprot_to_nodes[uniprot_id].add(node)

    print("\nUniProt-to-node mapping summary:")
    print(f"  UniProt IDs with AGI mapping: {len(uniprot_to_nodes):,}")
    print(
        "  Total UniProt-to-AGI pairs: "
        f"{sum(len(v) for v in uniprot_to_nodes.values()):,}"
    )

    return uniprot_to_nodes


def parse_complex_participant(entry: str) -> tuple[str | None, str]:
    """
    Parse one Complex Portal participant entry.

    Example:
        Q9XYZ1(2) -> ("Q9XYZ1", "2")

    Stoichiometry of 0 means unknown in Complex Portal-style data.
    """
    if entry is None or pd.isna(entry):
        return None, "0"

    entry = str(entry).strip()

    if not entry:
        return None, "0"

    if "(" not in entry or ")" not in entry:
        return entry.strip().upper(), "0"

    protein_id = entry.split("(", 1)[0].strip().upper()
    stoich = entry.split("(", 1)[1].split(")", 1)[0].strip()

    if not protein_id:
        return None, stoich or "0"

    return protein_id, stoich or "0"


def compute_oligomer_size(stoichiometry: str, unknown_stoichiometry: int) -> float:
    """
    Sum stoichiometry values when all stoichiometries are known.

    Returns NaN for unknown stoichiometry complexes.
    """
    if unknown_stoichiometry == 1:
        return np.nan

    if stoichiometry is None or pd.isna(stoichiometry):
        return np.nan

    values = []

    for value in str(stoichiometry).split(";"):
        value = value.strip()

        if not value:
            continue

        try:
            values.append(int(value))
        except ValueError:
            return np.nan

    if not values:
        return np.nan

    return float(sum(values))


def build_processed_complex_table(
    oligomer_data: str,
    uniprot_to_nodes: dict[str, set[str]],
) -> pd.DataFrame:
    """
    Parse the Complex Portal oligomer table into AGI-level complexes.

    Output columns:
        complex_string
        stoichiometry
        unknown_stoichiometry
        oligomer_size
        homomer
    """
    olig_df = pd.read_csv(oligomer_data, sep="\t", dtype=str)

    required_col = "Expanded participant list"

    if required_col not in olig_df.columns:
        raise ValueError(
            f"Oligomer data is missing required column '{required_col}'. "
            f"Observed columns: {list(olig_df.columns)}"
        )

    rows = []

    n_complexes_raw = len(olig_df)
    n_complexes_with_any_mapped_node = 0
    n_unmapped_participants = 0

    for _, r in olig_df.iterrows():
        participant_field = r[required_col]

        if participant_field is None or pd.isna(participant_field):
            continue

        participant_entries = str(participant_field).split("|")

        mapped_nodes = []
        stoichiometry = []
        unknown = False
        seen_nodes = set()

        for entry in participant_entries:
            uniprot_id, stoich = parse_complex_participant(entry)

            if uniprot_id is None:
                continue

            if stoich == "0":
                unknown = True

            node_ids = uniprot_to_nodes.get(uniprot_id, set())

            if not node_ids:
                n_unmapped_participants += 1
                continue

            for node_id in sorted(node_ids):
                if node_id not in seen_nodes:
                    mapped_nodes.append(node_id)
                    stoichiometry.append(stoich)
                    seen_nodes.add(node_id)

        if not mapped_nodes:
            continue

        n_complexes_with_any_mapped_node += 1

        complex_string = ";".join(mapped_nodes)
        stoich_string = ";".join(stoichiometry)
        unknown_int = int(unknown)

        oligomer_size = compute_oligomer_size(
            stoichiometry=stoich_string,
            unknown_stoichiometry=unknown_int,
        )

        homomer = int(len(mapped_nodes) == 1)

        rows.append(
            {
                "complex_string": complex_string,
                "stoichiometry": stoich_string,
                "unknown_stoichiometry": unknown_int,
                "oligomer_size": oligomer_size,
                "homomer": homomer,
            }
        )

    complex_df = pd.DataFrame(rows)

    if complex_df.empty:
        raise ValueError(
            "No complexes contained participants that could be mapped to AGI nodes. "
            "Check the UniProt mapping file and Complex Portal participant IDs."
        )

    print("\nComplex parsing summary:")
    print(f"  Raw complexes: {n_complexes_raw:,}")
    print(f"  Complexes with at least one mapped AGI: {n_complexes_with_any_mapped_node:,}")
    print(f"  Unmapped participant entries: {n_unmapped_participants:,}")
    print(f"  Processed AGI-level complexes: {len(complex_df):,}")
    print(
        "  Processed complexes with known stoichiometry: "
        f"{complex_df['unknown_stoichiometry'].eq(0).sum():,}"
    )
    print(
        "  Processed complexes with unknown stoichiometry: "
        f"{complex_df['unknown_stoichiometry'].eq(1).sum():,}"
    )

    return complex_df


def split_complex_string(value) -> list[str]:
    """
    Split a semicolon-delimited complex_string into AGI node IDs.
    """
    if value is None or pd.isna(value):
        return []

    out = []

    for part in str(value).split(";"):
        node = normalize_node_id(part)

        if node:
            out.append(node)

    return out


def build_node_level_oligomer_features(
    nodes: pd.DataFrame,
    complex_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build node-level oligomer annotations using exact AGI membership.
    """
    nodes = nodes.copy()
    nodes[NODE_COL] = nodes[NODE_COL].apply(normalize_node_id)

    nodes = nodes[nodes[NODE_COL].notna()].copy()

    # Initialize all binary oligomer-class columns.
    for col in HOMOMER_COLUMNS_BY_SIZE.values():
        nodes[col] = 0

    for col in HETEROMER_COLUMNS_BY_SIZE.values():
        nodes[col] = 0

    nodes["in_complex"] = 0
    nodes["complex_partners"] = ""
    nodes["complex_count_unknown"] = 0
    nodes["complex_count"] = 0

    node_to_complex_ids = defaultdict(set)
    node_to_partners = defaultdict(set)
    node_to_unknown_count = defaultdict(int)
    node_to_flags = defaultdict(lambda: defaultdict(int))

    for complex_id, row in complex_df.iterrows():
        proteins = split_complex_string(row["complex_string"])

        if not proteins:
            continue

        unknown = int(row.get("unknown_stoichiometry", 0)) == 1
        homomer = int(row.get("homomer", 0)) == 1
        oligomer_size = row.get("oligomer_size", np.nan)

        if pd.notna(oligomer_size):
            oligomer_size = int(oligomer_size)

        for protein in proteins:
            node_to_complex_ids[protein].add(complex_id)
            node_to_partners[protein].update(set(proteins) - {protein})

            if unknown:
                node_to_unknown_count[protein] += 1

            if pd.isna(row.get("oligomer_size", np.nan)):
                continue

            if homomer:
                col = HOMOMER_COLUMNS_BY_SIZE.get(oligomer_size)
            else:
                col = HETEROMER_COLUMNS_BY_SIZE.get(oligomer_size)

            if col:
                node_to_flags[protein][col] = 1

    # Apply collected features to nodes.
    for col in list(HOMOMER_COLUMNS_BY_SIZE.values()) + list(HETEROMER_COLUMNS_BY_SIZE.values()):
        nodes[col] = nodes[NODE_COL].apply(lambda x: node_to_flags[x].get(col, 0))

    nodes["in_complex"] = nodes[NODE_COL].apply(
        lambda x: int(len(node_to_complex_ids.get(x, set())) > 0)
    )

    nodes["complex_partners"] = nodes[NODE_COL].apply(
        lambda x: ";".join(sorted(node_to_partners.get(x, set())))
    )

    nodes["complex_count_unknown"] = nodes[NODE_COL].apply(
        lambda x: node_to_unknown_count.get(x, 0)
    )

    nodes["complex_count"] = nodes[NODE_COL].apply(
        lambda x: len(node_to_complex_ids.get(x, set()))
    )

    print("\nNode-level oligomer annotation summary:")
    print(f"  Nodes annotated: {len(nodes):,}")
    print(f"  Nodes in at least one complex: {nodes['in_complex'].sum():,}")
    print(f"  Nodes in complexes with unknown stoichiometry: {(nodes['complex_count_unknown'] > 0).sum():,}")

    for col in list(HOMOMER_COLUMNS_BY_SIZE.values()) + list(HETEROMER_COLUMNS_BY_SIZE.values()):
        print(f"  {col}: {nodes[col].sum():,}")

    return nodes


def main():
    parser = argparse.ArgumentParser(
        description="Process Complex Portal oligomer data into node-level annotations."
    )

    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--oligomer_data", required=True)
    parser.add_argument("--map_file", required=True)
    parser.add_argument("--nodes", required=True)

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read node table.
    nodes_df = pd.read_pickle(args.nodes)

    if NODE_COL not in nodes_df.columns:
        raise ValueError(f"Input nodes file must contain column '{NODE_COL}'.")

    nodes = nodes_df[[NODE_COL]].copy()
    nodes[NODE_COL] = nodes[NODE_COL].apply(normalize_node_id)

    if nodes[NODE_COL].duplicated().any():
        examples = nodes.loc[nodes[NODE_COL].duplicated(), NODE_COL].head(10).tolist()
        raise ValueError(
            "Input nodes contain duplicated node IDs after normalization. "
            f"Examples: {examples}"
        )

    # Build mapping and processed complex table.
    uniprot_to_nodes = build_uniprot_to_node_map(args.map_file)

    complex_df = build_processed_complex_table(
        oligomer_data=args.oligomer_data,
        uniprot_to_nodes=uniprot_to_nodes,
    )

    # Save uppercase AGI-level intermediate file for inspection.
    temp_out = (
        output_dir
        / f"{args.output_prefix}-{args.organism_tag}-processed-oligomers-temp.csv"
    )

    complex_df.to_csv(temp_out, index=False)

    print(f"\nWrote processed complex table: {temp_out}")

    # Build node-level output.
    node_features = build_node_level_oligomer_features(
        nodes=nodes,
        complex_df=complex_df,
    )

    per_node_out = (
        output_dir
        / f"{args.output_prefix}-{args.organism_tag}-processed-oligomers-per-node.csv"
    )

    node_features.to_csv(per_node_out, index=False)

    print(f"\nWrote node-level oligomer table: {per_node_out}")


if __name__ == "__main__":
    main()
