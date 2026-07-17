#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import pandas as pd
import py4cytoscape as p4c


def parse_args():
    parser = argparse.ArgumentParser(
        description="Test wk-shell decomposition in Cytoscape using an edge-list CSV."
    )

    parser.add_argument(
        "--edges",
        required=True,
        help="Input edge-list CSV file.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV file for node-level wk-shell results.",
    )

    parser.add_argument(
        "--source-col",
        default="source",
        help="Name of source node column in edge-list CSV. Default: source",
    )

    parser.add_argument(
        "--target-col",
        default="target",
        help="Name of target node column in edge-list CSV. Default: target",
    )

    parser.add_argument(
        "--network-title",
        default="wkshell_test_network",
        help="Cytoscape network title. Default: wkshell_test_network",
    )

    parser.add_argument(
        "--collection",
        default="wkshell_runs",
        help="Cytoscape collection name. Default: wkshell_runs",
    )

    parser.add_argument(
        "--drop-self-loops",
        action="store_true",
        help="Drop self-loop edges before importing into Cytoscape.",
    )

    return parser.parse_args()


def ensure_network_view(network_suid):
    """
    Ensure Cytoscape has a current network view for the selected network.

    The wk-shell app appears to require a current network view because it
    applies layout/coloring after calculating _wkshell.
    """
    print("Checking for existing network views...")

    try:
        views = p4c.get_network_views(network=network_suid)
        print(f"Existing network views: {views}")
    except Exception as e:
        print(f"Could not retrieve network views before creating one: {e}")
        views = []

    if not views:
        print("No network view found. Creating network view...")
        view_suid = p4c.create_view(network=network_suid)
        print(f"Created network view: {view_suid}")
    else:
        view_suid = views[0]
        print(f"Using existing network view: {view_suid}")

    print("Setting current network view...")
    p4c.set_current_view(network=network_suid)

    return view_suid


def main():
    args = parse_args()

    edge_file = Path(args.edges)
    output_file = Path(args.output)

    if not edge_file.exists():
        raise FileNotFoundError(f"Input edge file does not exist: {edge_file}")

    print("Connecting to Cytoscape...")
    p4c.cytoscape_ping()

    print(f"Reading edge list: {edge_file}")
    edges = pd.read_csv(edge_file)

    required_cols = [args.source_col, args.target_col]
    missing_cols = [col for col in required_cols if col not in edges.columns]

    if missing_cols:
        raise ValueError(
            f"Missing required edge-list column(s): {missing_cols}. "
            f"Available columns: {list(edges.columns)}"
        )

    edges = edges[required_cols].copy()
    edges = edges.dropna(subset=required_cols)

    edges[args.source_col] = edges[args.source_col].astype(str).str.strip()
    edges[args.target_col] = edges[args.target_col].astype(str).str.strip()

    edges = edges[
        (edges[args.source_col] != "")
        & (edges[args.target_col] != "")
    ].copy()

    if args.drop_self_loops:
        n_before = len(edges)
        edges = edges[edges[args.source_col] != edges[args.target_col]].copy()
        print(f"Dropped self-loops: {n_before - len(edges):,}")

    # Rename to py4cytoscape's default expected edge-table columns.
    # create_network_from_data_frames expects source, target, and interaction
    # unless overridden with source_id_list, target_id_list, and interaction_type_list.
    edges = edges.rename(
        columns={
            args.source_col: "source",
            args.target_col: "target",
        }
    )

    edges["interaction"] = "interacts_with"

    print(f"Edges after filtering: {len(edges):,}")

    if len(edges) == 0:
        raise ValueError("No edges remain after filtering.")

    print("Creating Cytoscape network...")
    network_suid = p4c.create_network_from_data_frames(
        nodes=None,
        edges=edges,
        title=args.network_title,
        collection=args.collection,
    )

    print(f"Created network SUID: {network_suid}")

    print("Setting current network...")
    p4c.set_current_network(network_suid)

    ensure_network_view(network_suid)

    print("Running wk-shell decomposition...")
    try:
        result = p4c.commands_post("wkshell decompose")
        print(result)
    except Exception as e:
        print("Warning: wkshell command reported an error.")
        print("This may occur after the node-table columns have already been written.")
        print(f"Command error was: {e}")

    print("Checking whether wk-shell columns were written...")
    available_columns = p4c.get_table_column_names(table="node")
    print("Node-table columns containing 'wks':")
    print([col for col in available_columns if "wks" in col.lower()])

    required_output_columns = ["name", "_wkshell", "_wks_percentile_bucket"]
    missing_output_columns = [
        col for col in required_output_columns
        if col not in available_columns
    ]

    if missing_output_columns:
        raise RuntimeError(
            "wk-shell columns were not written. "
            f"Missing columns: {missing_output_columns}"
        )

    print("Retrieving node table...")
    node_table = p4c.get_table_columns(
        table="node",
        columns=required_output_columns,
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    node_table.to_csv(output_file, index=False)

    print(f"Wrote wk-shell node table to: {output_file}")
    print()
    print(node_table.head())


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
