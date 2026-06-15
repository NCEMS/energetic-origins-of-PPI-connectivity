import argparse
import sys
from typing import List, Optional

import numpy as np
import pandas as pd


def get_base_locus_id(protein_id: str) -> str:
    """
    Strip isoform suffix from a protein ID.

    Examples:
        AT1G01010.1 -> AT1G01010
        AT1G01020.3 -> AT1G01020
        AT1G01030   -> AT1G01030
    """
    protein_id = str(protein_id).strip()
    return protein_id.rsplit(".", 1)[0] if "." in protein_id else protein_id


def is_isoform_1(protein_id: str) -> bool:
    """
    Return True if protein_id is explicitly isoform .1.
    """
    return str(protein_id).strip().endswith(".1")


def add_DeepTMHMM(nodes_df: pd.DataFrame, path_to_3line_file: str) -> pd.DataFrame:
    """
    Add DeepTMHMM information to nodes_df.

    nodes_df["node"] is expected to contain base locus IDs with isoform suffixes
    stripped, e.g. AT1G01010. DeepTMHMM output may contain isoform-specific IDs,
    e.g. AT1G01010.1. This function maps AT1G01010.1 back to AT1G01010 and
    uses isoform .1 preferentially.
    """

    with open(path_to_3line_file) as f:
        lines = f.readlines()

    model_output_dict = parse_model_output(lines)

    nodes_df = apply_model_results(nodes_df, model_output_dict)

    return nodes_df


def parse_model_output(lines: List[str]) -> dict[str, dict[str, str]]:
    """
    Parse DeepTMHMM 3-line output.

    The returned dictionary is keyed by base locus ID, not isoform-specific ID.

    Example:
        DeepTMHMM header: >AT1G01010.1 | GLOB
        returned key:     AT1G01010
    """

    parsed = {}
    i = 0

    expected_labels = {"TM", "SP", "GLOB", "SP+TM", "BETA"}

    while i < len(lines):
        header = lines[i].strip()

        if not header:
            i += 1
            continue

        if not header.startswith(">"):
            raise ValueError(
                f"Expected DeepTMHMM header line starting with '>' at line {i + 1}, "
                f"but found: {header}"
            )

        if i + 2 >= len(lines):
            raise ValueError(
                f"Incomplete DeepTMHMM 3-line record starting at line {i + 1}: {header}"
            )

        seq = lines[i + 1].strip()
        mask = lines[i + 2].strip()
        i += 3

        # Extract protein ID and classification
        parts = header[1:].split("|")
        protein_id = parts[0].strip()
        base_locus_id = get_base_locus_id(protein_id)

        class_label = parts[1].strip() if len(parts) > 1 else "UNKNOWN"

        if class_label not in expected_labels:
            print(
                "The classification of this protein took on an unexpected label:",
                protein_id,
                class_label,
            )
            sys.exit(1)

        record = {
            "protein_id": protein_id,
            "class": class_label,
            "sequence": seq,
            "mask": mask,
        }

        # Prefer isoform .1 when multiple isoforms exist for the same base locus.
        if base_locus_id not in parsed:
            parsed[base_locus_id] = record

        elif is_isoform_1(protein_id):
            parsed[base_locus_id] = record

        elif not is_isoform_1(parsed[base_locus_id]["protein_id"]):
            # If neither the existing nor current record is .1, keep the first one.
            # This avoids arbitrary replacement among .2, .3, etc.
            pass

    return parsed


def apply_model_results(
    df: pd.DataFrame,
    model_output_dict: dict[str, dict[str, str]],
) -> pd.DataFrame:
    """
    Add DeepTMHMM mask, trimmed sequence, and class to nodes_df.
    """

    df = df.copy()

    if "node" not in df.columns:
        raise ValueError("Input DataFrame must contain a 'node' column.")

    # Ensure nodes are also base locus IDs
    df["node"] = df["node"].apply(get_base_locus_id)

    def extract_mask(gene_id: str) -> str:
        return model_output_dict.get(gene_id, {}).get("mask", np.nan)

    def trim_sequence(gene_id: str) -> Optional[str]:
        data = model_output_dict.get(gene_id)

        if not data:
            return np.nan

        sequence = data["sequence"]
        mask = data["mask"]
        class_label = data["class"]

        if class_label in ["SP", "SP+TM"]:
            trimmed = "".join(
                aa for aa, m in zip(sequence, mask)
                if m != "S"
            )
        else:
            trimmed = sequence

        return trimmed.rstrip("*")

    def extract_class(gene_id: str) -> str:
        return model_output_dict.get(gene_id, {}).get("class", np.nan)

    def extract_deeptmhmm_protein_id(gene_id: str) -> str:
        return model_output_dict.get(gene_id, {}).get("protein_id", np.nan)

    df["DeepTMHMM_protein_id"] = df["node"].apply(extract_deeptmhmm_protein_id)
    df["mask"] = df["node"].apply(extract_mask)
    df["DeepTMHMM_trimmed_sequence"] = df["node"].apply(trim_sequence)
    df["DeepTMHMM_class"] = df["node"].apply(extract_class)

    n_total = len(df)
    n_annotated = df["DeepTMHMM_class"].notna().sum()

    print(f"Total nodes: {n_total:,}")
    print(f"Nodes annotated with DeepTMHMM_class: {n_annotated:,}")

    if n_annotated == 0:
        raise ValueError(
            "No DeepTMHMM annotations were added. Check whether node IDs and "
            "DeepTMHMM IDs use compatible naming."
        )

    print("\nDeepTMHMM class counts:")
    print(df["DeepTMHMM_class"].value_counts(dropna=False).to_string())

    return df


def main():
    parser = argparse.ArgumentParser(description="Add DeepTMHMM annotations to network")

    parser.add_argument(
        "--nodes",
        required=True,
        help="Path to the input nodes CSV file",
    )
    parser.add_argument(
        "--DeepTMHMM",
        required=True,
        help="Path to 3-line format prediction file from DeepTMHMM",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Path to output directory",
    )
    parser.add_argument(
        "--output_prefix",
        required=True,
        help="Prefix to be applied to output file",
    )
    parser.add_argument(
        "--output_suffix",
        required=True,
        help="Suffix to be applied to output file",
    )
    parser.add_argument(
        "--organism_tag",
        required=True,
        help="Tag to label the organism for this run",
    )

    args = parser.parse_args()

    nodes_df = pd.read_csv(args.nodes)

    nodes_df = add_DeepTMHMM(nodes_df, args.DeepTMHMM)

    nodes_df = nodes_df.drop(columns=["mask"])

    nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-DeepTMHMM-{args.output_suffix}.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
