#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


REQUIRED_LYU_COLUMNS = {"Protein", "Tm", "R2", "RSE", "gene_name"}
EXPERIMENTS = (
    ("lysate_cell_MS2", "lysate_ms2"),
    ("intact_cell_MS2", "intact_ms2"),
    ("intact_cell_MS3", "intact_ms3"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Preprocess the three Lyu et al. 2023 thermal-proteome-profiling "
            "datasets by mapping source UniProt accessions directly to the "
            "curated UniProtKB-AC column in an Arabidopsis node table."
        )
    )
    parser.add_argument("--nodes", required=True, help="Input node table (.pkl, .csv, or .tsv).")
    parser.add_argument("--lysate_ms2", required=True, help="Lyu lysate-cell MS2 CSV.")
    parser.add_argument("--intact_ms2", required=True, help="Lyu intact-cell MS2 CSV.")
    parser.add_argument("--intact_ms3", required=True, help="Lyu intact-cell MS3 CSV.")
    parser.add_argument("--output", required=True, help="Processed gene-level wide CSV.")
    parser.add_argument("--audit_output", required=True, help="All source rows with mapping results.")
    parser.add_argument("--unmapped_output", required=True, help="Rows lacking a usable node mapping.")
    parser.add_argument("--duplicate_output", required=True, help="Mapped duplicate gene/experiment rows.")
    parser.add_argument(
        "--ambiguous_output",
        required=True,
        help="Rows whose normalized UniProt accession occurs for multiple AGIs in the node table.",
    )
    parser.add_argument("--summary_output", required=True, help="Per-experiment summary TSV.")
    parser.add_argument("--node_col", default="node", help="AGI column in the node table.")
    parser.add_argument(
        "--uniprot_col",
        default="UniProtKB-AC",
        help="Curated UniProt accession column in the node table.",
    )
    return parser.parse_args()


def read_nodes(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in {".pkl", ".pickle"}:
        return pd.read_pickle(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")

    raise ValueError(
        f"Unsupported node-table format for {path}. Expected .pkl, .pickle, .csv, .tsv, or .txt."
    )


def normalize_agi(value: object) -> str | pd.NA:
    if pd.isna(value):
        return pd.NA

    text = str(value).strip().upper()
    if not text:
        return pd.NA

    return re.sub(r"\.\d+$", "", text)


def normalize_uniprot(value: object) -> str | pd.NA:
    """Return a canonical UniProt accession suitable for direct matching."""
    if pd.isna(value):
        return pd.NA

    text = str(value).strip().upper()
    if not text:
        return pd.NA

    # Accept FASTA-style values such as sp|Q9XXXX|ENTRY_ARATH.
    if "|" in text:
        tokens = [token.strip() for token in text.split("|") if token.strip()]
        if len(tokens) >= 2 and tokens[0] in {"SP", "TR"}:
            text = tokens[1]

    # The curated node table stores the selected protein accession. Treat a
    # source isoform suffix as equivalent to its canonical accession.
    return re.sub(r"-\d+$", "", text)


def split_uniprot_values(value: object) -> list[object]:
    """Defensively split a node field containing comma/semicolon-separated IDs."""
    if isinstance(value, (list, tuple, set, np.ndarray, pd.Series)):
        return list(value)
    if pd.isna(value):
        return []

    text = str(value).strip()
    if not text:
        return []

    return [token.strip() for token in re.split(r"[;,]", text) if token.strip()]


def build_curated_mapping(
    nodes_df: pd.DataFrame,
    node_col: str,
    uniprot_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = [column for column in (node_col, uniprot_col) if column not in nodes_df.columns]
    if missing:
        raise KeyError(
            f"Node table is missing required column(s): {missing}. "
            f"Available columns: {list(nodes_df.columns)}"
        )

    mapping = nodes_df[[node_col, uniprot_col]].copy()
    mapping["gene"] = mapping[node_col].map(normalize_agi)
    mapping["node_uniprot_raw"] = mapping[uniprot_col]
    mapping["uniprot"] = mapping[uniprot_col].map(split_uniprot_values)
    mapping = mapping.explode("uniprot", ignore_index=True)
    mapping["uniprot"] = mapping["uniprot"].map(normalize_uniprot)
    mapping = mapping.dropna(subset=["gene", "uniprot"])
    mapping = mapping[["gene", "uniprot", "node_uniprot_raw"]].drop_duplicates(
        subset=["gene", "uniprot"]
    )

    counts = mapping.groupby("uniprot")["gene"].nunique().rename("agi_count_for_uniprot")
    mapping = mapping.merge(counts, on="uniprot", how="left", validate="many_to_one")

    ambiguous_mapping = mapping.loc[mapping["agi_count_for_uniprot"].gt(1)].copy()
    unique_mapping = mapping.loc[mapping["agi_count_for_uniprot"].eq(1)].copy()
    unique_mapping = unique_mapping[["uniprot", "gene"]].drop_duplicates("uniprot")

    return (
        unique_mapping.sort_values("uniprot", kind="stable").reset_index(drop=True),
        ambiguous_mapping.sort_values(["uniprot", "gene"], kind="stable").reset_index(drop=True),
    )


def validate_lyu_columns(df: pd.DataFrame, path: Path) -> None:
    missing = sorted(REQUIRED_LYU_COLUMNS.difference(df.columns))
    if missing:
        raise KeyError(
            f"{path} is missing required column(s): {missing}. "
            f"Available columns: {list(df.columns)}"
        )


def process_experiment(
    path: str | Path,
    experiment_label: str,
    unique_mapping: pd.DataFrame,
    ambiguous_accessions: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    path = Path(path)
    source = pd.read_csv(path)
    validate_lyu_columns(source, path)

    source = source.copy()
    source.insert(0, "source_row", np.arange(1, len(source) + 1, dtype=int))
    source.insert(0, "source_file", path.name)
    source.insert(0, "experiment", experiment_label)

    source["source_gene_name"] = source["gene_name"]
    source["source_protein"] = source["Protein"]
    source["mapping_accession"] = source["gene_name"].map(normalize_uniprot)
    source["protein_accession"] = source["Protein"].map(normalize_uniprot)
    source["protein_gene_name_match"] = (
        source["mapping_accession"].notna()
        & source["protein_accession"].notna()
        & source["mapping_accession"].eq(source["protein_accession"])
    )

    for column in ("Tm", "R2", "RSE"):
        source[column] = pd.to_numeric(source[column], errors="coerce")

    merged = source.merge(
        unique_mapping,
        left_on="mapping_accession",
        right_on="uniprot",
        how="left",
        validate="many_to_one",
        sort=False,
    ).drop(columns=["uniprot"])

    merged["mapping_status"] = np.select(
        [
            merged["mapping_accession"].isna(),
            merged["mapping_accession"].isin(ambiguous_accessions),
            merged["gene"].isna(),
        ],
        [
            "missing_uniprot",
            "ambiguous_node_mapping",
            "unmapped_to_nodes",
        ],
        default="mapped_unique",
    )

    # Rank multiple source entries mapping to the same AGI within one
    # experiment. Prefer highest R2, then lowest RSE, then the earliest row.
    merged["selection_rank"] = pd.Series(pd.NA, index=merged.index, dtype="Int64")
    merged["selected_for_wide_table"] = False

    mapped_mask = merged["mapping_status"].eq("mapped_unique")
    mapped = merged.loc[mapped_mask].copy()
    mapped = mapped.sort_values(
        ["gene", "R2", "RSE", "source_row"],
        ascending=[True, False, True, True],
        na_position="last",
        kind="stable",
    )
    mapped["selection_rank"] = mapped.groupby("gene").cumcount() + 1
    mapped["selected_for_wide_table"] = mapped["selection_rank"].eq(1)

    merged.loc[mapped.index, "selection_rank"] = mapped["selection_rank"].astype("Int64")
    merged.loc[mapped.index, "selected_for_wide_table"] = mapped[
        "selected_for_wide_table"
    ]

    duplicate_rows = merged.loc[
        mapped_mask & merged["gene"].duplicated(keep=False)
    ].copy()

    selected = merged.loc[
        merged["selected_for_wide_table"], ["gene", "Tm", "R2", "RSE"]
    ].copy()
    selected = selected.rename(
        columns={
            "Tm": f"Lyu_{experiment_label}_Tm",
            "R2": f"Lyu_{experiment_label}_R2",
            "RSE": f"Lyu_{experiment_label}_RSE",
        }
    )

    audit_columns = [
        "experiment",
        "source_file",
        "source_row",
        "source_protein",
        "source_gene_name",
        "mapping_accession",
        "protein_accession",
        "protein_gene_name_match",
        "gene",
        "mapping_status",
        "selection_rank",
        "selected_for_wide_table",
        "Tm",
        "R2",
        "RSE",
    ]
    audit = merged[audit_columns].copy()

    summary = {
        "experiment": experiment_label,
        "source_file": path.name,
        "input_rows": len(source),
        "rows_with_mapping_accession": int(source["mapping_accession"].notna().sum()),
        "protein_gene_name_mismatches": int(
            (
                source["mapping_accession"].notna()
                & source["protein_accession"].notna()
                & ~source["protein_gene_name_match"]
            ).sum()
        ),
        "mapped_unique_rows": int(mapped_mask.sum()),
        "unmapped_rows": int(merged["mapping_status"].eq("unmapped_to_nodes").sum()),
        "missing_uniprot_rows": int(merged["mapping_status"].eq("missing_uniprot").sum()),
        "ambiguous_node_mapping_rows": int(
            merged["mapping_status"].eq("ambiguous_node_mapping").sum()
        ),
        "duplicate_mapped_rows": int(len(duplicate_rows)),
        "duplicate_gene_groups": int(
            duplicate_rows["gene"].nunique() if not duplicate_rows.empty else 0
        ),
        "selected_gene_rows": int(len(selected)),
    }

    return selected, audit, duplicate_rows, summary


def ensure_parent_dirs(paths: Iterable[str | Path]) -> None:
    for path in paths:
        Path(path).parent.mkdir(parents=True, exist_ok=True)


def main() -> None:
    args = parse_args()

    output_paths = [
        args.output,
        args.audit_output,
        args.unmapped_output,
        args.duplicate_output,
        args.ambiguous_output,
        args.summary_output,
    ]
    ensure_parent_dirs(output_paths)

    nodes_df = read_nodes(args.nodes)
    unique_mapping, ambiguous_node_mapping = build_curated_mapping(
        nodes_df,
        args.node_col,
        args.uniprot_col,
    )
    ambiguous_accessions = set(ambiguous_node_mapping["uniprot"].dropna())

    experiment_frames: list[pd.DataFrame] = []
    audit_frames: list[pd.DataFrame] = []
    duplicate_frames: list[pd.DataFrame] = []
    summaries: list[dict[str, object]] = []

    for experiment_label, argument_name in EXPERIMENTS:
        experiment_df, audit_df, duplicate_df, summary = process_experiment(
            path=getattr(args, argument_name),
            experiment_label=experiment_label,
            unique_mapping=unique_mapping,
            ambiguous_accessions=ambiguous_accessions,
        )
        experiment_frames.append(experiment_df)
        audit_frames.append(audit_df)
        duplicate_frames.append(duplicate_df)
        summaries.append(summary)

    audit = pd.concat(audit_frames, ignore_index=True)
    audit.to_csv(args.audit_output, index=False)

    unmapped = audit.loc[
        audit["mapping_status"].isin(["missing_uniprot", "unmapped_to_nodes"])
    ].copy()
    unmapped.to_csv(args.unmapped_output, index=False)

    ambiguous_source_rows = audit.loc[
        audit["mapping_status"].eq("ambiguous_node_mapping")
    ].copy()
    ambiguous_source_rows.to_csv(args.ambiguous_output, index=False)

    duplicates = pd.concat(duplicate_frames, ignore_index=True)
    duplicates.to_csv(args.duplicate_output, index=False)

    summary_df = pd.DataFrame(summaries)
    summary_df["node_mapping_unique_accessions"] = len(unique_mapping)
    summary_df["node_mapping_ambiguous_accessions"] = len(ambiguous_accessions)
    summary_df.to_csv(args.summary_output, sep="\t", index=False)

    harmonized = experiment_frames[0]
    for frame in experiment_frames[1:]:
        harmonized = harmonized.merge(
            frame,
            on="gene",
            how="outer",
            validate="one_to_one",
            sort=False,
        )

    harmonized = harmonized.sort_values("gene", kind="stable").reset_index(drop=True)
    harmonized.to_csv(args.output, index=False)

    print("Lyu et al. 2023 preprocessing complete")
    print(f"  Curated unique UniProt-to-AGI accessions: {len(unique_mapping):,}")
    print(f"  Ambiguous node-table UniProt accessions: {len(ambiguous_accessions):,}")
    print(f"  Output genes: {len(harmonized):,}")
    print(f"  Processed output: {args.output}")
    print(f"  Mapping audit: {args.audit_output}")
    print(f"  Unmapped source rows: {len(unmapped):,}")
    print(f"  Ambiguous source rows: {len(ambiguous_source_rows):,}")
    print(f"  Duplicate mapped rows: {len(duplicates):,}")
    print(f"  Summary: {args.summary_output}")


if __name__ == "__main__":
    main()
