#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


AGI_PATTERN = re.compile(r"AT[1-5CM]G\d{5}", flags=re.IGNORECASE)

ID_COLUMN_CANDIDATES = (
    "protein_id",
    "protein",
    "gene_name",
    "uniprotkb_ac",
    "uniprotkb_id",
    "uniprot_id",
    "uniprot",
    "accession",
    "entry",
    "identifier",
    "gene",
)
TM_COLUMN_CANDIDATES = (
    "melting_point",
    "meltome_melting_point",
    "meltpoint",
    "melt_point",
    "tm",
    "temperature",
)
R2_COLUMN_CANDIDATES = (
    "fit_r2",
    "meltome_fit_r2",
    "r2",
    "r_squared",
    "rsquared",
)
ORGANISM_COLUMN_CANDIDATES = (
    "organism",
    "organism_name",
    "species",
    "species_name",
    "taxon",
    "taxon_id",
    "taxonomy_id",
)
NODE_UNIPROT_CANDIDATES = (
    "UniProtKB-ID",
    "UniProtKB-AC",
    "uniprotkb_id",
    "uniprotkb_ac",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Preprocess an Arabidopsis Meltome Atlas table, map source identifiers "
            "to locus-level AGIs, and retain a compact Tm/R2 representation."
        )
    )
    parser.add_argument("--input", required=True, help="Input Meltome Atlas CSV.")
    parser.add_argument("--nodes", required=True, help="Input node table (.pkl, .csv, or .tsv).")
    parser.add_argument(
        "--idmapping",
        required=True,
        help="UniProt ARATH_3702_idmapping.dat file.",
    )
    parser.add_argument("--output", required=True, help="Processed gene-level CSV.")
    parser.add_argument("--audit_output", required=True, help="All source rows with mapping results.")
    parser.add_argument("--unmapped_output", required=True, help="Selected rows that could not be mapped.")
    parser.add_argument("--ambiguous_output", required=True, help="Rows mapping to multiple AGIs.")
    parser.add_argument("--duplicate_output", required=True, help="Duplicate mapped AGI rows.")
    parser.add_argument("--summary_output", required=True, help="Preprocessing summary TSV.")
    parser.add_argument("--node_col", default="node", help="AGI column in the node table.")
    parser.add_argument(
        "--uniprot_col",
        default="auto",
        help="UniProt column in the node table, or 'auto'.",
    )
    parser.add_argument(
        "--id_col",
        default="auto",
        help="Source identifier column in the Meltome CSV, or 'auto'.",
    )
    parser.add_argument(
        "--tm_col",
        default="auto",
        help="Source melting-temperature column, or 'auto'.",
    )
    parser.add_argument(
        "--r2_col",
        default="auto",
        help="Source fit-quality R2 column, or 'auto'.",
    )
    parser.add_argument(
        "--organism_col",
        default="auto",
        help="Optional organism column, 'auto', or 'none'.",
    )
    parser.add_argument(
        "--organism_pattern",
        default=r"(?i)(?:arabidopsis[ _.-]*thaliana|a[._ -]*thaliana|arath|3702)",
        help="Regex used to retain Arabidopsis rows when an organism column exists.",
    )
    return parser.parse_args()


def ensure_parent_dirs(paths: Iterable[str | Path]) -> None:
    for path in paths:
        Path(path).parent.mkdir(parents=True, exist_ok=True)


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


def normalized_column_name(name: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def choose_column(
    columns: list[str],
    requested: str,
    candidates: tuple[str, ...],
    label: str,
    optional: bool = False,
) -> str | None:
    if requested.lower() == "none":
        if optional:
            return None
        raise ValueError(f"{label} cannot be set to 'none'.")

    if requested.lower() != "auto":
        if requested not in columns:
            raise KeyError(
                f"Requested {label} column '{requested}' was not found. "
                f"Available columns: {columns}"
            )
        return requested

    normalized_to_originals: dict[str, list[str]] = defaultdict(list)
    for column in columns:
        normalized_to_originals[normalized_column_name(column)].append(column)

    for candidate in candidates:
        matches = normalized_to_originals.get(normalized_column_name(candidate), [])
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValueError(
                f"Automatic {label} detection is ambiguous for candidate '{candidate}': "
                f"{matches}. Pass the desired column explicitly."
            )

    if optional:
        return None

    raise KeyError(
        f"Could not automatically identify the {label} column. Tried: {list(candidates)}. "
        f"Available columns: {columns}. Pass the column explicitly."
    )


def normalize_agi(value: object) -> str | None:
    if pd.isna(value):
        return None
    matches = AGI_PATTERN.findall(str(value).upper())
    unique = sorted(set(match.upper() for match in matches))
    return unique[0] if len(unique) == 1 else None


def extract_agis(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {match.upper() for match in AGI_PATTERN.findall(str(value).upper())}


def normalize_alias(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().upper()
    if not text:
        return None
    text = text.strip('"\'')
    return text


def strip_uniprot_isoform(alias: str) -> str:
    return re.sub(r"-\d+$", "", alias)


def alias_candidates(value: object) -> set[str]:
    """Generate aliases from accessions, entry names, and UniProt pipe records."""
    normalized = normalize_alias(value)
    if normalized is None:
        return set()

    aliases: set[str] = {normalized}
    first_token = normalized.split()[0]
    aliases.add(first_token)

    for piece in re.split(r"[|;,]", normalized):
        piece = piece.strip()
        if piece:
            aliases.add(piece)
            aliases.add(piece.split()[0])

    expanded: set[str] = set()
    for alias in aliases:
        if not alias:
            continue
        expanded.add(alias)
        expanded.add(strip_uniprot_isoform(alias))
        expanded.add(re.sub(r"\.\d+$", "", alias))

    return {alias for alias in expanded if alias}


def split_node_uniprot_values(value: object) -> list[object]:
    if isinstance(value, (list, tuple, set, np.ndarray, pd.Series)):
        return list(value)
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [token.strip() for token in re.split(r"[;,]", text) if token.strip()]


def resolve_node_uniprot_column(nodes: pd.DataFrame, requested: str) -> str:
    columns = list(nodes.columns)
    if requested.lower() != "auto":
        if requested in columns:
            return requested
        # Preserve compatibility with pipelines that use AC rather than ID.
        fallback = choose_column(
            columns,
            "auto",
            tuple(normalized_column_name(x) for x in NODE_UNIPROT_CANDIDATES),
            "node UniProt",
        )
        print(
            f"WARNING: requested node UniProt column '{requested}' was absent; "
            f"using detected column '{fallback}'."
        )
        return str(fallback)

    detected = choose_column(
        columns,
        "auto",
        tuple(normalized_column_name(x) for x in NODE_UNIPROT_CANDIDATES),
        "node UniProt",
    )
    return str(detected)


def add_alias(mapping: dict[str, set[str]], alias: object, agis: Iterable[str]) -> None:
    agi_set = {agi for agi in agis if agi}
    if not agi_set:
        return
    for candidate in alias_candidates(alias):
        mapping[candidate].update(agi_set)


def build_node_alias_map(
    nodes: pd.DataFrame,
    node_col: str,
    uniprot_col: str,
) -> dict[str, set[str]]:
    missing = [column for column in (node_col, uniprot_col) if column not in nodes.columns]
    if missing:
        raise KeyError(
            f"Node table is missing required column(s): {missing}. "
            f"Available columns: {list(nodes.columns)}"
        )

    mapping: dict[str, set[str]] = defaultdict(set)
    for _, row in nodes[[node_col, uniprot_col]].iterrows():
        gene = normalize_agi(row[node_col])
        if gene is None:
            continue
        add_alias(mapping, gene, [gene])
        for value in split_node_uniprot_values(row[uniprot_col]):
            add_alias(mapping, value, [gene])
    return mapping


def build_idmapping_alias_map(path: str | Path) -> dict[str, set[str]]:
    path = Path(path)
    cross = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["UniProtKB-AC", "ID_type", "ID"],
        dtype="string",
        usecols=[0, 1, 2],
    )
    cross = cross.dropna(subset=["UniProtKB-AC", "ID_type", "ID"])

    accession_to_agis: dict[str, set[str]] = defaultdict(set)
    for _, row in cross.iterrows():
        accession = normalize_alias(row["UniProtKB-AC"])
        if accession is None:
            continue
        id_type = str(row["ID_type"]).strip()
        agis = extract_agis(row["ID"])
        if id_type == "Gene_OrderedLocusName" or agis:
            accession_to_agis[strip_uniprot_isoform(accession)].update(agis)

    mapping: dict[str, set[str]] = defaultdict(set)
    for _, row in cross.iterrows():
        accession = normalize_alias(row["UniProtKB-AC"])
        if accession is None:
            continue
        accession = strip_uniprot_isoform(accession)
        agis = accession_to_agis.get(accession, set())
        if not agis:
            continue
        add_alias(mapping, accession, agis)
        add_alias(mapping, row["ID"], agis)
        for agi in agis:
            add_alias(mapping, agi, [agi])

    return mapping


def map_source_identifier(
    value: object,
    node_alias_map: dict[str, set[str]],
    idmapping_alias_map: dict[str, set[str]],
) -> tuple[object, str, str, str]:
    direct_agis = extract_agis(value)
    if len(direct_agis) == 1:
        gene = next(iter(direct_agis))
        return gene, "mapped_unique", "direct_agi", gene
    if len(direct_agis) > 1:
        return pd.NA, "ambiguous_identifier", "direct_agi", ";".join(sorted(direct_agis))

    aliases = sorted(alias_candidates(value))
    if not aliases:
        return pd.NA, "missing_identifier", "none", ""

    node_hits: set[str] = set()
    idmapping_hits: set[str] = set()
    hit_aliases: list[str] = []

    for alias in aliases:
        node_values = node_alias_map.get(alias, set())
        idmapping_values = idmapping_alias_map.get(alias, set())
        if node_values or idmapping_values:
            hit_aliases.append(alias)
        node_hits.update(node_values)
        idmapping_hits.update(idmapping_values)

    all_hits = node_hits | idmapping_hits
    if not all_hits:
        return pd.NA, "unmapped", "none", ""

    if node_hits and idmapping_hits:
        method = "node_and_idmapping"
    elif node_hits:
        method = "node_uniprot"
    else:
        method = "idmapping_alias"

    if len(all_hits) == 1:
        return next(iter(all_hits)), "mapped_unique", method, ";".join(hit_aliases)

    return (
        pd.NA,
        "ambiguous_identifier",
        method,
        ";".join(sorted(all_hits)),
    )


def main() -> None:
    args = parse_args()
    output_paths = [
        args.output,
        args.audit_output,
        args.unmapped_output,
        args.ambiguous_output,
        args.duplicate_output,
        args.summary_output,
    ]
    ensure_parent_dirs(output_paths)
    Path(args.output).unlink(missing_ok=True)

    input_path = Path(args.input)
    header = pd.read_csv(input_path, nrows=0)
    columns = list(header.columns)

    id_col = choose_column(columns, args.id_col, ID_COLUMN_CANDIDATES, "source identifier")
    tm_col = choose_column(columns, args.tm_col, TM_COLUMN_CANDIDATES, "melting temperature")
    r2_col = choose_column(columns, args.r2_col, R2_COLUMN_CANDIDATES, "fit R2")
    organism_col = choose_column(
        columns,
        args.organism_col,
        ORGANISM_COLUMN_CANDIDATES,
        "organism",
        optional=True,
    )

    selected_columns = [str(id_col), str(tm_col), str(r2_col)]
    if organism_col is not None:
        selected_columns.append(organism_col)
    selected_columns = list(dict.fromkeys(selected_columns))

    source = pd.read_csv(input_path, usecols=selected_columns)
    source.insert(0, "source_row", np.arange(1, len(source) + 1, dtype=int))
    source["source_identifier"] = source[str(id_col)]
    source["MeltomeAtlas_Tm"] = pd.to_numeric(source[str(tm_col)], errors="coerce")
    source["MeltomeAtlas_R2"] = pd.to_numeric(source[str(r2_col)], errors="coerce")

    if organism_col is None:
        source["source_organism"] = pd.NA
        source["organism_selected"] = True
    else:
        source["source_organism"] = source[organism_col]
        source["organism_selected"] = (
            source[organism_col]
            .astype("string")
            .str.contains(args.organism_pattern, regex=True, na=False)
        )

    nodes = read_nodes(args.nodes)
    node_uniprot_col = resolve_node_uniprot_column(nodes, args.uniprot_col)
    node_alias_map = build_node_alias_map(nodes, args.node_col, node_uniprot_col)
    idmapping_alias_map = build_idmapping_alias_map(args.idmapping)

    mapping_results = source["source_identifier"].map(
        lambda value: map_source_identifier(value, node_alias_map, idmapping_alias_map)
    )
    source["gene"] = [result[0] for result in mapping_results]
    source["mapping_status"] = [result[1] for result in mapping_results]
    source["mapping_method"] = [result[2] for result in mapping_results]
    source["mapping_detail"] = [result[3] for result in mapping_results]

    source.loc[~source["organism_selected"], "mapping_status"] = "other_organism"
    source.loc[~source["organism_selected"], "gene"] = pd.NA

    source["value_status"] = np.where(
        source["MeltomeAtlas_Tm"].isna(), "missing_tm", "usable"
    )

    audit_columns = [
        "source_row",
        "source_identifier",
        "source_organism",
        "organism_selected",
        "gene",
        "mapping_status",
        "mapping_method",
        "mapping_detail",
        "value_status",
        "MeltomeAtlas_Tm",
        "MeltomeAtlas_R2",
    ]
    audit = source[audit_columns].copy()
    audit.to_csv(args.audit_output, index=False)

    unmapped = audit.loc[
        audit["organism_selected"]
        & audit["mapping_status"].isin(["missing_identifier", "unmapped"])
    ].copy()
    unmapped.to_csv(args.unmapped_output, index=False)

    ambiguous = audit.loc[
        audit["organism_selected"]
        & audit["mapping_status"].eq("ambiguous_identifier")
    ].copy()
    ambiguous.to_csv(args.ambiguous_output, index=False)

    usable = audit.loc[
        audit["organism_selected"]
        & audit["mapping_status"].eq("mapped_unique")
        & audit["value_status"].eq("usable"),
        ["source_row", "gene", "MeltomeAtlas_Tm", "MeltomeAtlas_R2"],
    ].copy()

    duplicate_mask = usable.duplicated(subset=["gene"], keep=False)
    duplicates = usable.loc[duplicate_mask].sort_values(
        ["gene", "source_row"], kind="stable"
    )
    duplicates.to_csv(args.duplicate_output, index=False)

    summary = pd.DataFrame(
        [
            {
                "source_file": input_path.name,
                "input_rows": len(source),
                "organism_column": organism_col if organism_col is not None else "none",
                "selected_organism_rows": int(source["organism_selected"].sum()),
                "identifier_column": id_col,
                "temperature_column": tm_col,
                "r2_column": r2_col,
                "node_uniprot_column": node_uniprot_col,
                "mapped_unique_rows": int(
                    (
                        source["organism_selected"]
                        & source["mapping_status"].eq("mapped_unique")
                    ).sum()
                ),
                "unmapped_rows": len(unmapped),
                "ambiguous_rows": len(ambiguous),
                "rows_with_numeric_tm": int(
                    (source["organism_selected"] & source["MeltomeAtlas_Tm"].notna()).sum()
                ),
                "usable_rows": len(usable),
                "unique_usable_genes": int(usable["gene"].nunique()),
                "duplicate_rows": len(duplicates),
                "duplicate_genes": int(duplicates["gene"].nunique()),
                "node_aliases": len(node_alias_map),
                "idmapping_aliases": len(idmapping_alias_map),
            }
        ]
    )
    summary.to_csv(args.summary_output, sep="\t", index=False)

    # Ambiguous mappings are deliberately excluded from `usable` above and are
    # retained in the ambiguity audit. They do not block preprocessing.
    if not ambiguous.empty:
        print(
            "WARNING: "
            f"{len(ambiguous)} selected row(s) map to multiple AGIs and were "
            f"excluded; review {args.ambiguous_output}"
        )

    # Duplicate mapped AGIs remain a hard failure because selecting or averaging
    # among multiple measurements would require an explicit study-level policy.
    if not duplicates.empty:
        raise ValueError(
            "Meltome Atlas preprocessing requires review: "
            f"{len(duplicates)} rows represent "
            f"{duplicates['gene'].nunique()} duplicated AGI(s); review "
            f"{args.duplicate_output}. No duplicate values were averaged or selected."
        )

    processed = usable.drop(columns="source_row")
    processed = processed.sort_values("gene", kind="stable").reset_index(drop=True)
    processed.to_csv(args.output, index=False)

    print("Meltome Atlas preprocessing complete")
    print(f"  Identifier column: {id_col}")
    print(f"  Temperature column: {tm_col}")
    print(f"  R2 column: {r2_col}")
    print(f"  Organism column: {organism_col if organism_col else 'none'}")
    print(f"  Node UniProt column: {node_uniprot_col}")
    print(f"  Output genes: {len(processed):,}")
    print(f"  Unmapped selected rows: {len(unmapped):,}")
    print(f"  Processed output: {args.output}")
    print(f"  Audit output: {args.audit_output}")
    print(f"  Summary output: {args.summary_output}")


if __name__ == "__main__":
    main()
