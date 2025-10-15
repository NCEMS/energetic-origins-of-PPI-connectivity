import os
import sys
import argparse
import pandas as pd
from typing import Dict, Iterable, Any

def get_domain_source(domain_id: str, prefix_to_source: Dict[str, str]) -> str:
    """
    Determine the domain source based on prefix.

    Args:
        domain_id: The domain annotation ID (e.g., "PF00001").
        prefix_to_source: Mapping of annotation ID prefixes to source names.

    Returns:
        The source name corresponding to the prefix, or "Other" if not matched.
    """
    for prefix in prefix_to_source:
        if domain_id.startswith(prefix):
            return prefix_to_source[prefix]
    return "Other"

def build_domain_summary(group: pd.DataFrame, sources" Iterable[str]) -> pd.Series:
    """
    Build a summary of domains for a given UniProtKB-AC group.

    Args:
        group: Subset of the domain dataframe for a single UniProtKB-AC.
        sources: Iterable of source names to summarize (e.g., ["Pfam", "SMART"]).

    Returns:
        A pandas Series with counts and domain ranges for each source.
    """
    result = {}
    for source in sources:
        sub = group[group["Source"] == source]
        if sub.empty:
            continue
        domain_ranges = {
            str(i + 1): f"{row['Start']}-{row['End']}"
            for i, row in sub.reset_index(drop=True).iterrows()
        }
        result[f"{source}_Ndomains"] = len(domain_ranges)
        result[f"{source}_domains"] = domain_ranges
    return pd.Series(result)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--input_data", required=True)
    args = parser.parse_args()

    # read and format input data
    domain_df = pd.read_csv(args.input_data, sep="\t")
    domain_df.columns = [
        "UniProtKB-AC",
        "InterPro_ID",
        "Description",
        "Domain_Annotation_ID",
        "Start",
        "End"
    ]

    # map annotation prefixes to source names
    prefix_to_source = {
        "PF": "Pfam",
        "PR": "PROSITE_PROFILE",
        "PS": "PROSITE_PATTERN",
        "SM": "SMART",
        "SSF": "SUPERFAMILY",
        "TIGR": "TIGRFAM",
        "G3DSA": "Gene3D",
        "PTHR": "PANTHER"
    }

    # annotate each row with the domain source
    domain_df["Source"] = domain_df["Domain_Annotation_ID"].apply(
        lambda x: get_domain_source(x, prefix_to_source)
    )

    # compute domain summary: one row per UniProt_ID
    summary_rows = []
    for uniprot_id, group in domain_df.groupby("UniProtKB-AC"):
        summary = build_domain_summary(group, prefix_to_source.values())
        summary["UniProtKB-AC"] = uniprot_id
        summary_rows.append(summary)

    summary_df = pd.DataFrame(summary_rows)

    # reorder columns so UniProt_ID is first
    cols = ["UniProtKB-AC"] + [col for col in summary_df.columns if col != "UniProtKB-AC"]
    summary_df = summary_df[cols]

    summary_df.to_csv(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-domain-annotations-processed.csv", index=False)
    summary_df.to_pickle(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-domain-annotations-processed.pkl")

if __name__ == "__main__":
    main()
