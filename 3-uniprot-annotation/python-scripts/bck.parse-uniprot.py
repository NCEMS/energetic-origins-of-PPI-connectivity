from Bio import SwissProt
import pandas as pd
import typing
from typing import List
import argparse
import re


def parse_localization(localizations: List[str]) -> List[str]:
    """
    Given a list of localization strings (e.g. from extract_block),
    removes extra metadata (e.g. curly-brace content, trailing punctuation)
    and entirely removes any text starting with "Note=",
    then splits the string into a list of core localization keywords.

    For example, given an input string like:
      "Endoplasmic reticulum membrane {ECO:0000269|PubMed:14562095, ECO:0000269|Ref.8};
       Multi-pass membrane protein {ECO:0000269|PubMed:14562095, ECO:0000269|Ref.8}.
       Note=This protein lacks a nuclear localization signal, it may interact..."
    this function will return:
      ["Endoplasmic reticulum membrane", "Multi-pass membrane protein"]

    It also leaves isoform tags intact. For instance, if one part is:
      "SUBCELLULAR LOCATION: [Isoform Mitochondrial]: Mitochondrion"
    it will become:
      "[Isoform Mitochondrial]: Mitochondrion"
    """
    keywords = []
    for loc in localizations:
        # Remove the entire note block (everything from "Note=" to the end)
        loc = re.sub(r"\s*Note=.*$", "", loc)
        # Choose a delimiter: if semicolon is present, use it; otherwise, split on period.
        if ";" in loc:
            parts = [p.strip() for p in loc.split(";") if p.strip()]
        else:
            # Split on period that is followed by a space (to avoid splitting "Ref.8")
            parts = [p.strip() for p in re.split(r"\.\s+", loc) if p.strip()]
        for part in parts:
            # Remove any curly brace content
            part = re.sub(r"\s*\{.*?\}", "", part).strip()
            # If the part starts with "SUBCELLULAR LOCATION:", remove that prefix.
            if part.startswith("SUBCELLULAR LOCATION:"):
                part = part[len("SUBCELLULAR LOCATION:") :].strip()
            # Remove any trailing punctuation (commas, periods)
            part = part.rstrip(",").rstrip(".")
            if part:
                keywords.append(part)
    return keywords


def parse_functions(raw_functions: List[str]) -> List[str]:
    """
    Given a list of raw function annotation strings (each possibly containing multiple
    function annotations separated by "FUNCTION:" and extra metadata in curly braces),
    this function splits them and returns a list of cleaned function strings.

     Example:
        Input: ["The B regulatory subunit might modulate substrate selectivity and catalytic activity,
               and also might direct the localization of the catalytic enzyme to a particular subcellular compartment.
               {ECO:0000269|PubMed:9154823}. FUNCTION: Multicopy suppressor of ROX3 and HSP60.
               {ECO:0000269|PubMed:9154823}."]
        Output: ["The B regulatory subunit might modulate substrate selectivity and catalytic activity,
                and also might direct the localization of the catalytic enzyme to a particular subcellular compartment",
                "Multicopy suppressor of ROX3 and HSP60"]
    """
    cleaned_functions = []
    for raw in raw_functions:
        # Split on "FUNCTION:" with optional surrounding whitespace
        parts = re.split(r"\s*FUNCTION:\s*", raw)
        for part in parts:
            # Remove anything inside curly braces (non-greedy)
            cleaned = re.sub(r"\{.*?\}", "", part)
            # Remove extra whitespace and any trailing punctuation
            cleaned = cleaned.strip().rstrip(".")
            if cleaned:
                cleaned_functions.append(cleaned)
    return cleaned_functions


def parse_ptm(ptm_strings: List[str]) -> List[str]:
    """
    Given a list of raw PTM comment strings, this function:
      1. Removes citation info in curly braces (e.g., {ECO:...}).
      2. Splits the string on "PTM:".
      3. Retains the text before the first "PTM:" (if any) as a separate PTM entry.
      4. Strips whitespace and trailing periods.

    Example:
      Input: [
        "Extensively O-glycosylated by PMT1 and PMT2. {ECO:...}. PTM: The GPI-anchor is attached ..."
      ]
      Output: [
        "Extensively O-glycosylated by PMT1 and PMT2",
        "The GPI-anchor is attached ..."
      ]
    """
    ptms = []
    for raw in ptm_strings:
        # Remove everything in curly braces
        raw_no_citations = re.sub(r"\{.*?\}", "", raw)

        # Split on "PTM:"
        parts = raw_no_citations.split("PTM:")

        # The first part is text before the first PTM (if any)
        # Keep it if it's not empty after stripping
        first_part = parts[0].strip().rstrip(".")
        if first_part:
            ptms.append(first_part)

        # Each subsequent part is text after a PTM:
        for p in parts[1:]:
            cleaned = p.strip().rstrip(".")
            if cleaned:
                ptms.append(cleaned)

    return ptms


def extract_ft_ptms(features: List[tuple]) -> List[str]:
    ptm_keywords = {"MOD_RES", "LIPID", "DISULFID", "CARBOHYD", "CROSSLNK", "GLYCOSYLATION"}
    ptms = []
    for feature in features:
        key, start, end, description = feature
        if key in ptm_keywords:
            ptms.append(f"{key} at {start}: {description}")
    return ptms


def extract_block(comments: List[str], header: str) -> List[str]:
    """
    Extracts comment blocks that start with the specified header.
    It collects the text following the header and any continuation lines
    until a new block (one of the known headers) is encountered.
    """
    blocks = []
    # Define known headers that mark the start of a new block.
    known_headers = [
        "FUNCTION:",
        "SUBCELLULAR LOCATION:",
        "SUBUNIT:",
        "MISCELLANEOUS:",
        "SIMILARITY:",
        "INDUCTION:",
        "DOMAIN:",
        "PTM:",
        "ALTERNATIVE PRODUCTS:",
        "DISRUPTION PHENOTYPE:",
        "DEVELOPMENTAL STAGE:",
        "MASS SPECTROMETRY:",
        "SEQUENCE CAUTION:",
        "ALTERNATIVE PRODUCTS",
        "BIOTECHNOLOGY",
        "CAUTION",
        "BIOPHYSICOCHEMICAL PROPERTIES:",
        "COFACTOR:",
        "CATALYTIC ACTIVITY:",
    ]
    i = 0
    while i < len(comments):
        if comments[i].startswith(header):
            block_text = comments[i][len(header) :].strip()
            i += 1
            # Collect continuation lines until we hit a new block that starts with a different header.
            while i < len(comments) and not any(
                comments[i].startswith(h) for h in known_headers if h != header
            ):
                block_text += " " + comments[i].strip()
                i += 1
            blocks.append(block_text.strip())
        else:
            i += 1
    return blocks


def extract_data(input_file: str, organism: str) -> List[dict]:
    """
    Reads in uniprot_sprot.dat and outputs a list of dictionaries with information for the specified organism

    Args:
        input_file (str): file path to the input uniprot_sprot.dat file
        organism (str): name of the organism whose information we want to parse

    Returns:
        list[str]
    """

    records_data = []

    with open(input_file) as handle:
        for record in SwissProt.parse(handle):
            if record.organism == organism:

                # print ("DEBUG comments:", record.comments)
                # break

                # extract localization information
                localization = extract_block(record.comments, "SUBCELLULAR LOCATION:")
                localization_keywords = parse_localization(localization)

                # extract function information
                functions = extract_block(record.comments, "FUNCTION:")
                parsed_functions = parse_functions(functions)

                # extract post-translational modification information
                ptms = extract_block(record.comments, "PTM:")
                ptms_ft = extract_ft_ptms(record.features)
                parsed_ptms = parse_ptm(ptms) + ptms_ft

                # save info for this entry to records_data
                records_data.append(
                    {
                        "EntryName": record.entry_name,
                        "PrimaryAccession": record.accessions[0],
                        "AllAccessions": ";".join(record.accessions),
                        "GeneName": record.gene_name,
                        "ProteinName": record.description,
                        "Organism": record.organism,
                        "SequenceLength": record.sequence_length,
                        "Sequence": record.sequence,
                        "GO_terms": ";".join(
                            ref[2] for ref in record.cross_references if ref[0] == "GO"
                        ).split(";"),
                        "CrossReferences": str(record.cross_references),
                        "Localization": localization,
                        "localization_keywords": localization_keywords,
                        "Function": functions,
                        "parsed_functions": parsed_functions,
                        "PTMs": ptms,
                        "parsed_PTMs": parsed_ptms,
                    }
                )

    return records_data


def main():

    parser = argparse.ArgumentParser(
        description="Process UniProt database to extract entries matching organism name."
    )
    parser.add_argument(
        "--input_file",
        default="data-files/uniprot_sprot.dat",
        type=str,
        help="Path to input UniProt.dat file",
    )
    parser.add_argument(
        "--output_file",
        default="processed-data/uniprot_sprot.csv",
        type=str,
        help="Path to the output file with UniProt information in a .csv",
    )
    parser.add_argument(
        "--organism",
        default="Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast).",
        type=str,
        help="Path to the directory containing AF2 structures for structure predictions",
    )
    args = parser.parse_args()

    # extract desired UniProt entries as a dictionary
    uniprot_entries = extract_data(args.input_file, args.organism)

    # convert to a pd.DataFrame
    df = pd.DataFrame(uniprot_entries)

    # save to file for later use in annotations
    df.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main()
