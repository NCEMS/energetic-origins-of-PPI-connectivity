from Bio import SwissProt
import pandas as pd
import typing
from typing import List
import argparse


def extract_data(input_file: str, organism: str) -> list[dict]:
    """
    Reads in uniprot_sprot.dat and outputs a list of lines with information for organism

    Args:
        input_file (str): file path to the input uniprot_sprot.dat file
        organism (str): name of the organism whose information we want to parse

    Returns:
        list[str]
    """

    def extract_functions(comments: str) -> list[str]:
        functions = []
        current = []

        for comment in comments:
            if comment.startswith("-!- FUNCTION:"):
                if current:
                    functions.append(" ".join(current).strip())
                    current = []
                current.append(comment[len("-!- FUNCTION:") :].strip())
            elif comment.startswith("CC       ") and current:
                current.append(comment.strip()[len("CC       ") :])
            else:
                if current:
                    functions.append(" ".join(current).strip())
                    current = []
        if current:
            functions.append(" ".join(current).strip())

        return functions

    def extract_localization(comments: str) -> list[str]:
        locations = []
        current = []

        for comment in comments:
            if comment.startswith("-!- SUBCELLULAR LOCATION:"):
                if current:
                    locations.append(" ".join(current).strip())
                    current = []
                current.append(comment[len("-!- SUBCELLULAR LOCATION:"):].strip())
            elif comment.startswith("CC       ") and current:
                current.append(comment.strip()[len("CC       "):])
            else:
                if current:
                    locations.append(" ".join(current).strip())
                    current = []

        if current:
            locations.append(" ".join(current).strip())

        return locations

    records_data = []

    with open(input_file) as handle:
        for record in SwissProt.parse(handle):
            if record.organism == organism:

                # extract localization information
                localization = extract_localization(record.comments)

                # extract function information
                functions = extract_functions(record.comments)

                # save info for this entry to records_data
                """
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
                        ),
                        "CrossReferences": str(record.cross_references),
                        "Localization": localization,
                        "Comments": " ".join(record.comments),
                        "Features": str(record.features),
                    }
                )
                """
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
                        ),
                        "CrossReferences": str(record.cross_references),
                        "Localization": localization,
                        "Function": functions
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

    uniprot_entries = extract_data(args.input_file, args.organism)

    df = pd.DataFrame(uniprot_entries)

    df.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main()
