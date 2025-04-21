import pandas as pd
import typing
from typing import List
from typing import Optional
import argparse


def add_DeepTMHMM(nodes_df: pd.DataFrame, path_to_3line_file: str) -> pd.DataFrame:
    """
    Carries out steps to add DeepTMHMM information to df

    Args:
        df (pd.DataFrame): input nodes_df to be annotated
        path_to_3line_file (str): file path to DeepTMHMM output

    Returns:
        pd.DataFrame
    """
    # read in 3line file as a list of lines
    with open(path_to_3line_file) as f:
        lines = f.readlines()

    # parse the list of lines
    model_output_dict = parse_model_output(lines)

    # add information to the DataFrame
    nodes_df = apply_model_results(nodes_df, model_output_dict)

    return nodes_df


def parse_model_output(lines: List[str]) -> dict[str, dict[str, str, str]]:
    """
    Takes list of lines output by DeepTMHMM and converts it to a dictionary

    Args:
        lines (List[str]): list of lines read in from a file

    Returns:
        dict[str, dict[str, str, str]]: dictionary; keys are gene IDs, values are dictionaries containing class_label, sequence, and a mask based on the DeepTMHMM output
    """

    parsed = {}
    i = 0

    # list of possible output labels from DeepTMHMM; if one of these is not found,
    # throws an error and terminates the run
    expected_labels = ["TM", "SP", "GLOB", "SP+TM", "BETA"]

    while i < len(lines):

        header = lines[i].strip()
        seq = lines[i + 1].strip()
        mask = lines[i + 2].strip()
        i += 3

        # Extract gene_id and classification (e.g., SP, TM, GLOB)
        parts = header[1:].split("|")
        gene_id = parts[0].strip()
        class_label = parts[1].strip() if len(parts) > 1 else "UNKNOWN"
        if class_label not in expected_labels:
            print(
                "The classification of this protein took on an unexpected label:",
                gene_id,
                class_label,
            )
            sys.exit()

        parsed[gene_id] = {"class": class_label, "sequence": seq, "mask": mask}

    return parsed


def apply_model_results(
    df: pd.DataFrame, model_output_dict: dict[str, dict[str, str, str]]
) -> pd.DataFrame:
    """
    Adds DeepTMHMM information to nodes_df

    Args:
        df (pd.DataFrame): input nodes_df to be annotated
        model_output_dict (dict[str, dict[str, str, str]]): dictionary; output from parse_model_output

    Returns:
        pd.DataFrame: df annotated with DeepTMHMM information
    """

    def extract_mask(gene_id: str) -> dict[str, str]:
        return model_output_dict.get(gene_id, {}).get("mask", None)

    def trim_sequence(gene_id: str) -> Optional[str]:

        data = model_output_dict.get(gene_id)
        if not data:
            return None

        sequence = data["sequence"]
        mask = data["mask"]
        class_label = data["class"]

        if class_label == "SP":
            trimmed = "".join([aa for aa, m in zip(sequence, mask) if m != "S"])
        else:
            trimmed = sequence

        return trimmed.rstrip("*")

    def extract_class(gene_id: str) -> dict[str, str]:
        return model_output_dict.get(gene_id, {}).get("class", None)

    # use apply with functions to update the input df
    df["mask"] = df["node"].apply(extract_mask)
    df["DeepTMHMM_trimmed_sequence"] = df["node"].apply(trim_sequence)
    df["DeepTMHMM_class"] = df["node"].apply(extract_class)

    return df


def main():

    parser = argparse.ArgumentParser(description="Add DeepTMHMM annotations to network")
    parser.add_argument(
        "--nodes",
        default="processed-data/0_nodes-centrality-seqs.csv",
        help="Path to the input nodes CSV file",
    )
    parser.add_argument(
        "--DeepTMHMM",
        default="DeepTMHMM-runs/s288c-results/all-predictions-s288c.3line",
        help="Path to 3line format prediction file from DeepTMHMM",
    )
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Path to output directory",
    )
    parser.add_argument(
        "--output_prefix",
        default="0_",
        help="Prefix to be applied to output file",
    )
    parser.add_argument(
        "--organism_tag", default="s288c", help="Tag to label the organism for this run"
    )
    args = parser.parse_args()

    # load the nodes csv file
    nodes_df = pd.read_csv(args.nodes)

    # insert DeepTMHMM annotations into df
    nodes_df = add_DeepTMHMM(nodes_df, args.DeepTMHMM)

    # replace missing values/nan with "None"
    nodes_df = nodes_df.replace("", "None")
    nodes_df = nodes_df.fillna("None")

    # save the updated df to file
    nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM.csv",
        index=False,
    )


if __name__ == "__main__":

    main()
