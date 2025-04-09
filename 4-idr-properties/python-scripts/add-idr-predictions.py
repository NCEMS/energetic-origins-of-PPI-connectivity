import metapredict as meta
import typing
from typing import Optional
from typing import Union
from typing import Dict
import pandas as pd
import numpy as np
import argparse


def predict_disorder(
    nodes_df: pd.DataFrame, disorder_threshold: float = 0.5
) -> pd.DataFrame:
    """
    Predict IDRs with metapredict and add output to nodes_df pd.DataFrame

    Args:
        nodes_df (pd.DataFrame): input DataFrame to be annotated with IDR information
        disorder_threshold (float): cutoff for when a PROTEIN is considered to be disordered

    Returns:
        pd.DataFrame
    """

    # create dictionary in format needed by metapredict
    map_nodes_to_seq = {
        k: v
        for k, v in zip(nodes_df["node"], nodes_df["trimmed_sequence"])
        if pd.notnull(v) and v.strip() != ""
    }

    # run metapredict
    disorder_predictions = meta.predict_disorder(map_nodes_to_seq)

    # add dictionary information to DataFrame
    nodes_df["disorder_predictions"] = nodes_df["node"].map(
        lambda x: disorder_predictions[x][1] if x in disorder_predictions else None
    )

    # add column with protein's fraction of disordered residues
    nodes_df["disorder_fraction"] = nodes_df["disorder_predictions"].apply(
        fraction_disordered
    )

    # add binary classification of protein as disordered/not disordered
    nodes_df["is_disordered"] = nodes_df["disorder_fraction"].apply(
        lambda x: 1 if x > disorder_threshold else 0
    )

    # return the updated DataFrame
    return nodes_df


def fraction_disordered(predictions: np.ndarray) -> Optional[float]:
    """
    Returns the fraction of residues above the disorder threshold of 0.5

    Args:
        predictions (np.ndarray): Array of prediction values, one per residue in the protein

    Returns:
        float or None: proportion of residues in the protein that are predicted to be disordered or None if no prediction was made for this protein
    """

    # this is the threshold mentioned in https://www.biorxiv.org/content/10.1101/2024.11.05.622168v1
    per_residue_disorder_cutoff = 0.5

    if predictions is None:
        return None

    return (predictions > per_residue_disorder_cutoff).sum() / len(predictions)


def count_IDRs(
    arr: Union[np.ndarray[np.float64], float, None],
    threshold: float = 0.5,
    min_length: int = 30,
) -> int:
    """
    Counts regions in an array where values are > threshold
    for at least min_length consecutive positions

    Args:
        arr (np.ndarray or float or None): The input array or a NaN placeholder.
        threshold (float): The threshold value to define a RESIDUE as disordered
        min_length (int): Minimum number of residues a region must be to count towards total

    Returns:
        int: Number of regions satisfying the threshold and length criteria
    """

    # catch instances of NaN in the array/values
    if arr is None or isinstance(arr, float) and np.isnan(arr):
        return 0

    # create a boolean array: True where value > threshold
    mask = arr > threshold
    count = 0
    current_run = 0

    for val in mask:
        if val:
            current_run += 1
        else:
            if current_run >= min_length:
                count += 1
            current_run = 0

    # check if the last run reached the threshold
    if current_run >= min_length:
        count += 1

    return count


def extract_IDR_seqs(
    row: pd.Series, threshold: float = 0.5, min_length: int = 30
) -> Optional[Dict[int, str]]:
    """
    Adds IDR sequences to input nodes_df as a Dict with keys as integers and values as sequence strings

    Args:
        row (pd.Series): input row from nodes_df
        threshold (float): cutoff above which a residue is considered to be disordered

    Returns:
        Updated nodes_df (pd.DataFrame) with IDR sequences in a dictionary
    """

    sequence = row["trimmed_sequence"]

    if sequence is None or pd.isna(sequence):
        return None

    scores = row["disorder_predictions"]
    regions = {}
    current_region = []
    region_start = None
    region_id = 1

    for i, score in enumerate(scores):
        if score > threshold:
            if region_start is None:
                region_start = i
            current_region.append(sequence[i])
        else:
            if region_start is not None and len(current_region) >= min_length:
                regions[region_id] = "".join(current_region)
                region_id += 1
            current_region = []
            region_start = None

    if region_start is not None and len(current_region) >= min_length:
        regions[region_id] = "".join(current_region)

    return regions


def main():

    parser = argparse.ArgumentParser(description="Add sequences to PPI network.")
    parser.add_argument(
        "--nodes",
        default="../1-network-centrality/processed-data/0_nodes-centrality.csv",
        help="Path to the nodes CSV file",
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
    parser.add_argument(
        "--min_idr_length",
        default=30,
        type=int,
        help="The minimum length of a disordered region for it to be counted",
    )
    parser.add_argument(
        "--disorder_threshold_aa",
        default=0.5,
        type=float,
        help="The cutoff above which a RESIDUE is considered to be disordered",
    )
    parser.add_argument(
        "--disorder_threshold_prot",
        default=0.5,
        type=float,
        help="The cutoff above which a PROTEIN is considered to be disordered",
    )
    args = parser.parse_args()

    # read in the previous step's nodes_df
    nodes_df = pd.read_csv(args.nodes)

    # use metapredict to predict IDRs
    nodes_df = predict_disorder(
        nodes_df, disorder_threshold=args.disorder_threshold_prot
    )

    # count the number of IDRs in each protein
    nodes_df["IDR_count"] = nodes_df["disorder_predictions"].apply(
        lambda x: count_IDRs(
            x, threshold=args.disorder_threshold_aa, min_length=args.min_idr_length
        )
    )

    # extract_IDR_seqs(row: pd.Series, threshold: float = 0.5, min_length: int = 30 )
    nodes_df["IDR_sequences"] = nodes_df.apply(
        lambda row: extract_IDR_seqs(
            row, threshold=args.disorder_threshold_aa, min_length=args.min_idr_length
        ),
        axis=1,
    )

    # write the output file
    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-UniProt-IDRs.pkl"
    )


if __name__ == "__main__":

    main()
