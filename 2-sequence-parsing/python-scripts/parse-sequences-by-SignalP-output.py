import os, sys
import typing
import argparse
import pandas as pd

def add_signalP_seq(nodes_df: pd.DataFrame, signalP_predictions: str) -> pd.DataFrame:

    """
    Add a column named signalP_trimmed_sequence that has the sequence of the protein after signal sequences as predicted by signalP have been cleaved off

    Args:

        nodes_df (pd.DataFrame): the current nodes_df to be updated
        signalP_predictions (str): path to the file prediction_results.txt output by signalP

    Returns:

        pd.DataFrame: the updated nodes_df
    """

    # load the signalP predictions and remove the # from the column names
    preds = pd.read_csv(signalP_predictions, sep="\t", skiprows=1)
    preds.columns = [col.lstrip("# ").strip() for col in preds.columns]

    preds["cleavage_site_start"] = preds["CS Position"].str.extract(r'CS pos: (\d+)-')[0].astype("Int64")

    nodes_df["cleavage_site_start"] = nodes_df["node"].map(preds.set_index("ID")["cleavage_site_start"])


    nodes_df["signalP_trimmed_sequence"] = nodes_df.apply(
        lambda row: row["sequence"][row["cleavage_site_start"]:] if pd.notna(row["cleavage_site_start"]) else pd.NA,
        axis=1
    )

    return nodes_df

def main():

    parser = argparse.ArgumentParser(description="Add sequences after signal sequence cleavage with SignalP")
    parser.add_argument(
        "--nodes",
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
        "--signalP_predictions", default="processed-data/signalP/prediction_results.txt", help="Path to file containing SignalP output"
    )
    args = parser.parse_args()

    # load the previous nodes file
    nodes_df = pd.read_csv(args.nodes)

    # add signalP_trimmed_sequence as a column; will be NaN if there was no signalP prediction for this protein sequence
    nodes_df = add_signalP_seq(nodes_df, args.signalP_predictions)

    # write the updated nodes_df to file
    nodes_df.to_csv(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP.csv", index=False)

if __name__ == "__main__":

    main()
