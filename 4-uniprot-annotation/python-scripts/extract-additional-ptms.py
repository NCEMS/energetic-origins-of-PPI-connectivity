import pandas as pd
import argparse


def main():

    # load arguments
    parser = argparse.ArgumentParser(
        description="Extract PTM information from proteomeXchange file"
    )
    parser.add_argument(
        "--input", required=True, help="Path to Yeast_GSB_phospho_all_prots_0125.csv"
    )
    parser.add_argument(
        "--output", required=True, help="Path to write UniProt-to-PTMs CSV"
    )
    args = parser.parse_args()

    # load the raw Phospho data
    df = pd.read_csv(args.input)

    # extract UniProt accession and add as a column to df
    df["UniProtKB-AC"] = df["Protein"].str.extract(r"\|([A-Z0-9]+)\|")

    # combine residue, position, and FLR category
    df["site"] = (
        df["PTM_residue"].astype(str)
        + df["Protein_pos"].astype(str)
        + "("
        + df["PTM_FLR_category"].astype(str)
        + ")"
    )

    # aggregate PTM sites per UniProt accession
    agg_df = (
        df.groupby("UniProtKB-AC")["site"]
        .apply(lambda x: ";".join(sorted(set(x))))
        .reset_index()
        .rename(columns={"site": "additional_PTMs"})
    )

    # save to file
    agg_df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
