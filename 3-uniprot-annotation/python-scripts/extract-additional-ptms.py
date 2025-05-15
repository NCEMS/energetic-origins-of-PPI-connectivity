import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to Yeast_GSB_phospho_all_prots_0125.csv")
    parser.add_argument("--output", required=True, help="Path to write UniProt-to-PTMs CSV")
    args = parser.parse_args()

    df = pd.read_csv(args.input)

    # Extract UniProt accession
    df["UniProtKB-AC"] = df["Protein"].str.extract(r"\|([A-Z0-9]+)\|")

    # Combine residue, position, and FLR category
    df["site"] = df["PTM_residue"].astype(str) + df["Protein_pos"].astype(str) + "(" + df["PTM_FLR_category"].astype(str) + ")"

    # Aggregate PTM sites per UniProt accession
    agg_df = (
        df.groupby("UniProtKB-AC")["site"]
        .apply(lambda x: ";".join(sorted(set(x))))
        .reset_index()
        .rename(columns={"site": "additional_PTMs"})
    )

    agg_df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
