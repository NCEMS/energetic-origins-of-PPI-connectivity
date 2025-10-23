import pandas as pd
import pint
import pint_pandas
import argparse


def main():

    parser = argparse.ArgumentParser(description="Merge meltome atlas data in")
    parser.add_argument("--nodes")
    parser.add_argument("--output_dir")
    parser.add_argument("--output_prefix")
    parser.add_argument("--output_suffix")
    parser.add_argument("--organism_tag")
    parser.add_argument("--meltome_data")
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)

    melt_df = pd.read_csv(args.meltome_data, sep="\t")
    melt_df.rename(
        columns={
            "Protein ID": "meltome-id",
            "Melting point [°C]": "meltome-melting-point",
            "Protein stability class": "meltome-stability-class",
            "AUC": "meltome-AUC",
            "Curve fit converged": "meltome-fit-converged",
            "R2 of curve fit": "meltome-fit-R2",
        },
        inplace=True,
    )

    melt_df["UniProtKB-AC"] = melt_df["meltome-id"].str.split("_").str[0]

    melt_df.drop(columns=["meltome-id"], inplace=True)

    nodes_df = nodes_df.merge(melt_df, on="UniProtKB-AC", how="left")

    output_path = f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-{args.output_suffix}.pkl"

    nodes_df.to_pickle(output_path)


if __name__ == "__main__":
    main()
