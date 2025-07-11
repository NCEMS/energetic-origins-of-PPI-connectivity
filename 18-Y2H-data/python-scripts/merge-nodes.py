import pandas as pd
import argparse

def prefix_y2h_columns(df):
    df = df.rename(columns={
        col: ("Y2H" + col if col.startswith("_") else "Y2H_" + col)
        for col in df.columns if col != "node"
    })
    return df

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Merge Y2H nodes into dataset")
    parser.add_argument("--input_nodes")
    parser.add_argument("--Y2H_centralities")
    parser.add_argument("--output_dir")
    parser.add_argument("--output_prefix")
    parser.add_argument("--organism_tag")
    args = parser.parse_args()

    # load the DataFrames
    df1 = pd.read_pickle(args.input_nodes)
    df2 = pd.read_csv(args.Y2H_centralities)

    # rename columns in df2
    df2 = prefix_y2h_columns(df2)

    # merge Y2H centrality information into the overall network
    df1 = df1.merge(df2, on="node", how="left")

    # write output
    output_path = f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-SignalP-UniProt-IDRs-albatross-cider-GhoshDill-Cagiada-Rosetta-FoldX-halflife-expr-speed-PTMGPT2-LiPMS-entanglement-chaperones-oligomers-domains-essential-Y2H.pkl"
    df1.to_pickle(output_path)

if __name__ == "__main__":
    main()
