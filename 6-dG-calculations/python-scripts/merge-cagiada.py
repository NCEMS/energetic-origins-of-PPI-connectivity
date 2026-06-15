import pandas as pd
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--cagiada", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)
    cagiada_df = pd.read_csv(args.cagiada)

    merged = pd.merge(nodes_df, cagiada_df, on="node", how="left")
    merged.to_pickle(args.output)


if __name__ == "__main__":
    main()
