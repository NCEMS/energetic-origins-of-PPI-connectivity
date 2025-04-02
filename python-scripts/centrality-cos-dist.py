import sys
sys.path.append("CentralityCosDist")
from centralitycosdist import CentralityCosDist
import pandas as pd
import argparse


def main():

    parser = argparse.ArgumentParser(
        description="Compute CentralityCosDist metric on network."
    )

    # CentralityCosDist/data/Seeds.tsv
    parser.add_argument("--seeds", help="Path to the seeds TSV file", required=True)

    # CentralityCosDist/data/Network_Centrality.csv
    parser.add_argument(
        "--centrality_metrics",
        help="Path to the centrality metrics CSV file",
        required=True,
    )

    args = parser.parse_args()

    # read in seeds information
    Seeds = set(open(args.seeds).read().splitlines()[1:])  # [1:] to remove header

    # read in centrality metric information
    df_centralites = pd.read_csv(args.centrality_metrics)

    # filter out seed needs for which we do not have centrality metrics
    Nodes = set(df_centralites.ID.to_list())
    Seeds = list(Nodes.intersection(Seeds))

    algorithm = CentralityCosDist(Centrality_file=args.centrality_metrics)

    algorithm.run(seed_nodes=Seeds)

    df_rank = algorithm.rank
    print(df_rank.head(10))

    print(algorithm.similarity_score.head(10))

    print(df_rank.loc[list(Seeds)])


if __name__ == "__main__":

    main()
