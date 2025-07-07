import pandas as pd
import argparse

def extract_edge():

    return

def main():

    # parse arguments
    parser = argparse.ArgumentParser(description="Reformat edge data")
    parser.add_argument("--input_edges")
    parser.add_argument("--from")
    parser.add_argument("--to")
    parser.add_argument("--column_name")
    parser.add_argument("--output_dir")
    parser.add_argument("--output_prefix")
    parser.add_argument("--organism_tag")
    args = parser.parse_args()

    df = pd.read_csv(args.input_edges)

    df = df.join(df["interaction"].apply(lambda x: pd.Series(extract_from_to(x, args.from, args.to))))

    df.to_csv(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-edges-reformatted.csv")

if __name__ == "__main__":
    main()

