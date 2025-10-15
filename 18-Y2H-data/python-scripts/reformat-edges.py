import pandas as pd
import argparse
import re
import typing


def extract_edge(df: pd.DataFrame, col_name: str, out1: str, out2: str) -> pd.DataFrame:
    """
    Args
        df (pd.DataFrame): input edge dataframe
        col_name (str): initial column name from which from and to nodes will be extracted
        out1 (str): output name for the "from" column
        out2 (str): output name for the "to" column

    Returns
        Updated pd.DataFrame object with "from" and "to" for each node
    """

    def extract_ids(s):
        # Extract all alphanumeric strings starting with Y and followed by 6 characters (e.g., YJL092W)
        temp = s.split()
        ids = [temp[0], temp[-1]]
        return ids if len(ids) == 2 else [None, None]

    df[[out1, out2]] = df[col_name].apply(lambda x: pd.Series(extract_ids(x)))
    return df


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Reformat edge data")
    parser.add_argument("--input_edges")
    parser.add_argument("--from_name")
    parser.add_argument("--to_name")
    parser.add_argument("--column_name")
    parser.add_argument("--output_dir")
    parser.add_argument("--output_prefix")
    parser.add_argument("--organism_tag")
    args = parser.parse_args()

    # Load and process the DataFrame
    df = pd.read_csv(args.input_edges)
    df = extract_edge(df, args.column_name, args.from_name, args.to_name)

    # Write output
    output_path = f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-edges-reformatted.csv"
    df.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
