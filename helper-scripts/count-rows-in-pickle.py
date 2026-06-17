#!/usr/bin/env python3

import argparse
import pandas as pd
import pint
import pint_pandas

def main():
    parser = argparse.ArgumentParser(
        description="Print the number of rows in a pandas pickle file."
    )
    parser.add_argument("pickle_file", help="Path to input .pkl file")
    args = parser.parse_args()

    df = pd.read_pickle(args.pickle_file)

    print(f"File: {args.pickle_file}")
    print(f"Rows: {len(df):,}")
    print(f"Columns: {len(df.columns):,}")


if __name__ == "__main__":
    main()
