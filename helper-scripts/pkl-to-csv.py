#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Convert a pandas pickle file to CSV."
    )
    parser.add_argument("input", type=Path, help="Input .pkl file")
    parser.add_argument(
        "-o", "--output", type=Path,
        help="Output CSV path (default: input filename with .csv extension)",
    )
    parser.add_argument(
        "--index", action="store_true",
        help="Include the DataFrame index in the CSV",
    )
    args = parser.parse_args()

    output = args.output or args.input.with_suffix(".csv")
    data = pd.read_pickle(args.input)

    if not isinstance(data, (pd.DataFrame, pd.Series)):
        parser.error("The pickle must contain a pandas DataFrame or Series.")

    data.to_csv(output, index=args.index)
    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
