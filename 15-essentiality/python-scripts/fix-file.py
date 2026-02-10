import argparse
import sys


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_file", required=True)
    args = parser.parse_args()

    # number of columns expected per row of the file; hard-coded based on SGD file format
    expected_cols = 10

    with open(args.input_file, "r", encoding="utf-8") as infile, open(
        args.output_file, "w", encoding="utf-8"
    ) as outfile:
        for i, line in enumerate(infile):
            parts = line.rstrip("\n").split("\t")
            if len(parts) == expected_cols:
                outfile.write("\t".join(parts) + "\n")
            elif len(parts) > expected_cols:
                fixed = parts[: expected_cols - 1] + [
                    " ".join(parts[expected_cols - 1 :])
                ]
                outfile.write("\t".join(fixed) + "\n")
            else:
                print(f"Skipping malformed line {i + 1}: only {len(parts)} fields")


if __name__ == "__main__":
    main()
