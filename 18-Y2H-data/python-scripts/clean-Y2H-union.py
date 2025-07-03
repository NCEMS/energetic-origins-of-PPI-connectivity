import pandas as pd
import sys
import argparse

def fix_sgd_name(name):
    if name.startswith("Y") and len(name) == 8:
        # insert hyphen before last character
        return name[:-1] + "-" + name[-1]
    else:
        return name

def main():

    parser = argparse.ArgumentParser(description="Clean CCSB data")
    parser.add_argument("--input_edges")
    parser.add_argument("--output_edges")
    args = parser.parse_args()

    # assume the input file is tab delimited
    df = pd.read_csv(args.input_edges, sep="\t", names=["source", "target"])

    all_proteins = pd.concat([df["source"], df["target"]]).unique()

    non_sgd = [
        p for p in all_proteins
        if not (
            p.startswith("Y")
            and len(p) == 7
            and (p.endswith("W") or p.endswith("C"))
        )
    ]

    print("Proteins with non-SGD style names:")
    for p in non_sgd:
        print(p)

    """
    Other than the proteins with an ordered locus name that is just missing a hyphen,
    we have these six proteins that do not match the standard naming scheme:

    MEL1   - SGD gives YSC0019 when you search MEL1
    Q0085  - mitochondrial gene, expected
    TORF1  - no results for yeast on SGD or UniProt
    TORF19 - no results for yeast on SGD or UniProt
    TORF21 - no results for yeast on SGD or UniProt
    TORF47 - no results for yeast on SGD or UniProt
    """

    # MEL1 > YSC0019 (looked up on SGD)
    df["source"] = df["source"].replace("MEL1", "YSC0019")
    df["target"] = df["target"].replace("MEL1", "YSC0019")

    # remove rows containing TORF proteins in source or target
    #torfs = {"TORF1", "TORF19", "TORF21", "TORF47"}
    #df = df[
    #    ~(
    #        df["source"].isin(torfs) |
    #        df["target"].isin(torfs)
    #    )
    #]


    # "fix" the ordered locus names that are missing a hyphen
    df["source_fixed"] = df["source"].apply(fix_sgd_name)
    df["target_fixed"] = df["target"].apply(fix_sgd_name)

    df = df.rename(columns={"source":"old_source", 
                            "target":"old_target",
                            "source_fixed":"source",
                            "target_fixed":"target"}
                  )

    # drop the old columnets
    df = df.drop(columns=["old_source", "old_target"])

    # save the resulting modified df to file as a csv
    df.to_csv(args.output_edges, index=False)

if __name__ == "__main__":
    main()
