import typing
from typing import Optional
from typing import Union
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np
import argparse
import re

HEADER_RE = re.compile(r"^>(\S+).*?\bIDR_START=(\d+)\b.*?\bIDR_END=(\d+)\b")


def parse_metapredict_fasta(path: str, min_len: int) -> pd.DataFrame:
    """
    Parse metapredict FASTA-like file and return a DataFrame with columns:
      - node: str (token after '>' up to first space)
      - IDR_count: int (number of IDRs with length >= min_len)  [threshold applies ONLY here]
      - IDR_sequences: dict[str, str] mapping "1","2",... -> IDR sequences (no newlines)
      - IDR_ranges:    dict[str, str] mapping "1","2",... -> "start-end" (1-based, inclusive),
                       where start = IDR_START + 1 (systematic off-by-one fix), end = IDR_END
      - N_aa_disordered: int (sum of lengths of all IDR sequences; no threshold)
    """
    # store per-node list of (sequence, start_corrected, end)
    segments: Dict[str, List[Tuple[str, int, int]]] = {}

    current_node: str | None = None
    current_start: int | None = None
    current_end: int | None = None
    current_seq_parts: List[str] = []

    def _flush_current() -> None:
        nonlocal current_node, current_start, current_end, current_seq_parts
        if current_node is not None and current_seq_parts:
            seq = "".join(current_seq_parts).replace("\n", "").replace("\r", "").replace(" ", "")
            seq = seq.replace("*", "").upper()
            # apply off-by-one fix to START here (1-based inclusive range)
            start_fixed = (current_start or 0) + 1
            end = current_end or 0
            segments.setdefault(current_node, []).append((seq, start_fixed, end))
        # reset
        current_node = None
        current_start = None
        current_end = None
        current_seq_parts = []

    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                _flush_current()
                m = HEADER_RE.match(line)
                if not m:
                    raise ValueError(f"Header line not recognized: {line}")
                current_node = m.group(1)
                current_start = int(m.group(2))
                current_end = int(m.group(3))
            else:
                if current_node is None:
                    raise ValueError("Found sequence line before any header.")
                current_seq_parts.append(line)
    _flush_current()

    rows = []
    for node, segs in segments.items():
        # keep all sequences/ranges; threshold only affects the idr_count
        idr_sequences = {str(i + 1): seq for i, (seq, _, _) in enumerate(segs)}
        idr_ranges    = {str(i + 1): f"{start}-{end}" for i, (_, start, end) in enumerate(segs)}
        idr_count     = sum(1 for (seq, _, _) in segs if len(seq) >= min_len)
        n_disordered  = sum(len(seq) for (seq, _, _) in segs)

        rows.append({
            "node": node,
            "IDR_count": idr_count,
            "IDR_sequences": idr_sequences,
            "IDR_ranges": idr_ranges,
            "N_aa_disordered": n_disordered,
        })

    return pd.DataFrame(
        rows,
        columns=["node", "IDR_count", "IDR_sequences", "IDR_ranges", "N_aa_disordered"]
    )


def get_disprot_thresholds(
    disprot: str, percentiles: List[float], organism_string: str
) -> pd.DataFrame:
    """
    Extracts thresholds for determining when a protein is disordered based on percentiles of the fraction disordered data within DisProt

    Args:
        disprot (str): path to the DisProt database file
        percentlies (List(float)]: list of percentile cutoffs to be computed
        organism_string (str): string that is used to label the organism of interest for this run in DisProt

    Returns:
        Dict[float, float] with keys as percentiles and values as the disorder fraction at each percentile
    """

    # parse DisProt to select organism of interest and unique entries by accession code
    dp = pd.read_csv(disprot, sep="\t")

    dp_unique = (
        dp.loc[dp["organism"] == organism_string]
          .drop_duplicates(subset="acc", keep="first")
          .copy()
    )

    dp_unique["disorder_content"] = (
        pd.to_numeric(dp_unique["disorder_content"], errors="coerce") / 100.0
    )

    return {round(p, 2): dp_unique["disorder_content"].quantile(p) for p in percentiles}


def main():

    parser = argparse.ArgumentParser(description="Add IDR information to PPI network")
    parser.add_argument(
        "--nodes",
        help="Path to the nodes CSV file",
    )
    parser.add_argument(
        "--idr_predictions",
        help="Path to the metapredict output file",
    )
    parser.add_argument(
        "--output_dir",
        help="Path to output directory",
    )
    parser.add_argument(
        "--output_prefix",
        help="Prefix to be applied to output file",
    )
    parser.add_argument(
        "--organism_tag", default="s288c", help="Tag to label the organism for this run"
    )
    parser.add_argument(
        "--min_idr_length",
        type=int,
        help="The minimum length of a disordered region for it to be counted",
    )
    parser.add_argument(
        "--seq_column",
        help="Name of the column from which sequence information should be used"
    )
    parser.add_argument("--disprot", help="Path to DisProt database file")
    parser.add_argument("--disprot_organism_name", help="Exact string of the organism name to match within DisProt")
    args = parser.parse_args()

    # read in the previous step's nodes_df
    nodes_df = pd.read_csv(args.nodes)

    # parse the metapredict output to get IDR sequences and residue ranges
    idr_df = parse_metapredict_fasta(args.idr_predictions, args.min_idr_length)

    # merge into nodes_df and compute fraction disordered
    nodes_df = nodes_df.merge(idr_df, on="node", how="left")
    nodes_df["N_aa_disordered"] = nodes_df["N_aa_disordered"].fillna(0).astype(int)

    # backfill IDR_count with 0's where appropriate
    nodes_df["IDR_count"] = nodes_df.get("IDR_count", 0).fillna(0).astype(int)

    # ensure empty dict when no IDRs
    for col in ["IDR_sequences", "IDR_ranges"]:
        if col in nodes_df.columns:
            nodes_df[col] = nodes_df[col].apply(lambda x: x if isinstance(x, dict) else {})

    # compute disorder fraction
    seq_len = nodes_df[args.seq_column].str.len()
    nodes_df["disorder_fraction"] = np.where(
        (seq_len > 0),
        nodes_df["N_aa_disordered"] / seq_len,
        np.nan
    ).astype("float64")

    # calculate thresholds for determining what is and is not an IDR based on DisProt database
    percentiles = np.arange(0.05, 1.05, 0.05)
    disprot_thresholds = get_disprot_thresholds(
        args.disprot, percentiles, args.disprot_organism_name
    )
    print(
        "Whether or not a protein is disordered will be determined using the following percentile thresholds from DisProt:\n",
        disprot_thresholds,
    )

    # add binary classifications of protein as disordered/not disordered based on DisProt percentiles
    for percentile, cutoff in disprot_thresholds.items():
        colname = f"is_disordered_{percentile:.2f}"
        nodes_df[colname] = nodes_df["disorder_fraction"].apply(
            lambda x: 1 if pd.notnull(x) and x >= cutoff else 0
        )

    # write the output file
    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-IDRs.pkl"
    )


if __name__ == "__main__":

    main()
