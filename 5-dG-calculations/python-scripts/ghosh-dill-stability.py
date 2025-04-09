import pandas as pd
import pint
import pint_pandas
import argparse
import numpy as np

def compute_Ghosh_Dill_dG(nodes_df: pd.DataFrame, T: float):
    """
    Compute the stability of a protein from chain length and temperature alone using Eq. 1 from Ghosh and Dill Biophys. J. 2010 equation

    Args:
        nodes_df (pd.DataFrame): output from network-analysis.py loaded as a pd.DataFrame
        T (float): temperature, units of Kelvin

    Returns:
        pd.DataFrame
    """

    # setup Pint and register with pandas
    ureg = pint.UnitRegistry()
    Q_ = ureg.Quantity
    pint_pandas.PintType.ureg = ureg

    # add units to input temperature
    T = Q_(T, "kelvin")

    # get required information
    dG_df = nodes_df[["node", "trimmed_sequence"]].copy()

    # add protein length
    dG_df["L"] = dG_df["trimmed_sequence"].str.len()

    # define constants
    T_h = Q_(373.5, "kelvin")
    T_s = Q_(385.0, "kelvin")

    # compute ΔH(L)
    dG_df["dH"] = dG_df["L"].apply(lambda L: Q_(-5.03 * L - 41.6, "kilojoule / mole"))

    # compute ΔC_p(L)*[T - T_h]
    dG_df["dCp"] = dG_df["L"].apply(
        lambda L: Q_(-0.062 * L + 0.53, "kilojoule / mole / kelvin")
    )

    # compute ΔS(L)
    dG_df["dS"] = dG_df["L"].apply(
        lambda L: Q_(-16.8 * L - 85.0, "joule / mole / kelvin").to(
            "kilojoule / mole / kelvin"
        )
    )

    # compute ΔG
    dG = dG_df.apply(
        lambda row: (
            row["dH"]
            + row["dCp"] * (T - T_h)
            - T * row["dS"]
            - T * row["dCp"] * np.log(T.magnitude / T_s.magnitude)
        ).to("kilocalorie / mole"),
        axis=1,
    )
    dG_df["dG"] = dG

    nodes_df = pd.merge(nodes_df, dG_df, on="node", how="left")

    return nodes_df

def main():

    parser = argparse.ArgumentParser(
        description="Run stability predictions using ESM inverse folding and/or Ghosh & Dill 2010 Eq. 1."
    )
    parser.add_argument(
        "--nodes", required=True, help="Output from network-analysis.py"
    )
    parser.add_argument(
        "--output_dir",
        default="processed-data",
        help="Directory where results will be saved (default: 'outputs')",
    )
    parser.add_argument("--output_prefix", default="0_", help="Prefix for output files")
    parser.add_argument(
        "--temperature",
        default=303.15,
        type=float,
        help="Temperature in kelvin for Ghosh and Dill equation",
    )
    parser.add_argument(
        "--organism_tag",
        help="Organism label for this run")
    args = parser.parse_args()

    nodes_df = pd.read_pickle(args.nodes)

    nodes_df = compute_Ghosh_Dill_dG(nodes_df, args.temperature)

    nodes_df.to_pickle(
        f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-nodes-centrality-seqs-DeepTMHMM-UniProt-IDRs-albatross-cider-GhoshDill.pkl"
    )

if __name__ == "__main__":
    main()
