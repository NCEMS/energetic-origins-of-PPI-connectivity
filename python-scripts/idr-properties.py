from sparrow import Protein
from sparrow.predictors import batch_predict
from Bio import SeqIO
import argparse
import typing
import pandas as pd


def run_albatross(fasta_file_path, nodes_df):
    """
    Run Albatross on a set of sequences found in a FASTA-format file and save the output to pd.DataFrame

    Args:
        fasta_file_path (str): path to the input fasta file
        nodes_df (pd.DataFrame): DataFrame to which information will be added

    Returns:
        Updated nodes_df with IDR information from Albatross inserted into new columns
    """

    # Write the header
    #out.write(
    #    "ID\tAsphericity\tRadius_of_gyration\tRadius_of_gyration_scaled\tEnd_to_end_distance\tEnd_to_end_distance_scaled\tScaling_exponent\tPrefactor\n"
    #)

    seqs = []
    nodes = []
    for rec in SeqIO.parse(fasta_file_path, "fasta"):
        seq = [x for x in str(rec.seq)]
        seq = [x for x in seq if x != "*"]
        seq = "".join(seq)
        P = Protein(seq)
        nodes.append(rec.id)
        seqs.append(P)

    # compute end to end distance
    end_to_end_distance = batch_predict.batch_predict(seqs, network="re")
    end_to_end_distance = [
        round(end_to_end_distance[i][1], 3) for i in range(len(nodes))
    ]

    # compute end to end distance scaled
    end_to_end_distance_scaled = batch_predict.batch_predict(seqs, network="scaled_re")
    end_to_end_distance_scaled = [
        round(end_to_end_distance_scaled[i][1], 3) for i in range(len(nodes))
    ]

    # compute radius of gyration
    radius_of_gyration = batch_predict.batch_predict(seqs, network="rg")
    radius_of_gyration = [round(radius_of_gyration[i][1], 3) for i in range(len(nodes))]

    # compute radius of gyration scaled
    radius_of_gyration_scaled = batch_predict.batch_predict(seqs, network="scaled_rg")
    radius_of_gyration_scaled = [
        round(radius_of_gyration_scaled[i][1], 3) for i in range(len(nodes))
    ]

    # compute asphericity
    asphericity = batch_predict.batch_predict(seqs, network="asphericity")
    asphericity = [round(asphericity[i][1], 3) for i in range(len(nodes))]

    # compute scaling exponent
    scaling_exponent = batch_predict.batch_predict(seqs, network="scaling_exponent")
    scaling_exponent = [round(scaling_exponent[i][1], 3) for i in range(len(nodes))]

    # compute prefactor
    prefactor = batch_predict.batch_predict(seqs, network="prefactor")
    prefactor = [round(prefactor[i][1], 3) for i in range(len(nodes))]

    # insert IDR information into a pd.DataFrame object
    idr_data = {
        "node": nodes,
        "idr_asphericity": asphericity,
        "idr_rg": radius_of_gyration,
        "idr_rg_scaled": radius_of_gyration_scaled,
        "idr_end_to_end": end_to_end_distance,
        "idr_end_to_end_scaled": end_to_end_distance_scaled,
        "idr_scaling_exponent": scaling_exponent,
        "idr_prefactor": prefactor,
    }

    print(idr_data)

    # convert to a DataFrame
    idr_df = pd.DataFrame(idr_data)

    # merge IDR information into nodes_df
    nodes_df = pd.merge(nodes_df, idr_df, left_on="node", right_on="node", how="left")

    return nodes_df


def run_cider():

    # Output file
    out = open(args.output, "w")
    out.write(
        "ID\tKappa\tLength\tFCR\tNCPR\tSHD\tSCD\t" "Molecular_weight\tF_Neg\tF_Pos\t"
    )

    # Add columns for each amino acid (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    for aa in amino_acids:
        out.write(f"Amino_Acid_Fraction_{aa}\t")

    # Add the remaining columns
    out.write(
        "Hydrophobicity\tFraction_aromatic\tFraction_aliphatic\tFraction_polar\tComplexity\tMolecular_weight\n"
    )

    # Parse the fasta file and get the IDR properties
    for rec in tqdm(SeqIO.parse(args.input, "fasta")):
        seq = str(rec.seq)
        if seq[-1] == "*":
            seq = seq[:-1]
        SeqOb = Protein(seq)

        # Retrieve each of the requested properties
        kappa = SeqOb.kappa
        length = len(SeqOb)
        fcr = SeqOb.FCR
        ncpr = SeqOb.NCPR
        scd = SeqOb.SCD
        shd = SeqOb.SHD
        molecular_weight = SeqOb.molecular_weight
        f_neg = SeqOb.fraction_negative
        f_pos = SeqOb.fraction_positive
        amino_acid_fractions = SeqOb.amino_acid_fractions

        hydrophobicity = SeqOb.hydrophobicity
        aromatic = SeqOb.fraction_aromatic
        aliphatic = SeqOb.fraction_aliphatic
        polar = SeqOb.fraction_polar
        complexity = SeqOb.complexity
        mw = SeqOb.molecular_weight

        # Writing the extracted data to the output file
        out.write(
            f"{rec.id}\t{kappa}\t{length}\t{fcr}\t{ncpr}\t{shd}\t{scd}\t{molecular_weight}\t{f_neg}\t{f_pos}\t"
        )

        # Write each amino acid fraction to the file
        for aa in amino_acids:
            out.write(
                f"{amino_acid_fractions.get(aa, 0):.4f}\t"
            )  # Use 0 if the amino acid is not present

        # Writing the remaining data to the file
        out.write(
            f"{hydrophobicity}\t{aromatic}\t{aliphatic}\t{polar}\t{complexity}\t{mw}\n"
        )

    # Close the output file
    out.close()

    return


def main():

    parser = argparse.ArgumentParser(
        description="Annotate proteins with IDR properties"
    )
    parser.add_argument("--input_fasta", help="Path to input fasta file")
    parser.add_argument(
        "--input_node_file",
        help="Path to the node file generated by previous pipeline step",
    )
    parser.add_argument("--output_dir", help="Path to output directory")
    parser.add_argument("--output_prefix", help="Prefix to be applied to output file")
    args = parser.parse_args()

    # read in the nodes_df generated so far during the pipeline
    nodes_df = pd.read_csv(args.input_node_file)

    # run albatross
    run_albatross(args.input_fasta, nodes_df)

    # run CIDER via SPARROW
    # run_cider()

    nodes_df.to_csv(
        f"{args.output_dir}/{args.output_prefix}_network_nodes_with_IDR_properties.csv",
        index=False,
    )


if __name__ == "__main__":

    main()
