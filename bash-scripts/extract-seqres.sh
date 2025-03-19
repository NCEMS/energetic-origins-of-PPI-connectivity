#!/bin/bash

# check if the user provided a directory argument
if [[ -z "$1" ]]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# assign the first argument to DIR
DIR="$1"

# define a mapping function for 3-letter to 1-letter amino acid codes
declare -A aa_map=(
    [ALA]=A [ARG]=R [ASN]=N [ASP]=D [CYS]=C
    [GLN]=Q [GLU]=E [GLY]=G [HIS]=H [ILE]=I
    [LEU]=L [LYS]=K [MET]=M [PHE]=F [PRO]=P
    [SER]=S [THR]=T [TRP]=W [TYR]=Y [VAL]=V
)

# check if the directory exists
if [[ ! -d "$DIR" ]]; then
    echo "Error: Directory '$DIR' does not exist."
    exit 1
fi

# iterate over all PDB files in the specified directory
for pdb_file in "$DIR"/*.pdb; do
    [[ -f "$pdb_file" ]] || continue  # Skip if no PDB files found
    fasta_file="${pdb_file%.pdb}.fasta"  # Change extension to .fasta
    seq=""

    # extract SEQRES lines and process amino acid sequences
    while read -r line; do
        if [[ $line == SEQRES* ]]; then
            # extract amino acid codes (positions 5+ in SEQRES lines)
            aa_list=($(echo "$line" | awk '{for (i=5; i<=NF; i++) print $i}'))

            # convert to one-letter codes and append to sequence
            for aa in "${aa_list[@]}"; do
                seq+="${aa_map[$aa]:-X}"  # Use 'X' for unknown residues
            done
        fi
    done < "$pdb_file"

    # write the sequence to the corresponding FASTA file
    echo ">$(basename "$pdb_file" .pdb)" > "$fasta_file"
    echo "$seq" >> "$fasta_file"
done

# create a file indicating all FASTA files have been created
touch "$DIR/.all_fasta_created"

echo "FASTA files created successfully in $DIR"
