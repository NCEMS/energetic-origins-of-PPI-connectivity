#!/bin/bash

CONFIG_FILE=$1

#STAGES=("0-download-inputs"
#         "1-network-centrality"
#         "2-sequence-parsing"
#         "3-uniprot-annotation"
#         "4-idr-properties"
#         "5-dG-calculations"
#         "6-protein-half-life"
#         "7-protein-expression"
#         "flatten"
#       )

#STAGES=("1-network-centrality")

#STAGES=("0-download-inputs"
#         "1-network-centrality"
#         "2-sequence-parsing"
#         "3-uniprot-annotation"
#         "4-idr-properties"
#         "5-dG-calculations"
#       )

#STAGES=("6-Rosetta-scoring")

#STAGES=("8-protein-half-life"
#        "9-protein-expression"
#        "flatten")
#STAGES=("flatten")

STAGES=("1-network-centrality")

for STAGE in "${STAGES[@]}"; do
    echo -e "\nRunning stage $STAGE with config $CONFIG_FILE"
    snakemake --snakefile "$STAGE/Snakefile" --configfile "$CONFIG_FILE" --cores 50 --use-conda --conda-frontend conda || exit 1
done
