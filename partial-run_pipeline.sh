#!/bin/bash

CONFIG_FILE=$1

#STAGES=("0-download-inputs")
#STAGES=("1-network-centrality")
#STAGES=("2-sequence-parsing")
#STAGES=("3-uniprot-annotation")
#STAGES=("4-idr-properties")
#STAGES=("5-dG-calculations")
#STAGES=("6-Rosetta-scoring")
#STAGES=("7-FoldX-scoring")
#STAGES=("8-protein-half-life")
#STAGES=("9-protein-expression")
#STAGES=("flatten")

STAGES=("6-Rosetta-scoring" "7-FoldX-scoring" "8-protein-half-life" "9-protein-expression" "flatten")

for STAGE in "${STAGES[@]}"; do
    echo -e "\nRunning stage $STAGE with config $CONFIG_FILE"
    snakemake --snakefile "$STAGE/Snakefile" --configfile "$CONFIG_FILE" --cores 50 --use-conda --conda-frontend conda || exit 1
done
