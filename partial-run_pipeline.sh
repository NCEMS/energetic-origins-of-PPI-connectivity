#!/bin/bash

CONFIG_FILE=$1

STAGES=("0-download-inputs"
        "1-network-centrality"
        "2-sequence-parsing"
        "3-uniprot-annotation"
        "4-idr-properties"
       )


STAGES=("5-dG-calculations"
        "6-Rosetta-scoring"
        "7-FoldX-scoring"
        "8-protein-half-life"
        "9-protein-expression"
        "flatten"
        )

STAGES=("10-translation-speed" "flatten")

STAGES=("11-predict-PTMs" "flatten")

STAGES=("11-predict-PTMs")

for STAGE in "${STAGES[@]}"; do
    echo -e "\nRunning stage $STAGE with config $CONFIG_FILE"
    snakemake --snakefile "$STAGE/Snakefile" --configfile "$CONFIG_FILE" --cores 4 --use-conda --conda-frontend conda || exit 1
done
