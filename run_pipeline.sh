#!/bin/bash

CONFIG_FILE=$1

STAGES=("0-download-inputs"
        "1-network-centrality"
        "2-sequence-parsing"
        "3-uniprot-annotation"
        "4-idr-properties"
        "5-dG-calculations"
        "6-Rosetta-scoring"
        "7-FoldX-scoring"
        "8-protein-half-life"
        "9-protein-expression"
        "10-translation-speed"
        "11-predict-PTMs"
        "flatten"
       )

for STAGE in "${STAGES[@]}"; do
    echo -e "\nRunning stage $STAGE with config $CONFIG_FILE"

    if [ "$STAGE" == "5-dG-calculations" ]; then
        snakemake --snakefile "$STAGE/Snakefile" --configfile "$CONFIG_FILE" -j 2 --cores 2 --use-conda --conda-frontend conda || exit 1
    else
        snakemake --snakefile "$STAGE/Snakefile" --configfile "$CONFIG_FILE" --cores 4 --use-conda --conda-frontend conda || exit 1
    fi
done
