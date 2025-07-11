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
        "flatten"
       )


STAGES=("1-network-centrality"
        "2-sequence-parsing"
        "3-uniprot-annotation"
         )

#STAGES=("4-idr-properties")

#STAGES=("flatten")
#STAGES=("3-uniprot-annotation")

#STAGES=("5-dG-calculations" "6-Rosetta-scoring" "7-FoldX-scoring" "8-protein-half-life" "9-protein-expression")

#STAGES=("11-predict-PTMs")

#STAGES=("10-translation-speed")

#STAGES=("flatten")

#STAGES=("12-LiP-MS" "flatten")

#STAGES=("12-LiP-MS")

#STAGES=("12-LiP-MS" "13-entanglement" "flatten")

#STAGES=("14-chaperones" "flatten")

#STAGES=("15-oligomers")

#STAGES=("16-domain-annotations")

#STAGES=("17-essentiality")

#STAGES=("18-Y2H-data")

STAGES=("flatten")

for STAGE in "${STAGES[@]}"; do
    echo -e "\nRunning stage $STAGE with config $CONFIG_FILE"
    snakemake --snakefile "$STAGE/Snakefile" --configfile "$CONFIG_FILE" --cores 4 --use-conda --conda-frontend conda || exit 1
done
