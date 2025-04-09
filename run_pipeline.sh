#!/bin/bash

CONFIG_FILE=$1

#STAGES=("0-download-inputs" 
#        "1-network-centrality" 
#        "2-sequence-parsing" 
#        "3-uniprot-annotation" 
#        "4-idr-properties" 
#        "5-dG-calculations")

STAGES=("0-download-inputs"
         "1-network-centrality"
         "2-sequence-parsing"
         "3-uniprot-annotation"
         "4-idr-properties"
         "5-dG-calculations"
       )

#STAGES=("0-download-inputs")

for STAGE in "${STAGES[@]}"; do
    echo -e "\nRunning stage $STAGE with config $CONFIG_FILE"
    snakemake --snakefile "$STAGE/Snakefile" --configfile "$CONFIG_FILE" --cores 4 --use-conda || exit 1
done
