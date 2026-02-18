#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${1:-}"

if [[ -z "$CONFIG_FILE" ]]; then
  echo "Usage: $0 <configfile>"
  exit 2
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "ERROR: config file not found: $CONFIG_FILE"
  exit 2
fi

STAGES=(
  "0-download-inputs"
  "1-network-centrality"
  "2-sequence-parsing"
  "3-uniprot-annotation"
  "4-idr-properties"
  "5-dG-calculations"
  "6-Rosetta-scoring"
  "7-protein-half-life"
  "8-protein-expression"
  "9-translation-efficiency"
  "10-predict-PTMs"
  "11-entanglement"
  "12-chaperones"
  "13-oligomers"
  "14-domain-annotations"
  "15-essentiality"
  "16-Y2H-data"
  "17-meltome-atlas"
  "18-finalize"
)

# choose cores explicitly (optional)
CORES="${CORES:-all}"

for STAGE in "${STAGES[@]}"; do
  echo -e "\nRunning stage $STAGE with config $CONFIG_FILE"

  SNAKEFILE="$STAGE/Snakefile"
  if [[ ! -f "$SNAKEFILE" ]]; then
    echo "ERROR: missing Snakefile: $SNAKEFILE"
    exit 2
  fi

  if [[ "$STAGE" == "2-sequence-parsing" ]]; then
    # do not run SignalP
    snakemake -c "$CORES" --use-conda --conda-frontend conda \
      --snakefile "$SNAKEFILE" --configfile "$CONFIG_FILE" \
      --omit-from run_SignalP

  elif [[ "$STAGE" == "5-dG-calculations" ]]; then
    # do not run Cagiada dG predictions
    snakemake -c "$CORES" --use-conda --conda-frontend conda \
      --snakefile "$SNAKEFILE" --configfile "$CONFIG_FILE" \
      --omit-from cagiada_stability

  elif [[ "$STAGE" == "6-Rosetta-scoring" ]]; then
  # do not run Rosetta
    snakemake -c "$CORES" --use-conda --conda-frontend conda \
      --snakefile "$SNAKEFILE" --configfile "$CONFIG_FILE" \
      add_scores_to_nodes

  elif [[ "$STAGE" == "10-predict-PTMs" ]]; then
  # do not run PTMGPT2
    snakemake -c "$CORES" --use-conda --conda-frontend conda \
      --snakefile "$SNAKEFILE" --configfile "$CONFIG_FILE" \
      process_PTMs

  else
  # standard command for all other steps
    snakemake -c "$CORES" --use-conda --conda-frontend conda \
      --snakefile "$SNAKEFILE" --configfile "$CONFIG_FILE"
  fi
done
