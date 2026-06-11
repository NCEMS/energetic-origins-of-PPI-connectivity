#!/usr/bin/env bash
set -euo pipefail

python curate-athaliana-ppi-evidence.py \
  --biogrid "data-files/BIOGRID-ARABIDOPSIS-THALIANA-5.0.258.tab3.xlsx" \
  --table-s1 "data-files/Table S1.xlsx" \
  --in-vivo-output data-files/ppi_in_vivo_edges.csv \
  --in-vitro-output data-files/ppi_in_vitro_edges.csv \
  --unknown-methods-output data-files/ppi_unknown_methods.csv
