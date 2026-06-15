#!/usr/bin/env bash
set -euo pipefail

python build-athaliana-ppi-networks.py \
  --in-vivo-input data-files/ppi_in_vivo_edges.csv \
  --in-vitro-input data-files/ppi_in_vitro_edges.csv \
  --output-dir data-files
