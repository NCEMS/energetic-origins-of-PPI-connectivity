#!/usr/bin/env bash
set -euo pipefail

DIR="${1:-}"

if [[ -z "${DIR}" ]]; then
  echo "Usage: $0 <TARGET_DIR>" >&2
  exit 1
fi

if [[ ! -d "${DIR}" ]]; then
  echo "Error: ${DIR} is not a directory." >&2
  exit 1
fi

echo "Scanning ${DIR} for compressed files..."

# .tar.gz and .tgz
shopt -s nullglob
for f in "${DIR}"/*.tar.gz "${DIR}"/*.tgz; do
  echo "Extracting tar.gz archive: ${f}"
  tar -xzf "${f}" -C "${DIR}"
  rm -f "${f}"
done

# .tar
for f in "${DIR}"/*.tar; do
  echo "Extracting tar archive: ${f}"
  tar -xf "${f}" -C "${DIR}"
  rm -f "${f}"
done

# .gz (not .tar.gz)
for f in "${DIR}"/*.gz; do
  [[ "${f}" == *.tar.gz ]] && continue
  echo "Decompressing gzip file: ${f}"
  gunzip -f "${f}"
done
shopt -u nullglob

echo "All compressed files unpacked in ${DIR}"
