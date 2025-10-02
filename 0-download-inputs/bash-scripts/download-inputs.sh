#!/usr/bin/env bash
set -euo pipefail

DIR="${1:-}"
MANIFEST="${2:-}"

if [[ -z "${DIR}" || -z "${MANIFEST}" ]]; then
  echo "Usage: $0 <TARGET_DIR> <MANIFEST_TSV>" >&2
  exit 1
fi

mkdir -p "${DIR}"

have_cmd() { command -v "$1" >/dev/null 2>&1; }

download_file() {
  local url="$1" out="$2"
  if have_cmd curl; then
    curl -L --retry 5 --retry-delay 10 -o "$out" "$url"
  elif have_cmd wget; then
    wget --tries=5 --wait=10 --retry-connrefused -O "$out" "$url"
  else
    echo "Error: Neither curl nor wget are available." >&2
    exit 1
  fi
}

extract_file() {
  local file="$1" mode="${2:-auto}"

  # If user forces 'none', skip.
  if [[ "${mode}" == "none" ]]; then
    echo "No extraction (forced): ${file}"
    return 0
  fi

  # Infer mode if 'auto'
  if [[ "${mode}" == "auto" ]]; then
    if [[ "$file" == *.tar.gz || "$file" == *.tgz ]]; then
      mode="tar.gz"
    elif [[ "$file" == *.tar ]]; then
      mode="tar"
    elif [[ "$file" == *.gz ]]; then
      mode="gz"
    else
      mode="none"
    fi
  fi

  case "${mode}" in
    tar.gz)
      echo "Extracting TAR.GZ: ${file}"
      tar -xzf "${file}" -C "${DIR}" && rm -f "${file}"
      ;;
    tar)
      echo "Extracting TAR: ${file}"
      tar -xf "${file}" -C "${DIR}" && rm -f "${file}"
      ;;
    gz)
      echo "Extracting GZ: ${file}"
      gunzip -f "${file}"
      ;;
    none)
      echo "No extraction needed: ${file}"
      ;;
    *)
      echo "Unknown extract mode '${mode}' for file '${file}'" >&2
      exit 1
      ;;
  esac
}

# read manifest
while IFS=$'\t' read -r url filename extract || [[ -n "${url}${filename}${extract}" ]]; do
  # skip blank lines & comment
  [[ -z "${url:-}" ]] && continue
  [[ "${url:0:1}" == "#" ]] && continue

  # default extract
  extract="${extract:-auto}"

  out="${DIR}/${filename}"
  echo "Downloading: ${url} -> ${out}"
  download_file "${url}" "${out}"
  extract_file "${out}" "${extract}"
done < "${MANIFEST}"

shopt -s nullglob
for gz in "${DIR}"/*.gz; do
  echo "Post-pass gunzip: ${gz}"
  gunzip -f "${gz}" || true
done
shopt -u nullglob

# clean-up of *.cif from AF2 downloads
find "${DIR}" -maxdepth 1 -type f -name '*.cif' -print -delete

echo "All files processed into ${DIR}"
