#!/bin/bash

# NOTE WELL - the Yeast Interactome dataset must be downloaded in .cys format and then exported as an edge graph
#             for analysis in Python with NetworkX. The xgmml file format was not working for me. 

# set the download directory
DIR="data-files"

# make sure it exists
mkdir -p "$DIR"

# files to download
URLS=(
	"http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz"
	"http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"
	"https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002311_559292_YEAST_v4.tar"
)

# function to download a file
download_file() {

	local url="$1"
	local output="$DIR/$(basename "$url")"

	# Use curl if available, otherwise use wget
	if command -v curl &>/dev/null; then
		curl -L --retry 5 --retry-delay 10 -o "$output" "$url"
	elif command -v wget &>/dev/null; then
		wget --tries=5 --wait=10 --retry-connrefused -O "$output" "$url"
	else
		echo "Error: Neither curl nor wget are available." >&2
		exit 1
	fi
}

# function to extract compressed files
extract_file() {
    local file="$1"

    if [[ "$file" == *.tar.gz || "$file" == *.tgz ]]; then
        echo "Extracting TAR.GZ: $file"
        tar -xzf "$file" -C "$DIR" && rm "$file"
    elif [[ "$file" == *.tar ]]; then
        echo "Extracting TAR: $file"
        tar -xf "$file" -C "$DIR" && rm "$file"
    elif [[ "$file" == *.gz && "$file" != *.tar.gz ]]; then
        echo "Extracting GZ: $file"
        gunzip "$file"
    else
        echo "No extraction needed: $file"
    fi
}


# loop through URLs and download each file
for url in "${URLS[@]}"; do
	filename="$DIR/$(basename "$url")"
	download_file "$url"
	extract_file "$filename"
done

echo "All files downloaded to $DIR."
