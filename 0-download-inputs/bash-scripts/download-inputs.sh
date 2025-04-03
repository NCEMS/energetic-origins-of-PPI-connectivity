#!/bin/bash

# NOTE WELL - the Yeast Interactome dataset must be downloaded in .cys format and then exported as an edge graph
#             for analysis in Python with NetworkX. The xgmml file format was not working for me.

# set the download directory
DIR=$1

# make sure it exists
mkdir -p "$DIR"

# also create directory for processed data (not being used right now!)
#mkdir -p processed-data

# files to download
URLS=(
	"http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz"
	"http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"
	"https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002311_559292_YEAST_v4.tar"
	"https://sid.erda.dk/share_redirect/eIZVVNEd8B"
	"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/YEAST_559292_idmapping.dat.gz"
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz"
)

# define local files names (must match order of URLs)
FILENAMES=(
	"orf_trans.fasta.gz"
	"SGD_features.tab"
	"UP000002311_559292_YEAST_v4.tar"
	"esm_if1_gvp4_t16_142M_UR50.pt"
	"YEAST_559292_idmapping.dat.gz"
        "uniprot_sprot.dat.gz"
)

# function to download a file
download_file() {
	local url="$1"
	local output="$2"

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

# ensure URLS and FILENAMES arrays have the same length
if [[ ${#URLS[@]} -ne ${#FILENAMES[@]} ]]; then
	echo "Error: The number of URLs does not match the number of filenames." >&2
	exit 1
fi

# loop through URLs and download each file with the specified name
for i in "${!URLS[@]}"; do
	url="${URLS[$i]}"
	filename="$DIR/${FILENAMES[$i]}"
	download_file "$url" "$filename"
	extract_file "$filename"
done

# unpack the PDB files from AF2
gunzip data-files/*gz

# remove the .cif files, we will not need them
rm $DIR/*cif

touch "$DIR/.all_files_downloaded"

echo "All files downloaded and saved to $DIR with custom filenames."
