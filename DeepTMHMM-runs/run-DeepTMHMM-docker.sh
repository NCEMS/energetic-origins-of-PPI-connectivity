#!/bin/bash

# Record the start time
start_time=$(date +%s)

# Ensure correct usage
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <FASTA_FILE> <OUTPUT_DIR> <N_CHUNKS>"
    exit 1
fi

FASTA=$1
OUTDIR=$2
N_CHUNKS=$3

# Create the output directory and convert it to an absolute path
mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")

# Create a shared temporary directory for embeddings and probabilities
TMPDIR="$OUTDIR/tmp"
mkdir -p "$TMPDIR/embeddings" "$TMPDIR/probabilities"

# Step 1: Split the FASTA file
echo "Splitting FASTA into $N_CHUNKS chunks..."
if ! command -v seqkit &> /dev/null; then
    echo "Error: seqkit is required but not installed. Please install it using 'conda install -c bioconda seqkit' or 'brew install seqkit'."
    exit 1
fi

SPLIT_DIR="$OUTDIR/split_fasta"
mkdir -p "$SPLIT_DIR"
seqkit split --by-part "$N_CHUNKS" "$FASTA" -O "$SPLIT_DIR"

# Step 2: Process each chunk in series with resume capability
for chunk in "$SPLIT_DIR"/*.fasta; do
    # Ensure we have a valid chunk file
    [ -e "$chunk" ] || { echo "Error: No split FASTA files found!"; exit 1; }
    
    CHUNK_NAME=$(basename "$chunk" .fasta)
    CHUNK_OUTDIR="$OUTDIR/$CHUNK_NAME"
    mkdir -p "$CHUNK_OUTDIR"

    # Check if this chunk has already been processed
    if [ -f "$CHUNK_OUTDIR/done.txt" ]; then
        echo "Skipping chunk: $CHUNK_NAME (already processed)"
        continue
    fi

    echo "Processing chunk: $CHUNK_NAME..."

    # **Cleanup:** Clear the temporary directories before processing this chunk
    rm -rf "$TMPDIR/embeddings"/* "$TMPDIR/probabilities"/*

    # Remove any existing files in the chunk's directory (to ensure a clean state)
    rm -rf "$CHUNK_OUTDIR"/*

    # Create expected output files in the chunk's directory
    touch "$CHUNK_OUTDIR/predicted_topologies.3line" \
          "$CHUNK_OUTDIR/TMRs.gff3" \
          "$CHUNK_OUTDIR/deeptmhmm_results.md" \
          "$CHUNK_OUTDIR/plot.png"

    # Run the Docker container with appropriate mounts
    docker run --platform linux/amd64 \
        -v "$(realpath "$chunk")":/openprotein/prot_seqs.fasta \
        -v "$TMPDIR/embeddings":/openprotein/embeddings \
        -v "$TMPDIR/probabilities":/openprotein/probabilities \
        -v "$CHUNK_OUTDIR/predicted_topologies.3line":/openprotein/predicted_topologies.3line \
        -v "$CHUNK_OUTDIR/TMRs.gff3":/openprotein/TMRs.gff3 \
        -v "$CHUNK_OUTDIR/deeptmhmm_results.md":/deeptmhmm_results.md \
        -v "$CHUNK_OUTDIR/plot.png":/openprotein/plot.png \
        dtu/deeptmhmm:1.0.24 \
        python3 predict.py --fasta /openprotein/prot_seqs.fasta

    # Check if the Docker command succeeded; if yes, mark the chunk as done.
    if [ "$?" -eq 0 ]; then
        echo "done" > "$CHUNK_OUTDIR/done.txt"
        echo "Chunk $CHUNK_NAME processed successfully."
    else
        echo "Error processing chunk $CHUNK_NAME. Please check logs."
    fi

    # Clear the temporary directories between chunks
    rm -rf "$TMPDIR/embeddings"/* "$TMPDIR/probabilities"/*
done

# Record the end time and calculate elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Processing complete. Results are in $OUTDIR."
echo "Total elapsed time: ${elapsed_time} seconds"
