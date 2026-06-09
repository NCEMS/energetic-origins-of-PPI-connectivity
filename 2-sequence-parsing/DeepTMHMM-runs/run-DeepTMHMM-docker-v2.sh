#!/bin/bash

set -u

# record the start time
start_time=$(date +%s)

# ensure correct usage
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <FASTA_FILE> <OUTPUT_DIR> <N_CHUNKS>"
    exit 1
fi

# get command-line arguments
FASTA=$1
OUTDIR=$2
N_CHUNKS=$3

# parallelization settings
N_PARALLEL=8
CPUS_PER_CONTAINER=12
TOTAL_CPUS_REQUESTED=$((N_PARALLEL * CPUS_PER_CONTAINER))

echo "Parallel DeepTMHMM settings:"
echo "  Chunks in parallel:       $N_PARALLEL"
echo "  CPUs per Docker chunk:    $CPUS_PER_CONTAINER"
echo "  Total CPUs requested:     $TOTAL_CPUS_REQUESTED"

# create the output directory and convert it to an absolute path
mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")

# split the FASTA file
echo "Splitting FASTA into $N_CHUNKS chunks..."

if ! command -v seqkit &> /dev/null; then
    echo "Error: seqkit is required but not installed."
    echo "Install with: conda install -c bioconda seqkit"
    exit 1
fi

SPLIT_DIR="$OUTDIR/split_fasta"
mkdir -p "$SPLIT_DIR"

# Use --force so reruns do not fail or warn because split_fasta already exists.
seqkit split --by-part "$N_CHUNKS" "$FASTA" -O "$SPLIT_DIR" --force

# Function to process one FASTA chunk
process_chunk() {
    chunk="$1"

    [ -e "$chunk" ] || {
        echo "Error: chunk file does not exist: $chunk"
        return 1
    }

    CHUNK_NAME=$(basename "$chunk" .fasta)
    CHUNK_OUTDIR="$OUTDIR/$CHUNK_NAME"
    CHUNK_TMPDIR="$OUTDIR/tmp/$CHUNK_NAME"

    mkdir -p "$CHUNK_OUTDIR"
    mkdir -p "$CHUNK_TMPDIR/embeddings" "$CHUNK_TMPDIR/probabilities"

    # Resume behavior: skip chunks already marked complete
    if [ -f "$CHUNK_OUTDIR/done.txt" ]; then
        echo "Skipping chunk: $CHUNK_NAME (already processed)"
        return 0
    fi

    echo "Processing chunk: $CHUNK_NAME"

    # Clean any partial outputs from a previous failed attempt
    rm -rf "$CHUNK_OUTDIR"/*
    rm -rf "$CHUNK_TMPDIR/embeddings"/* "$CHUNK_TMPDIR/probabilities"/*

    # Create expected output files for Docker bind mounts
    touch "$CHUNK_OUTDIR/predicted_topologies.3line" \
          "$CHUNK_OUTDIR/TMRs.gff3" \
          "$CHUNK_OUTDIR/deeptmhmm_results.md" \
          "$CHUNK_OUTDIR/plot.png"

    # Run DeepTMHMM in Docker
    docker run --rm --platform linux/amd64 \
        --cpus "$CPUS_PER_CONTAINER" \
        -e OMP_NUM_THREADS="$CPUS_PER_CONTAINER" \
        -e OPENBLAS_NUM_THREADS="$CPUS_PER_CONTAINER" \
        -e MKL_NUM_THREADS="$CPUS_PER_CONTAINER" \
        -e NUMEXPR_NUM_THREADS="$CPUS_PER_CONTAINER" \
        -v "$(realpath "$chunk")":/openprotein/prot_seqs.fasta \
        -v "$(realpath "$CHUNK_TMPDIR/embeddings")":/openprotein/embeddings \
        -v "$(realpath "$CHUNK_TMPDIR/probabilities")":/openprotein/probabilities \
        -v "$(realpath "$CHUNK_OUTDIR/predicted_topologies.3line")":/openprotein/predicted_topologies.3line \
        -v "$(realpath "$CHUNK_OUTDIR/TMRs.gff3")":/openprotein/TMRs.gff3 \
        -v "$(realpath "$CHUNK_OUTDIR/deeptmhmm_results.md")":/deeptmhmm_results.md \
        -v "$(realpath "$CHUNK_OUTDIR/plot.png")":/openprotein/plot.png \
        dtu/deeptmhmm:1.0.24 \
        python3 predict.py --fasta /openprotein/prot_seqs.fasta

    exit_code=$?

    # Mark complete only if Docker succeeded and the main output exists/non-empty
    if [ "$exit_code" -eq 0 ] && [ -s "$CHUNK_OUTDIR/predicted_topologies.3line" ]; then
        echo "done" > "$CHUNK_OUTDIR/done.txt"
        echo "Chunk $CHUNK_NAME processed successfully."
        rm -rf "$CHUNK_TMPDIR"
        return 0
    else
        echo "Error processing chunk $CHUNK_NAME."
        echo "Docker exit code: $exit_code"
        rm -f "$CHUNK_OUTDIR/done.txt"
        return 1
    fi
}

export -f process_chunk
export OUTDIR
export CPUS_PER_CONTAINER

# Process chunks in parallel
echo "Starting parallel DeepTMHMM processing..."

find "$SPLIT_DIR" -maxdepth 1 -name "*.fasta" | sort | \
    xargs -n 1 -P "$N_PARALLEL" bash -c 'process_chunk "$0"'

xargs_exit_code=$?

# record the end time and calculate elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo
echo "Processing finished."
echo "Results are in: $OUTDIR"
echo "Total elapsed time: ${elapsed_time} seconds"

if [ "$xargs_exit_code" -ne 0 ]; then
    echo "Warning: one or more chunks failed. Rerun the script to retry incomplete chunks."
    exit "$xargs_exit_code"
fi
