#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="/home/jovyan/work"
EXAMPLES_DIR="/home/jovyan/examples"

# If /work is empty, seed it with examples (so users edit copies)
if [ -z "$(ls -A "$WORK_DIR" 2>/dev/null || true)" ]; then
  mkdir -p "$WORK_DIR"
  cp -r "$EXAMPLES_DIR/"* "$WORK_DIR/" 2>/dev/null || true
fi

# Hand off to the stock Jupyter starter
exec start-notebook.sh "$@"
