#!/usr/bin/env bash

set -euo pipefail

for file in union-*; do
    # Skip if no union-* files exist
    [[ -e "$file" ]] || continue

    base="${file#union-}"

    cp "$file" "in-vivo-${base}"
    cp "$file" "in-vitro-${base}"
    cp "$file" "intersection-in-vivo-${base}"
    cp "$file" "intersection-in-vitro-${base}"
done
