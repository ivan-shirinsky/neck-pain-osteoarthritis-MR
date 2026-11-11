#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

# Primary MR
for file in data/processed/primary/*.tsv; do
    ln -sfn "../$file" "results/$(basename "$file")"
done

# Replication
for dir in data/processed/replication/*; do
    name=$(basename "$dir")
    ln -sfn "../$dir" "results/$name"
done

# MVMR
mkdir -p results/mvmr
for file in data/processed/mvmr/*; do
    ln -sfn "../../$file" \
        "results/mvmr/$(basename "$file")"
done

echo "Prepared input links created."
