#!/usr/bin/env bash
set -euo pipefail

# Extract outcome associations for six bidirectional MR analyses.
#
# Inputs:
#   GWAS_DIR/{neck,knee,hand,hip}.gz
#   INSTRUMENT_DIR/{neck,knee,hand,hip}_exposure.tsv
#
# Outputs:
#   OUT_DIR/<outcome>_outcome_by_<exposure>.tsv
#
# GWAS rsid column: 16.
# Reproduces the original extraction procedure.

GWAS_DIR="${GWAS_DIR:?Set GWAS_DIR}"
INSTRUMENT_DIR="${INSTRUMENT_DIR:-results/full_pipeline}"
OUT_DIR="${OUT_DIR:-results/full_pipeline}"

mkdir -p "$OUT_DIR/intermediate"

pairs=(
    "neck knee"
    "neck hand"
    "neck hip"
    "knee neck"
    "hand neck"
    "hip neck"
)

for pair in "${pairs[@]}"; do
    read -r exposure outcome <<< "$pair"

    exposure_file="$INSTRUMENT_DIR/${exposure}_exposure.tsv"
    outcome_file="$GWAS_DIR/${outcome}.gz"
    rsids_file="$OUT_DIR/intermediate/${exposure}_rsids.txt"
    output="$OUT_DIR/${outcome}_outcome_by_${exposure}.tsv"

    echo "Processing ${exposure} -> ${outcome}"

    if [[ ! -f "$exposure_file" || ! -f "$outcome_file" ]]; then
        echo "Missing input for ${exposure} -> ${outcome}" >&2
        exit 1
    fi

    awk -F '\t' 'NR > 1 {print $16}' \
        "$exposure_file" > "$rsids_file"

    gzip -cd "$outcome_file" |
        awk -F '\t' '
            NR==FNR {snp[$1]; next}
            FNR==1 || ($16 in snp)
        ' "$rsids_file" - > "$output"

    echo "Saved: $output"
done

echo "Outcome extraction completed."
