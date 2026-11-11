#!/usr/bin/env bash
set -euo pipefail

# Reproduce instrument selection for the primary MR analysis.
# GWAS Catalog: GCST90480571, GCST90566800,
#               GCST90566797, GCST90566798.
#
# Requirements: PLINK 1.9, EUR LD reference panel.
# Parameters: P < 5e-8, r2 = 0.001, window = 10000 kb.
#
# Input files: neck.gz, knee.gz, hand.gz, hip.gz
# Expected columns: p_value (8), rsid (16).

GWAS_DIR="${GWAS_DIR:?Set GWAS_DIR}"
LD_PREFIX="${LD_PREFIX:?Set LD_PREFIX}"
OUT_DIR="${OUT_DIR:-results/full_pipeline}"
PLINK_BIN="${PLINK_BIN:-plink}"

mkdir -p "$OUT_DIR/intermediate"

for trait in neck knee hand hip; do

    echo "======================================"
    echo "Processing: $trait"
    echo "======================================"

    input="$GWAS_DIR/$trait.gz"
    filtered="$OUT_DIR/intermediate/${trait}_filtered.tsv"
    formatted="$OUT_DIR/intermediate/${trait}_formatted.tsv"
    clumped="$OUT_DIR/intermediate/${trait}_clumped"
    snps="$OUT_DIR/intermediate/${trait}_snps.txt"
    output="$OUT_DIR/${trait}_exposure.tsv"

    if [[ ! -f "$input" ]]; then
        echo "Missing GWAS: $input" >&2
        exit 1
    fi

    gzip -cd "$input" |
        awk -F '\t' '
            NR == 1 || ($8+0 < 5e-8 && $16 ~ /^rs/)
        ' OFS='\t' > "$filtered"

    awk -F '\t' '
        BEGIN {OFS="\t"}
        NR == 1 {print "SNP", "P"}
        NR > 1 {print $16, $8}
    ' "$filtered" > "$formatted"

    "$PLINK_BIN" \
        --bfile "$LD_PREFIX" \
        --clump "$formatted" \
        --clump-p1 5e-8 \
        --clump-r2 0.001 \
        --clump-kb 10000 \
        --out "$clumped"

    awk 'NR>1 && $3 ~ /^rs/ {print $3}' \
        "${clumped}.clumped" > "$snps"

    head -n 1 "$filtered" > "$output"

    awk -F '\t' '
        NR==FNR {keep[$1]; next}
        ($16 in keep)
    ' "$snps" "$filtered" >> "$output"

    echo "Selected instruments:"
    awk 'NR>1 {n++} END {print n+0}' "$output"
    echo "Saved: $output"

done

echo "Instrument selection completed."
