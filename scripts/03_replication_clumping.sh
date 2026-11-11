#!/usr/bin/env bash
set -euo pipefail

# Independent instrument selection for replication GWAS:
# UK Biobank neck/shoulder pain (Meng et al., 2020)
# Tao et al. (2025)
#
# Thresholds: P < 5e-8 and P < 5e-7
# LD: r2 = 0.001, window = 10000 kb

REPLICATION_GWAS_DIR="${REPLICATION_GWAS_DIR:?Set REPLICATION_GWAS_DIR}"
LD_PREFIX="${LD_PREFIX:?Set LD_PREFIX}"
OUT_DIR="${OUT_DIR:-results/replication_clumping}"

mkdir -p "$OUT_DIR"

for study in UKB Tao; do

    if [[ "$study" == "UKB" ]]; then
        input="$REPLICATION_GWAS_DIR/neckshoulderpain2_f6159_v3.bgeniemodify.txt.gz"
    else
        input="$REPLICATION_GWAS_DIR/Tao_2025_GWAS/3vsno3_2.txt.gz"
    fi

    [[ -f "$input" ]] || {
        echo "Missing GWAS: $input" >&2
        exit 1
    }

    for threshold in 5e-8 5e-7 1e-6; do

        [[ "$study" == "UKB" && "$threshold" == "1e-6" ]] && continue

        prefix="$OUT_DIR/${study}_${threshold}"

        echo "Processing $study, P < $threshold"

        gzip -cd "$input" |
        awk -v threshold="$threshold" '
            BEGIN {OFS="\t"; print "SNP","P"}
            NR>1 &&
            $2 ~ /^rs[0-9]+$/ &&
            $10 != "NA" &&
            ($10+0)>0 &&
            ($10+0)<threshold {
                print $2,$10
            }
        ' > "${prefix}_input.tsv"

        plink \
            --bfile "$LD_PREFIX" \
            --clump "${prefix}_input.tsv" \
            --clump-p1 "$threshold" \
            --clump-r2 0.001 \
            --clump-kb 10000 \
            --out "$prefix"

        awk 'NR>1 && $3 ~ /^rs/ {print $3}' \
            "${prefix}.clumped" |
            sort > "${prefix}_snps.txt"

        echo -n "Independent SNPs: "
        wc -l < "${prefix}_snps.txt"

    done
done

echo "Replication instrument selection completed."
