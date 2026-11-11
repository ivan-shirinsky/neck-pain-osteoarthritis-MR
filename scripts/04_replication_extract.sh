#!/usr/bin/env bash
set -euo pipefail

# Reproduce replication exposure and outcome datasets.
# Outcome GWAS: GO1 without UK Biobank, 16 September 2024.

GWAS_DIR="${REPLICATION_GWAS_DIR:?Set REPLICATION_GWAS_DIR}"
SNP_DIR="${SNP_DIR:-results/replication_clumping}"
OUT_DIR="${OUT_DIR:-results/replication_inputs}"

mkdir -p "$OUT_DIR"

specs=(
    "UKB 5e-8"
    "UKB 5e-7"
    "Tao 5e-8"
    "Tao 5e-7"
    "Tao 1e-6"
)

for spec in "${specs[@]}"; do
    read -r study threshold <<< "$spec"

    name="${study}_${threshold}"
    out="$OUT_DIR/$name"
    mkdir -p "$out"

    snps="$SNP_DIR/${name}_snps.txt"

    [[ -s "$snps" ]] || {
        echo "Missing SNP list: $snps" >&2
        exit 1
    }

    echo "=== $name ==="

    if [[ "$study" == "UKB" ]]; then

        input="$GWAS_DIR/neckshoulderpain2_f6159_v3.bgeniemodify.txt.gz"

        gzip -cd "$input" |
        awk -v list="$snps" '
            BEGIN {
                while ((getline s < list) > 0) keep[s]=1
            }
            NR==1 || ($2 in keep) {print}
        ' > "$out/exposure.tsv"

    else

        input="$GWAS_DIR/Tao_2025_GWAS/3vsno3_2.txt.gz"

        gzip -cd "$input" |
        awk -F '\t' -v list="$snps" '
            BEGIN {
                OFS="\t"
                while ((getline s < list) > 0) keep[s]=1
                print "SNP","CHR","POS","A1","A2",
                      "AF1","BETA","SE","P","N","INFO","F"
            }
            $2 in keep {
                F=($8/$9)^2
                print $2,$1,$3,$4,$5,
                      $7,$8,$9,$10,$6,$11,F
            }
        ' > "$out/exposure.tsv"

    fi

    expected=$(wc -l < "$snps")
    actual=$(awk 'END {print NR-1}' "$out/exposure.tsv")

    echo "Exposure SNPs: $actual / $expected"

    [[ "$actual" -eq "$expected" ]] || {
        echo "ERROR: Exposure SNP count mismatch" >&2
        exit 1
    }

    for phenotype in KneeOA HandOA HipOA; do

        input="$GWAS_DIR/GO1SummaryStats_noUKBB/GO1.${phenotype}.withoutUKBB.16Sept2024.txt"
        output="$out/${phenotype}_outcome.tsv"

        [[ -s "$input" ]] || {
            echo "Missing outcome GWAS: $input" >&2
            exit 1
        }

        awk -v study="$study" '
            BEGIN {OFS="\t"}

            FNR==NR {
                if (FNR>1) {
                    if (study=="UKB") {
                        key=$1 ":" $3
                        snp[key]=$2
                    } else {
                        key=$2 ":" $3
                        snp[key]=$1
                    }
                    a0[key]=toupper($4)
                    a1[key]=toupper($5)
                }
                next
            }

            FNR==1 {
                print "SNP","CPTID","EA","NEA",
                      "EAF","BETA","SE","P","N"
                next
            }

            {
                split($1, parts, /[:_]/)
                key=parts[1] ":" parts[2]

                if (key in snp) {
                    ea=toupper($2)
                    nea=toupper($3)

                    if ((ea==a0[key] && nea==a1[key]) ||
                        (ea==a1[key] && nea==a0[key])) {
                        print snp[key],$1,ea,nea,
                              $4,$6,$7,$8,$11
                    }
                }
            }
        ' "$out/exposure.tsv" "$input" > "$output"

        n=$(awk 'END {print NR-1}' "$output")
        echo "$phenotype: $n matched SNPs"
    done
done

echo "Replication data extraction completed."
