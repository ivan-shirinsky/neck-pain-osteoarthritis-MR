# Neck pain and osteoarthritis: Mendelian randomization analyses

Code and analysis-ready data for a study of the bidirectional
relationships between neck pain and osteoarthritis of the knee,
hand and hip.

## Contents

- `R/` — R scripts for MR, sensitivity analyses and MVMR.
- `scripts/` — scripts for preparing inputs from GWAS summary statistics.
- `data/processed/primary/` — analysis-ready primary MR data.
- `data/processed/replication/` — replication datasets.
- `data/processed/mvmr/` — harmonised MVMR data and reference results.
- `reference/` — GWAS study characteristics.
- `provenance/` — mapping of analytical code to original notebooks.

Large GWAS summary statistics and LD reference panels are not
distributed in this repository.

## Reproducibility status

### Primary bidirectional MR

The primary analysis was reproduced from the original GWAS
summary statistics, including instrument selection, LD clumping,
outcome extraction, harmonisation and MR estimation.

### Replication analyses

Replication analyses using the GWAS reported by Meng et al.
and Tao et al., with GO1 non-UK Biobank OA outcomes, were
reproduced from source GWAS data.

All 12 IVW estimates (beta, SE and P-value) matched the
previously obtained results.

### Multivariable MR

The MVMR analyses were reproduced using the included
harmonised input datasets.

These analyses include BMI adjustment, conditional instrument
strength, heterogeneity statistics and jackknife diagnostics.

The repository does not currently provide a validated
end-to-end MVMR preprocessing workflow starting from raw GWAS.

## Requirements

R and the packages required by the individual analysis scripts.

The primary and replication workflows also use the
TwoSampleMR package.

Full GWAS preprocessing requires PLINK 1.9 and a suitable
European LD reference panel.

## Running analyses from prepared data

Clone or download this repository and open a terminal in
its root directory.

First create links to the included analysis-ready datasets:

    bash scripts/00_link_prepared_inputs.sh

Then check the available inputs:

    Rscript R/00_check_inputs.R "$(pwd)"

Primary MR:

    Rscript R/02_primary_MR.R "$(pwd)"
    Rscript R/02b_primary_effects_for_figures.R "$(pwd)"
    Rscript R/03_primary_diagnostics.R "$(pwd)"

Replication:

    Rscript R/04_replication_from_harmonised.R "$(pwd)"

MVMR:

    Rscript R/05a_MVMR_Qhet.R "$(pwd)"
    Rscript R/05b_MVMR_Qhet_jackknife.R "$(pwd)"
    Rscript R/06_MVMR_tables_from_results.R "$(pwd)"

Additional sensitivity analyses:

    Rscript R/07_MR_PRESSO_Lasso.R "$(pwd)"

Individual scripts may require additional R packages.
See their headers and console output.

## Effect scales

Primary MR: odds ratios are expressed per doubling of
genetically predicted odds of the exposure.

Replication: odds ratios are expressed per 10-percentage-point
increase in genetically predicted probability of neck or
shoulder pain.

The primary MR uses a Bonferroni significance threshold
of 0.05 / 6.

## Data sources

The original GWAS summary statistics are not redistributed.

See `data/README.md` and
`reference/S1_GWAS_study_characteristics.csv` for
information on the source datasets.

## Reproducibility scope

The included analysis-ready datasets permit reproduction of
statistical analyses without downloading the full GWAS files.

The original GWAS preprocessing steps were also checked
locally for the primary MR and replication analyses.

The MVMR preprocessing stage is documented separately and
is not claimed to be independently reproduced from raw GWAS
by this public workflow.
