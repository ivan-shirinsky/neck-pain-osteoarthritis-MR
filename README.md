# Neck pain and osteoarthritis: Mendelian randomization analyses

This repository contains the R and Bash code and analysis-ready datasets used to investigate the bidirectional relationships between neck pain and osteoarthritis (OA) of the knee, hand, and hip using Mendelian randomization (MR).

The analyses include:

- Primary bidirectional two-sample MR.
- Replication analyses using alternative GWAS datasets of neck or shoulder pain.
- Multivariable MR (MVMR) adjusting for genetically predicted body mass index (BMI).
- Sensitivity analyses and instrument-strength diagnostics.

The repository is intended to support transparency and reproducibility of the statistical analyses.

## Repository structure

```text
neck-pain-osteoarthritis-MR/
├── R/
│   ├── 00_check_inputs.R
│   ├── 02_primary_MR.R
│   ├── 02b_primary_effects_for_figures.R
│   ├── 03_primary_diagnostics.R
│   ├── 04_replication_from_harmonised.R
│   ├── 04a_replication_harmonise.R
│   ├── 05a_MVMR_Qhet.R
│   ├── 05b_MVMR_Qhet_jackknife.R
│   ├── 06_MVMR_tables_from_results.R
│   ├── 07_MR_PRESSO_Lasso.R
│   └── 08_MVMR_from_GWAS.R
├── scripts/
│   ├── 00_link_prepared_inputs.sh
│   ├── 01_primary_select_instruments.sh
│   ├── 02_primary_extract_associations.sh
│   ├── 03_replication_select_instruments.sh
│   └── 04_replication_extract_associations.sh
├── data/
│   ├── README.md
│   └── processed/
│       ├── primary/
│       ├── replication/
│       └── mvmr/
├── reference/
├── results/
├── .gitignore
└── README.md
```

### `scripts/` — GWAS preprocessing

This directory contains Bash scripts for preprocessing large GWAS summary-statistics datasets used in the **primary bidirectional MR and replication analyses**.

Preprocessing for multivariable MR (MVMR) is implemented separately in `R/08_MVMR_from_GWAS.R`.

**Input preparation**

- `00_link_prepared_inputs.sh`: creates symbolic links to the included analysis-ready datasets.

**Primary MR**

- `01_primary_select_instruments.sh`: selects genome-wide significant genetic instruments and performs LD clumping.
- `02_primary_extract_associations.sh`: extracts the corresponding SNP associations from outcome GWAS datasets.

**Replication MR**

- `03_replication_select_instruments.sh`: selects instruments from alternative neck or shoulder pain GWAS datasets (Meng et al. and Tao et al.).
- `04_replication_extract_associations.sh`: extracts the corresponding exposure and outcome associations, including OA associations from GWAS datasets excluding UK Biobank participants.

Full GWAS preprocessing requires the original summary-statistics files, PLINK 1.9, and an appropriate European ancestry LD reference panel.

### `R/` — Statistical analyses

The R scripts implement the statistical analyses and diagnostics.

| Script | Purpose |
|---|---|
| **Primary MR** | |
| `02_primary_MR.R` | Primary bidirectional MR |
| `02b_primary_effects_for_figures.R` | Effect estimates for main figures |
| `03_primary_diagnostics.R` | Primary MR sensitivity diagnostics |
| `07_MR_PRESSO_Lasso.R` | Additional sensitivity analyses |
| **Replication MR** | |
| `04a_replication_harmonise.R` | Harmonisation of replication datasets |
| `04_replication_from_harmonised.R` | Replication and exploratory MR analyses |
| **Multivariable MR** | |
| `08_MVMR_from_GWAS.R` | Main MVMR pipeline: GWAS preprocessing, instrument selection, harmonisation and estimation |
| `05a_MVMR_Qhet.R` | Weak-instrument-robust MVMR heterogeneity analysis |
| `05b_MVMR_Qhet_jackknife.R` | Jackknife confidence intervals |
| `06_MVMR_tables_from_results.R` | MVMR results tables |

The full MVMR script includes GWAS preprocessing, instrument selection, harmonisation, MVMR estimation, and diagnostic analyses.

Unlike the Bash preprocessing scripts for primary and replication MR, this implementation uses R for processing large GWAS files and PLINK for LD clumping.

The full MVMR workflow requires substantially more computational resources than analyses using the included harmonised datasets.

## GWAS data sources

The primary analyses use the following GWAS datasets:

| Phenotype | GWAS accession |
|---|---|
| Neck pain | GCST90480571 |
| Knee OA | GCST90566800 |
| Hand OA | GCST90566797 |
| Hip OA | GCST90566798 |

The MVMR analyses additionally use BMI associations from the GIANT consortium GWAS reported by Locke et al. (2015), available through the OpenGWAS dataset `ieu-a-835`.

Replication analyses use alternative GWAS datasets of neck or shoulder pain reported by Meng et al. (2020) and Tao et al. (2025), together with GO Consortium OA GWAS datasets excluding UK Biobank participants.

Further information is provided in `data/README.md` and `reference/S1_GWAS_study_characteristics.csv`.

## Requirements

The analyses were developed and tested in an R environment on Linux.

The main R dependencies include:

- TwoSampleMR
- ieugwasr
- MVMR
- MendelianRandomization
- data.table
- dplyr
- readr
- tibble
- purrr
- tidyr

Additional packages are required for some sensitivity analyses.

PLINK 1.9 and a European LD reference panel are required for GWAS preprocessing.

Package requirements and computational settings are documented in the individual scripts.

## Reproducing analyses from prepared data

The repository includes compact, analysis-ready datasets, allowing statistical analyses to be reproduced without downloading the complete GWAS summary statistics.

Clone the repository:

```bash
git clone https://github.com/ivan-shirinsky/neck-pain-osteoarthritis-MR.git
cd neck-pain-osteoarthritis-MR
```

Prepare input paths:

```bash
bash scripts/00_link_prepared_inputs.sh
```

Check the available datasets and R packages:

```bash
Rscript R/00_check_inputs.R "$(pwd)"
```

### Primary MR

```bash
Rscript R/02_primary_MR.R "$(pwd)"
Rscript R/02b_primary_effects_for_figures.R "$(pwd)"
Rscript R/03_primary_diagnostics.R "$(pwd)"
```

### Replication MR

```bash
Rscript R/04_replication_from_harmonised.R "$(pwd)"
```

### Multivariable MR

```bash
Rscript R/05a_MVMR_Qhet.R "$(pwd)"
Rscript R/05b_MVMR_Qhet_jackknife.R "$(pwd)"
Rscript R/06_MVMR_tables_from_results.R "$(pwd)"
```

### Additional sensitivity analyses

```bash
Rscript R/07_MR_PRESSO_Lasso.R "$(pwd)"
```

Generated outputs are written to the `results/` directory, which is excluded from version control.

## Full MVMR workflow from GWAS summary statistics

The full MVMR implementation is provided in:

```text
R/08_MVMR_from_GWAS.R
```

This script requires five GWAS summary-statistics datasets:

- Neck pain
- Knee OA
- Hand OA
- Hip OA
- BMI

The script uses environment variables to specify the location of the GWAS files, output directory, PLINK executable, and European LD reference panel.

For example:

```bash
MVMR_DATA_DIR="/path/to/gwas" \
MVMR_RESULTS_DIR="results/mvmr_from_raw" \
EUR_BFILE="/path/to/EUR_reference" \
Rscript R/08_MVMR_from_GWAS.R
```

The supplied GWAS directory must contain the files in the formats expected by the script.

**Verification status:** The full MVMR implementation has been recovered from the original analysis workflow. 

## Effect scales and statistical interpretation

**Primary MR**

Effect estimates are expressed as odds ratios per doubling of genetically predicted odds of the exposure.

The primary analyses use a Bonferroni-corrected significance threshold of 0.05/6.

**Replication MR**

The alternative neck or shoulder pain GWAS datasets provide effect estimates on the linear probability scale.

Replication MR odds ratios are therefore expressed per 10-percentage-point increase in genetically predicted probability of neck or shoulder pain.

**Multivariable MR**

MVMR estimates account for genetically predicted BMI.

Instrument strength is assessed using conditional F-statistics. Additional analyses include heterogeneity and weak-instrument-robust diagnostics.

Interpretation of the MVMR results should take account of conditional instrument strength.

## Data availability

Large GWAS summary-statistics datasets and LD reference panels are not redistributed in this repository.

The included compact datasets contain GWAS-derived summary associations and harmonised analysis inputs.

No individual-level participant data are included.

The original GWAS datasets should be obtained from their respective sources in accordance with applicable access and usage conditions.

## Repository status



The analytical code and prepared datasets are available for reproducibility purposes. 
