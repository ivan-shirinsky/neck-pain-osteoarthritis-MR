# Input data

This repository includes compact, analysis-ready datasets
for primary MR, replication MR, and multivariable MR (MVMR).

Full GWAS summary statistics are not redistributed.

## Primary MR

Prepared exposure and outcome association datasets are
available in `data/processed/primary/`.

The original GWAS preprocessing workflow is implemented in:

- `scripts/01_primary_select_instruments.sh`
- `scripts/02_primary_extract_associations.sh`

Instrument selection uses genome-wide significance
(P < 5e-8) and LD clumping (r² < 0.001, 10 Mb)
with the 1000 Genomes European reference panel.

## Replication MR

Prepared replication datasets are available in
`data/processed/replication/`.

The preprocessing scripts are:

- `scripts/03_replication_select_instruments.sh`
- `scripts/04_replication_extract_associations.sh`

The analyses use alternative neck or shoulder pain GWAS
from Meng et al. (2020) and Tao et al. (2025), with
OA outcome datasets excluding UK Biobank participants.

## Multivariable MR

Prepared harmonised MVMR datasets and results are available
in `data/processed/mvmr/`.

The full MVMR pipeline is implemented in:

`R/08_MVMR_from_GWAS.R`

It includes GWAS preprocessing, instrument selection,
LD clumping, harmonisation, MVMR estimation and diagnostics.

Full preprocessing requires the original GWAS summary
statistics, PLINK and the European LD reference panel.

## Data sources

Primary GWAS:
- Neck pain: GCST90480571
- Knee OA: GCST90566800
- Hand OA: GCST90566797
- Hip OA: GCST90566798

BMI:
- GIANT consortium, Locke et al. (2015)
- IEU OpenGWAS: ieu-a-835

Replication:
- Meng et al. (2020)
- Tao et al. (2025)
- GO Consortium non-UK Biobank OA datasets

## Usage

From the repository root:

    bash scripts/00_link_prepared_inputs.sh

See the main README for instructions on running the analyses.

No individual-level participant data are included.
