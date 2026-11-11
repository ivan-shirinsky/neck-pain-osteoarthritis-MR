# Input data (not included)

Raw full GWAS summary statistics are not redistributed in this repository.

Primary MR requires exposure and outcome tables in `results/`, created from source GWAS with PLINK clumping using 1000 Genomes EUR panel (r² < 0.001, 10 Mb, P < 5e-8). The historical preparation shell fragment is in `scripts/01_exposure_clumping_ORIGINAL_REVIEW.sh` but needs manual refactoring.

Replication expects:
- `results/sensitivity_noUKB/{KneeOA,HandOA,HipOA}_harmonised.csv`
- `results/sensitivity_Meng2020_5e-7/{KneeOA,HandOA,HipOA}_harmonised.csv`
- `results/sensitivity_Tao2025/{KneeOA,HandOA,HipOA}_harmonised.csv`
- `results/Tao_25SNP_GO1/{KneeOA,HandOA,HipOA}_harmonised.tsv`
- `results/Tao_25SNP_GO1/exposure.tsv`

MVMR expects previously generated `results/mvmr_denis_ubuntu_validation` inputs (`mvmr_results.csv`, harmonised `.rds` files) and summary inputs created by the external Denis validation script, which was referenced in the notebooks but NOT attached. The full MVMR pipeline cannot yet be rerun without that script.

Study datasets: primary neck pain MVP GCST90480571; OA GO2 knee GCST90566800, hand GCST90566797, hip GCST90566798; BMI GIANT Locke 2015, IEU OpenGWAS `ieu-a-835`. Replication neck/shoulder pain Meng 2020 (Figshare 7699583), Tao 2025 (Figshare 31362121), GO1 OA outcomes excluding UK Biobank from Musculoskeletal Knowledge Portal (2024-09-16 snapshot). See manuscript for complete references, license and usage restrictions.
