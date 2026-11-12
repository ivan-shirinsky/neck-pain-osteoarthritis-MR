# primary_MR.R
# Extracted/adapted from user's Org-mode research notebooks.
# primary_MR.org lines 358-439
# Requires derived data in results/: <trait>_exposure.tsv and <outcome>_outcome_by_<exposure>.tsv
# TwoSampleMR::mr() runs all configured methods; compare selected IVW-MRE results with manuscript.


args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else getwd()
if (!file.exists(file.path(project, "README.md"))) stop("Run with project root as argument or working directory")
setwd(project)

primary_input_dir <- Sys.getenv(
  "PRIMARY_INPUT_DIR",
  unset = "data/processed/primary"
)

  library(TwoSampleMR)
  library(dplyr)
  library(readr)
  library(purrr)
  library(stats)

  # ---- Helper: Compute SE from OR and P-value ----
  compute_se_from_or_and_p <- function(or, pval) {
    z <- abs(qnorm(pval / 2, lower.tail = FALSE))
    se <- abs(log(or)) / z
    se[is.nan(se) | is.infinite(se)] <- NA
    return(se)
  }

  # ---- MR Pairs ----
  mr_pairs <- list(
    c("neck", "knee"),
    c("neck", "hand"),
    c("neck", "hip"),
    c("knee", "neck"),
    c("hand", "neck"),
    c("hip", "neck")
  )

  # ---- Batch MR Function ----
  run_mr_pair <- function(exposure, outcome) {
    exposure_file <- paste0(paste0(primary_input_dir, "/"), exposure, "_exposure", ".tsv")
    outcome_file <- paste0(paste0(primary_input_dir, "/"), outcome, "_outcome_by_", exposure, ".tsv")

    message("Running MR: ", exposure, " → ", outcome)

    exp <- read_tsv(exposure_file, show_col_types = FALSE) %>%
      mutate(
        beta = log(odds_ratio),
        se = if_else(is.na(standard_error), compute_se_from_or_and_p(odds_ratio, p_value), standard_error),
        SNP = rsid
      )

    out <- read_tsv(outcome_file, show_col_types = FALSE) %>%
      mutate(
        beta = log(odds_ratio),
        se = if_else(is.na(standard_error), compute_se_from_or_and_p(odds_ratio, p_value), standard_error),
        SNP = rsid
      ) %>%
      filter(SNP %in% exp$SNP)

    if (nrow(out) == 0 || nrow(exp) == 0) {
      warning("No SNPs available for pair ", exposure, " → ", outcome)
      return(NULL)
    }

    exposure_dat <- format_data(exp, type = "exposure", snp_col = "SNP", beta_col = "beta",
                                se_col = "se", effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele", eaf_col = "effect_allele_frequency",
                                pval_col = "p_value")

    outcome_dat <- format_data(out, type = "outcome", snp_col = "SNP", beta_col = "beta",
                               se_col = "se", effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele", eaf_col = "effect_allele_frequency",
                               pval_col = "p_value")

    dat <- harmonise_data(exposure_dat, outcome_dat)
saveRDS(dat, file = paste0("results/harmonised_", exposure, "_", outcome, ".rds"))
res <- mr(dat)
    res$exposure <- exposure
    res$outcome <- outcome
    
# Reorder columns to put exposure first
res <- res %>%
  select(exposure, outcome, everything())
    res
  }

  # ---- Run all pairs ----
  results <- map_dfr(mr_pairs, ~ run_mr_pair(.x[1], .x[2]))

  # ---- Save results ----
  write_csv(results, "results/mr_results_all_pairs.csv")
saveRDS (results, "results/main_mr_results.rds")

