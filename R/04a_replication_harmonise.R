#!/usr/bin/env Rscript

# Reproduce harmonisation of Meng and Tao replication datasets.
# Method: TwoSampleMR::harmonise_data(action = 2).
# Does not overwrite original processed datasets.

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(readr)
  library(dplyr)
})

INPUT <- Sys.getenv(
  "REPLICATION_INPUT_DIR",
  unset = "results/replication_inputs"
)

OUTPUT <- Sys.getenv(
  "REPLICATION_HARMONISED_DIR",
  unset = "results/replication_harmonised"
)

specs <- data.frame(
  study = c("UKB", "UKB", "Tao", "Tao"),
  threshold = c("5e-8", "5e-7", "5e-8", "1e-6"),
  directory = c(
    "sensitivity_noUKB",
    "sensitivity_Meng2020_5e-7",
    "sensitivity_Tao2025",
    "Tao_25SNP_GO1"
  )
)

for (i in seq_len(nrow(specs))) {

  study <- specs$study[i]
  threshold <- specs$threshold[i]

  source_dir <- file.path(INPUT, paste0(study, "_", threshold))
  target_dir <- file.path(OUTPUT, specs$directory[i])

  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n===", study, threshold, "===\n")

  exposure_file <- file.path(source_dir, "exposure.tsv")

  if (study == "UKB") {

    raw <- read.table(
      exposure_file,
      header = TRUE,
      sep = "",
      check.names = FALSE
    )

    exp <- format_data(
      raw,
      type = "exposure",
      snp_col = "rsid",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "a_1",
      other_allele_col = "a_0",
      eaf_col = "af",
      pval_col = "p",
      samplesize_col = "N"
    )

  } else {

    raw <- read_tsv(exposure_file, show_col_types = FALSE)

    exp <- format_data(
      raw,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      eaf_col = "AF1",
      pval_col = "P",
      samplesize_col = "N",
      chr_col = "CHR",
      pos_col = "POS"
    )

  }

  if (study == "Tao" && threshold == "1e-6") {
    ok <- file.copy(
      exposure_file,
      file.path(target_dir, "exposure.tsv"),
      overwrite = TRUE
    )
    if (!ok) stop("Could not copy Tao exposure file")
  }

  for (phenotype in c("KneeOA", "HandOA", "HipOA")) {

    outcome_file <- file.path(
      source_dir,
      paste0(phenotype, "_outcome.tsv")
    )

    out <- read_outcome_data(
      filename = outcome_file,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "EA",
      other_allele_col = "NEA",
      eaf_col = "EAF",
      pval_col = "P",
      samplesize_col = "N"
    )

    dat <- harmonise_data(exp, out, action = 2)

    # Original Meng and Tao strict datasets retained mr_keep SNPs.
    # Tao exploratory retains the full 25-SNP harmonised dataset.

    if (study == "UKB" || threshold == "5e-8") {
      dat <- dat[dat$mr_keep %in% TRUE, ]
    }

    extension <- if (
      study == "UKB" || threshold == "5e-8"
    ) "csv" else "tsv"

    output_file <- file.path(
      target_dir,
      paste0(phenotype, "_harmonised.", extension)
    )

    if (extension == "csv") {
      write_csv(dat, output_file)
    } else {
      write_tsv(dat, output_file)
    }

    cat(
      phenotype,
      ": harmonised =", nrow(dat),
      "; mr_keep =", sum(dat$mr_keep %in% TRUE),
      "\n"
    )

    if (study == "Tao" && threshold == "1e-6") {
      selected <- dat |>
        filter(mr_keep %in% TRUE, pval.exposure < 5e-7)

      cat(
        "Tao exploratory P < 5e-7:",
        nrow(selected), "\n"
      )
    }

  }
}

cat("\nReplication harmonisation completed.\n")
