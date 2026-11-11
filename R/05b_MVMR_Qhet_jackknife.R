# 05b_MVMR_Qhet_jackknife.R
# Extracted/adapted from user's Org-mode research notebooks.
# STATUS: REVIEW REQUIRED. Not independently executed in this environment.
# MVMR.org lines 2221-2305
# Requires results/mvmr inputs (.rds); computes leave-one-variant-out jackknife.


args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else getwd()
if (!file.exists(file.path(project, "README.md"))) stop("Run with project root as argument or working directory")
setwd(project)

library(MVMR)

project <- normalizePath(project, mustWork = TRUE)

input_dir <- file.path(
  project,
  "data/processed/mvmr"
)

output_dir <- file.path(
  project,
  "results/mvmr_weak_instrument_robust"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

directions <- c(
  "neck_plus_BMI_to_knee",
  "neck_plus_BMI_to_hand",
  "neck_plus_BMI_to_hip",
  "knee_plus_BMI_to_neck",
  "hip_plus_BMI_to_neck",
  "hand_plus_BMI_to_neck"
)

results <- list()

for (d in directions) {

  cat("\n====================================\n")
  cat(d, "\n")
  cat("====================================\n")

  outfile <- file.path(
    output_dir,
    paste0("qhet_jackknife_", d, ".rds")
  )

  if (file.exists(outfile)) {
    cat("Loading previously calculated result\n")
    res <- readRDS(outfile)

  } else {

    x <- readRDS(file.path(
      input_dir,
      paste0("mvmr_harmonised_", d, ".rds")
    ))

    cat("SNPs:", nrow(x$mvmr_input), "\n")

    res <- tryCatch(
      MVMR::qhet_mvmr(
        r_input = x$mvmr_input,
        pcor = diag(2),
        CI = TRUE,
        CI_method = "jackknife",
        ncores = 2
      ),
      error = function(e) {
        message("ERROR: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(res)) {
      saveRDS(res, outfile)
    }
  }

  print(res)

  results[[d]] <- res
}

saveRDS(
  results,
  file.path(output_dir, "qhet_jackknife_all.rds")
)

cat("\n=== COMPLETED ===\n")
