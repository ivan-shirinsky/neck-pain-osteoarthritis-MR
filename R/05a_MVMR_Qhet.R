# 05a_MVMR_Qhet.R
# Extracted/adapted from user's Org-mode research notebooks.
# STATUS: REVIEW REQUIRED. Not independently executed in this environment.
# MVMR.org lines 1906-2033
# Requires results/mvmr_denis_ubuntu_validation inputs (.rds); calculates Q-min estimates.


args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else getwd()
if (!file.exists(file.path(project, "README.md"))) stop("Run with project root as argument or working directory")
setwd(project)

library(MVMR)

# Directories
project <- normalizePath(project, mustWork = TRUE)

input_dir <- file.path(
  project,
  "data/processed/mvmr"
)

output_dir <- file.path(
  project,
  "results/mvmr_weak_instrument_robust"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Six directions
directions <- c(
  "neck_plus_BMI_to_knee",
  "neck_plus_BMI_to_hand",
  "neck_plus_BMI_to_hip",
  "knee_plus_BMI_to_neck",
  "hip_plus_BMI_to_neck",
  "hand_plus_BMI_to_neck"
)

results <- list()

cat("MVMR version:",
    as.character(packageVersion("MVMR")), "\n")

for (d in directions) {

  cat("\n====================================\n")
  cat(d, "\n")
  cat("====================================\n")

  # Load previously harmonised data
  x <- readRDS(
    file.path(
      input_dir,
      paste0("mvmr_harmonised_", d, ".rds")
    )
  )

  dat <- x$mvmr_input

  cat("Number of SNPs:", nrow(dat), "\n")

  # Existing IVW results
  cat("\n=== Standard MVMR-IVW ===\n")
  print(x$ivw)

  cat("\n=== Conditional F ===\n")
  print(x$cond_f)

  # Weak-instrument robust estimation
  # Assumes zero covariance between exposure GWAS
  cat("\n=== Q-het MVMR ===\n")

  qhet <- tryCatch(
    MVMR::qhet_mvmr(
      r_input = dat,
      pcor = diag(2),
      CI = FALSE,
      ncores = 1
    ),
    error = function(e) {
      message("Q-het error: ", conditionMessage(e))
      NULL
    }
  )

  print(qhet)

  # Weak-instrument robust heterogeneity test
  cat("\n=== Q-min heterogeneity test ===\n")

  qmin <- tryCatch(
    MVMR::pleiotropy_mvmr(
      r_input = dat,
      gencov = 0,
      estimator = "qmin"
    ),
    error = function(e) {
      message("Q-min error: ", conditionMessage(e))
      NULL
    }
  )

  print(qmin)

  results[[d]] <- list(
    n_snp = nrow(dat),
    ivw = x$ivw,
    conditional_F = x$cond_f,
    qhet = qhet,
    qmin = qmin
  )

  saveRDS(
    results[[d]],
    file.path(
      output_dir,
      paste0("qhet_", d, ".rds")
    )
  )
}

# Save complete results
saveRDS(
  results,
  file.path(output_dir, "qhet_all_results.rds")
)

cat("\n====================================\n")
cat("Q-HET ANALYSIS COMPLETED\n")
cat("Results saved to:\n")
cat(output_dir, "\n")
