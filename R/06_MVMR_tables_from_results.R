# 06_MVMR_tables_from_results.R
# Extracted/adapted from user's Org-mode research notebooks.
# MVMR.org lines 2364-2853
# Converts numerically calculated IVW and Q-het results to machine-readable CSV.
# DOES NOT recalculate the IVW or Q-het fit; run prerequisite analyses first.
# Inspect formatting and labels before public release.


args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else getwd()
if (!file.exists(file.path(project, "README.md"))) stop("Run with project root as argument or working directory")
setwd(project)



# ============================================================
# 1. Directories
# ============================================================

project <- normalizePath(project, mustWork = TRUE)

ivw_dir <- file.path(
  project,
  "data/processed/mvmr"
)

qhet_dir <- file.path(
  project,
  "results/mvmr_weak_instrument_robust"
)

out_dir <- file.path(
  project,
  "results/tables"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. Directions
# ============================================================

directions <- c(
  "neck_plus_BMI_to_knee",
  "neck_plus_BMI_to_hand",
  "neck_plus_BMI_to_hip",
  "knee_plus_BMI_to_neck",
  "hip_plus_BMI_to_neck",
  "hand_plus_BMI_to_neck"
)

labels <- c(
  "Neck pain -> Knee OA",
  "Neck pain -> Hand OA",
  "Neck pain -> Hip OA",
  "Knee OA -> Neck pain",
  "Hip OA -> Neck pain",
  "Hand OA -> Neck pain"
)

# ============================================================
# 3. Helper functions
# ============================================================

fmt <- function(x, digits = 3) {
  sprintf(paste0("%.", digits, "f"), x)
}

fmt_p <- function(x) {

  ifelse(
    is.na(x),
    "NA",
    ifelse(
      x < 0.001,
      sprintf("%.2e", x),
      sprintf("%.3f", x)
    )
  )
}

# OR per doubling of exposure odds
or_doubling <- function(beta) {
  exp(beta * log(2))
}

ci_doubling <- function(beta, se) {

  z <- qnorm(0.975)

  lower <- exp((beta - z * se) * log(2))
  upper <- exp((beta + z * se) * log(2))

  paste0(
    fmt(lower),
    "-",
    fmt(upper)
  )
}

# ============================================================
# 4. Read existing MVMR-IVW results
# ============================================================

ivw <- read.csv(
  file.path(ivw_dir, "mvmr_results.csv"),
  stringsAsFactors = FALSE
)

stopifnot(
  all(directions %in% ivw$direction)
)

# Primary exposure rows
main <- ivw[ivw$exposure != "BMI", ]

main <- main[
  match(directions, main$direction),
]

# BMI rows
bmi <- ivw[ivw$exposure == "BMI", ]

bmi <- bmi[
  match(directions, bmi$direction),
]

stopifnot(
  nrow(main) == 6,
  nrow(bmi) == 6,
  identical(main$direction, directions),
  identical(bmi$direction, directions)
)

# ============================================================
# 5. Supplementary Table S6
#    Standard multivariable IVW
# ============================================================

s4 <- data.frame(
  Direction = labels,

  SNPs = main$n_snp,

  OR = fmt(or_doubling(main$beta)),

  `95% CI` = mapply(
    ci_doubling,
    main$beta,
    main$se
  ),

  P = fmt_p(main$p_value),

  `Conditional F (primary)` = fmt(
    main$F_conditional,
    2
  ),

  `Conditional F (BMI)` = fmt(
    bmi$F_conditional,
    2
  ),

  `BMI beta` = fmt(bmi$beta),

  `BMI SE` = fmt(bmi$se),

  `BMI P` = fmt_p(bmi$p_value),

  Q = fmt(main$Q_stat, 2),

  `Q P` = fmt_p(main$Q_pval),

  check.names = FALSE
)

# ============================================================
# 6. Supplementary Table S7
#    Weak-instrument robust Q-het MVMR
# ============================================================

qhet_rows <- lapply(seq_along(directions), function(i) {

  d <- directions[i]

  # Existing jackknife results
  f <- file.path(
    qhet_dir,
    paste0("qhet_jackknife_", d, ".rds")
  )

  stopifnot(file.exists(f))

  x <- readRDS(f)

  stopifnot(
    nrow(x) == 2,
    all(c(
      "Effect Estimates",
      "Std. Error"
    ) %in% names(x))
  )

  # Exposure 1: primary phenotype
  beta <- as.numeric(
    x[1, "Effect Estimates"]
  )

  se <- as.numeric(
    x[1, "Std. Error"]
  )

  # Exposure 2: BMI
  bmi_beta <- as.numeric(
    x[2, "Effect Estimates"]
  )

  bmi_se <- as.numeric(
    x[2, "Std. Error"]
  )

  # Approximate two-sided Wald P from jackknife SE
  p <- 2 * pnorm(
    -abs(beta / se)
  )

  # Previously calculated Q-min test
  qfile <- file.path(
    qhet_dir,
    paste0("qhet_", d, ".rds")
  )

  stopifnot(file.exists(qfile))

  original <- readRDS(qfile)

  qmin <- original$qmin

  stopifnot(
    !is.null(qmin$Qstat),
    !is.null(qmin$Qpval)
  )

  data.frame(
    Direction = labels[i],

    SNPs = main$n_snp[i],

    OR = fmt(or_doubling(beta)),

    `95% CI` = ci_doubling(
      beta, se
    ),

    `P (approx.)` = fmt_p(p),

    `Conditional F (primary)` = fmt(
      main$F_conditional[i], 2
    ),

    `Conditional F (BMI)` = fmt(
      bmi$F_conditional[i], 2
    ),

    `BMI beta` = fmt(bmi_beta),

    `BMI SE` = fmt(bmi_se),

    `Q-min` = fmt(qmin$Qstat, 2),

    `Q-min P` = fmt_p(qmin$Qpval),

    check.names = FALSE
  )
})

s5 <- do.call(
  rbind,
  qhet_rows
)

rownames(s5) <- NULL

# ============================================================
# 7. Save CSV files
# ============================================================

write.csv(
  s4,
  file.path(
    out_dir,
    "Supplementary_Table_S6_MVMR_IVW.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  s5,
  file.path(
    out_dir,
    "Supplementary_Table_S7_MVMR_Qhet.csv"
  ),
  row.names = FALSE,
  na = ""
)


cat("\n=== S6: MVMR IVW ===\n"); print(s4, row.names=FALSE)
cat("\n=== S7: Q-het MVMR ===\n"); print(s5, row.names=FALSE)
