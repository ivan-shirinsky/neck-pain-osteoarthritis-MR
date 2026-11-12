#!/usr/bin/env Rscript
# =============================================================================
#  Multivariable Mendelian randomization of neck pain, osteoarthritis of the
#  knee / hip / hand, and body mass index
#
#  Repository : https://github.com/ivan-shirinsky/neck-pain-osteoarthritis-MR
#  Zenodo DOI : not assigned
#  Maintainer : Ivan Shirinsky
#  License    : MIT
#  Version    : 1.0.0
#  Date       : 2026-11-01
#
#  SUMMARY
#  -------
#  Two-sample multivariable Mendelian randomization (MVMR) analysis of six
#  directional models (three exposures × two outcomes), each adjusted for
#  genetically predicted body mass index (BMI):
#
#    1. neck pain + BMI  ->  knee OA
#    2. neck pain + BMI  ->  hand OA
#    3. neck pain + BMI  ->  hip OA
#    4. knee OA   + BMI  ->  neck pain
#    5. hip OA    + BMI  ->  neck pain
#    6. hand OA   + BMI  ->  neck pain
#
#  Primary estimator : MVMR-IVW  (MVMR package, v0.4.9)
#  Sensitivity       : MVMR-Egger, MVMR-Median (MendelianRandomization pkg)
#  Diagnostics       : conditional F-statistic, Cochran's Q (pleiotropy_mvmr)
#
#  DATA SOURCES (full summary statistics; NOT pre-filtered to instruments)
#  ----------------------------------------------------------------------
#    neck pain   : MVP                 GCST90480571
#    knee OA     : GO Consortium 2025  GCST90566800
#    hip OA      : GO Consortium 2025  GCST90566798
#    hand OA     : GO Consortium 2025  GCST90566797
#    BMI         : GIANT Locke 2015    ieu-a-835 (n = 322,154)
#
#  REQUIRED ENVIRONMENT
#  --------------------
#    R            >= 4.2      
#    data.table   >= 1.14
#    readr        >= 2.0
#    dplyr        >= 1.0
#    tibble       >= 3.0
#    TwoSampleMR  >= 0.5
#    MVMR         >= 0.4.9    (for ivw_mvmr / strength_mvmr / pleiotropy_mvmr)
#    ieugwasr     >= 1.0
#    MendelianRandomization >= 0.6  (for mr_mvegger / mr_mvmedian)
#
#  OPTIONAL (for local LD clumping instead of OpenGWAS API)
#  ------------------------------------------------------
#    PLINK 1.9                 (plink_win64_*.zip / plink_linux_*)
#    1000 Genomes EUR panel    (http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz)
#    Set USE_LOCAL_PLINK = TRUE and point PLINK_BIN / EUR_BFILE below.
#
#  OPTIONAL (for API clumping — default)
#  -------------------------------------
#    A valid OpenGWAS JWT saved in ~/.Renviron as
#       OPENGWAS_JWT=eyJhbGci...
#    Register at https://openGWAS.org/profile to generate a token.
#
#  INPUT LAYOUT
#  ------------
#    <DATA_DIR>/neck.gz    or any of the following column-name conventions:
#      rsid/SNP, A1/effect_allele, A2/other_allele, beta/b/OR,
#      se/standard_error, p/pval/p_value, EAF/effect_allele_frequency
#    The script auto-detects column names via detect_columns().
#
#  PIPELINE
#  --------
#    PHASE 0  Pre-slim each full GWAS to 7 canonical columns (gzip).
#             One-time operation; results cached in <CACHE_DIR>/slim/.
#    PHASE 1  For each trait: stream slim file, keep p < 5e-8,
#             LD-clump (r2 < 0.001, 10,000 kb), cache instruments.
#    PHASE 2  For each of the 6 directions: pool instruments, extract
#             SNP subsets from slim, harmonise, run MVMR-IVW + sensitivity.
#
#  OUTPUTS (under <RESULTS_DIR>)
#  -----------------------------
#    mvmr_results.csv                 summary table, 12 rows (6 × 2 exposures)
#    mvmr_results.rds                 same, R-serialised
#    mvmr_harmonised_<direction>.rds  harmonised data + MVMR inputs + results
#    mvmr_sessionInfo.txt             reproducibility metadata
#    cache/slim/<trait>.tsv.gz        slim GWAS (intermediate, regenerable)
#    cache/instr_<trait>.rds          cached instrument sets
#
#  HOW TO CITE
#  -----------
#    Shirinsky et al. (2026). Multivariable Mendelian randomization of neck
#    pain, osteoarthritis and body mass index. Zenodo.
#    https://doi.org/not assigned
#
#  =============================================================================


# ---------------------------------------------------------------------------
# 0.  USER CONFIGURATION
# ---------------------------------------------------------------------------
# Edit these paths to match your local environment.
# Relative paths are resolved against the working directory.

USE_LOCAL_PLINK <- TRUE      # local PLINK with European LD reference

# Required external inputs are supplied through environment variables.
# Results are saved inside the current repository by default.

DATA_DIR <- path.expand(Sys.getenv("MVMR_DATA_DIR", ""))
if (!nzchar(DATA_DIR))
  stop("Set MVMR_DATA_DIR to the directory containing five GWAS files")

RESULTS_DIR <- path.expand(Sys.getenv(
  "MVMR_RESULTS_DIR", "results/mvmr_from_raw"
))

PLINK_BIN <- Sys.getenv("PLINK_BIN", Sys.which("plink"))
EUR_BFILE <- path.expand(Sys.getenv("EUR_BFILE", ""))
if (!nzchar(EUR_BFILE))
  stop("Set EUR_BFILE to the PLINK European reference prefix")

# Analysis parameters (defaults correspond to the published analysis)
P_THRESH  <- 5e-8
CLUMP_R2  <- 0.001
CLUMP_KB  <- 10000
N_EXPOSURES <- 2
GCOV_APPROX <- "zero"        # MVMR 0.4: gencov is a scalar; 0 = independence
CHUNK_SIZE  <- 500000L


# ---------------------------------------------------------------------------
# 1.  PACKAGES
# ---------------------------------------------------------------------------

required_pkgs <- c(
  "data.table", "readr", "dplyr", "tibble",
  "TwoSampleMR", "MVMR", "ieugwasr", "MendelianRandomization"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tibble)
  library(TwoSampleMR)
  library(MVMR)
  library(ieugwasr)
})


# ---------------------------------------------------------------------------
# 2.  PATHS AND ENVIRONMENT CHECK
# ---------------------------------------------------------------------------

CACHE_DIR <- file.path(RESULTS_DIR, "cache")
SLIM_DIR  <- file.path(CACHE_DIR, "slim")
for (d in c(RESULTS_DIR, CACHE_DIR, SLIM_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

gwas_paths <- list(
  neck = file.path(DATA_DIR, "neck.gz"),
  knee = file.path(DATA_DIR, "knee.gz"),
  hip  = file.path(DATA_DIR, "hip.gz"),
  hand = file.path(DATA_DIR, "hand.gz"),
  bmi  = file.path(DATA_DIR, "bmi.gz")
)

gwas_origins <- list(
  neck = "MVP GCST90480571",
  knee = "GO Consortium 2025 GCST90566800",
  hip  = "GO Consortium 2025 GCST90566798",
  hand = "GO Consortium 2025 GCST90566797",
  bmi  = "GIANT Locke 2015 (ieu-a-835, n = 322,154)"
)

message("== Environment check ==")

miss_files <- gwas_paths[!file.exists(unlist(gwas_paths))]
if (length(miss_files) > 0)
  stop("Missing GWAS files:\n  ", paste(unlist(miss_files), collapse = "\n  "),
       "\nExpected under: ", DATA_DIR)
message("  GWAS files   : OK (", length(gwas_paths), ")")

if (isTRUE(USE_LOCAL_PLINK)) {
  if (!file.exists(PLINK_BIN))
    stop("PLINK not found: ", PLINK_BIN)
  candidates <- unique(c(
    EUR_BFILE,
    file.path(dirname(PLINK_BIN), "1kg.v3", "EUR"),
    file.path(dirname(PLINK_BIN), "EUR"),
    "C:/tools/1kg.v3/EUR", "C:/data/1kg.v3/EUR"
  ))
  candidates <- candidates[nzchar(candidates)]
  found <- NULL
  for (bf in candidates) {
    if (all(file.exists(paste0(bf, c(".bed", ".bim", ".fam"))))) {
      found <- bf; break
    }
  }
  if (is.null(found))
    stop("EUR reference not found. Tried:\n  ",
         paste(candidates, collapse = "\n  "))
  EUR_BFILE <- found
  message("  PLINK        : ", PLINK_BIN)
  message("  EUR bfile    : ", EUR_BFILE)
} else {
  if (Sys.getenv("OPENGWAS_JWT") == "")
    stop("OPENGWAS_JWT is not set. Either set it in ~/.Renviron ",
         "or switch USE_LOCAL_PLINK <- TRUE.")
  message("  Clumping mode: OpenGWAS API")
}

message("  Data dir     : ", DATA_DIR)
message("  Results dir  : ", RESULTS_DIR)
message("  Cache dir    : ", CACHE_DIR)
message("  MVMR         : ", as.character(packageVersion("MVMR")))
message("  TwoSampleMR  : ", as.character(packageVersion("TwoSampleMR")))
message("  MR pkg       : ",
        as.character(packageVersion("MendelianRandomization")))
message("== check passed ==\n")


# ---------------------------------------------------------------------------
# 3.  GENERIC HELPERS
# ---------------------------------------------------------------------------

# Recover a standard error from an OR and a two-sided p-value.
compute_se_from_or_and_p <- function(or, pval) {
  z  <- abs(qnorm(pval / 2, lower.tail = FALSE))
  se <- abs(log(or)) / z
  se[!is.finite(se)] <- NA_real_
  se
}

# Thousands-separated number for logging.
fmt_num <- function(x) {
  tryCatch(
    base::format(as.numeric(x), big.mark = ",", scientific = FALSE,
                 trim = TRUE),
    error = function(e) as.character(x)
  )
}

# Parse p-value strings that may contain "<" or ">".
clean_p_vec <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(gsub("^[<>]=?", "", trimws(as.character(x)))))
}

# Read the header of a (possibly gzipped) GWAS and map column names to
# canonical ones. Accepts most common conventions.
detect_columns <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  header <- readLines(con, n = 1L, warn = FALSE)
  nms <- strsplit(header, "\t", fixed = TRUE)[[1]]
  nm_low <- tolower(nms)
  pick <- function(cands) {
    for (c in cands) {
      idx <- which(nm_low == tolower(c))
      if (length(idx)) return(nms[idx[1]])
    }
    NULL
  }
  list(
    snp  = pick(c("snp","rsid","rs_id","markername","variant_id")),
    ea   = pick(c("effect_allele","a1","ea","alt")),
    oa   = pick(c("other_allele","a2","oa","ref")),
    beta = pick(c("beta","b","effect","es")),
    or   = pick(c("odds_ratio","or","oddsratio")),
    se   = pick(c("standard_error","se","stderr","sebeta")),
    p    = pick(c("p_value","p","pval","pvalue","p_val")),
    eaf  = pick(c("effect_allele_frequency","eaf","freq1.hapmap",
                  "freq","af","frq"))
  )
}

# Verify a slim file exists and is readable (regenerate if not).
is_valid_slim <- function(path) {
  if (!file.exists(path) || file.size(path) < 100) return(FALSE)
  tryCatch({
    df <- readr::read_tsv(path, n_max = 3, show_col_types = FALSE,
                          progress = FALSE)
    nrow(df) >= 1 && all(c("SNP","p_value") %in% names(df))
  }, error = function(e) FALSE)
}

# Safe numeric extraction from an S4 object (MendelianRandomization results).
slot_num <- function(obj, candidates, idx = 1L) {
  if (is.null(obj)) return(NA_real_)
  for (nm in candidates) {
    val <- tryCatch(methods::slot(obj, nm), error = function(e) NULL)
    if (!is.null(val) && length(val) >= idx)
      return(suppressWarnings(as.numeric(val[idx])))
  }
  for (nm in candidates) {
    val <- tryCatch(obj[[nm]], error = function(e) NULL)
    if (!is.null(val) && length(val) >= idx)
      return(suppressWarnings(as.numeric(val[idx])))
  }
  NA_real_
}


# ---------------------------------------------------------------------------
# 4.  PHASE 0 — SLIM GWAS FILES
# ---------------------------------------------------------------------------
# Reduces each full GWAS to seven canonical columns:
#   SNP, effect_allele, other_allele, effect_allele_frequency,
#   beta, standard_error, p_value
# Written as gzip. Skipped if the slim file already exists and is valid.

prepare_slim_files <- function(force = FALSE) {
  for (trait in names(gwas_paths)) {
    slim <- file.path(SLIM_DIR, paste0(trait, ".tsv.gz"))
    if (!force && is_valid_slim(slim)) {
      message("  [slim ", trait, "] valid, ",
              round(file.size(slim)/1e6, 1), " MB")
      next
    }
    
    orig <- gwas_paths[[trait]]
    message("  [slim ", trait, "] reading ", basename(orig), " ...")
    cols <- detect_columns(orig)
    
    if (is.null(cols$snp))  stop("[", trait, "] no SNP column")
    if (is.null(cols$ea) || is.null(cols$oa))
      stop("[", trait, "] no allele columns")
    if (is.null(cols$beta) && is.null(cols$or))
      stop("[", trait, "] no beta or OR column")
    if (is.null(cols$p))    stop("[", trait, "] no p-value column")
    
    wanted <- unique(c(cols$snp, cols$ea, cols$oa, cols$beta, cols$or,
                       cols$se, cols$p, cols$eaf))
    wanted <- wanted[!is.na(wanted)]
    
    dt <- tryCatch(
      data.table::fread(orig, select = wanted, showProgress = FALSE,
                        na.strings = c("NA","NaN",".","","-","--","N/A","na")),
      error = function(e) {
        warning("[", trait, "] selective read failed (", conditionMessage(e),
                "); falling back to full read")
        data.table::fread(orig, showProgress = FALSE,
                          na.strings = c("NA","NaN",".","","-","--","N/A","na"))
      }
    )
    
    message("    read ", fmt_num(nrow(dt)), " rows x ", ncol(dt), " cols")
    
    if (cols$snp != "SNP")           setnames(dt, cols$snp, "SNP")
    if (cols$ea  != "effect_allele") setnames(dt, cols$ea,  "effect_allele")
    if (cols$oa  != "other_allele")  setnames(dt, cols$oa,  "other_allele")
    if (cols$p   != "p_value")       setnames(dt, cols$p,   "p_value")
    if (!is.null(cols$eaf) && cols$eaf != "effect_allele_frequency")
      setnames(dt, cols$eaf, "effect_allele_frequency")
    if (!is.null(cols$beta) && cols$beta != "beta")
      setnames(dt, cols$beta, "beta")
    if (!is.null(cols$or) && cols$or != "odds_ratio")
      setnames(dt, cols$or, "odds_ratio")
    if (!is.null(cols$se) && cols$se != "standard_error")
      setnames(dt, cols$se, "standard_error")
    
    # p-value cleanup (underflow to .Machine$double.xmin)
    dt[, p_value := clean_p_vec(p_value)]
    dt[p_value == 0 & !is.na(p_value), p_value := .Machine$double.xmin]
    
    # beta
    if (!"beta" %in% names(dt)) {
      if (!"odds_ratio" %in% names(dt)) stop("[", trait, "] no beta/OR")
      dt[, beta := log(as.numeric(odds_ratio))]
    } else {
      dt[, beta := as.numeric(beta)]
    }
    
    # SE
    if (!"standard_error" %in% names(dt)) {
      dt[, standard_error := NA_real_]
    } else {
      dt[, standard_error := as.numeric(standard_error)]
    }
    if (any(is.na(dt$standard_error)) && "odds_ratio" %in% names(dt)) {
      na_idx <- is.na(dt$standard_error)
      dt[na_idx, standard_error := compute_se_from_or_and_p(
        as.numeric(odds_ratio[na_idx]), as.numeric(p_value[na_idx]))]
    }
    
    if (!"effect_allele_frequency" %in% names(dt))
      dt[, effect_allele_frequency := NA_real_]
    dt[, effect_allele_frequency := as.numeric(effect_allele_frequency)]
    
    dt[, effect_allele := toupper(as.character(effect_allele))]
    dt[, other_allele  := toupper(as.character(other_allele))]
    
    n_sig <- sum(!is.na(dt$p_value) & dt$p_value < P_THRESH)
    message("    p < 5e-8 in source: ", fmt_num(n_sig))
    
    dt <- dt[
      !is.na(beta) & !is.na(standard_error) & standard_error > 0 &
        !is.na(p_value) & p_value > 0 & p_value <= 1 &
        !is.na(effect_allele) & !is.na(other_allele) &
        nchar(effect_allele) == 1 & nchar(other_allele) == 1 &
        effect_allele != other_allele
    ]
    
    dt <- dt[, .(SNP, effect_allele, other_allele, effect_allele_frequency,
                 beta, standard_error, p_value)]
    dt <- unique(dt, by = "SNP")
    message("    retained: ", fmt_num(nrow(dt)), " unique SNPs")
    
    if (nrow(dt) == 0)
      stop("[", trait, "] 0 rows after cleaning — check column mapping")
    
    data.table::fwrite(dt, slim, sep = "\t", quote = FALSE, na = "NA")
    rm(dt); invisible(gc(verbose = FALSE))
    message("    slim OK (", round(file.size(slim)/1e6, 1), " MB)")
  }
  invisible(gc(verbose = FALSE))
}


# ---------------------------------------------------------------------------
# 5.  STREAMING READER FOR SLIM FILES
# ---------------------------------------------------------------------------
# Reads the slim file in chunks and applies SNP and/or p-value filters.
# Memory footprint scales with chunk size (default 500k rows), not with the
# full GWAS size.

read_slim_filtered <- function(path, snp_filter = NULL, p_filter = NULL) {
  if (!file.exists(path)) stop("Slim file not found: ", path)
  
  results <- list()
  callback <- function(chunk, pos) {
    if (!is.null(snp_filter))
      chunk <- chunk[chunk$SNP %in% snp_filter, , drop = FALSE]
    if (!is.null(p_filter))
      chunk <- chunk[!is.na(chunk$p_value) & chunk$p_value < p_filter, ,
                     drop = FALSE]
    if (nrow(chunk) > 0)
      results[[length(results) + 1]] <<- chunk
    NULL
  }
  
  readr::read_tsv_chunked(
    path,
    callback = readr::SideEffectChunkCallback$new(callback),
    chunk_size = CHUNK_SIZE,
    col_types = readr::cols(
      SNP                     = readr::col_character(),
      effect_allele           = readr::col_character(),
      other_allele            = readr::col_character(),
      effect_allele_frequency = readr::col_double(),
      beta                    = readr::col_double(),
      standard_error          = readr::col_double(),
      p_value                 = readr::col_double()
    ),
    progress = FALSE
  )
  
  if (length(results) == 0) return(tibble())
  dplyr::bind_rows(results)
}


# ---------------------------------------------------------------------------
# 6.  LD CLUMPING
# ---------------------------------------------------------------------------

clump_snps <- function(df) {
  if (nrow(df) == 0) return(df)
  args <- list(
    dat      = df %>% dplyr::rename(rsid = SNP, pval = p_value),
    clump_kb = CLUMP_KB,
    clump_r2 = CLUMP_R2,
    clump_p  = P_THRESH
  )
  if (isTRUE(USE_LOCAL_PLINK)) {
    args$plink_bin <- PLINK_BIN
    args$bfile     <- EUR_BFILE
  }
  out <- tryCatch(
    do.call(ieugwasr::ld_clump, args),
    error = function(e)
      stop("ld_clump failed (mode = ",
           if (USE_LOCAL_PLINK) "LOCAL PLINK" else "OpenGWAS API",
           "): ", conditionMessage(e), call. = FALSE)
  )
  as.data.frame(out) %>% dplyr::rename(SNP = rsid, p_value = pval)
}


# ---------------------------------------------------------------------------
# 7.  PHASE 1 — INSTRUMENT CACHE
# ---------------------------------------------------------------------------

get_instruments_cached <- function(trait) {
  cache_file <- file.path(CACHE_DIR, paste0("instr_", trait, ".rds"))
  if (file.exists(cache_file)) {
    ins <- readRDS(cache_file)
    message("  [", trait, "] from cache: ", nrow(ins), " SNPs")
    return(ins)
  }
  
  message("  [", trait, "] streaming + filter p<5e-8 + clumping ...")
  sig <- read_slim_filtered(file.path(SLIM_DIR, paste0(trait, ".tsv.gz")),
                            p_filter = P_THRESH)
  message("    significant SNPs: ", fmt_num(nrow(sig)))
  if (nrow(sig) == 0) { saveRDS(sig, cache_file); return(sig) }
  
  clumped <- clump_snps(sig)
  rm(sig); invisible(gc(verbose = FALSE))
  message("    after clumping: ", nrow(clumped), " SNPs")
  saveRDS(clumped, cache_file)
  clumped
}


# ---------------------------------------------------------------------------
# 8.  FORMAT HELPERS
# ---------------------------------------------------------------------------

format_exposure <- function(df, id_label) {
  format_data(
    df, type = "exposure",
    snp_col = "SNP", beta_col = "beta", se_col = "standard_error",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency", pval_col = "p_value"
  ) %>% mutate(id.exposure = id_label)
}

format_outcome <- function(df, id_label) {
  format_data(
    df, type = "outcome",
    snp_col = "SNP", beta_col = "beta", se_col = "standard_error",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency", pval_col = "p_value"
  ) %>% mutate(id.outcome = id_label)
}


# ---------------------------------------------------------------------------
# 9.  PHASE 2 — MVMR FOR ONE DIRECTION
# ---------------------------------------------------------------------------

run_mvmr_direction <- function(exposure_name, outcome_name) {
  direction_label <- paste0(exposure_name, "_plus_BMI_to_", outcome_name)
  message("\n====================================================")
  message("MVMR: ", direction_label)
  message("====================================================")
  
  exp_instr <- readRDS(file.path(CACHE_DIR,
                                 paste0("instr_", exposure_name, ".rds")))
  bmi_instr <- readRDS(file.path(CACHE_DIR, "instr_bmi.rds"))
  if (nrow(exp_instr) == 0 || nrow(bmi_instr) == 0) {
    warning("No instruments for ", direction_label); return(NULL)
  }
  
  union_snps <- unique(c(exp_instr$SNP, bmi_instr$SNP))
  message("Union SNPs: ", length(union_snps))
  
  exp_sub <- read_slim_filtered(
    file.path(SLIM_DIR, paste0(exposure_name, ".tsv.gz")),
    snp_filter = union_snps)
  bmi_sub <- read_slim_filtered(
    file.path(SLIM_DIR, "bmi.tsv.gz"), snp_filter = union_snps)
  out_sub <- read_slim_filtered(
    file.path(SLIM_DIR, paste0(outcome_name, ".tsv.gz")),
    snp_filter = union_snps)
  message("  subsets: ", nrow(exp_sub), " / ", nrow(bmi_sub),
          " / ", nrow(out_sub), " SNPs")
  
  if (nrow(exp_sub) < 3 || nrow(bmi_sub) < 3 || nrow(out_sub) < 3) {
    warning("Too few SNPs for ", direction_label)
    rm(exp_sub, bmi_sub, out_sub); invisible(gc(verbose = FALSE))
    return(NULL)
  }
  
  out_dat  <- format_outcome(out_sub, outcome_name)
  harm_exp <- harmonise_data(format_exposure(exp_sub, exposure_name), out_dat)
  harm_bmi <- harmonise_data(format_exposure(bmi_sub, "BMI"),           out_dat)
  rm(exp_sub, bmi_sub, out_sub, out_dat); invisible(gc(verbose = FALSE))
  
  harm_exp <- harm_exp[harm_exp$mr_keep, ]
  harm_bmi <- harm_bmi[harm_bmi$mr_keep, ]
  
  common <- intersect(harm_exp$SNP, harm_bmi$SNP)
  if (length(common) < 3) {
    warning("Fewer than 3 common SNPs for ", direction_label); return(NULL)
  }
  harm_exp <- harm_exp[match(common, harm_exp$SNP), ]
  harm_bmi <- harm_bmi[match(common, harm_bmi$SNP), ]
  
  BXGs   <- cbind(harm_exp$beta.exposure, harm_bmi$beta.exposure)
  seBXGs <- cbind(harm_exp$se.exposure,   harm_bmi$se.exposure)
  BR     <- harm_exp$beta.outcome
  seBR   <- harm_exp$se.outcome
  
  # --- MVMR (package MVMR): IVW + conditional F + Q ------------------------
  mvmr_input <- format_mvmr(BXGs   = BXGs,
                            BYG    = BR,
                            seBXGs = seBXGs,
                            seBYG  = seBR,
                            RSID   = common)
  
  gcov <- switch(GCOV_APPROX, zero = 0, identity = 0, 0)
  
  ivw_res <- ivw_mvmr(mvmr_input)
  cond_f  <- tryCatch(suppressWarnings(strength_mvmr(mvmr_input, gencov = gcov)),
                      error = function(e) { warning("strength_mvmr: ",
                                                    conditionMessage(e)); NULL })
  q_res   <- tryCatch(suppressWarnings(pleiotropy_mvmr(mvmr_input, gencov = gcov)),
                      error = function(e) { warning("pleiotropy_mvmr: ",
                                                    conditionMessage(e)); NULL })
  
  # --- MVMR-Egger and MVMR-Median (MendelianRandomization package) ---------
  egger_res  <- NULL
  median_res <- NULL
  mv_input <- tryCatch(
    MendelianRandomization::mr_mvinput(bx = BXGs, bxse = seBXGs,
                                       by = BR,   byse = seBR),
    error = function(e) { warning("mr_mvinput: ", conditionMessage(e)); NULL }
  )
  if (!is.null(mv_input)) {
    egger_res <- tryCatch(
      MendelianRandomization::mr_mvegger(mv_input),
      error = function(e) { warning("mr_mvegger: ", conditionMessage(e)); NULL }
    )
    median_res <- tryCatch(
      MendelianRandomization::mr_mvmedian(mv_input),
      error = function(e) { warning("mr_mvmedian: ", conditionMessage(e)); NULL }
    )
  }
  
  # --- save harmonised data and all results --------------------------------
  saveRDS(list(
    direction  = direction_label,
    exposure   = harm_exp,
    bmi        = harm_bmi,
    common     = common,
    mvmr_input = mvmr_input,
    ivw        = ivw_res,
    cond_f     = cond_f,
    q          = q_res,
    egger      = egger_res,
    median     = median_res,
    meta       = list(
      n_snp        = length(common),
      gcov_approx  = GCOV_APPROX,
      mvmr_version = as.character(packageVersion("MVMR")),
      clump_mode   = if (USE_LOCAL_PLINK) "local_plink" else "opengwas_api",
      timestamp    = Sys.time()
    )
  ),
  file = file.path(RESULTS_DIR,
                   paste0("mvmr_harmonised_", direction_label, ".rds")))
  
  # --- build the two summary rows (exposure + BMI) -------------------------
  build_row <- function(expo_label, i, is_binary) {
    est  <- ivw_res[i, "Estimate"]
    se   <- ivw_res[i, "Std. Error"]
    pval <- ivw_res[i, "Pr(>|t|)"]
    
    f_i <- if (is.null(cond_f)) {
      NA_real_
    } else if (is.data.frame(cond_f) || is.matrix(cond_f)) {
      if (nrow(cond_f) >= 1 && ncol(cond_f) >= i)
        suppressWarnings(as.numeric(cond_f[1, i]))
      else if (nrow(cond_f) >= i && ncol(cond_f) >= 1)
        suppressWarnings(as.numeric(cond_f[i, 1]))
      else NA_real_
    } else if (is.atomic(cond_f) && length(cond_f) >= i) {
      as.numeric(cond_f[i])
    } else NA_real_
    
    egger_int_est <- slot_num(egger_res, c("Intercept"))
    egger_int_se  <- slot_num(egger_res, c("StdError.Int", "StdError.Intercept"))
    egger_int_p   <- slot_num(egger_res, c("Pvalue.Int",   "Pvalue.Intercept"))
    
    median_beta <- slot_num(median_res, c("Estimate", "Est"))
    median_p    <- slot_num(median_res, c("Pvalue", "Pvalue.Est"), idx = i)
    if (!is.na(median_beta) && i == 2) {
      # median_res is per-exposure; second exposure uses index 2
      median_beta <- slot_num(median_res, c("Estimate", "Est"), idx = 2)
      median_p    <- slot_num(median_res, c("Pvalue", "Pvalue.Est"), idx = 2)
    }
    
    tibble(
      direction                 = direction_label,
      exposure                  = expo_label,
      outcome                   = outcome_name,
      n_snp                     = length(common),
      beta                      = est,
      se                        = se,
      p_value                   = pval,
      OR_per_unit               = exp(est),
      OR_per_unit_lci           = exp(est - 1.96 * se),
      OR_per_unit_uci           = exp(est + 1.96 * se),
      OR_per_doubling           = if (is_binary) exp(est * log(2)) else NA_real_,
      OR_per_doubling_lci       = if (is_binary) exp((est - 1.96*se) * log(2)) else NA_real_,
      OR_per_doubling_uci       = if (is_binary) exp((est + 1.96*se) * log(2)) else NA_real_,
      F_conditional             = f_i,
      Q_stat                    = if (!is.null(q_res) && !is.null(q_res$Qstat))
        as.numeric(q_res$Qstat) else NA_real_,
      Q_df                      = as.integer(length(common) - N_EXPOSURES),
      Q_pval                    = if (!is.null(q_res) && !is.null(q_res$Qpval))
        as.numeric(q_res$Qpval) else NA_real_,
      egger_intercept           = egger_int_est,
      egger_intercept_se        = egger_int_se,
      egger_intercept_p         = egger_int_p,
      median_beta               = median_beta,
      median_p_value            = median_p,
      median_OR_per_unit        = exp(median_beta),
      median_OR_per_doubling    = if (is_binary && !is.na(median_beta))
        exp(median_beta * log(2)) else NA_real_,
      gwas_exposure             = gwas_origins[[exposure_name]],
      gwas_covariate_BMI        = gwas_origins[["bmi"]],
      gwas_outcome              = gwas_origins[[outcome_name]],
      gcov_approx               = GCOV_APPROX,
      clump_mode                = if (USE_LOCAL_PLINK) "local_plink" else "opengwas_api",
      mvmr_version              = as.character(packageVersion("MVMR"))
    )
  }
  
  res <- bind_rows(
    build_row(exposure_name, 1, is_binary = TRUE),
    build_row("BMI",         2, is_binary = FALSE)
  )
  
  rm(exp_instr, bmi_instr, harm_exp, harm_bmi, mvmr_input,
     ivw_res, cond_f, q_res, egger_res, median_res)
  invisible(gc(verbose = FALSE))
  res
}


# ---------------------------------------------------------------------------
# 10.  MAIN
# ---------------------------------------------------------------------------

message("=== PHASE 0: pre-slim GWAS files ===")
prepare_slim_files()
message("=== PHASE 0 done ===\n")

message("=== PHASE 1: instrument cache ===")
for (trait in names(gwas_paths)) get_instruments_cached(trait)
invisible(gc(verbose = FALSE))
message("=== PHASE 1 done ===\n")

message("=== PHASE 2: MVMR per direction ===")
directions <- list(
  c("neck", "knee"),
  c("neck", "hand"),
  c("neck", "hip"),
  c("knee", "neck"),
  c("hip",  "neck"),
  c("hand", "neck")
)

results_list <- list()
for (i in seq_along(directions)) {
  d <- directions[[i]]
  row <- tryCatch(
    run_mvmr_direction(d[1], d[2]),
    error = function(e) {
      message("!! ", paste(d, collapse = " -> "), " failed: ",
              conditionMessage(e))
      NULL
    }
  )
  if (!is.null(row)) results_list[[length(results_list) + 1]] <- row
  if (exists("row")) rm(row)
  invisible(gc(verbose = FALSE))
}

mvmr_results <- bind_rows(results_list)
if (nrow(mvmr_results) == 0) stop("MVMR results table is empty")

write_csv(mvmr_results, file.path(RESULTS_DIR, "mvmr_results.csv"))
saveRDS(mvmr_results,   file.path(RESULTS_DIR, "mvmr_results.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(RESULTS_DIR, "mvmr_sessionInfo.txt"))

message("\n== Done ==")
message("Results -> ", file.path(RESULTS_DIR, "mvmr_results.csv"))
message("Session info -> ", file.path(RESULTS_DIR, "mvmr_sessionInfo.txt"))
message("Clumping: ", if (USE_LOCAL_PLINK)
  paste0("LOCAL PLINK (", PLINK_BIN, ")") else "OpenGWAS API")

# ---------------------------------------------------------------------------
# End of script
# ---------------------------------------------------------------------------