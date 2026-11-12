# replication_MR.R
# Extracted/adapted from user's Org-mode research notebooks.
# replication.org final tables block, approx. lines 3514-3771
# Requires 12 harmonised phenotype datasets plus Tao exploratory exposure.tsv.
# Calculates replication results and saves machine-readable CSV, no Word.


args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else getwd()
if (!file.exists(file.path(project, "README.md"))) stop("Run with project root as argument or working directory")
setwd(project)

library(TwoSampleMR)
library(dplyr)
library(readr)
library(tidyr)

ROOT <- Sys.getenv(
  'REPLICATION_HARMONISED_DIR',
  unset = file.path(project, 'data', 'processed', 'replication')
)
OUT <- file.path(project, 'results', 'tables')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
phenotypes <- c('KneeOA', 'HandOA', 'HipOA')
analysis_order <- c('Meng strict', 'Meng exploratory', 'Tao strict', 'Tao exploratory')

spec <- tibble::tribble(
  ~Analysis,           ~Directory,                     ~Extension, ~Threshold, ~Expected,
  'Meng strict',       'sensitivity_noUKB',             'csv',      '5e-8',      3L,
  'Meng exploratory',  'sensitivity_Meng2020_5e-7',     'csv',      '5e-7',     10L,
  'Tao strict',        'sensitivity_Tao2025',           'csv',      '5e-8',      3L,
  'Tao exploratory',   'Tao_25SNP_GO1',                 'tsv',      '5e-7',     NA_integer_
)

read_group <- function(i, ph) {
  s <- spec[i, ]
  path <- file.path(ROOT, s$Directory, paste0(ph, '_harmonised.', s$Extension))
  if (!file.exists(path)) stop('Missing file: ', path)
  d <- if (s$Extension == 'csv') read_csv(path, show_col_types = FALSE) else
    read_tsv(path, show_col_types = FALSE)
  needed <- c('SNP', 'mr_keep', 'pval.exposure', 'beta.exposure',
              'se.exposure', 'effect_allele.exposure', 'other_allele.exposure')
  if (!all(needed %in% names(d))) stop('Missing columns: ', path)
  d <- d |> filter(mr_keep %in% TRUE)
  if (s$Analysis == 'Tao exploratory') d <- d |> filter(pval.exposure < 5e-7)
  if (anyDuplicated(d$SNP)) stop('Duplicate SNPs in ', path)
  expected <- if (s$Analysis == 'Tao exploratory') {
    c(KneeOA = 16L, HandOA = 14L, HipOA = 14L)[[ph]]
  } else s$Expected
  if (nrow(d) != expected) {
    stop(s$Analysis, ' / ', ph, ': expected ', expected, ' SNPs; found ', nrow(d))
  }
  d$Analysis <- s$Analysis
  d$Threshold <- s$Threshold
  d$Outcome <- ph
  d
}

all_data <- bind_rows(lapply(seq_len(nrow(spec)), function(i) {
  bind_rows(lapply(phenotypes, function(ph) read_group(i, ph)))
}))

# S4: instrument rows for each analysis, including all 17 selected Tao exploratory SNPs.
# The Tao exposure file is used because 1–3 SNPs were absent after harmonisation.
s6_from_harmonised <- all_data |>
  filter(Analysis != 'Tao exploratory') |>
  distinct(Analysis, SNP, .keep_all = TRUE) |>
  transmute(Analysis, Threshold, SNP,
            Effect_allele = effect_allele.exposure,
            Other_allele = other_allele.exposure,
            Beta = beta.exposure, SE = se.exposure,
            P = pval.exposure, F = (beta.exposure / se.exposure)^2)

tao_path <- file.path(ROOT, 'Tao_25SNP_GO1', 'exposure.tsv')
if (!file.exists(tao_path)) stop('Missing Tao exposure file: ', tao_path)
tao_exp <- read_tsv(tao_path, show_col_types = FALSE) |> filter(P < 5e-7)
if (n_distinct(tao_exp$SNP) != 17L || nrow(tao_exp) != 17L)
  stop('Expected exactly 17 Tao exploratory exposure SNPs')
if (!all(c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P') %in% names(tao_exp)))
  stop('Unexpected Tao exposure column names')
tao_s6 <- tao_exp |>
  transmute(Analysis = 'Tao exploratory', Threshold = '5e-7', SNP,
            Effect_allele = A1, Other_allele = A2,
            Beta = BETA, SE, P, F = (BETA / SE)^2)

s6 <- bind_rows(s6_from_harmonised, tao_s6) |>
  mutate(Analysis = factor(Analysis, levels = analysis_order)) |>
  arrange(Analysis, P) |>
  mutate(Analysis = as.character(Analysis))
if (nrow(s6) != 33L || anyDuplicated(s6[c('Analysis', 'SNP')]))
  stop('S4 must have 33 unique analysis–SNP rows')

used <- all_data |> distinct(Analysis, Outcome, SNP) |>
  mutate(Used = 'Yes') |>
  pivot_wider(names_from = Outcome, values_from = Used, values_fill = 'No')
s6 <- s6 |> left_join(used, by = c('Analysis', 'SNP')) |>
  mutate(across(all_of(phenotypes), ~replace_na(.x, 'No')))

# S5: calculate exactly the same MR methods for all four analyses.
# mr_ivw_mre = multiplicative random-effects IVW; not mr_ivw.
mr_one <- function(d) {
  r <- TwoSampleMR::mr(d, method_list = c('mr_ivw_mre', 'mr_weighted_median'))
  if (nrow(r) != 2L) stop('Unexpected number of MR estimates')
  r |>
    transmute(Analysis = d$Analysis[1], Threshold = d$Threshold[1],
              Outcome = d$Outcome[1],
              Method = case_when(
                method == 'Inverse variance weighted (multiplicative random effects)' ~ 'IVW',
                method == 'Weighted median' ~ 'Weighted median',
                TRUE ~ NA_character_),
              SNPs = nsnp, Beta = b, SE = se,
              OR_10pp = exp(0.1 * b),
              CI_lower = exp(0.1 * (b - 1.96 * se)),
              CI_upper = exp(0.1 * (b + 1.96 * se)), P = pval)
}

groups <- split(all_data, interaction(all_data$Analysis, all_data$Outcome, drop = TRUE))
s7 <- bind_rows(lapply(groups, mr_one)) |>
  mutate(Analysis = factor(Analysis, levels = analysis_order),
         Outcome = factor(Outcome, levels = phenotypes),
         Method = factor(Method, levels = c('IVW', 'Weighted median'))) |>
  arrange(Analysis, Outcome, Method) |>
  mutate(across(c(Analysis, Outcome, Method), as.character))
if (nrow(s7) != 24L || anyNA(s7$Method)) stop('S5 validation failed')

actual <- all_data |> count(Analysis, Outcome, name = 'Actual')
check <- s7 |> distinct(Analysis, Outcome, SNPs) |>
  left_join(actual, by = c('Analysis', 'Outcome'))
if (anyNA(check$Actual) || any(check$SNPs != check$Actual))
  stop('MR SNP counts do not match harmonised data')

# Save complete precision, without Word formatting.
write_csv(s6, file.path(OUT, 'S4_replication_instruments.csv'))
write_csv(s7, file.path(OUT, 'S5_replication_MR.csv'))
cat("\n=== S4: INSTRUMENTS ===\n"); print(s6, n = Inf)
cat("\n=== S5: IVW AND WEIGHTED MEDIAN ===\n"); print(tibble::as_tibble(s7), n = Inf, width = Inf)
cat("\nS4 rows:", nrow(s6), "S5 rows:", nrow(s7), "\n")
