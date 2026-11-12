# Reproduce MR-PRESSO and MR-Lasso from harmonised primary MR data.
# Run from project root: Rscript R/primary_PRESSO_Lasso.R
required <- c("MRPRESSO", "MendelianRandomization")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

pairs <- data.frame(exposure = c("neck","neck","neck","knee","hand","hip"),
                    outcome  = c("knee","hand","hip","neck","neck","neck"),
                    stringsAsFactors = FALSE)
dir.create("results/reproduced", showWarnings = FALSE, recursive = TRUE)

presso_rows <- list()
lasso_rows <- list()
full_presso <- list()

for (i in seq_len(nrow(pairs))) {
  ex <- pairs$exposure[i]; out <- pairs$outcome[i]
  key <- paste(ex, out, sep = "_to_")
  path <- sprintf("results/harmonised_%s_%s.rds", ex, out)
  cat("\n===", key, "===\n")
  if (!file.exists(path)) stop("Missing input: ", path)
  dat <- readRDS(path)
  needed <- c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
  if (!all(needed %in% names(dat))) stop("Missing required columns in ", path)
  dat <- dat[stats::complete.cases(dat[, needed]), , drop = FALSE]
  if (anyDuplicated(dat$SNP)) stop("Duplicate SNPs in ", path)
  cat("SNPs:", nrow(dat), "\n")
  if (nrow(dat) < 4L) stop("Too few instruments in ", path)

  # Fixed seed makes the Monte Carlo calculation repeatable under the same
  # R/package/RNG versions; it is NOT guaranteed to match the original run.
  set.seed(20261111 + i)
  press <- tryCatch(
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat,
      NbDistribution = 10000, SignifThreshold = 0.05
    ), error = function(e) e)

  if (inherits(press, "error")) {
    warning("MR-PRESSO ", key, ": ", conditionMessage(press))
    presso_rows[[key]] <- data.frame(Direction = key, Type = NA_character_,
      Beta = NA_real_, P = NA_real_, SNPs = nrow(dat),
      Error = conditionMessage(press))
  } else {
    full_presso[[key]] <- press
    main <- as.data.frame(press$`Main MR results`)
    typ <- if ("MR Analysis" %in% names(main)) as.character(main$`MR Analysis`) else rownames(main)
    if (is.null(typ) || any(!nzchar(typ))) typ <- paste0("row_", seq_len(nrow(main)))
    presso_rows[[key]] <- data.frame(Direction = key, Type = typ,
      Beta = as.numeric(main$`Causal Estimate`),
      P = as.numeric(main$`P-value`), SNPs = nrow(dat), Error = NA_character_)
    cat("MR-PRESSO:\n")
    print(presso_rows[[key]], row.names = FALSE)
  }

  lass <- tryCatch({
    inp <- MendelianRandomization::mr_input(
      bx = dat$beta.exposure, bxse = dat$se.exposure,
      by = dat$beta.outcome, byse = dat$se.outcome, snps = dat$SNP)
    MendelianRandomization::mr_lasso(inp)
  }, error = function(e) e)
  if (inherits(lass, "error")) {
    warning("MR-Lasso ", key, ": ", conditionMessage(lass))
    lasso_rows[[key]] <- data.frame(Direction = key, Beta = NA_real_,
       P = NA_real_, SNPs = nrow(dat), Error = conditionMessage(lass))
  } else {
    lasso_rows[[key]] <- data.frame(Direction = key,
      Beta = as.numeric(lass@Estimate), P = as.numeric(lass@Pvalue),
      SNPs = nrow(dat), Error = NA_character_)
    cat("MR-Lasso:\n"); print(lasso_rows[[key]], row.names = FALSE)
  }
}

presso <- do.call(rbind, presso_rows); rownames(presso) <- NULL
lasso <- do.call(rbind, lasso_rows); rownames(lasso) <- NULL
write.csv(presso, "results/reproduced/MR_PRESSO_all_estimates.csv", row.names = FALSE, na = "")
write.csv(lasso, "results/reproduced/MR_Lasso_estimates.csv", row.names = FALSE, na = "")
saveRDS(full_presso, "results/reproduced/MR_PRESSO_full_results.rds")
cat("\nSaved to results/reproduced/\n")
