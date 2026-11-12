# primary_diagnostics.R: machine-readable outputs supporting manuscript figures
# Requires harmonised RDS saved by primary_MR.R. Reruns diagnostic calculations.
args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else getwd()
setwd(project)
suppressPackageStartupMessages(library(TwoSampleMR))
pairs <- list(c('neck','knee'), c('neck','hand'), c('neck','hip'),
              c('knee','neck'), c('hand','neck'), c('hip','neck'))
out <- file.path('results','tables'); dir.create(out, recursive=TRUE, showWarnings=FALSE)
concat <- function(xs) do.call(rbind, xs)
f_list <- list(); f_sum <- list(); diag_list <- list(); loo_list <- list(); scatter_list <- list()
for (pair in pairs) {
  exposure <- pair[1]; outcome <- pair[2]
  fn <- file.path('results', paste0('harmonised_',exposure,'_',outcome,'.rds'))
  if (!file.exists(fn)) stop('Missing harmonised dataset: ',fn)
  d <- readRDS(fn)
  d <- d[!is.na(d$mr_keep) & d$mr_keep &
         is.finite(d$beta.exposure) & is.finite(d$se.exposure) &
         is.finite(d$beta.outcome) & is.finite(d$se.outcome), ,drop=FALSE]
  if (!nrow(d)) stop('No variants for ',exposure,' -> ',outcome)
  pair_id <- paste(exposure,outcome,sep='_to_')
  f <- (d$beta.exposure/d$se.exposure)^2
  f_list[[pair_id]] <- data.frame(Exposure=exposure, SNP=d$SNP,
      beta.exposure=d$beta.exposure,se.exposure=d$se.exposure,
      F_statistic=f,p=d$pval.exposure)
  f_sum[[pair_id]] <- data.frame(Exposure=exposure,Outcome=outcome,SNPs=length(f),
      Mean_F=mean(f),Median_F=median(f),Min_F=min(f),Max_F=max(f))
  hetero <- TwoSampleMR::mr_heterogeneity(d)
  egger <- TwoSampleMR::mr_pleiotropy_test(d)
  hetero$Exposure <- exposure; hetero$Outcome <- outcome
  # Retain entire tests (including IVW and Egger heterogeneity), rather than only row 1.
  egger$Exposure <- exposure; egger$Outcome <- outcome
  diag_list[[pair_id]] <- list(heterogeneity=hetero, intercept=egger)
  loo <- tryCatch(TwoSampleMR::mr_leaveoneout(d),error=function(e) {
       warning('Leave-one-out unavailable for ',pair_id,': ',conditionMessage(e)); NULL })
  if (!is.null(loo)) {loo$Exposure <- exposure;loo$Outcome <- outcome; loo_list[[pair_id]] <- loo}
  scatter_list[[pair_id]] <- data.frame(Exposure=exposure,Outcome=outcome,SNP=d$SNP,
       beta_exposure=d$beta.exposure,se_exposure=d$se.exposure,
       beta_outcome=d$beta.outcome,se_outcome=d$se.outcome)
}
write.csv(concat(f_list),file.path(out,'S2_per_SNP_F.csv'),row.names=FALSE)
write.csv(concat(f_sum),file.path(out,'S3_F_summary.csv'),row.names=FALSE)
write.csv(concat(lapply(diag_list,`[[`,'heterogeneity')),
          file.path(out,'primary_heterogeneity.csv'),row.names=FALSE)
write.csv(concat(lapply(diag_list,`[[`,'intercept')),
          file.path(out,'primary_Egger_intercepts.csv'),row.names=FALSE)
if(length(loo_list)) write.csv(concat(loo_list),file.path(out,'primary_leave_one_out.csv'),row.names=FALSE)
write.csv(concat(scatter_list),file.path(out,'primary_scatter_SNP_data.csv'),row.names=FALSE)
cat('\n=== Supplementary S2 per-SNP F ===\n'); print(concat(f_list))
cat('\n=== Supplementary S3 summary F ===\n'); print(concat(f_sum))
cat('\n=== Heterogeneity ===\n'); print(concat(lapply(diag_list,`[[`,'heterogeneity')))
cat('\n=== MR Egger intercepts ===\n'); print(concat(lapply(diag_list,`[[`,'intercept')))
cat('\nLOO and scatter inputs written to CSV.\n')
