# Check which inputs are present before attempting the reproducibility runs.
args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else getwd()
if (!file.exists(file.path(project, 'README.md'))) stop('Pass the repository root as argument.')
traits <- c('neck', 'knee', 'hand', 'hip')
pairs <- list(c('neck','knee'), c('neck','hand'), c('neck','hip'),
              c('knee','neck'), c('hand','neck'), c('hip','neck'))
primary <- unique(c(file.path('data/processed/primary',paste0(traits,'_exposure.tsv')),
                    vapply(pairs, function(p) file.path('data/processed/primary',paste0(p[2],'_outcome_by_',p[1],'.tsv')),character(1))))
rep_dirs <- c('sensitivity_noUKB', 'sensitivity_Meng2020_5e-7', 'sensitivity_Tao2025', 'Tao_25SNP_GO1')
rep_ext <- c('csv','csv','csv','tsv')
replication <- unlist(Map(function(dir,ext) file.path('data/processed/replication', dir, paste0(c('KneeOA','HandOA','HipOA'),'_harmonised.',ext)),rep_dirs,rep_ext))
replication <- c(replication, 'data/processed/replication/Tao_25SNP_GO1/exposure.tsv')
dirs <- c('neck_plus_BMI_to_knee','neck_plus_BMI_to_hand','neck_plus_BMI_to_hip',
          'knee_plus_BMI_to_neck','hip_plus_BMI_to_neck','hand_plus_BMI_to_neck')
mvmr <- c(file.path('data/processed/mvmr',paste0('mvmr_harmonised_',dirs,'.rds')),
          'data/processed/mvmr/mvmr_results.csv')
expected <- data.frame(stage=c(rep('primary',length(primary)),rep('replication',length(replication)),rep('MVMR',length(mvmr))),
                       file=c(primary,replication,mvmr), stringsAsFactors=FALSE)
expected$present <- file.exists(file.path(project,expected$file))
cat('\n=== Required input availability ===\n')
print(aggregate(present~stage,expected,function(x) sprintf('%d/%d',sum(x),length(x))),row.names=FALSE)
missing <- subset(expected,!present)
if (nrow(missing)) {
  cat('\n=== Missing files ===\n')
  print(missing[,c('stage','file')],row.names=FALSE)
} else cat('\nAll expected input files are present.\n')
cat('\n=== R package availability ===\n')
packages <- c('TwoSampleMR','dplyr','readr','purrr','tidyr','MVMR')
print(data.frame(package=packages,installed=vapply(packages,requireNamespace,logical(1),quietly=TRUE)),row.names=FALSE)
outfile <- file.path(project,'input_check.csv')
write.csv(expected,outfile,row.names=FALSE)
cat('\nReport: ',outfile,'\n',sep='')
cat('Note: successful input check is NOT evidence of numerical reproducibility.\n')
