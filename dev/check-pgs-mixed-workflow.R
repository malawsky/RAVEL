# Fresh synthetic mixed fit plus the same reporting/export path used by the script.
pkgload::load_all('.',quiet=TRUE)
RAVEL:::.ensure_lavaan_runtime()
source('tests/testthat/helper-mixed-validation.R')
options(ravel.workflow_functions_only=TRUE)
e<-new.env(parent=globalenv());sys.source('analysis/parental_deviation/ravel_parental_deviation.R',e)
d<-mixed_validation_data(N_families=300L,family_variance=.5,sex_effect=0,age_effect=0,seed=715L)
roles<-c('child','mother','father')
pgs<-setNames(unlist(lapply(roles,function(r)paste0('PGS',1:2,'_',r))),unlist(lapply(roles,function(r)paste0(r,1:2))))
z<-e$ravel_parental_deviation(d,outcome='Y',pgs_columns=pgs,subject_id='IID',maternal_id='M_ID',
  family_id='FID',random_family=TRUE,pi=.4,pi_se=NULL,fit_constraints='estimated_state',include_comparison=TRUE)
stopifnot(any(grepl('Plain PGS1',z$analysis)),any(grepl('Plain PGS2',z$analysis)))
comp<-attr(z,'pgs_comparators');stopifnot(length(comp$fits)==2L,is.null(comp$error),all(is.finite(comp$table$estimate)))
r<-attr(z,'fit_records');for(i in seq_along(r)) {
  r[[i]]$model<-'Estimated state';r[[i]]$outcome<-'Y';r[[i]]$parental_sibs<-'none';r[[i]]$genetic_roles<-roles
  r[[i]]$settings<-attr(z,'settings');r[[i]]$diagnostic<-attr(z,'one_factor');r[[i]]$pgs_comparators<-comp
}
attr(z,'report_fit_records')<-r
report<-e$ravel_reporting_tables(z)
stopifnot(length(unique(report$fit_statistics$source_fit[report$fit_statistics$engine=='pgs_comparator']))==2L)
stopifnot(all(c('pgs_input_covariance','pgs_correlation','observed_PGS_covariance','observed_PGS_correlation')%in%report$genetic_structure$block))
aux<-report$fit_statistics[report$fit_statistics$engine=='lavaan_auxiliary',]
stopifnot(any(aux$statistic=='AIC' & is.finite(aux$value)),any(aux$statistic=='rmsea' & is.finite(aux$value)))
pgsstats<-report$fit_statistics[report$fit_statistics$engine=='pgs_comparator',]
stopifnot(all(is.finite(pgsstats$value[pgsstats$statistic=='AIC'])))
folder<-tempfile('ravel-pgs-mixed-export-',tmpdir='/tmp');dir.create(folder)
paths<-e$write_ravel_workflow_outputs(z,file.path(folder,'result'))
stopifnot(length(paths)==7L,all(file.exists(unlist(paths))))
saveRDS(z,'/tmp/ravel-pgs-mixed-workflow.rds')
cat('Fresh mixed workflow: both vanilla regressions, auxiliary fit indices, PGS moments and Excel exports passed.\n')
cat(folder,'\n')
