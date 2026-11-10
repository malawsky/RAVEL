pkgload::load_all('.',quiet=TRUE)
RAVEL:::.ensure_lavaan_runtime()
source('tests/testthat/helper-family-validation.R')
options(ravel.workflow_functions_only=TRUE)
e<-new.env(parent=globalenv());sys.source('analysis/parental_deviation/ravel_parental_deviation.R',e)
d<-family_validation_data(family_validation_population(),N=500L,seed=420L)
d$IID<-d$FID
roles<-c('child','mother','father')
pgs<-setNames(unlist(lapply(roles,function(r)paste0('PGS',1:2,'_',r))),unlist(lapply(roles,function(r)paste0(r,1:2))))
for(k in 1:2)for(r in roles)d[[paste0('PGS',k,'_',r)]]<-d[[paste0('PGS',k,'_',r)]]*10^(-4-k)
z<-e$ravel_parental_deviation(d,outcome='Y',pgs_columns=pgs,subject_id='IID',pi=.4,pi_se=NULL,fit_constraints='estimated_state',include_comparison=TRUE)
st<-attr(z,'settings');fit<-attr(z,'fit');stopifnot(fit$converged,nrow(st$score_scaling)==6L)
tab<-e$.ravel_report_pgs_structure(fit,st)
S<-lavaan::lavInspect(fit$fit,'cov.ov');v<-st$score_scaling$scale;names(v)<-st$score_scaling$canonical
rows<-tab[tab$block=='pgs_input_covariance',]
expected<-S[cbind(rows$lhs,rows$rhs)]*v[rows$lhs]*v[rows$rhs]
stopifnot(isTRUE(all.equal(rows$estimate,unname(expected),check.attributes=FALSE)),all(rows$scale=='Original input PGS units'))
stopifnot(any(grepl('Plain PGS1',z$analysis)),any(grepl('Plain PGS2',z$analysis)))
saveRDS(z,'/tmp/ravel-fresh-pgs-workflow.rds')
cat('Fresh raw-PGS workflow, scaling metadata, both vanilla regressions: passed\n')
