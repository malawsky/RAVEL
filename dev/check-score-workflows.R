# Run from the repository root: Rscript dev/check-score-workflows.R
# Uses generated data and the repository HEAD for a primary-preprocessing
# regression check; it never runs user SETTINGS or writes fitted RDS objects.
pkgload::load_all('.',quiet=TRUE)
library(testthat)
options(ravel.workflow_functions_only=TRUE)
load_workflow <- function(path) { e<-new.env(parent=globalenv());sys.source(path,envir=e);e }
candidate <- load_workflow('analysis/parental_deviation/ravel_parental_deviation.R')
legacy_source <- system2('git',c('show','HEAD:analysis/parental_deviation/ravel_parental_deviation.R'),stdout=TRUE)
if (!is.null(attr(legacy_source,'status'))) stop('Cannot read the repository HEAD workflow for the regression comparison.')
legacy <- new.env(parent=globalenv());eval(parse(text=legacy_source),legacy)
cols <- c(child1='PGS1_child',child2='PGS2_child',mother1='PGS1_mother',mother2='PGS2_mother',father1='PGS1_father',father2='PGS2_father')
fixture <- function(n=1500L) {
 set.seed(9102L)
 C<-matrix(c(1,.65,.65,.65,1,.3,.65,.3,1),3);K<-matrix(c(1,.5,.5,.5,1,0,.5,0,1),3)
 G<-matrix(rnorm(n*3),n,3)%*%chol(C);E1<-matrix(rnorm(n*3),n,3)%*%chol(K);E2<-matrix(rnorm(n*3),n,3)%*%chol(K)
 d<-data.frame(Y=drop(G%*%c(.32,.13,.08))+rnorm(n),sex=rbinom(n,1,.5),IID=sprintf('I%05d',seq_len(n)),M_ID=sprintf('M%05d',seq_len(n)),FID=sprintf('F%05d',seq_len(n)))
 for(j in 1:3){role<-c('child','mother','father')[j];d[[paste0('PGS1_',role)]]<-G[,j]+sqrt(.4)*E1[,j];d[[paste0('PGS2_',role)]]<-G[,j]+sqrt(1.3)*E2[,j]}
 d
}
record_main <- function(e) {e$main_calls<-list();e$ravel_fit<-function(...) {e$main_calls[[length(e$main_calls)+1L]]<-list(...);stop('Controlled primary failure for integration testing')};invisible(e)}
record_diagnostic <- function(e) {e$diagnostic_calls<-list();e$ravel_score_diagnostic<-function(...) {e$diagnostic_calls[[length(e$diagnostic_calls)+1L]]<-list(...);list(status='fixture_status',reason='Routing fixture',N=NA_integer_)};invisible(e)}
run <- function(e,d,...) suppressMessages(e$ravel_parental_deviation(d,pgs_columns=cols,pi=.4,pi_se=NULL,fit_constraints='equilibrium',pi_grid=NULL,include_comparison=FALSE,...))
record_main(candidate);record_main(legacy)

test_that('auxiliary patch preserves primary preprocessing and table columns',{
 d<-fixture(100);record_diagnostic(candidate)
 old<-run(legacy,d,covariates='sex');new<-run(candidate,d,covariates='sex')
 expect_identical(candidate$main_calls[[1]],legacy$main_calls[[1]])
 expect_identical(names(new),names(old))
 expect_true(any(new$parameter=='one_factor_test',na.rm=TRUE))
 a<-candidate$diagnostic_calls[[1]]
 expect_identical(a$dat[unname(cols)],candidate$main_calls[[1]][[1]][unname(cols)])
 expect_identical(a$covariates,'sex');expect_identical(a$cluster,'FID');expect_identical(a$subject_id,'IID')
})

test_that('original outer FID and missing records survive legacy preprocessing',{
 d<-fixture(280);d$FID<-rep(sprintf('F%02d',1:4),each=70);d$Y[7]<-NA;d$FID[8]<-NA
 record_main(candidate);record_diagnostic(candidate)
 z<-run(candidate,d,family_id='M_ID')
 a<-candidate$diagnostic_calls[[1]]
 expect_identical(a$cluster,'FID');expect_identical(a$dat$FID,d$FID)
 expect_equal(nrow(a$dat),280);expect_true(is.na(a$dat$Y[7]))
 expect_false('FID'%in%names(candidate$main_calls[[1]][[1]]))
 expect_identical(attr(z,'one_factor')$workflow$primary_excluded_rows,7L)
})

test_that('large real families keep primary clustering and a single family errors explicitly',{
 d<-fixture(280);d$FID<-rep(sprintf('F%02d',1:4),each=70)
 record_main(candidate);record_diagnostic(candidate)
 z<-run(candidate,d)
 expect_identical(candidate$main_calls[[1]]$cluster,'FID')
 expect_true('FID'%in%names(candidate$main_calls[[1]][[1]]))
 expect_true(attr(z,'settings')$clustered_cross)
 d$FID[8]<-NA_character_
 expect_error(run(candidate,d),'missing IDs cannot be treated as independent')
 d$FID<-'one_family'
 expect_error(run(candidate,d),'at least two distinct IDs')
})

test_that('explicit independent diagnostic and legacy binary limitation remain explicit',{
 d<-fixture(100);record_main(candidate);record_diagnostic(candidate)
 run(candidate,d,family_id=NULL,diagnostic_cluster=NULL)
 a<-candidate$diagnostic_calls[[1]]
 expect_true('cluster'%in%names(a));expect_null(a$cluster);expect_null(a$family_id)
 d$Y<-as.numeric(d$Y>0);d$FID<-rep(sprintf('F%02d',1:2),each=50)
 record_main(candidate);record_diagnostic(candidate)
 z<-run(candidate,d,binary_outcome=TRUE)
 expect_null(candidate$main_calls[[1]]$cluster)
 expect_match(attr(z,'settings')$family_note,'SEs treat rows as independent')
 expect_identical(candidate$diagnostic_calls[[1]]$cluster,'FID')
 expect_identical(candidate$diagnostic_calls[[1]]$outcome,'binary_liability')
})

test_that('diagnostic target age stays in input units when age is also a covariate',{
 d<-fixture(100);d<-d[rep(seq_len(nrow(d)),each=2),];d$age<-rep(c(8,10),100)
 record_main(candidate);record_diagnostic(candidate)
 run(candidate,d,design='longitudinal',covariates=c('sex','age'),age_center=NULL)
 a<-candidate$diagnostic_calls[[1]]
 expect_equal(a$target_age,9);expect_identical(a$dat$age,d$age)
 record_diagnostic(candidate)
 run(candidate,d,design='longitudinal',covariates=c('sex','age'),age_center=10)
 expect_equal(candidate$diagnostic_calls[[1]]$target_age,10)
})

test_that('actual unadjusted diagnostic runs despite covariates and primary failure',{
 rm('ravel_score_diagnostic',envir=candidate)
 d<-fixture();z<-run(candidate,d,covariates='sex')
 diag<-attr(z,'one_factor');row<-z[z$parameter%in%'one_factor_test',]
 expect_identical(diag$status,'ok',info=diag$reason)
 expect_equal(row$diagnostic_N,1500)
 expect_match(row$note,'Covariates were not used to compute the one-factor test',fixed=TRUE)
 expect_false(diag$workflow$primary_fit_converged)
 expect_equal(length(diag$selection$rows_used),1500)
})

test_that('actual longitudinal diagnostic selects one age per IID before missingness',{
 d<-fixture();d<-d[rep(seq_len(nrow(d)),each=2),];d$age<-rep(c(8,10),1500)
 # The missing target-age record is selected and dropped, never replaced by age8.
 d$Y[2]<-NA_real_
 z<-run(candidate,d,design='longitudinal',covariates='sex',age_center=10)
 diag<-attr(z,'one_factor');row<-z[z$parameter%in%'one_factor_test',]
 expect_equal(diag$selection$N_selected,1500)
 expect_equal(diag$selection$N_incomplete_selected,1)
 expect_equal(diag$N,1499)
 expect_true(all(d$age[diag$selection$rows_used]==10))
 expect_false('I00001'%in%diag$selection$individuals_used)
 expect_match(row$note,'Genetic interactions were not used',fixed=TRUE)
 expect_identical(diag$status,'ok',info=diag$reason)
})

test_that('parent phenotype wrapper forwards the focal adult ID and diagnostic cluster',{
 parent<-load_workflow('analysis/parental_deviation/ravel_parental_phenotypes.R')
 captured<-NULL
 engine<-new.env();engine$ravel_parental_deviation<-function(...){captured<<-list(...);stop('Controlled source stop')}
 parent$.ravel_parental_engine<-function()engine
 d<-fixture(100);d$adult_id<-paste0('A',seq_len(nrow(d)))
 suppressMessages(parent$ravel_parental_phenotypes(d,outcome='Y',pgs_columns=cols,sib_columns=NULL,
   subject_id='adult_id',diagnostic_cluster='FID',models=list(Equilibrium=list(fit_constraints='equilibrium')),
   outcome_models=list(own='mother'),include_sem_audit=FALSE))
 expect_identical(captured$subject_id,'adult_id');expect_identical(captured$diagnostic_cluster,'FID')
})

cat('Score workflow integration checks passed.\n')
