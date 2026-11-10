# Synthetic numerical validation; lme4/numDeriv are test-only oracles.
source("R/ravel_pgs_comparators.R")
set.seed(6271)
outdir <- "dev/pgs-comparator-validation"; dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
checks <- list()
check <- function(label,condition,detail="") {
  checks[[length(checks)+1L]] <<- data.frame(check=label,passed=isTRUE(condition),detail=as.character(detail))
  if(!isTRUE(condition))stop(label,": ",detail)
}
near <- function(label,a,b,tol=1e-6)check(label,max(abs(a-b))<tol,max(abs(a-b)))
roles <- c("child","mother","father"); pgs <- setNames(paste0(rep(roles,2),rep(1:2,each=3)),paste0(rep(roles,2),rep(1:2,each=3)))
make_data <- function(nf=220,nc=2,nt=1) {
  d<-expand.grid(visit=seq_len(nt),kid=seq_len(nc),family=seq_len(nf)); d$M_ID<-factor(d$family)
  d$IID<-factor((d$family-1)*nc+d$kid);d$FID<-factor(ceiling(d$family/2))
  d$age<-if(nt==1)rnorm(nrow(d),8,1)else 6+2*d$visit
  d$sex<-rep(rbinom(nf*nc,1,.5),each=nt)
  for(k in 1:2) {
    d[[paste0("mother",k)]]<-rnorm(nf)[d$family]
    d[[paste0("father",k)]]<-rnorm(nf)[d$family]
    d[[paste0("child",k)]]<-rnorm(nf*nc)[as.integer(d$IID)]
  }
  d
}
design <- function(d,k=1,interactions=character()) {
  Z<-sapply(roles,function(role){id<-if(role=="child")d$IID else d$M_ID;x<-d[[paste0(role,k)]];v<-x[!duplicated(id)];(x-mean(v))/sd(v)})
  colnames(Z)<-c("beta_D","beta_M","beta_F")
  X<-cbind(intercept=1,Z,fixed_age=d$age-8,fixed_sex=d$sex-.5)
  for(mod in interactions)X<-cbind(X,setNames(as.data.frame(Z*(d[[mod]]-if(mod=="age")8 else .5)),paste0(colnames(Z),":",mod)))
  as.matrix(X)
}
run <- function(d,...) .ravel_pgs_comparators(d,"Y",pgs,covariates=c("age","sex"),interaction_center=c(age=8,sex=.5),...)
d<-make_data(nf=160);X<-design(d);d$Y<-drop(X%*%c(.1,.4,.2,.1,.05,.1))+rnorm(nrow(d))
a<-run(d)
for(k in 1:2) {
  f<-a$fits[[k]];Xk<-design(d,k);ols<-lm.fit(Xk,d$Y);B<-solve(crossprod(Xk));U<-rowsum(Xk*ols$residuals,d$FID,reorder=FALSE)
  near(paste("OLS coefficients",k),f$parameters,ols$coefficients)
  near(paste("OLS CR0 covariance",k),f$coefficient_vcov,B%*%crossprod(U)%*%B)
  expectedll<-sum(dnorm(ols$residuals,sd=sqrt(mean(ols$residuals^2)),log=TRUE))
  near(paste("OLS log likelihood",k),f$fit_statistics[["logLik"]],expectedll)
  tab<-a$table[a$table$regression==paste0("PGS",k),]
  near(paste("Parental contrast",k),tab$estimate[tab$parameter=="delta_M_mid"],f$parameters[["beta_M"]]-f$parameters[["beta_F"]])
}
small<-d;small[unname(pgs)]<-lapply(small[unname(pgs)],function(x)x*1e-9)
as<-run(small)
near("Raw PGS scale invariance coefficients",a$table$estimate,as$table$estimate)
near("Raw PGS scale invariance SE",a$table$se,as$table$se)
for(k in 1:2){at<-a$pgs_moments$regression==paste0("PGS",k)&a$pgs_moments$block=="observed_PGS_correlation";near(paste("Tiny raw PGS correlation",k),a$pgs_moments$estimate[at],as$pgs_moments$estimate[at]);near(paste("Tiny raw PGS correlation SE",k),a$pgs_moments$se[at],as$pgs_moments$se[at])}
# Family ML against independent established mixed-model implementation.
d$Y<-drop(X%*%c(.1,.4,.2,.1,.05,.1))+rnorm(160,sd=.8)[d$family]+rnorm(nrow(d),sd=.7)
m<-run(d,family_random=TRUE);f<-m$fits$PGS1
oracledata<-data.frame(Y=d$Y,X[,-1],family=d$M_ID,check.names=FALSE)
formula<-as.formula(paste("Y~",paste(colnames(X)[-1],collapse="+"),"+(1|family)"))
oracle<-lme4::lmer(formula,oracledata,REML=FALSE,control=lme4::lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
near("Family ML coefficients",f$parameters,lme4::fixef(oracle),1e-5)
near("Family ML likelihood",f$fit_statistics[["logLik"]],as.numeric(logLik(oracle)),1e-6)
check("Family ML inference status",f$status=="ok",f$status)
# Full finite-difference likelihood Hessian and family score sandwich oracle.
par<-f$internal_parameters;p<-length(f$par);q<-ncol(X);Xs<-sweep(X,2,f$design_scaling,"/")
familyrows<-split(seq_len(nrow(d)),d$M_ID)
nll<-function(theta,which=seq_along(familyrows)) {
  sum(vapply(familyrows[which],function(at){S<-diag(exp(2*theta[1]),length(at))+exp(2*theta[2]);r<-d$Y[at]-drop(Xs[at,,drop=FALSE]%*%theta[p+seq_len(q)]);.5*(as.numeric(determinant(S,log=TRUE)$modulus)+sum(r*solve(S,r))+length(at)*log(2*pi))},numeric(1)))
}
H<-numDeriv::hessian(nll,par);U<-t(vapply(seq_along(familyrows),function(i)numDeriv::grad(function(z)-nll(z,i),par),numeric(length(par))))
U<-rowsum(U,as.character(d$FID[vapply(familyrows,`[`,integer(1),1)]),reorder=FALSE)
B<-solve(H);V<-B%*%crossprod(U)%*%B
near("Full joint sandwich including variance nuisance",f$vcov,V,2e-6)
# Longitudinal intercept/slope, family effect, and PGS age/sex interactions.
dl<-make_data(nf=180,nc=2,nt=3);XL<-design(dl,interactions=c("age","sex"));beta<-c(.1,.3,.2,.1,.1,.1,.04,-.02,.03,.1,-.1,.05)
uf<-rnorm(180,.0,.65);ui<-rnorm(360,0,.5);us<-.1*ui+rnorm(360,0,.15)
dl$Y<-drop(XL%*%beta)+uf[dl$family]+ui[as.integer(dl$IID)]+us[as.integer(dl$IID)]*(dl$age-8)+rnorm(nrow(dl),0,.7)
ml<-run(dl,family_random=TRUE,individual_random=TRUE,random_slope="age",interactions=c("age","sex"));fl<-ml$fits$PGS1
oldata<-data.frame(Y=dl$Y,XL[,-1],fam=dl$M_ID,iid=dl$IID,age_c=dl$age-8,check.names=FALSE)
lf<-as.formula(paste("Y~",paste(sprintf("`%s`",colnames(XL)[-1]),collapse="+"),"+(1|fam)+(1+age_c|iid)"))
lo<-lme4::lmer(lf,oldata,REML=FALSE,control=lme4::lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=200000)))
near("Longitudinal ML coefficients and interactions",fl$parameters,lme4::fixef(lo),2e-5)
near("Longitudinal ML likelihood",fl$fit_statistics[["logLik"]],as.numeric(logLik(lo)),2e-5)
check("Longitudinal inference status",fl$status=="ok",fl$status)
sv<-as.data.frame(lme4::VarCorr(lo))
near("Longitudinal family variance",fl$variance_components$estimate[fl$variance_components$parameter=="family_variance"],sv$vcov[sv$grp=="fam"],2e-5)
near("Longitudinal slope variance",fl$variance_components$estimate[fl$variance_components$parameter=="slope_variance"],sv$vcov[sv$grp=="iid"&!is.na(sv$var1)&sv$var1=="age_c"&is.na(sv$var2)],2e-5)
# Conventional binary comparator uses probit GLM, despite requested random effects.
d$Y<-rbinom(nrow(d),1,pnorm(drop(X%*%c(-.1,.4,.2,.1,.05,.1))))
b<-run(d,binary_outcome=TRUE,family_random=TRUE,individual_random=TRUE)
bf<-b$fits$PGS1;go<-glm.fit(X,d$Y,family=binomial("probit"))
near("Binary probit coefficients",bf$parameters,go$coefficients)
near("Binary probit AIC",bf$fit_statistics[["AIC"]],go$aic)
check("Binary actual random effects metadata",!bf$family_random&&!bf$individual_random&&is.null(bf$random_slope))
probitnll<-function(z,idx=seq_len(nrow(X))) -sum(pnorm((2*d$Y[idx]-1)*drop(X[idx,,drop=FALSE]%*%z),log.p=TRUE))
Hb<-numDeriv::hessian(probitnll,bf$parameters);ub<-t(vapply(split(seq_len(nrow(d)),d$FID),function(i)numDeriv::grad(function(z)-probitnll(z,i),bf$parameters),numeric(ncol(X))))
Bb<-solve(Hb);near("Binary observed Hessian FID sandwich",bf$coefficient_vcov,Bb%*%crossprod(ub)%*%Bb,1e-6)
check("Binary finite robust inference",all(is.finite(b$table$se))&&all(b$table$status=="ok"))
# One SE cluster: points survive and uncertainty is marked unavailable.
d$FID<-factor(1);one<-run(d,binary_outcome=TRUE)
check("Single cluster retains points",all(is.finite(one$table$estimate)))
check("Single cluster suppresses invalid SE",all(is.na(one$table$se))&&all(one$table$status=="inference_unavailable"))
# Sparse sibling design, similar to a cohort with only five percent sibling pairs.
ds<-make_data(nf=500);ds<-ds[ds$kid==1|ds$family<=25,];XS<-design(ds)
ds$Y<-drop(XS%*%c(.1,.3,.2,.1,.1,.1))+rnorm(500,sd=1.1)[ds$family]+rnorm(nrow(ds),sd=.6)
sp<-run(ds,family_random=TRUE);sf<-sp$fits$PGS1
sdata<-data.frame(Y=ds$Y,XS[,-1],family=ds$M_ID,check.names=FALSE)
so<-lme4::lmer(formula,sdata,REML=FALSE,control=lme4::lmerControl(optimizer="bobyqa"))
near("Sparse sibling ML coefficients",sf$parameters,lme4::fixef(so),2e-5)
near("Sparse sibling ML likelihood",sf$fit_statistics[["logLik"]],as.numeric(logLik(so)),2e-5)
check("Sparse sibling full sandwich finite",sf$status=="ok"&&all(is.finite(sf$coefficient_vcov)),sf$status)
# Actual zero-variance boundary: retain point estimates but avoid regular Wald SEs.
dz<-make_data(nf=120);zv<-rnorm(120);dz$Y<-zv[dz$family]*ifelse(dz$kid==1,1,-1)
zero<-run(dz,family_random=TRUE)
check("Boundary retains coefficients",all(is.finite(zero$table$estimate)))
check("Boundary reports status and no regular SE",all(zero$table$status=="variance_boundary")&&all(is.na(zero$table$se)),paste(unique(zero$table$status),collapse=","))
# Moment inference independently from the full means/raw-products delta method.
r<-d[unname(pgs[1:3])];names(r)<-roles;d$FID<-factor(ceiling(d$family/2))
mm<-.ravel_pgs_moments(r,d$IID,d$FID,"PGS1")
rr<-as.matrix(r);ii<-which(lower.tri(matrix(0,3,3),diag=TRUE),arr.ind=TRUE)
raw<-cbind(rr,sapply(seq_len(nrow(ii)),function(j)rr[,ii[j,1]]*rr[,ii[j,2]]));meanraw<-colMeans(raw)
um<-rowsum(sweep(raw,2,meanraw,"-"),d$FID,reorder=FALSE);vm<-crossprod(um)/nrow(raw)^2
map<-function(t){S<-matrix(0,3,3);S[ii]<-t[4:9];S[cbind(ii[,2],ii[,1])]<-t[4:9];S<-S-tcrossprod(t[1:3]);c(S[ii],cov2cor(S)[ii])}
jm<-numDeriv::jacobian(map,meanraw)
near("PGS moment covariance/correlation estimates",mm$estimate,map(meanraw))
near("PGS moment FID delta SE",mm$se,sqrt(pmax(0,diag(jm%*%vm%*%t(jm)))),1e-7)
# Force separation; a numerical solution must not masquerade as regular inference.
d$Y<-as.numeric(X[,"beta_D"]>0);sep<-run(d,binary_outcome=TRUE)
check("Separation marks nonregular fit",sep$fits$PGS1$status%in%c("separation","not_converged"))
check("Separation withholds regular SE",all(is.na(sep$table$se[sep$table$regression=="PGS1"])))
write.csv(do.call(rbind,checks),file.path(outdir,"checks.csv"),row.names=FALSE)
saveRDS(list(ordinary=a,family=m,longitudinal=ml,binary=b),"/tmp/ravel-pgs-comparator-validation-fits.rds")
writeLines(capture.output(sessionInfo()),file.path(outdir,"session.txt"))
cat(length(checks),"checks passed\n")
