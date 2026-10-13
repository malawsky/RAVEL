score_signal_population <- function(parental_sibs = 'none', mode = 'equilibrium',
                                    n_gen = 3L, pi = .4, rho = .35,
                                    noise = c(.4, 1.3), beta = NULL,
                                    alternative = c('null','unequal','proportional'),
                                    departure = .10) {
  alternative <- match.arg(alternative)
  roles <- switch(parental_sibs, none=c('child','mother','father'),
    maternal=c('child','mother','father','mat_sib'),
    paternal=c('child','mother','father','pat_sib'),
    both=c('child','mother','father','mat_sib','pat_sib'))
  if (is.null(roles)) stop('Unknown sibling mode')
  V <- if (mode == 'equilibrium') 1/(1-rho) else 1
  if (mode != 'equilibrium' && n_gen > 1L) for (j in seq_len(n_gen-1L)) V <- .5*(1+rho)*V + .5
  vp <- 1+(1-pi)*(V-1); k <- (1-pi)*rho*V
  own <- vp-.5; across <- k*(1-.5/V); between <- k*(1-.5/V)^2
  child_parent <- .5*(vp+k); child_sib <- .5*(own+across)
  G <- matrix(c(.5*(vp+k)+.5, child_parent,child_parent,child_sib,child_sib,
    child_parent,vp,k,own,across, child_parent,k,vp,across,own,
    child_sib,own,across,vp,between, child_sib,across,own,between,vp),5,5)
  K <- matrix(c(1,.5,.5,.25,.25, .5,1,0,.5,0, .5,0,1,0,.5,
    .25,.5,0,1,0, .25,0,.5,0,1),5,5)
  ids <- match(roles,c('child','mother','father','mat_sib','pat_sib'))
  G <- G[ids,ids,drop=FALSE]; K <- K[ids,ids,drop=FALSE]
  nr <- length(roles)
  S <- rbind(cbind(G+noise[1]*K,G),cbind(G,G+noise[2]*K))
  L <- rbind(diag(sqrt(diag(G))),diag(sqrt(diag(G))))
  C <- cov2cor(G)
  if (is.null(beta)) beta <- c(.32,.13,.08,.05,-.04)[ids]
  c_null <- as.vector(L %*% C %*% beta)
  a <- numeric(nr); a[1] <- departure
  added <- if (alternative=='unequal') c(a,rep(0,nr)) else if(alternative=='proportional') as.vector(L %*% a) else rep(0,2*nr)
  cy <- c_null+added
  # Add enough outcome variance for the error-associated alternative; fixed
  # positive Schur margin avoids accidental boundary alternatives.
  vy <- sum(beta*(C %*% beta))+.8
  vy <- max(vy, as.numeric(crossprod(cy,solve(S,cy)))+.6)
  full <- rbind(cbind(S,cy),c(cy,vy))
  nm <- c(paste0('PGS1_',roles),paste0('PGS2_',roles),'Y')
  dimnames(full) <- list(nm,nm)
  list(covariance=full, score_covariance=S, phenotype_covariance=cy,
    score_loading=L, genetic_covariance=C, beta=beta, roles=roles,
    alternative=alternative, absorbed_beta=beta+solve(C,a),
    mode=mode, pi=pi,n_gen=n_gen, parental_sibs=parental_sibs)
}
score_signal_data <- function(pop, N=2000L, seed=82031L, empirical=FALSE) {
  set.seed(seed)
  out <- as.data.frame(MASS::mvrnorm(N,rep(0,ncol(pop$covariance)),pop$covariance,empirical=empirical))
  names(out) <- colnames(pop$covariance)
  out
}
score_signal_scientific_fit <- function(dat,pop,...) {
  RAVEL::ravel_fit(dat,fit_constraints=pop$mode,parental_sibs=pop$parental_sibs,
    pi=pop$pi,am_n_generations=pop$n_gen,...)
}

# These exact-moment fixtures independently construct the family covariance.
# They do not use the package's measurement or genetic constraint builders.
score_signal_cache <- new.env(parent=emptyenv())
score_signal_fixture <- function(roles='none',mode='equilibrium',seed=82031L,
                                 empirical=TRUE,...) {
  key <- paste(roles,mode,seed,empirical,sep=':')
  if (!exists(key,score_signal_cache,inherits=FALSE)) {
    pop <- score_signal_population(roles,mode=mode)
    dat <- score_signal_data(pop,seed=seed,empirical=empirical)
    value <- list(pop=pop,dat=dat,fit=score_signal_scientific_fit(dat,pop,...))
    assign(key,value,score_signal_cache)
  }
  get(key,score_signal_cache,inherits=FALSE)
}


test_that('constructive observable nesting holds for every supported family size', {
  for (roles in c('none','maternal','paternal','both')) {
    p <- score_signal_population(roles)
    S <- p$score_covariance; cy <- p$phenotype_covariance
    v <- p$covariance[nrow(p$covariance),nrow(p$covariance)]
    q <- solve(S,cy); ey <- v-as.numeric(crossprod(cy,q))
    expect_gt(ey,0)
    expect_equal(as.vector(S%*%q),cy,tolerance=1e-12)
    expect_equal(as.numeric(crossprod(q,S%*%q))+ey,v,tolerance=1e-12)
    expect_equal(as.vector(p$score_loading%*%p$genetic_covariance%*%p$beta),cy,
      tolerance=1e-12)
    # Arbitrary affine units preserve the normalized covariance restriction.
    z <- seq(.3,3,length.out=2*length(p$roles))
    SS <- S*outer(z,z); cc <- cy*z
    expect_equal(as.vector(SS%*%solve(SS,cc)),cc,tolerance=1e-12)
  }
})


test_that('proportional score-specific signal is a genuine diagnostic blind spot', {
  for (roles in c('none','maternal','paternal','both')) {
    p <- score_signal_population(roles,alternative='proportional',departure=.12)
    absorbed <- as.vector(p$score_loading%*%p$genetic_covariance%*%p$absorbed_beta)
    expect_equal(absorbed,p$phenotype_covariance,tolerance=1e-12)
    expect_gt(max(abs(p$beta-p$absorbed_beta)),.05)
    q <- score_signal_population(roles,alternative='unequal',departure=.12)
    n <- length(q$roles)
    c1 <- q$phenotype_covariance[seq_len(n)]/diag(q$score_loading[seq_len(n),])
    c2 <- q$phenotype_covariance[n+seq_len(n)]/diag(q$score_loading[n+seq_len(n),])
    expect_gt(max(abs(c1-c2)),.05)
  }
})


test_that('actual constrained fits preserve every score restriction and the role-specific df', {
  for (roles in c('none','maternal','paternal','both')) {
    x <- score_signal_fixture(roles)
    z <- ravel_test_score_signal(x$fit)
    expect_identical(z$status,'ok',info=paste(roles,z$reason))
    if (!identical(z$status,'ok')) next
    n <- length(x$pop$roles)
    expect_equal(z$df,n)
    expect_equal(z$diagnostics$rank_difference,n)
    expect_equal(z$diagnostics$H0$rank$covariance_rank,4*n+1)
    expect_equal(z$diagnostics$score$rank$covariance_rank,3*n)
    expect_equal(z$diagnostics$H1$rank$covariance_rank,5*n+1)
    expect_lt(z$diagnostics$nesting$maximum_scaled_moment_error,1e-9)
    expect_gt(z$diagnostics$nesting$alternative_outcome_residual,0)
    for (fitname in c('H0','score','H1')) {
      r <- z$diagnostics[[fitname]]$rank
      expect_true(r$identified)
      expect_equal(r$rank,r$tangent_dimension)
      expect_true(all(r$ranks_at_tolerances==r$rank))
    }
    pt <- lavaan::parameterTable(z$fits$H1)
    yp <- pt[pt$lhs=='Y' & pt$op=='~',]
    expect_setequal(yp$rhs,setdiff(names(x$dat),'Y'))
    expect_equal(nrow(yp),2*n)
    expect_true(all(yp$free>0))
    yc <- pt[pt$op=='~~' & xor(pt$lhs=='Y',pt$rhs=='Y'),]
    expect_true(all(yc$free==0 & yc$est==0))
    expect_false(any(pt$op==':=' & grepl('\\b(bD|bm|bf|bms|bfs)\\b',pt$rhs)))
    restrictions <- function(f) {
      p <- lavaan::parameterTable(f)
      sort(paste(p$lhs[p$op=='=='],p$rhs[p$op=='=='],sep=' == '))
    }
    expect_identical(restrictions(z$fits$H0),restrictions(z$fits$H1))
    expect_identical(restrictions(z$fits$H0),restrictions(z$fits$score))
    expect_equal(as.integer(lavaan::lavInspect(z$fits$H1,'nobs')),nrow(x$dat))
    expect_equal(z$raw_statistic,unname(2*diff(z$log_likelihood)),tolerance=1e-10)
    expect_match(z$test,'satorra.2000',fixed=TRUE)
    expect_equal(z$adjustment$jacobian,'delta')
  }
})


test_that('finite and estimated-state history restrictions retain the expected tangent dimension', {
  for (mode in c('n_gen','estimated_state')) for (roles in c('none','maternal','paternal','both')) {
    x <- score_signal_fixture(roles,mode)
    z <- ravel_test_score_signal(x$fit)
    expect_identical(z$status,'ok',info=paste(mode,roles,z$reason))
    if (!identical(z$status,'ok')) next
    n <- length(x$pop$roles)
    expect_equal(z$df,n)
    expect_equal(z$diagnostics$rank_difference,n)
    expect_identical(z$diagnostics$conditional_on_pi,mode=='n_gen')
    expect_equal(z$external_pi,.4)
    expect_lt(z$diagnostics$nesting$maximum_scaled_moment_error,1e-9)
  }
})


test_that('diagnostic statistics are invariant to positive score-unit changes and post-fit correction', {
  x <- score_signal_fixture('both',empirical=FALSE)
  z <- ravel_test_score_signal(x$fit)
  d <- x$dat
  scale <- c(.3,2,1.4,3,.8,1.6,.25,4,1.7,.65)
  d[,seq_len(10)] <- sweep(d[,seq_len(10)],2,scale,'*')
  d[,seq_len(10)] <- sweep(d[,seq_len(10)],2,seq(-3,3,length.out=10),'+')
  f <- score_signal_scientific_fit(d,x$pop)
  zz <- ravel_test_score_signal(f)
  expect_identical(z$status,'ok',info=z$reason)
  expect_identical(zz$status,'ok',info=zz$reason)
  expect_equal(z$raw_statistic,zz$raw_statistic,tolerance=1e-5)
  expect_equal(z$statistic,zz$statistic,tolerance=2e-4)
  expect_equal(z$p_value,zz$p_value,tolerance=1e-5)
  changed <- x$fit
  changed$pi_correction <- ravel_correct_pi(changed,pi=.45,pi_se=.02)
  changed$beta_total_assumed <- changed$beta_total_assumed+10
  after <- ravel_test_score_signal(changed)
  expect_identical(z$status,after$status)
  expect_equal(z$statistic,after$statistic,tolerance=1e-10)
  expect_equal(z$p_value,after$p_value,tolerance=1e-10)
  expect_equal(z$external_pi,after$external_pi)
  expect_identical(x$fit$beta_vis,changed$beta_vis)
})


# Independent constrained-information calculation of the robust difference.
score_signal_reference_lrt <- function(f0,f1,df) {
  nullspace <- function(A) {
    if(!nrow(A)) return(diag(ncol(A)))
    z <- svd(A,nu=0,nv=ncol(A))
    rank <- sum(z$d>1e-8*max(z$d))
    z$v[,seq.int(rank+1L,ncol(A)),drop=FALSE]
  }
  tangent <- lapply(list(f0,f1),function(f) {
    A <- f@Model@con.jac
    nullspace(A[attr(A,'ceq.idx'),,drop=FALSE])
  })
  J <- lapply(list(f0,f1),function(f)as.matrix(lavaan::lavInspect(f,'delta')))
  H <- MASS::ginv(J[[2]]%*%tangent[[2]])%*%(J[[1]]%*%tangent[[1]])
  A <- t(nullspace(t(H)))%*%t(tangent[[2]])
  P <- lavaan::lavTech(f1,'information')
  N <- tangent[[2]]
  Pn <- N%*%solve(crossprod(N,P%*%N))%*%t(N)
  Q <- Pn%*%t(A)%*%solve(A%*%Pn%*%t(A))%*%A%*%Pn
  W <- lavaan::lavTech(f1,'WLS.V')[[1]]
  Gamma <- lavaan::lavTech(f1,'Gamma')[[1]]
  UG <- W%*%J[[2]]%*%Q%*%t(J[[2]])%*%W%*%Gamma
  scale <- sqrt(df/sum(diag(UG%*%UG)))
  shift <- df-scale*sum(diag(UG))
  raw <- 2*(lavaan::fitMeasures(f1,'logl')-lavaan::fitMeasures(f0,'logl'))
  c(statistic=unname(raw*scale+shift),scale=scale,shift=shift)
}


test_that('Gaussian ML and robust MLR use their actual nested comparison', {
  x <- score_signal_fixture('none',empirical=FALSE)
  before <- serialize(x$fit,NULL)
  robust <- ravel_test_score_signal(x$fit)
  expect_identical(serialize(x$fit,NULL),before)
  f <- x$fit
  f$fit <- lavaan::sem(f$model,data=x$dat,std.lv=TRUE,fixed.x=FALSE,
    conditional.x=FALSE,estimator='ML',missing='fiml')
  gaussian <- ravel_test_score_signal(f)
  expect_identical(robust$status,'ok',info=robust$reason)
  expect_identical(gaussian$status,'ok',info=gaussian$reason)
  expect_equal(gaussian$statistic,gaussian$raw_statistic,tolerance=1e-10)
  expect_equal(gaussian$p_value,pchisq(gaussian$statistic,gaussian$df,lower.tail=FALSE))
  expect_equal(robust$raw_statistic,gaussian$raw_statistic,tolerance=1e-5)
  copies <- lapply(c('H0','H1'),function(nm)
    .ravel_score_signal_lrt_fit(robust$fits[[nm]],robust$diagnostics[[nm]]$rank))
  args <- c(copies,list(method='satorra.2000'))
  available <- names(formals(lavaan::lavTestLRT))
  args[[if('a_method' %in% available)'a_method' else 'A.method']] <- 'delta'
  args[[if('scaled_shifted' %in% available)'scaled_shifted' else 'scaled.shifted']] <- TRUE
  table <- do.call(lavaan::lavTestLRT,args)
  row <- which(is.finite(table[['Df diff']]))
  expect_equal(robust$statistic,table[['Chisq diff']][row],tolerance=1e-10)
  expect_equal(robust$p_value,table[['Pr(>Chisq)']][row],tolerance=1e-10)
  expect_true(is.finite(robust$adjustment$scale))
  independent <- score_signal_reference_lrt(robust$fits$H0,robust$fits$H1,robust$df)
  expect_equal(robust$statistic,unname(independent['statistic']),tolerance=1e-7)
  expect_equal(robust$adjustment$scale,unname(independent['scale']),tolerance=1e-7)
  expect_equal(robust$adjustment$shift,unname(independent['shift']),tolerance=1e-7)
})


test_that('unsupported designs and data mismatches return missing inference', {
  x <- score_signal_fixture()
  check <- function(f,status,data=NULL) {
    z <- ravel_test_score_signal(f,data=data)
    expect_identical(z$status,status,info=z$reason)
    expect_true(all(is.na(c(z$statistic,z$df,z$p_value))))
    expect_true(nzchar(z$reason))
  }
  f <- x$fit; f$fit <- NULL; check(f,'needs_fit')
  f <- x$fit; f$outcome <- 'binary_liability'; check(f,'unsupported_outcome')
  f <- x$fit; f$fixed_effects <- 'age'; check(f,'unsupported_covariates')
  f <- x$fit; f$imputed <- TRUE; check(f,'unsupported_imputed')
  f <- x$fit; f$score_signal_input$dependent_families <- TRUE; check(f,'unsupported_design')
  f <- x$fit; f$score_signal_input$incomplete <- TRUE; check(f,'unsupported_missingness')
  f <- x$fit; f$constraint_spec$mode <- 'none'; f$fixed_crossperson <- f$fixed_crossperson_corr <- NULL
  check(f,'unsupported_free_measurement')
  d <- x$dat; d$Y[1] <- d$Y[1]+.1; check(x$fit,'data_mismatch',d)
  d <- x$dat; d$FID <- rep(seq_len(nrow(d)/2),each=2); check(x$fit,'unsupported_design',d)
  d <- x$dat; attr(d,'imputed') <- TRUE; check(x$fit,'unsupported_imputed',d)
  check(x$fit,'data_mismatch',x$dat[-1,])
  expect_identical(ravel_test_score_signal(x$fit,data=x$dat)$status,'ok')
  f <- score_signal_scientific_fit(x$dat,x$pop,fix_parent_paths=TRUE)
  check(f,'unsupported_fixed_paths')
})


test_that('singular and boundary null models do not receive ordinary chi-square inference', {
  for (scenario in c('singular','boundary')) {
    p <- score_signal_population(rho=if(scenario=='singular') 0 else .35,
      noise=if(scenario=='boundary') c(0,1.3) else c(.4,1.3))
    d <- score_signal_data(p,empirical=TRUE)
    f <- suppressWarnings(score_signal_scientific_fit(d,p))
    z <- ravel_test_score_signal(f)
    expect_false(identical(z$status,'ok'),info=scenario)
    expect_true(all(is.na(c(z$statistic,z$df,z$p_value))))
    expect_match(z$reason,'singular|weak|boundary|converge|positive|fit',ignore.case=TRUE)
  }
})


test_that('the fitted test detects unequal signal but cannot detect proportional contamination', {
  for (roles in c('none','both')) {
    p <- score_signal_population(roles,alternative='unequal',departure=.18)
    d <- score_signal_data(p,empirical=TRUE)
    f <- score_signal_scientific_fit(d,p)
    z <- ravel_test_score_signal(f)
    expect_identical(z$status,'ok',info=paste(roles,z$reason))
    expect_lt(z$p_value,.001)
    p <- score_signal_population(roles,alternative='proportional',departure=.12)
    d <- score_signal_data(p,empirical=TRUE)
    f <- score_signal_scientific_fit(d,p)
    z <- ravel_test_score_signal(f)
    expect_identical(z$status,'ok',info=paste(roles,z$reason))
    expect_gt(z$p_value,.9)
    pe <- lavaan::parameterEstimates(f$fit)
    effect <- pe$est[match(c('bD','bm','bf',if(roles=='both')c('bms','bfs')),pe$label)]
    expect_equal(effect,as.vector(p$absorbed_beta)*sqrt((nrow(d)-1)/nrow(d)),tolerance=1e-4)
    expect_gt(max(abs(effect-p$beta)),.05)
  }
})


test_that('the fixed-residual benchmark supports all roles and preserves covariance units', {
  for (roles in c('none','maternal','paternal','both')) {
    x <- score_signal_fixture(roles)
    sr <- switch(roles,none=character(0),maternal='ms',paternal='fs',both=c('ms','fs'))
    corr <- c(e1_mf=0,e1_mc=.5,e1_fc=.5,e2_mf=0,e2_mc=.5,e2_fc=.5,
      .parental_sib_residual_benchmarks(sr))
    initial <- .ravel_score_signal_capture(
      ravel_fit(x$dat,parental_sibs=roles,fixed_crossperson_corr=corr))
    expect_null(initial$error)
    expect_true(all(initial$warnings %in% c('NaNs produced','NA/NaN function evaluation')))
    f <- initial$value
    z <- ravel_test_score_signal(f)
    expect_identical(z$status,'ok',info=paste(roles,z$reason))
    if (!identical(z$status,'ok')) next
    n <- length(x$pop$roles)
    expect_equal(z$df,n)
    expect_equal(z$diagnostics$score$rank$covariance_rank,4*n+n*(n-1)/2)
  }
  x <- score_signal_fixture('none',empirical=FALSE)
  fixed <- c(e1_mf=0,e1_mc=.2,e1_fc=.2,e2_mf=0,e2_mc=.65,e2_fc=.65)
  f <- ravel_fit(x$dat,fixed_crossperson=fixed)
  z <- ravel_test_score_signal(f)
  d <- x$dat
  u <- setNames(c(.3,2,1.4,3,.8,1.6),names(d)[1:6])
  d[,1:6] <- sweep(d[,1:6],2,u,'*')
  scaled <- fixed
  scaled['e1_mc'] <- fixed['e1_mc']*u['PGS1_mother']*u['PGS1_child']
  scaled['e1_fc'] <- fixed['e1_fc']*u['PGS1_father']*u['PGS1_child']
  scaled['e2_mc'] <- fixed['e2_mc']*u['PGS2_mother']*u['PGS2_child']
  scaled['e2_fc'] <- fixed['e2_fc']*u['PGS2_father']*u['PGS2_child']
  g <- ravel_fit(d,fixed_crossperson=scaled)
  zz <- ravel_test_score_signal(g)
  expect_identical(z$status,'ok',info=z$reason)
  expect_identical(zz$status,'ok',info=zz$reason)
  expect_equal(z$raw_statistic,zz$raw_statistic,tolerance=2e-5)
  expect_equal(z$statistic,zz$statistic,tolerance=2e-5)
  expect_equal(z$p_value,zz$p_value,tolerance=2e-6)
})
