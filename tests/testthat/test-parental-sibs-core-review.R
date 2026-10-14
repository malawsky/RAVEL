# Independent generative construction: choose parents, then conditional
# grandparental midparents and independent segregation innovations.
review_sib_population <- function(pi, rho, V, roles = c('ms','fs')) {
  w <- c(sqrt(1-pi),sqrt(pi)); I <- diag(2)
  S <- I + (V-1)*tcrossprod(w); A <- S - I/2
  K <- rho*V*tcrossprod(w)
  H <- A %*% solve(S); Q <- A - H %*% A
  # primitive vector: mother, father, lineage innovations, sibling segregation,
  # and child segregation. Only parent blocks are correlated a priori.
  P <- matrix(0,14,14)
  P[1:2,1:2] <- P[3:4,3:4] <- S
  P[1:2,3:4] <- P[3:4,1:2] <- K
  P[5:6,5:6] <- P[7:8,7:8] <- Q
  P[9:10,9:10] <- P[11:12,11:12] <- P[13:14,13:14] <- I/2
  W <- matrix(0,14,14)
  # output vector: child, mother, father, maternal sibling, paternal sibling,
  # maternal lineage mean, paternal lineage mean, each with two components.
  W[1:2,1:2] <- W[1:2,3:4] <- I/2
  W[1:2,13:14] <- I
  W[3:4,1:2] <- W[5:6,3:4] <- I
  W[11:12,1:2] <- W[13:14,3:4] <- H
  W[11:12,5:6] <- W[13:14,7:8] <- I
  W[7:8,] <- W[11:12,]; W[7:8,9:10] <- I
  W[9:10,] <- W[13:14,]; W[9:10,11:12] <- I
  SZ <- W %*% P %*% t(W)
  ww <- kronecker(diag(7),matrix(w,1))
  ST <- ww %*% SZ %*% t(ww)
  # Both lineages affect Y; true specific maternal=.15, paternal=-.07.
  by <- c(.3,.15,-.07,0,0,.35,-.22)
  cy <- as.vector(ww %*% SZ %*% t(ww) %*% by)
  gy <- as.vector(SZ[seq(1,13,2),] %*% t(ww) %*% by)
  sel <- c(1:3, match(roles,c('ms','fs'))+3)
  sg <- SZ[seq(1,13,2)[sel],seq(1,13,2)[sel]]
  sdg <- sqrt(diag(sg)); R <- cov2cor(sg)
  bv <- solve(R,gy[sel]/sdg)
  bt <- solve(ST[sel,sel],cy[sel])
  names(bv) <- names(bt) <- c('bD','bm','bf',c(ms='bms',fs='bfs')[roles])
  dimnames(R) <- list(c('c','m','f',roles),c('c','m','f',roles))
  cv <- c(if ('ms' %in% roles) c(kap_m=R['m','ms'],chi_ms=R['f','ms'],eta_ms=R['c','ms']),
          if ('fs' %in% roles) c(kap_f=R['f','fs'],chi_fs=R['m','fs'],eta_fs=R['c','fs']),
          if (length(roles)==2) c(zeta=R['ms','fs']))
  list(bv=bv, bt=bt, R=R, varT=diag(ST)[sel],
       cyG=gy[sel]/sdg, varY=as.numeric(t(by)%*%ST%*%by)+.7,
       growth=diag(sg)[1]/diag(sg)[2],
       sib=list(roles=roles,b=bv[-(1:3)],cov=cv),
       GH=SZ[seq(1,13,2)[sel],seq(2,14,2)[sel]]/
         outer(sdg,sqrt(diag(SZ)[seq(2,14,2)[sel]])))
}


review_sib_scores <- function(N = 3000L, seed = 910L) {
  p <- review_sib_population(.4, .3, 1.15)
  # Two independent error processes with pedigree relatedness. Reliability
  # differs between score versions; observed units differ across family roles.
  kin <- matrix(c(1,.5,.5,.25,.25, .5,1,0,.5,0,
    .5,0,1,0,.5, .25,.5,0,1,0, .25,0,.5,0,1), 5, 5)
  loading <- c(rep(sqrt(.7),5),rep(sqrt(.4),5))
  factor_map <- rbind(diag(5),diag(5)) * loading
  score_cov <- factor_map %*% p$R %*% t(factor_map)
  score_cov[1:5,1:5] <- score_cov[1:5,1:5] + .3*kin
  score_cov[6:10,6:10] <- score_cov[6:10,6:10] + .6*kin
  cy <- as.vector(factor_map %*% p$R %*% p$bv)
  full <- rbind(cbind(score_cov,cy),c(cy,1+sum(p$bv*(p$R %*% p$bv))))
  set.seed(seed)
  d <- as.data.frame(MASS::mvrnorm(N,rep(0,11),full,empirical=TRUE))
  names(d) <- c(paste0('PGS1_',c('child','mother','father','mat_sib','pat_sib')),
    paste0('PGS2_',c('child','mother','father','mat_sib','pat_sib')),'Y')
  d
}

test_that("independent lineage construction recovers parent-specific total effects", {
  scenarios <- list(c(pi=.4,rho=.3,V=1.15,n=2),c(pi=.8,rho=-.2,V=.9,n=2),
    c(pi=.2,rho=.2,V=1.25,n=0),c(pi=.4,rho=.2,V=1.7,n=-1))
  for (s in scenarios) for (roles in list('ms','fs',c('ms','fs'))) {
    mode <- if(s[['n']]>0) 'n_gen' else if(s[['n']]==0) 'equilibrium' else 'estimated_state'
    p <- review_sib_population(s[['pi']],s[['rho']],s[['V']],roles)
    z <- .compute_scale_estimates(p$bv[1:3],p$R[1,2],p$R[1,3],p$R[2,3],
      s[['pi']],mode,max(1,s[['n']]),sib=p$sib)
    expect_equal(unname(z$sib_total$Sigma_xh),unname(p$GH),tolerance=1e-9)
    expect_equal(unname(z$beta_total_assumed),unname(p$bt[1:3]),tolerance=1e-9)
    if('ms' %in% roles) expect_equal(z$sib_total$raw[['delta_m']],.15,tolerance=1e-9)
    if('fs' %in% roles) expect_equal(z$sib_total$raw[['delta_f']],-.07,tolerance=1e-9)
  }
})

test_that("an impossible joint hidden covariance is not reported as a correction", {
  sib <- list(roles='ms',b=c(bms=.05),
    cov=c(kap_m=.999,chi_ms=.18,eta_ms=.565))
  # Visible covariance is positive definite, but almost-identical visible
  # mother/sibling liabilities cannot have the assumed different hidden links.
  S <- .sib_latent_matrix(.565,.565,.18,sib)
  expect_gt(min(eigen(S,symmetric=TRUE)$values),0)
  expect_warning(z <- .compute_scale_estimates(c(.3,.1,.1),.565,.565,.18,.4,
    'n_gen',1L,sib=sib),'parental-sibling moments are incompatible')
  expect_true(all(is.na(z$beta_total_assumed)))
  expect_true(all(is.na(z$sib_total$std)))
})

test_that("parent-specific contrasts and their SEs survive independent score rescaling", {
  d <- review_sib_scores()
  a <- ravel_fit(d,pi=.4,parental_sibs='both',am_n_generations=2L)
  d$PGS1_mat_sib <- 7*d$PGS1_mat_sib + 4
  d$PGS2_mat_sib <- .3*d$PGS2_mat_sib - 2
  b <- ravel_fit(d,pi=.4,parental_sibs='both',am_n_generations=2L)
  expect_true(a$converged && b$converged)
  expect_equal(a$parental_sibs$corrected$raw[['delta_m']],.15,tolerance=1e-4)
  expect_equal(a$parental_sibs$corrected$raw,b$parental_sibs$corrected$raw,tolerance=1e-4)
  expect_equal(a$lavaan_inference$sib_total_raw$se,b$lavaan_inference$sib_total_raw$se,tolerance=1e-4)
  # Parent-minus-sibling inference needs their covariance, not independent SEs.
  V <- lavaan::vcov(a$fit)
  se <- sqrt(V['bm','bm']+V['bms','bms']-2*V['bm','bms'])
  expect_equal(a$parental_sibs$lavaan_table$se[a$parental_sibs$lavaan_table$label=='delta_m'],
    unname(se),tolerance=1e-9)
})

test_that("binary sibling fits report the retained sample and covariate models use all scores", {
  d <- review_sib_scores()
  d$Z <- .2*d$PGS1_pat_sib + sin(seq_len(nrow(d)))
  d$Y <- as.numeric(d$Y+.1*d$Z > 0)
  d[1:120,c('PGS1_mat_sib','PGS2_mat_sib')] <- NA
  d[121:200,c('PGS1_pat_sib','PGS2_pat_sib')] <- NA
  expect_warning(fit <- ravel_fit(d,parental_sibs='both',outcome='binary_liability',
    fixed_effects='Z'),'listwise deletion')
  expect_true(fit$converged)
  expect_equal(fit$N,nrow(d)-200)
  expect_equal(as.numeric(lavaan::lavInspect(fit$fit,'nobs')),fit$N)
  expect_equal(unname(fit$parental_sibs$n_with_sib),rep(fit$N,2))
  expect_true(all(c('PGS1_mat_sib','PGS2_pat_sib','Z') %in% lavaan::lavNames(fit$fit,'ov')))
})

test_that("nonfinite sibling scores are rejected before fitting", {
  d <- review_sib_scores(1000)
  d$PGS1_mat_sib[1] <- Inf
  expect_error(ravel_fit(d,parental_sibs='maternal'),'finite or missing')
})


test_that("explicit sibling residual restrictions use correlation or covariance units", {
  d <- review_sib_scores()
  benchmarks <- .parental_sib_residual_benchmarks(c('ms','fs'))
  fixed <- ravel_fit(d,parental_sibs='both',fixed_crossperson_corr=benchmarks)
  expect_true(fixed$converged)
  expect_equal(fixed$parental_sibs$resid_corr[names(benchmarks)],benchmarks,tolerance=1e-6)
  units <- ravel_fit(d,parental_sibs='maternal',fixed_crossperson=c(e1_m_ms=.12))
  expect_true(units$converged)
  expect_equal(units$parental_sibs$resid_cov[['e1_m_ms']],.12,tolerance=1e-7)
  # Re-expressing one score changes the covariance but not its fixed correlation.
  d$PGS1_mat_sib <- 3*d$PGS1_mat_sib
  scaled <- ravel_fit(d,parental_sibs='both',fixed_crossperson_corr=benchmarks)
  expect_equal(scaled$parental_sibs$resid_corr[names(benchmarks)],benchmarks,tolerance=1e-6)
  expect_equal(scaled$parental_sibs$resid_cov[['e1_m_ms']],
    3*fixed$parental_sibs$resid_cov[['e1_m_ms']],tolerance=1e-5)
  expect_error(ravel_fit(d,parental_sibs='maternal',fixed_crossperson_corr=c(e1_f_fs=.5)),
    'Invalid fixed_crossperson_corr')
  expect_error(ravel_fit(d,parental_sibs='maternal',fixed_crossperson=c(e1_m_ms=.1),
    fixed_crossperson_corr=c(e1_m_ms=.5)),'both the covariance and correlation scales')
  expect_error(ravel_fit(d,parental_sibs='maternal',fit_constraints='equilibrium',
    fixed_crossperson_corr=c(e1_m_ms=.1)),'conflict')
  expect_error(ravel_fit(d,parental_sibs='maternal',fixed_crossperson_corr=c(e1_m_ms=1.2)),
    'lie in')
  expect_error(ravel_fit(d,parental_sibs='maternal',fixed_crossperson_corr=c(e1_m_ms=.5,e1_m_ms=.1)),
    'unique parameter names')
  expect_error(ravel_fit(d,parental_sibs='maternal',fit_constraints='equilibrium',
    fixed_crossperson=c(e1_f_ms=0),fixed_crossperson_corr=c(e1_f_ms=0)),
    'both the covariance and correlation scales')
  # Explicit matching values do not conflict with the constrained defaults.
  matching <- ravel_fit(d,parental_sibs='maternal',fit_constraints='n_gen',
    fixed_crossperson_corr=.parental_sib_residual_benchmarks('ms'))
  expect_true(matching$converged)
})


test_that("an incompatible sibling correction preserves visible single-call pi inference", {
  original <- .hidden_correction_sib_extend
  testthat::local_mocked_bindings(.hidden_correction_sib_extend = function(...) {
    out <- original(...)
    out$beta_total[] <- out$beta_total_std[] <- NA_real_
    out$invalid_reason <- 'incompatible_parental_sibling_covariance'
    out
  })
  d <- review_sib_scores(2000)
  expect_warning(fit <- ravel_fit(d,parental_sibs='maternal',pi=.4,pi_se=.05),
    'The visible fit is retained')
  expect_true(fit$converged)
  expect_true(all(is.finite(fit$beta_vis)))
  expect_match(fit$pi_correction_error,'incompatible')
  expect_null(fit$pi_correction)
  tab <- ravel_to_dataframe(fit,include_raw=TRUE,include_vanilla=FALSE)
  expect_output(print_simple(fit),'Corrected inference unavailable')
  expect_identical(attr(tab,'pi_correction_error'),fit$pi_correction_error)
  expect_true(all(is.finite(tab$se[tab$scale=='visible'])))
  expect_true(all(is.na(tab$se[tab$scale!='visible'])))
  expect_true(all(is.na(tab$estimate[tab$scale!='visible' & grepl('beta_|delta_',tab$parameter)])))
  # A malformed user joint covariance is still an error, not swallowed by the
  # retained-visible-fit path.
  cross <- setNames(rep(1,10),c(.report_core_labels,.sib_core_labels('ms')))
  expect_error(ravel_fit(d,parental_sibs='maternal',pi=.4,pi_se=.05,core_pi_cov=cross),
    'positive semidefinite')
})
# Body intended for tests/testthat/test-parental-sibs-core-review.R (after the
# existing "an incompatible sibling correction preserves visible single-call pi
# inference" test; it reuses review_sib_scores() defined in that file).

test_that("a failed pi correction does not leak fixed-pi parental-contrast intervals", {
  d <- review_sib_scores(2000)
  ref <- suppressWarnings(ravel_fit(d, parental_sibs = 'maternal', pi = .4))$beta_vis
  original <- .hidden_correction_sib_extend
  # Incompatible only at the full-sample estimate, so bootstrap replicates
  # still store finite fixed-pi total draws (a borderline sample).
  testthat::local_mocked_bindings(.hidden_correction_sib_extend = function(hc, beta_vis, Sigma_xx, h, roles) {
    out <- original(hc, beta_vis, Sigma_xx, h, roles)
    if (max(abs(beta_vis[c('bD','bm','bf')] - ref[c('bD','bm','bf')])) < 1e-12) {
      out$beta_total[] <- out$beta_total_std[] <- NA_real_
      out$invalid_reason <- 'incompatible_parental_sibling_covariance'
    }
    out
  })
  fit <- suppressWarnings(ravel_fit(d, parental_sibs = 'maternal', pi = .4, pi_se = .05,
    n_boot = 20, seed = 3))
  expect_match(fit$pi_correction_error, 'incompatible')
  expect_gt(sum(stats::complete.cases(fit$bootstrap$total_std_draws)), 1)
  for (u in c('bootstrap', 'lavaan')) {
    tab <- ravel_to_dataframe(fit, uncertainty = u, include_raw = TRUE, include_vanilla = FALSE,
      parent_contrasts = TRUE)
    pc <- tab[tab$parameter %in% c('beta_mid', 'delta_M_mid', 'delta_P_mid'), ]
    vis <- pc[pc$scale == 'visible', ]
    cor <- pc[pc$scale != 'visible', ]
    expect_equal(nrow(vis), 3L)
    expect_equal(nrow(cor), 6L)
    expect_true(all(is.finite(c(vis$se, vis$ci_low, vis$ci_high, vis$p_value))))
    expect_true(all(is.na(c(cor$se, cor$ci_low, cor$ci_high, cor$p_value))))
    expect_true(all(cor$uncertainty_source == 'unavailable (incompatible pi/assortment state)'))
    direct <- .ravel_parent_contrasts(fit, .resolve_uncertainty(fit, u), 'corrected_std')
    expect_true(all(is.na(c(direct$se, direct$ci_low, direct$ci_high))))
  }
  out <- capture.output(print_simple(fit, parent_contrasts = TRUE, uncertainty = 'bootstrap'))
  expect_false(any(grepl('_mid +NA \\[', out)))
})

test_that("a failed pi correction clears every corrected interval in the table", {
  d <- review_sib_scores(2000)
  ref <- suppressWarnings(ravel_fit(d, parental_sibs = 'maternal', pi = .4))$beta_vis
  original <- .hidden_correction_sib_extend
  testthat::local_mocked_bindings(.hidden_correction_sib_extend = function(hc, beta_vis, Sigma_xx, h, roles) {
    out <- original(hc, beta_vis, Sigma_xx, h, roles)
    if (max(abs(beta_vis[c('bD','bm','bf')] - ref[c('bD','bm','bf')])) < 1e-12) {
      out$beta_total[] <- out$beta_total_std[] <- NA_real_
      out$invalid_reason <- 'incompatible_parental_sibling_covariance'
    }
    out
  })
  fit <- suppressWarnings(ravel_fit(d, parental_sibs = 'maternal', pi = .4, pi_se = .05, n_boot = 20, seed = 3))
  for (u in c('bootstrap', 'lavaan')) {
    tab <- ravel_to_dataframe(fit, uncertainty = u, include_raw = TRUE, include_vanilla = FALSE)
    cor <- tab[tab$scale %in% c('corrected_std', 'corrected_raw'), ]
    expect_gt(nrow(cor), 0L)
    expect_true(all(is.na(c(cor$se, cor$ci_low, cor$ci_high, cor$p_value))))
  }
})
