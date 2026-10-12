# Independent moment fixture: the accumulated parental state and current
# assortment can differ from every constant-strength integer-generation path.
estimated_state_fixture <- function(parent_total_variance=1.3, rho_total=.35,
                                    h=.4, noise=c(.4,1.6), N=1000L,
                                    units=rep(1,6), beta=c(.4,.1,.05)) {
  q <- 1-h
  vp <- 1+q*(parent_total_variance-1)
  spouse <- q*rho_total*parent_total_variance
  vc <- (vp+spouse)/2+.5
  pc <- (vp+spouse)/2
  V <- matrix(c(vc,pc,pc,pc,vp,spouse,pc,spouse,vp),3)
  E <- matrix(c(1,.5,.5,.5,1,0,.5,0,1),3)
  S <- rbind(cbind(V+noise[1]*E,V),cbind(V,V+noise[2]*E))
  visible <- stats::cov2cor(V)
  loading <- rbind(diag(sqrt(diag(V))),diag(sqrt(diag(V))))
  sy <- loading %*% visible %*% beta
  S <- rbind(cbind(S,sy),c(t(sy),as.numeric(t(beta)%*%visible%*%beta)+.8))
  nm <- c(paste0("PGS1_",c("child","mother","father")),
          paste0("PGS2_",c("child","mother","father")),"Y")
  dimnames(S) <- list(nm,nm)
  i <- seq_len(N)
  Q <- qr.Q(qr(cbind(1,outer(i,seq_len(7L),function(i,j)
    sin(i*j/17)+cos(i*(j+1)/23)))))
  dat <- as.data.frame(sqrt(N)*Q[,-1,drop=FALSE]%*%chol(S))
  names(dat) <- nm
  # Mimic independently standardized supplied PGS weights, then allow arbitrary
  # observed units. Neither operation changes the latent reliability or state.
  dat[,1:6] <- sweep(scale(dat[,1:6]),2,units,"*")
  list(data=dat,visible=visible,beta=beta,vp=vp,growth=vc/vp,
    reliability=cbind(PGS1=diag(V)/(diag(V)+noise[1]),
                      PGS2=diag(V)/(diag(V)+noise[2])))
}

test_that("estimated state recovers equilibrium, finite, and changing histories", {
  r <- .35
  states <- list(
    equilibrium=c(1/(1-r),r),
    finite_three=c(.5*(1+r)*(.5*(1+r)+.5)+.5,r),
    changed_assortment=c(1.45,-.25))
  for (label in names(states)) {
    state <- states[[label]]
    x <- estimated_state_fixture(state[1],state[2])
    fit <- ravel_fit(x$data,h=.4,fit_constraints="estimated_state")
    expect_true(fit$converged,info=label)
    expect_equal(as.numeric(lavaan::fitMeasures(fit$fit,"df")),14,tolerance=1e-4)
    expect_equal(unname(fit$latent_cov_vis),
      unname(c(x$visible[2,3],x$visible[1,2],x$visible[1,3])),tolerance=3e-5)
    expect_equal(unname(fit$beta_vis),x$beta,tolerance=3e-5)
    pe <- lavaan::parameterEstimates(fit$fit)
    expect_equal(pe$est[pe$label=="am_growth"],x$growth,tolerance=3e-5)
    expect_equal(pe$est[pe$label=="am_vG_parent"],x$vp,tolerance=3e-5)
    for (j in 1:2) {
      rel <- fit$reliability[fit$reliability$PGS==paste0("PGS",j),]
      expect_equal(rel$reliability,
        x$reliability[match(rel$role,c("child","mother","father")),j],tolerance=3e-5)
    }
    expect_gt(abs(diff(fit$reliability$reliability[fit$reliability$role=="mother"])),.1)
    expect_equal(unname(fit$crossperson_resid_corr),c(0,.5,.5,0,.5,.5),tolerance=3e-5)
    expect_equal(fit$am_correction,"estimated_state")
    expect_false(fit$constraint_spec$depends_on_pi)
    expect_true(is.na(fit$constraint_spec$n_generations))
    expect_true(is.na(fit$am_n_generations))
  }
})

test_that("estimated-state inference refuses an unidentified reliability direction", {
  x <- estimated_state_fixture(parent_total_variance=1, rho_total=0,
    beta=rep(0,3), N=2000L)
  bad <- suppressWarnings(ravel_fit(x$data,h=.4,fit_constraints="estimated_state"))
  expect_false(bad$converged)
  expect_match(bad$error,"identification")
  expect_false(RAVEL:::.ravel_estimated_state_fit_valid(bad$fit))
  # A positive covariance from a generalized inverse must not hide the missing
  # moment direction: the tangent-rank check is needed independently of PSD.
  pseudoinverse <- bad$fit
  pseudoinverse@vcov$vcov <- diag(nrow(lavaan::vcov(pseudoinverse)))
  expect_false(RAVEL:::.ravel_estimated_state_fit_valid(pseudoinverse))
  expect_error(ravel_correct_pi(bad,.4,pi_se=.02),"converged, admissible")
  # A live object created before the guard must also be checked, even if it
  # claims convergence and a user supplies a seemingly valid core covariance.
  bad$converged <- TRUE
  bad$constraint_spec <- list(mode="estimated_state")
  expect_error(ravel_correct_pi(bad,.4,pi_se=.02),"locally identified")
  informative <- estimated_state_fixture(beta=rep(0,3))
  good <- ravel_fit(informative$data,h=.4,fit_constraints="estimated_state")
  expect_true(good$converged)
  expect_true(RAVEL:::.ravel_estimated_state_fit_valid(good$fit))
})

test_that("estimated-state identification uses actual fixed-path and covariate dimensions", {
  x <- estimated_state_fixture(beta=c(.4,0,0))
  fixed <- ravel_fit(x$data,h=.4,fit_constraints="estimated_state",fix_parent_paths=TRUE)
  expect_true(fixed$converged)
  expect_true(all(is.finite(ravel_correct_pi(fixed,.4,pi_se=.02)$se)))
  # Add an independent measured covariate with a nonzero outcome coefficient.
  d <- x$data
  z <- seq_len(nrow(d))
  z <- stats::residuals(stats::lm(sin(z/3)~as.matrix(d)))
  d$X <- as.numeric(scale(z))
  d$Y <- d$Y + .2*d$X
  covariate <- ravel_fit(d,h=.4,fit_constraints="estimated_state",fixed_effects="X")
  expect_true(covariate$converged)
  expect_true(all(is.finite(ravel_correct_pi(covariate,.4,pi_se=.02,
    covariate_mapping_confirmed=TRUE)$se)))
})

test_that("estimated-state fitting ignores pi and generations and preserves units", {
  x <- estimated_state_fixture()
  a <- ravel_fit(x$data,h=.2,fit_constraints="estimated_state",am_n_generations=2L)
  b <- ravel_fit(x$data,h=.6,fit_constraints="estimated_state",am_n_generations=99L)
  expect_identical(a$model,b$model)
  expect_equal(a$beta_vis,b$beta_vis,tolerance=1e-12)
  expect_equal(a$latent_cov_vis,b$latent_cov_vis,tolerance=1e-12)
  expect_equal(a$reliability$reliability,b$reliability$reliability,tolerance=1e-12)
  expect_equal(unname(RAVEL:::.ravel_constraint_core_h_sensitivity(a)),rep(0,6))
  y <- estimated_state_fixture(units=c(.3,2,1.4,3,.8,1.6))
  c <- ravel_fit(y$data,h=.4,fit_constraints="estimated_state")
  expect_true(c$converged)
  expect_equal(a$beta_vis,c$beta_vis,tolerance=3e-5)
  expect_equal(a$latent_cov_vis,c$latent_cov_vis,tolerance=3e-5)
  expect_equal(a$reliability$reliability,c$reliability$reliability,tolerance=3e-5)
})

test_that("estimated-state fit and correction modes cannot be mixed", {
  x <- estimated_state_fixture()
  expect_error(ravel_fit(x$data,h=.4,fit_constraints="estimated_state",
    am_correction="n_gen"),"must agree")
  expect_error(ravel_fit(x$data,h=.4,am_correction="estimated_state"),
    "requires fit_constraints")
  spec <- list(mode="estimated_state",depends_on_pi=FALSE,n_generations=NA_integer_)
  expect_true(RAVEL:::.ravel_check_constraint_correction(spec,.7,"estimated_state",100L))
  expect_error(RAVEL:::.ravel_check_constraint_correction(spec,.4,"equilibrium",1L),"regime")
})
