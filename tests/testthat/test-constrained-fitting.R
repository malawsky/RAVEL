# Deterministic Gaussian moment fixtures, derived from a raw rank-one AM
# recursion independently of the constraint builder and correction helpers.
constraint_fixture <- function(mode="equilibrium", n=3L, h=.4, rho_total=.35,
                               noise=c(.4,1.3), N=1000L) {
  q <- 1-h
  full_parent <- 1
  if (mode == "equilibrium") full_parent <- 1/(1-rho_total)
  else if (n > 1L) for (i in seq_len(n-1L))
    full_parent <- .5*(1+rho_total)*full_parent+.5
  full_child <- .5*(1+rho_total)*full_parent+.5
  vp <- 1+q*(full_parent-1)
  vc <- 1+q*(full_child-1)
  spouse_cov <- q*rho_total*full_parent
  pc_cov <- .5*(vp+spouse_cov)
  V <- matrix(c(vc,pc_cov,pc_cov,pc_cov,vp,spouse_cov,pc_cov,spouse_cov,vp),3)
  E <- matrix(c(1,.5,.5,.5,1,0,.5,0,1),3)
  S <- rbind(cbind(V+noise[1]*E,V),cbind(V,V+noise[2]*E))
  visible <- stats::cov2cor(V)
  beta <- c(.4,.1,.05)
  loading <- rbind(diag(sqrt(diag(V))),diag(sqrt(diag(V))))
  sy <- loading %*% visible %*% beta
  S <- rbind(cbind(S,sy), c(t(sy), as.numeric(t(beta)%*%visible%*%beta)+.8))
  names <- c(paste0("PGS1_",c("child","mother","father")),
             paste0("PGS2_",c("child","mother","father")),"Y")
  dimnames(S) <- list(names,names)
  # Orthogonal zero-mean columns give exactly the requested sample covariance.
  x <- seq_len(N)
  Q <- qr.Q(qr(cbind(1,outer(x,seq_len(7L),function(i,j) sin(i*j/17)+cos(i*(j+1)/23)))))
  dat <- as.data.frame(sqrt(N-1)*Q[,-1,drop=FALSE] %*% chol(S))
  names(dat) <- names
  list(data=dat, covariance=S, beta=beta, visible=visible,
       reliabilities=cbind(PGS1=diag(V)/(diag(V)+noise[1]),
                           PGS2=diag(V)/(diag(V)+noise[2])))
}

constraint_fit <- function(fixture, mode, h=.4, n=3L, ...) {
  ravel_fit(fixture$data, h=h, fit_constraints=mode, am_n_generations=n, ...)
}

test_that("equilibrium imposes 13 covariance parameters without tying PGS reliabilities", {
  x <- constraint_fixture()
  fit <- constraint_fit(x,"equilibrium")
  expect_true(fit$converged)
  expect_equal(as.numeric(lavaan::fitMeasures(fit$fit,"df")),15,tolerance=1e-5)
  expect_equal(unname(fit$latent_cov_vis),unname(c(x$visible[2,3],x$visible[1,2],x$visible[1,3])),tolerance=2e-5)
  expect_equal(unname(fit$beta_vis),x$beta*sqrt((nrow(x$data)-1)/nrow(x$data)),tolerance=2e-5)
  for(j in 1:2) {
    rel <- fit$reliability[fit$reliability$PGS==paste0("PGS",j),]
    expect_equal(unname(rel$reliability),unname(rep(x$reliabilities[1,j],3)),tolerance=2e-5)
  }
  expect_gt(abs(diff(fit$reliability$reliability[fit$reliability$role=="mother"])),.1)
  expect_equal(unname(fit$crossperson_resid_corr),c(0,.5,.5,0,.5,.5),tolerance=2e-5)
  expect_equal(fit$am_correction,"equilibrium")
  expect_false(fit$constraint_spec$depends_on_pi)
  expect_true(all(is.finite(fit$lavaan_inference$total_std$se)))
})

test_that("finite constraints reproduce known generation states and child reliability", {
  for(n in c(1L,3L,15L,60L)) {
    x <- constraint_fixture("n_gen",n)
    fit <- constraint_fit(x,"n_gen",n=n)
    expect_true(fit$converged,info=paste("generation",n))
    expect_equal(as.numeric(lavaan::fitMeasures(fit$fit,"df")),15,tolerance=1e-4)
    expect_equal(unname(fit$latent_cov_vis),unname(c(x$visible[2,3],x$visible[1,2],x$visible[1,3])),tolerance=3e-5)
    for(j in 1:2) {
      rel <- fit$reliability[fit$reliability$PGS==paste0("PGS",j),]
      expect_equal(rel$reliability,x$reliabilities[match(rel$role,c("child","mother","father")),j],tolerance=3e-5)
    }
    expect_equal(fit$constraint_spec$depends_on_pi,n>1)
  }
})

test_that("constrained fitting is invariant to separate score units", {
  x <- constraint_fixture("n_gen",3L)
  y <- x
  y$data[,1:6] <- sweep(y$data[,1:6],2,c(.3,2,1.4,3,.8,1.6),"*")
  a <- constraint_fit(x,"n_gen")
  b <- constraint_fit(y,"n_gen")
  expect_true(a$converged && b$converged)
  expect_equal(a$beta_vis,b$beta_vis,tolerance=3e-5)
  expect_equal(a$latent_cov_vis,b$latent_cov_vis,tolerance=3e-5)
  expect_equal(a$reliability$reliability,b$reliability$reliability,tolerance=3e-5)
})

test_that("constraint choices preserve legacy defaults and reject conflicts", {
  x <- constraint_fixture()
  a <- ravel_fit(x$data,h=.4,am_correction="equilibrium")
  b <- ravel_fit(x$data,h=.4,am_correction="equilibrium",fit_constraints="none")
  expect_identical(a$model,b$model)
  expect_equal(a$beta_vis,b$beta_vis,tolerance=1e-12)
  expect_error(constraint_fit(x,"equilibrium",am_correction="n_gen"),"must agree")
  expect_error(constraint_fit(x,"n_gen",n=1.5),"single integer")
  expect_error(constraint_fit(x,"n_gen",n=Inf),"single integer")
  expect_error(constraint_fit(x,"equilibrium",crossperson="equal_pc"),"crossperson")
  expect_error(constraint_fit(x,"equilibrium",fixed_crossperson=c(e1_mc=.5)),"conflict")
  expect_error(constraint_fit(x,"equilibrium",fixed_crossperson_corr=c(e1_mc=.4)),"conflict")
})

test_that("finite-fit pi sensitivity includes the fit dependence", {
  x <- constraint_fixture("n_gen",3L)
  fit <- constraint_fit(x,"n_gen")
  sensitivity <- RAVEL:::.ravel_constraint_core_h_sensitivity(fit)
  expect_named(sensitivity,c("bD","bm","bf","rhoMF","eta_m","eta_f"))
  expect_true(all(is.finite(sensitivity)))
  expect_gt(max(abs(sensitivity)),1e-4)
  f1 <- constraint_fit(x,"n_gen",h=.402)
  f0 <- constraint_fit(x,"n_gen",h=.398)
  check <- (c(f1$beta_vis,f1$latent_cov_vis)-c(f0$beta_vis,f0$latent_cov_vis))/.004
  expect_equal(sensitivity,check,tolerance=.015)
  expect_error(RAVEL:::.ravel_check_constraint_correction(fit$constraint_spec,.5,"n_gen",3L),"refit")
  expect_error(RAVEL:::.ravel_check_constraint_correction(fit$constraint_spec,.4,"equilibrium",3L),"regime")
})

test_that("sampled constrained fits support covariates and binary liability", {
  set.seed(8871)
  for(mode in c("equilibrium","n_gen")) {
    x <- constraint_fixture(mode,3L)
    x$data <- as.data.frame(MASS::mvrnorm(3000,rep(0,7),x$covariance))
    plain <- constraint_fit(x,mode)
    expect_true(plain$converged)
    expect_true(RAVEL:::.ravel_constraint_fit_valid(plain$fit))
    x$data$Z <- stats::rnorm(nrow(x$data))
    x$data$Y <- x$data$Y+.2*x$data$Z
    adjusted <- constraint_fit(x,mode,fixed_effects="Z")
    expect_true(adjusted$converged)
    expect_true(all(is.finite(adjusted$lavaan_inference$total_std$se)))
    x$data$Y <- as.integer(x$data$Y>0)
    binary <- constraint_fit(x,mode,outcome="binary_liability",fixed_effects="Z")
    expect_true(binary$converged)
    expect_true(all(is.finite(binary$lavaan_inference$total_std$se)))
  }
})
