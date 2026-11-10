test_that("cached nonlinear tangents do not veto a valid final score model", {
  set.seed(9102)
  n <- 1000L
  G <- MASS::mvrnorm(n, rep(0,3), matrix(c(1,.65,.65,.65,1,.3,.65,.3,1),3))
  K <- matrix(c(1,.5,.5,.5,1,0,.5,0,1),3)
  d <- data.frame(FID=seq_len(n))
  for(j in 1:2) d[paste0("PGS",j,"_",c("child","mother","father"))] <-
    G+MASS::mvrnorm(n,rep(0,3),c(.4,1.3)[j]*K)
  d$Y <- drop(G%*%c(.32,.13,.08))+rnorm(n)
  fit <- ravel_fit(d,fit_constraints="estimated_state",pi=.4)
  baseline <- ravel_test_score_signal(fit)
  expect_true(baseline$status %in% c("ok","ok_with_warnings"),info=baseline$reason)
  score <- baseline$fits$score
  before <- serialize(score,NULL)
  stale <- score
  p <- ncol(lavaan::lavInspect(stale,"delta"))
  stale@Model@eq.constraints.K <- diag(p)
  rank <- .ravel_score_signal_rank(stale)
  expect_false(rank$stored_tangent_match)
  expect_true(rank$tangent_match)
  expect_true(.ravel_score_signal_validate_fit(stale)$valid)
  fixed <- .ravel_score_signal_lrt_fit(stale,rank)
  expect_true(.ravel_score_signal_rank(fixed)$stored_tangent_match)
  expect_equal(lavaan::coef(fixed),lavaan::coef(score))
  expect_identical(serialize(score,NULL),before)
  # Exercise the complete test, with stale cached bases in every refitted model.
  original_fit <- .ravel_score_signal_fit
  testthat::local_mocked_bindings(.ravel_score_signal_fit=function(...) {
    ans <- original_fit(...)
    ans@Model@eq.constraints.K <- diag(ncol(lavaan::lavInspect(ans,"delta")))
    ans
  },.package="RAVEL")
  result <- ravel_test_score_signal(fit)
  expect_true(result$status %in% c("ok","ok_with_warnings"),info=result$reason)
  expect_equal(result$p_value,baseline$p_value,tolerance=1e-7)
  expect_equal(result$statistic,baseline$statistic,tolerance=1e-7)
  expect_false(result$diagnostics$score$rank$stored_tangent_match)
  for (bad_K in list(matrix(NA_real_,p,rank$tangent_dimension),
      matrix(0,p,rank$tangent_dimension),matrix(0,p,1L))) {
    altered <- score; altered@Model@eq.constraints.K <- bad_K
    rr <- .ravel_score_signal_rank(altered)
    expect_false(rr$stored_tangent_match)
    expect_true(rr$tangent_match)
    expect_true(.ravel_score_signal_validate_fit(altered)$valid)
  }
  # A nonfinite actual constraint Jacobian must still fail validation.
  invalid <- score; invalid@Model@con.jac[1,1] <- NA_real_
  expect_error(.ravel_score_signal_validate_fit(invalid),"constraint Jacobian")
})
