projection_test_source <- function(covariates = FALSE, pi = .3) {
  roles <- c("mother", "father")
  moments <- function(theta, pi) {
    GG <- matrix(c(1, theta["rho"], theta["rho"], 1), 2,
      dimnames = list(roles, roles))
    M <- diag(sqrt(1 - pi), 2)
    dimnames(M) <- list(roles, roles)
    TT <- (1 - pi) * GG + diag(pi, 2)
    b <- setNames(theta[c("bm", "bf")], roles)
    CW <- if (covariates) matrix(c(.3, -.1), 2, 1, dimnames = list(roles, "W")) else
      matrix(0, 2, 0, dimnames = list(roles, NULL))
    WW <- if (covariates) matrix(1, 1, 1, dimnames = list("W", "W")) else matrix(0, 0, 0)
    list(beta_visible = b, Sigma_GG = GG, beta_total_raw = b / sqrt(1 - pi),
      Sigma_TT = TT, M = M, Sigma_GW = CW, Sigma_WW = WW,
      outcome_variance = 2, outcome_covariate_cov = if (covariates) .4 else numeric(0))
  }
  theta <- c(bm = .3, bf = -.1, rho = .4)
  V <- diag(c(.003, .002, .001))
  dimnames(V) <- list(names(theta), names(theta))
  list(family_descriptor = list(theta = theta, vcov = V, roles = roles,
    pi = pi, pi_se = 0, moments = moments, has_covariates = covariates,
    has_interactions = FALSE, N = 1000, outcome_variable = "maternal_height"))
}

test_that("reduced corrected regressions use total covariance before projection", {
  f <- projection_test_source()
  r <- ravel_family_regression(f, "mother")
  m <- f$family_descriptor$moments(f$family_descriptor$theta, .3)
  expected <- drop(m$Sigma_TT["mother", , drop = FALSE] %*% m$beta_total_raw)
  expect_equal(unname(r$coefficients$total_raw), unname(expected))
  expect_false(isTRUE(all.equal(unname(r$coefficients$total_raw), unname(m$beta_total_raw["mother"]))))
  expect_named(r$coefficients$total_raw, "mother")
  expect_named(r$variance_explained$total_raw,
    c("genetic_variance", "genetic_R2", "conditional_genetic_variance", "partial_R2"))
  expect_true(all(is.finite(r$table$se)))
  expect_true(all(r$table$inference_status == "ok"))
})

test_that("covariate-adjusted projections use the Schur complement", {
  f <- projection_test_source(TRUE)
  blocked <- ravel_family_regression(f, "mother")
  expect_true(all(is.finite(blocked$coefficients$visible)))
  expect_true(all(is.na(blocked$coefficients$total_raw)))
  expect_match(blocked$correction_status, "Cov\\(T")
  r <- ravel_family_regression(f, "mother", covariate_mapping_confirmed = TRUE)
  m <- f$family_descriptor$moments(f$family_descriptor$theta, .3)
  TW <- t(m$M) %*% m$Sigma_GW
  conditional <- m$Sigma_TT - tcrossprod(TW)
  target <- drop(conditional["mother", , drop = FALSE] %*% m$beta_total_raw) /
    conditional["mother", "mother"]
  expect_equal(unname(r$coefficients$total_raw), unname(target))
  expect_true(all(is.finite(r$table$se)))
})

test_that("pi zero with covariates has finite uncertainty without extra assumptions", {
  f <- projection_test_source(TRUE, pi = 0)
  r <- ravel_family_regression(f, "mother")
  expect_equal(r$coefficients$total_std, r$coefficients$visible)
  expect_true(all(is.finite(r$table$se)))
})

test_that("external pi variance changes corrected but not visible uncertainty", {
  f <- projection_test_source()
  fixed <- ravel_family_regression(f, "mother")
  uncertain <- ravel_family_regression(f, "mother", pi_se = .03)
  expect_equal(fixed$table$se[fixed$table$scale == "visible"],
    uncertain$table$se[uncertain$table$scale == "visible"])
  expect_gt(uncertain$vcov["total_raw:mother", "total_raw:mother"],
    fixed$vcov["total_raw:mother", "total_raw:mother"])
  expect_match(uncertain$table$inference_method[uncertain$table$scale == "total_raw"][1], "pi_delta")
})

test_that("different outcome subsets reuse a validated per-fit pi sensitivity", {
  f <- projection_test_source()
  calls <- 0L
  th <- f$family_descriptor$theta
  f$family_descriptor$pi_dependent <- TRUE
  f$family_descriptor$pi_sensitivity <- .ravel_family_cached_sensitivity(function(pi) {
    calls <<- calls + 1L
    setNames(c(.1, -.03, .04), names(th))
  }, names(th))
  first <- ravel_family_regression(f, "mother", pi_se = .03)
  second <- ravel_family_regression(f, c("mother", "father"), pi_se = .03)
  expect_identical(calls, 1L)
  expect_true(all(is.finite(first$table$se)))
  expect_true(all(is.finite(second$table$se)))
  # The cached and uncached inference are identical.
  f$family_descriptor$pi_sensitivity <- function(pi) setNames(c(.1, -.03, .04), names(th))
  uncached <- ravel_family_regression(f, "mother", pi_se = .03)
  expect_equal(first$vcov, uncached$vcov, tolerance = 1e-12)
})

test_that("pi sensitivity caching is local and excludes invalid results", {
  calls <- 0L
  compute <- function(pi) {
    calls <<- calls + 1L
    c(b = if (calls == 1L) NA_real_ else pi)
  }
  cached <- .ravel_family_cached_sensitivity(compute, "b")
  expect_error(cached(.3), "must be finite")
  expect_equal(cached(.3), c(b = .3))
  expect_equal(cached(.3), c(b = .3))
  expect_identical(calls, 2L)
  expect_equal(cached(.4), c(b = .4))
  expect_identical(calls, 3L)
  other_fit <- .ravel_family_cached_sensitivity(compute, "b")
  expect_equal(other_fit(.3), c(b = .3))
  expect_identical(calls, 4L)
})

test_that("invalid conditional total covariance preserves visible results", {
  f <- projection_test_source(TRUE)
  old <- f$family_descriptor$moments
  f$family_descriptor$moments <- function(theta, pi) {
    m <- old(theta, pi)
    m$M <- m$M * 10
    m
  }
  r <- ravel_family_regression(f, "mother", covariate_mapping_confirmed = TRUE)
  expect_true(is.finite(r$coefficients$visible["mother"]))
  expect_true(is.na(r$coefficients$total_raw["mother"]))
  expect_match(r$correction_status, "not positive definite")
})

test_that("full parameter bootstrap retains percentile reporting", {
  f <- projection_test_source()
  d <- f$family_descriptor
  set.seed(329)
  z <- matrix(rnorm(30 * 3), 30, 3) * .01
  f$family_descriptor$bootstrap_draws <- sweep(z, 2, d$theta, "+")
  colnames(f$family_descriptor$bootstrap_draws) <- names(d$theta)
  r <- ravel_family_regression(f, "mother")
  expect_identical(r$uncertainty, "bootstrap")
  expect_true(all(r$table$inference_method == "bootstrap_percentile"))
  expect_true(all(is.finite(r$table$ci_low)))
})

test_that("visible interaction projections use full fitted product covariance", {
  roles <- c("mother", "father")
  dn <- c(roles, paste(roles, "age", sep = ":"))
  S <- matrix(c(1, .3, .2, .1, .3, 1, -.1, .3,
                .2, -.1, 1.5, .5, .1, .3, .5, 2), 4, 4)
  dimnames(S) <- list(dn, dn)
  b <- setNames(c(.3, -.1, .15, .04), dn)
  m <- list(beta_visible = b[roles], Sigma_GG = S[1:2, 1:2],
    beta_total_raw = b[roles], Sigma_TT = S[1:2, 1:2],
    M = matrix(c(1, 0, 0, 1), 2, dimnames = list(roles, roles)),
    design_cov_visible = S, design_beta_visible = b,
    design_covariate_cov_visible = matrix(0, 4, 0),
    design_names = dn, design_roles = rep(roles, 2),
    design_moderators = c("", "", "age", "age"), design_units = c(1, 1, 2, 2),
    Sigma_WW = matrix(0, 0, 0), outcome_variance = 2)
  out <- .ravel_family_projection_values(m, roles, "mother", 0, FALSE)$value
  idx <- c(1, 3)
  target <- drop(solve(S[idx, idx], S[idx, ] %*% b)) / c(1, 2)
  expect_equal(unname(out[c("visible:mother", "visible:mother:age")]), unname(target))
  expect_equal(unname(out[c("total_std:mother", "total_std:mother:age")]), unname(target))
})

test_that("correlated legacy pi cannot silently become independent", {
  f <- projection_test_source()
  f$family_descriptor$correlated_external_pi <- TRUE
  expect_error(ravel_family_regression(f, "mother", pi_se = .03), "parameter_pi_cov")
  cross <- setNames(rep(0, 3), names(f$family_descriptor$theta))
  expect_s3_class(ravel_family_regression(f, "mother", pi_se = .03,
    parameter_pi_cov = cross), "ravel_family_regression")
})

test_that("a corrected negative outcome residual is diagnosed", {
  f <- projection_test_source(pi = .98)
  r <- ravel_family_regression(f, "mother")
  expect_true(is.finite(r$coefficients$visible["mother"]))
  expect_true(is.na(r$coefficients$total_raw["mother"]))
  expect_match(r$correction_status, "negative outcome residual")
})

test_that("few bootstrap draws permit an explicit auto fallback", {
  f <- projection_test_source()
  th <- f$family_descriptor$theta
  f$family_descriptor$bootstrap_draws <- matrix(rep(th, each = 3), 3,
    dimnames = list(NULL, names(th)))
  expect_warning(r <- ravel_family_regression(f, "mother"), "auto uses joint lavaan")
  expect_identical(r$uncertainty, "lavaan")
  expect_error(ravel_family_regression(f, "mother", uncertainty = "bootstrap"), "Insufficient")
})

test_that("live binary and interaction projections retain the source scales", {
  p <- family_validation_population(pi = 0, rho = .2, parent_variance = 1.25,
    beta = c(child = .1, mother = .2, father = -.05, mat_sib = 0, pat_sib = 0),
    lineage = c(0, 0))
  d <- family_validation_data(p, N = 1500, exact = FALSE)
  d$Y <- as.integer(d$Y > 0)
  f <- ravel_fit(d, outcome = "binary_liability", pi = 0,
    fit_constraints = "equilibrium", am_correction = "equilibrium")
  expect_true(f$converged)
  r <- ravel_family_regression(f, c("child", "mother", "father"))
  expect_equal(unname(r$coefficients$visible), unname(f$beta_vis), tolerance = 1e-8)
  expect_true(all(is.finite(r$table$se)))

  d <- family_validation_data(p, N = 1500, exact = FALSE)
  set.seed(72)
  d$age <- runif(nrow(d), 8, 20)
  d$Y <- d$Y + .02 * (d$age - 10) * d$T_mother
  f <- ravel_fit(d, pi = 0, fit_constraints = "equilibrium", am_correction = "equilibrium",
    interactions = "age", interaction_roles = c("child", "mother", "father"),
    interaction_center = c(age = 10))
  expect_true(f$converged)
  r <- ravel_family_regression(f, c("child", "mother", "father"))
  expect_equal(unname(r$coefficients$visible[1:3]), unname(f$beta_vis), tolerance = 1e-8)
  expect_equal(unname(r$coefficients$visible[4:6]), f$interactions$table$estimate, tolerance = 1e-8)
  expect_true(all(is.finite(r$table$se)))
  reduced <- ravel_family_regression(f, "mother")
  expect_true(all(is.finite(reduced$table$se)))
})
