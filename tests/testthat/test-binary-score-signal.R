binary_signal_fixture <- local({
  cache <- new.env(parent = emptyenv())
  function(alternative = FALSE) {
    key <- as.character(alternative)
    if (exists(key, cache, inherits = FALSE)) return(cache[[key]])
    set.seed(58121L)
    n <- 2200L
    C <- matrix(c(1, .6, .6, .6, 1, .2, .6, .2, 1), 3L)
    K <- matrix(c(1, .5, .5, .5, 1, 0, .5, 0, 1), 3L)
    g <- matrix(rnorm(n * 3L), n, 3L) %*% chol(C)
    e1 <- matrix(rnorm(n * 3L), n, 3L) %*% chol(K)
    e2 <- matrix(rnorm(n * 3L), n, 3L) %*% chol(K)
    y <- .4 * g[, 1L] + .1 * g[, 2L] + .05 * g[, 3L] + rnorm(n)
    if (alternative) y <- y + .65 * e1[, 1L]
    d <- data.frame(Y = as.numeric(y > .6))
    roles <- c("child", "mother", "father")
    for (j in 1:3) {
      d[[paste0("PGS1_", roles[j])]] <- g[, j] + sqrt(.6) * e1[, j]
      d[[paste0("PGS2_", roles[j])]] <- g[, j] + sqrt(1.3) * e2[, j]
    }
    fit <- suppressWarnings(ravel_fit(d, outcome = "binary_liability", fit_constraints = "equilibrium", pi = 0))
    # Failed H0 fits retain their syntax but not all settings; supply the known
    # generating fit contract, as the automatic wrapper does.
    fit$outcome <- "binary_liability"
    fit$constraint_spec <- fit$constraint_spec %||% list(mode = "equilibrium", h = 0,
      n_generations = NA_integer_, depends_on_pi = FALSE)
    answer <- .ravel_test_binary_score_signal(fit, d)
    cache[[key]] <- list(data = d, fit = fit, answer = answer)
    cache[[key]]
  }
})

test_that("binary score diagnostic propagates measurement uncertainty and recognizes score-specific signal", {
  null <- binary_signal_fixture(FALSE)$answer
  alt <- binary_signal_fixture(TRUE)$answer
  expect_identical(null$status, "ok")
  expect_equal(null$df, 3L)
  expect_true(is.finite(null$p_value) && null$p_value >= 0 && null$p_value <= 1)
  expect_lt(null$diagnostics$measurement_vcov_relative_error, 1e-7)
  expect_gt(null$diagnostics$null_conditional_variance_slack, 0)
  expect_identical(alt$status, "ok")
  expect_lt(alt$p_value, .001)
  expect_match(null$test, "Wald")
})

test_that("binary score diagnostic preserves explicit limitations", {
  x <- binary_signal_fixture(FALSE)
  bad <- x$data; bad$Y[] <- 0
  expect_identical(.ravel_test_binary_score_signal(x$fit, bad)$status, "invalid_data")
  bad <- x$data; bad$PGS1_child[1L] <- NA_real_
  expect_identical(.ravel_test_binary_score_signal(x$fit, bad)$status, "unsupported_missingness")
  original <- x$fit
  adjusted <- original; adjusted$fixed_effects <- "sex"
  expect_identical(.ravel_test_binary_score_signal(adjusted, x$data)$status, "unsupported_design")
  expect_identical(original, x$fit)
})

test_that("binary score map agrees with the independently extracted covariance matrices", {
  sf <- binary_signal_fixture(FALSE)$answer$fits$score
  scores <- lavaan::lavNames(sf, "ov"); latents <- lavaan::lavNames(sf, "lv")
  map <- .ravel_binary_score_map(lavaan::coef(sf), lavaan::parameterTable(sf), scores, latents)
  m <- lavaan::lavInspect(sf, "est")
  expected <- solve(lavaan::lavInspect(sf, "cov.ov"), m$lambda %*% m$psi)
  expect_equal(unname(map$A), unname(expected), tolerance = 1e-10)
  expect_equal(as.numeric(map$conditional), as.numeric(m$psi - m$psi %*% t(m$lambda) %*% expected), tolerance = 1e-10)
})

test_that("binary diagnostic preserves fixed outcome paths and dependence restrictions", {
  x <- binary_signal_fixture(FALSE)
  constrained <- x$fit
  constrained$model <- sub("bm*Gm + bf*Gf", "0*Gm + 0*Gf", constrained$model, fixed = TRUE)
  expect_identical(.ravel_test_binary_score_signal(constrained, x$data)$status, "unsupported_fixed_paths")
  dependent <- x$fit
  dependent$score_signal_input$dependent_families <- TRUE
  expect_identical(.ravel_test_binary_score_signal(dependent, x$data)$status, "unsupported_design")
  d <- x$data; d$FID <- rep(seq_len(nrow(d)/2), each = 2)
  expect_identical(.ravel_test_binary_score_signal(x$fit, d)$status, "unsupported_design")
})

test_that("binary diagnostic is invariant to affine score units", {
  x <- binary_signal_fixture(FALSE)
  d <- x$data
  scores <- setdiff(names(d), "Y")
  # Stay interior to the original model's absolute 1e-8 variance bounds.
  scale <- c(1e-3, 1e3, .3, 1.8, 4, .03)
  shift <- c(.01, 100, -1, 7, -.2, 3)
  for (j in seq_along(scores)) d[[scores[j]]] <- d[[scores[j]]] * scale[j] + shift[j]
  z <- .ravel_test_binary_score_signal(x$fit, d)
  expect_identical(z$status, "ok", info = z$reason)
  expect_equal(z$statistic, x$answer$statistic, tolerance = 2e-5)
  expect_equal(z$p_value, x$answer$p_value, tolerance = 2e-5)
  expect_equal(z$df, x$answer$df)
})

test_that("binary diagnostic avoids nearly proportional score rows in its restriction basis", {
  # This independent null draw previously selected two nearly proportional
  # rows for the same role (condition number > 3e7), losing the Wald covariance.
  set.seed(620014L)
  n <- 1500L
  C <- matrix(c(1, .6, .6, .6, 1, .2, .6, .2, 1), 3L)
  K <- matrix(c(1, .5, .5, .5, 1, 0, .5, 0, 1), 3L)
  g <- matrix(rnorm(n * 3L), n, 3L) %*% chol(C)
  e1 <- matrix(rnorm(n * 3L), n, 3L) %*% chol(K)
  e2 <- matrix(rnorm(n * 3L), n, 3L) %*% chol(K)
  y <- drop(g %*% c(.4, .1, .05)) + rnorm(n)
  d <- data.frame(Y = as.numeric(y > .6))
  roles <- c("child", "mother", "father")
  d[paste0("PGS1_", roles)] <- g + sqrt(.6) * e1
  d[paste0("PGS2_", roles)] <- g + sqrt(1.3) * e2
  z <- .ravel_test_binary_score_signal(binary_signal_fixture(FALSE)$fit, d)
  expect_identical(z$status, "ok", info = z$reason)
  expect_lt(z$diagnostics$restriction_basis_condition, 10)
  expect_equal(length(z$diagnostics$restriction_basis_scores), 3L)
  eigenvalues <- eigen(z$diagnostics$restriction_vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(eigenvalues) / max(eigenvalues), .01)
})
