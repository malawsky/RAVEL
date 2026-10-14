# Integration coverage for the preferred one-call external-pi interface.
fit_pi_cache <- new.env(parent = emptyenv())
fit_pi_data <- function() {
  env <- new.env()
  data("ravel_example", package = "RAVEL", envir = env)
  dat <- env$ravel_example
  names(dat) <- sub("^pgs", "PGS", names(dat))
  dat
}
fit_pi_cached <- function(mode = "estimated_state", bootstrap = FALSE) {
  key <- paste(mode, bootstrap)
  if (!exists(key, fit_pi_cache, inherits = FALSE)) {
    fit_pi_cache[[key]] <- ravel_fit(fit_pi_data(), pi = .4, pi_se = .02,
      fit_constraints = mode, ci_level = .9, n_boot = if (bootstrap) 12L else 0L,
      seed = 9184L)
  }
  fit_pi_cache[[key]]
}
fit_pi_rows <- function(fit, source = "auto") {
  table <- ravel_to_dataframe(fit, uncertainty = source, include_raw = TRUE,
    include_vanilla = FALSE)
  table$parameter[table$parameter == "beta_F"] <- "beta_P"
  table[table$scale != "visible", ]
}

test_that("single-call pi uncertainty equals the two-step correction across fit modes", {
  for (mode in c("none", "equilibrium", "estimated_state")) {
    fitted <- fit_pi_cached(mode)
    expect_true(fitted$converged, info = mode)
    separate <- ravel_correct_pi(fitted, pi = .4, pi_se = .02, ci_level = .9)
    expect_equal(fitted$pi_correction, separate)
    exported <- fit_pi_rows(fitted)
    expect_equal(exported$estimate, separate$estimate, tolerance = 1e-9)
    expect_equal(exported$se, separate$se, tolerance = 1e-9)
    expect_equal(exported$ci_low, separate$ci_low, tolerance = 1e-9)
    expect_equal(exported$ci_high, separate$ci_high, tolerance = 1e-9)
    expect_equal(exported$p_value, separate$p_value, tolerance = 1e-9)
    expect_equal(unique(exported$pi_se), .02)
    expect_true(all(exported$inference == "lavaan_pi_delta_wald"))
    expect_equal(fit_pi_rows(fitted, "lavaan")$se, separate$se)
  }
})

test_that("the pi alias preserves legacy fitting and zero SE preserves fixed-pi inference", {
  legacy <- ravel_fit(fit_pi_data(), h = .4, fit_constraints = "estimated_state")
  alias <- ravel_fit(fit_pi_data(), pi = .4, fit_constraints = "estimated_state")
  certain <- ravel_fit(fit_pi_data(), pi = .4, pi_se = 0, fit_constraints = "estimated_state")
  expect_true(legacy$converged)
  expect_null(legacy$pi_correction)
  expect_null(alias$pi_correction)
  expect_equal(alias$beta_vis, legacy$beta_vis)
  expect_equal(alias$lavaan_inference, legacy$lavaan_inference)
  expect_equal(alias$model, legacy$model)
  expect_equal(certain$lavaan_inference, legacy$lavaan_inference)
  expect_equal(fit_pi_rows(certain)$se, fit_pi_rows(legacy)$se, tolerance = 1e-5)
  expect_equal(fit_pi_rows(certain)$estimate, fit_pi_rows(legacy)$estimate, tolerance = 1e-9)
  matching <- ravel_fit(fit_pi_data(), h = .4, pi = .4, pi_se = .02,
    fit_constraints = "estimated_state")
  expect_true(matching$converged)
  expect_equal(matching$pi_correction$estimate, certain$pi_correction$estimate)
})

test_that("pi inputs and mapping review fail before fitting when invalid", {
  empty <- data.frame()
  expect_error(ravel_fit(empty, h = .3, pi = .4), "h and pi")
  expect_error(ravel_fit(empty, pi_se = .02), "explicitly")
  for (bad in list(NA_real_, Inf, -.1, 1, c(.2, .3), "0.4", .4 + 1i))
    expect_error(ravel_fit(empty, pi = bad), "pi must")
  for (bad in list(NA_real_, Inf, -.1, c(.1, .2), "0.02", 1e200))
    expect_error(ravel_fit(empty, pi = .4, pi_se = bad), "pi_se")
  expect_error(ravel_fit(empty, pi = 0, pi_se = .02), "boundary")
  expect_error(ravel_fit(empty, pi = .4, core_pi_cov = 1:6), "requires pi_se")
  expect_error(ravel_fit(empty, pi = .4, pi_se = .02, core_pi_cov = 1:6), "six core")
})

test_that("printing and export carry external-pi inference and allow new pi corrections", {
  fitted <- fit_pi_cached()
  tab <- as.data.frame(fitted, include_raw = TRUE, include_vanilla = FALSE)
  expect_equal(tab, ravel_to_dataframe(fitted, include_raw = TRUE, include_vanilla = FALSE))
  expect_equal(attr(tab, "pi_uncertainty")$raw_scale_se, .02)
  expect_equal(attr(tab, "pi_correction"), fitted$pi_correction)
  expect_equal(attr(tab, "estimated_state"), attr(fitted$pi_correction, "estimated_state"))
  expect_equal(attr(tab, "ci_level"), .9)
  expect_true(all(grepl("conditional on pi", tab$uncertainty_source[tab$scale == "visible"])))
  for (printer in list(print_simple, print_ravel_results)) {
    text <- paste(capture.output(printer(fitted)), collapse = "\n")
    expect_match(text, "External pi = 0.4000; SE = 0.0200", fixed = TRUE)
    expect_match(text, "Stored visible inference is conditional on pi", fixed = TRUE)
    expect_match(text, "delta-Wald", fixed = TRUE)
    ignored <- capture.output(printed <- withVisible(printer(fitted, uncertainty = "none")))
    expect_false(printed$visible)
  }
  expect_error(ravel_to_dataframe(fitted, uncertainty = "none"), "not available")
  changed <- ravel_correct_pi(fitted, pi = .45, pi_se = .03, ci_level = .9)
  conditional <- fitted
  conditional$pi_correction <- NULL
  conditional$pi_se <- NULL
  expected <- ravel_correct_pi(conditional, pi = .45, pi_se = .03, ci_level = .9)
  expect_equal(changed, expected)
  expect_gt(max(abs(changed$estimate - fitted$pi_correction$estimate)), .001)
  expect_equal(fitted$pi_correction$pi, rep(.4, 8L))
})

test_that("bootstrap core inference stays conditional while total uncertainty includes pi", {
  fitted <- fit_pi_cached("estimated_state", bootstrap = TRUE)
  expect_true(fitted$converged)
  expect_gte(fitted$bootstrap$n_successful, 10L)
  expected <- ravel_correct_pi(fitted, pi = .4, pi_se = .02, ci_level = .9,
    uncertainty = "bootstrap")
  expect_equal(fitted$pi_correction, expected)
  tab <- as.data.frame(fitted, include_raw = TRUE, include_vanilla = FALSE)
  visible <- tab$scale == "visible"
  expect_true(all(tab$inference[visible] == "bootstrap_percentile_CI_normal_p"))
  expect_true(all(tab$inference[!visible] == "bootstrap_covariance_pi_delta_wald"))
  expect_equal(tab$ci_low[visible], unname(fitted$bootstrap$visible$ci["lower", ]))
  expect_equal(tab$se[!visible], expected$se)
  expect_equal(fit_pi_rows(fitted, "bootstrap")$ci_low, expected$ci_low)
  model <- ravel_correct_pi(fitted, pi = .4, pi_se = .02, ci_level = .9,
    uncertainty = "lavaan")
  expect_equal(fit_pi_rows(fitted, "lavaan")$se, model$se)
  # Asking for a different source must not change cached or conditional results.
  expect_equal(fitted$pi_correction, expected)
  expect_equal(fitted$bootstrap$visible$ci["lower", ],
    setNames(tab$ci_low[visible], colnames(fitted$bootstrap$visible$ci)))
  expect_gt(max(abs(fitted$pi_correction$se[1:4] - fitted$bootstrap$total_std$se)), 1e-6)
})

test_that("finite-generation single-call pi uncertainty includes fitting sensitivity", {
  fixture <- cov_test_fixture("continuous", correlated = FALSE,
    with_covariates = FALSE, n = 5000L)
  fitted <- ravel_fit(fixture$dat, pi = .4, pi_se = .02,
    fit_constraints = "n_gen", am_n_generations = 3L)
  expect_true(fitted$converged)
  sensitivity <- attr(fitted$pi_correction, "pi_uncertainty")$core_pi_sensitivity
  expect_gt(max(abs(sensitivity)), 1e-5)
  expect_null(fitted$constraint_refit$pi)
  expect_null(fitted$constraint_refit$pi_se)
  separate <- ravel_correct_pi(fitted, pi = .4, pi_se = .02,
    core_pi_sensitivity = sensitivity)
  expect_equal(separate, fitted$pi_correction)
  expect_equal(fit_pi_rows(fitted)$se, separate$se)
  expect_gt(max(abs(attr(separate, "jacobian")[, "pi"] -
    attr(separate, "partial_jacobian")[, "pi"])), 1e-5)
  expect_error(ravel_correct_pi(fitted, pi = .45, pi_se = .02), "Refit")
})


test_that("one-call covariance and covariate options preserve the explicit correction contract", {
  base <- fit_pi_cached("equilibrium")
  V <- attr(base$pi_correction, "joint_input_covariance")[1:6, 1:6]
  cross <- .02 * .1 * V[, "bD"] / sqrt(V["bD", "bD"])
  fitted <- ravel_fit(fit_pi_data(), h = .4, pi_se = .02,
    core_pi_cov = cross, fit_constraints = "equilibrium", ci_level = .9)
  expected <- ravel_correct_pi(base, pi = .4, pi_se = .02,
    core_pi_cov = cross, ci_level = .9)
  expect_equal(fitted$pi_correction, expected)
  expect_false(attr(fitted$pi_correction, "pi_uncertainty")$independent_of_core)
  expect_gt(max(abs(fitted$pi_correction$se - base$pi_correction$se)), 1e-6)
  fixture <- cov_test_fixture("continuous", correlated = FALSE,
    with_covariates = TRUE, n = 5000L)
  covariates <- ravel_fit(fixture$dat, pi = .4, pi_se = .02,
    fixed_effects = fixture$covariates, fit_constraints = "equilibrium",
    covariate_mapping_confirmed = TRUE)
  expect_true(covariates$converged)
  expect_equal(covariates$pi_correction,
    ravel_correct_pi(covariates, pi = .4, pi_se = .02,
      covariate_mapping_confirmed = TRUE))
  # The retired confirmation flag is ignored: same correction without it.
  expect_equal(ravel_correct_pi(covariates, pi = .4, pi_se = .02), covariates$pi_correction)
})
