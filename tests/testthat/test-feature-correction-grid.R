feature_grid_report <- function(mode = "none", generations = 3L, bootstrap = FALSE) {
  state <- RAVEL:::.am_state_components(.3, .4,
    if (mode == "equilibrium") "equilibrium" else "n_gen", generations)
  eta <- .5 * (1 + state$rho_G) * sqrt(state$vG / state$vGc)
  theta <- setNames(c(.3, .08, .04, state$rho_G, eta, eta), RAVEL:::.report_core_labels)
  A <- rbind(diag(5), c(0, 0, 0, 0, 1))
  V <- A %*% diag(1e-5, 5) %*% t(A)
  dimnames(V) <- list(names(theta), names(theta))
  packet <- list(schema = "ravel_visible_transfer_v1",
    diagnostics = list(admissible = TRUE, estimated_state_identified = TRUE),
    settings = list(uncertainty = if (bootstrap) "bootstrap" else "lavaan",
      residual_model = "fixed", n_complete = 10000L, ci_level = .95,
      constraint_spec = if (mode == "none") NULL else
        list(mode = mode, h = .4,
          n_generations = if (mode == "n_gen") generations else NA_integer_,
          depends_on_pi = mode == "n_gen" && generations > 1L)),
    theta = theta, vcov_core = V, outcome = list(name = "trait"),
    software = list(transform_fingerprint = RAVEL:::.report_transform_fingerprint()))
  if (bootstrap) {
    draws <- sweep(outer(seq(-.002, .002, length.out = 20), c(1, .3, -.2, .1, .15, .15)),
      2, theta, "+")
    colnames(draws) <- names(theta)
    packet$bootstrap <- list(core_draws = draws, minimum_success = 16L)
  }
  list(packets = list(fixed = packet),
    summary = data.frame(outcome_variable = "trait", regression = "RAVEL",
      residual_model = "fixed", outcome_scale = "prepared_phenotype_SD", N = 10000L),
    settings = list(ci_level = .95, fit_constraints = mode, h = .4,
      am_n_generations = generations))
}

test_that("omitted grid histories inherit constrained modes and preserve the legacy default", {
  for (mode in c("none", "n_gen", "equilibrium", "estimated_state")) {
    out <- feature_grid_report(mode)
    grid <- ravel_correct_grid(out, .4)
    expected <- ravel_correct_packet(out$packets$fixed, h = .4)
    expect_true(all(grid$correction_status == "ok"), info = mode)
    expect_equal(grid$estimate, expected$estimate, info = mode)
    expect_equal(grid$se, expected$se, info = mode)
    generations <- switch(mode, none = 15, n_gen = 3, equilibrium = Inf, estimated_state = NA_real_)
    expect_equal(unique(grid$n_generations), generations, info = mode)
    expect_identical(unique(grid$am_correction), if (mode == "none") "n_gen" else mode)
    expect_false("pi_se" %in% names(grid))
  }
})

test_that("estimated-state grids sweep pi without manufacturing generations or scalar attributes", {
  out <- feature_grid_report("estimated_state")
  before <- out
  grid <- ravel_correct_grid(out, c(.3, .4, .4))
  expect_identical(out, before)
  expect_equal(nrow(grid), 16L)
  expect_true(all(is.na(grid$n_generations)))
  expect_true(all(grid$am_correction == "estimated_state"))
  expect_equal(unique(grid$scenario_id), 1:2)
  for (h in c(.3, .4)) {
    expected <- ravel_correct_packet(out$packets$fixed, h = h)
    expect_equal(grid$estimate[grid$hidden_fraction == h], expected$estimate)
    expect_equal(grid$se[grid$hidden_fraction == h], expected$se)
  }
  expect_null(attr(grid, "estimated_state"))
  expect_error(ravel_correct_grid(out, .4, n_generations = 15), "Do not supply n_generations")
  expect_error(ravel_correct_grid(out, .4, n_generations = Inf), "Do not supply n_generations")
  out$packets$fixed$diagnostics$estimated_state_identified <- NULL
  expect_error(ravel_correct_grid(out, .4), "identification check")
})

test_that("explicit history and pi mismatches are retained as failed cells", {
  finite <- feature_grid_report("n_gen")
  expect_warning(grid <- ravel_correct_grid(finite, c(.4, .5), c(3, 15, Inf)), "unavailable")
  valid <- grid$hidden_fraction == .4 & grid$n_generations == 3
  expect_true(all(grid$correction_status[valid] == "ok"))
  expect_true(all(grid$correction_status[!valid] == "failed"))
  expect_true(all(is.na(grid$estimate[!valid])))
  expect_true(all(nzchar(grid$error[!valid])))
  equilibrium <- feature_grid_report("equilibrium")
  expect_warning(grid <- ravel_correct_grid(equilibrium, .4, c(3, Inf)), "unavailable")
  expect_true(all(grid$correction_status[is.infinite(grid$n_generations)] == "ok"))
  expect_true(all(grid$correction_status[grid$n_generations == 3] == "failed"))
  # n=1 fitting is pi-independent, so its default-history sweep is legitimate.
  first <- feature_grid_report("n_gen", generations = 1L)
  expect_true(all(ravel_correct_grid(first, c(.3, .4))$correction_status == "ok"))
})

test_that("grid pi SE matches scalar full covariance propagation in every compatible mode", {
  sensitivity <- setNames(c(.01, -.02, .01, .03, .01, .01), RAVEL:::.report_core_labels)
  for (mode in c("none", "equilibrium", "estimated_state", "n_gen")) {
    out <- feature_grid_report(mode)
    s <- if (mode == "n_gen") sensitivity else NULL
    grid <- ravel_correct_grid(out, .4, pi_se = .02, core_pi_sensitivity = s)
    expected <- ravel_correct_pi(out$packets$fixed, .4, pi_se = .02, core_pi_sensitivity = s)
    expect_equal(grid$estimate, expected$estimate, info = mode)
    expect_equal(grid$se, expected$se, info = mode)
    expect_equal(grid$ci_low, expected$ci_low, info = mode)
    expect_equal(grid$ci_high, expected$ci_high, info = mode)
    expect_true(all(grid$pi_se == .02))
    expect_true(all(grid$pi_uncertainty_source == "se"))
    expect_true(all(grid$inference == "lavaan_pi_delta_wald"))
    expect_null(attr(grid, "covariance"))
    zero <- ravel_correct_grid(out, .4, pi_se = 0)
    conditional <- ravel_correct_grid(out, .4)
    expect_equal(zero$estimate, conditional$estimate, tolerance = 1e-12)
    expect_equal(zero$se, conditional$se, tolerance = 1e-8)
  }
  out <- feature_grid_report("n_gen")
  expect_warning(missing <- ravel_correct_grid(out, .4, pi_se = .02), "unavailable")
  expect_true(all(missing$correction_status == "failed"))
  expect_true(all(grepl("core_pi_sensitivity", missing$error)))
  listed <- ravel_correct_grid(out, .4, pi_se = .02,
    core_pi_sensitivity = list(fixed = sensitivity))
  direct <- ravel_correct_grid(out, .4, pi_se = .02, core_pi_sensitivity = sensitivity)
  expect_equal(listed, direct)
})

test_that("pi-enabled grid bootstrap uses the documented covariance-delta route", {
  out <- feature_grid_report("estimated_state", bootstrap = TRUE)
  grid <- ravel_correct_grid(out, c(.3, .4), pi_se = .02)
  for (h in c(.3, .4)) {
    scalar <- ravel_correct_pi(out$packets$fixed, h, pi_se = .02)
    rows <- grid$hidden_fraction == h
    expect_equal(grid$estimate[rows], scalar$estimate)
    expect_equal(grid$se[rows], scalar$se)
    expect_equal(grid$ci_low[rows], scalar$ci_low)
  }
  expect_true(all(grid$inference == "bootstrap_covariance_pi_delta_wald"))
  expect_true(all(grid$n_boot_success == 20L))
  expect_identical(attr(grid, "settings")$uncertainty, "bootstrap")
  override <- ravel_correct_grid(out, .4, pi_se = .02, uncertainty = "lavaan")
  expect_true(all(override$inference == "lavaan_pi_delta_wald"))
})

test_that("grid uncertainty validation and missing-model history remain explicit", {
  out <- feature_grid_report("estimated_state")
  for (se in list(-.1, NA_real_, Inf, c(.01, .02), 1i, "0.02"))
    expect_error(ravel_correct_grid(out, .4, pi_se = se), "pi_se")
  s <- setNames(rep(0, 6), RAVEL:::.report_core_labels)
  expect_error(ravel_correct_grid(out, .4, core_pi_sensitivity = s), "requires pi_se")
  expect_error(ravel_correct_grid(out, .4, pi_se = .02,
    core_pi_sensitivity = list(free = s)), "selected residual models")
  expect_error(ravel_correct_grid(out, .4, pi_se = .02, core_pi_sensitivity = 1:6), "six core")
  expect_warning(boundary <- ravel_correct_grid(out, c(0, .4), pi_se = .02), "unavailable")
  expect_true(all(boundary$correction_status[boundary$hidden_fraction == 0] == "failed"))
  expect_true(all(boundary$correction_status[boundary$hidden_fraction == .4] == "ok"))
  out$packets["fixed"] <- list(NULL)
  expect_warning(absent <- ravel_correct_grid(out, .4), "unavailable")
  expect_true(all(absent$am_correction == "estimated_state"))
  expect_true(all(is.na(absent$n_generations)))
  expect_true(all(absent$correction_status == "unavailable_model"))
})

test_that("grid retains binary outcome and interval metadata for downstream plotting", {
  out <- feature_grid_report("estimated_state")
  out$packets$fixed$outcome$type <- "binary_liability"
  out$summary$outcome_scale <- "SEM_probit_liability"
  out$settings$outcome_type <- "binary_liability"
  grid <- ravel_correct_grid(out, .4, pi_se = .02, ci_level = .9)
  expect_true(all(grid$outcome_scale == "SEM_probit_liability"))
  expect_equal(attr(grid, "ci_level"), .9)
  expect_equal(attr(grid, "settings")$ci_level, .9)
  expect_identical(attr(grid, "settings")$outcome_type, "binary_liability")
})

test_that("mixed packet histories receive distinct scenario IDs without overriding either fit", {
  out <- feature_grid_report("equilibrium")
  estimated <- feature_grid_report("estimated_state")$packets$fixed
  estimated$settings$residual_model <- "free"
  out$packets$free <- estimated
  out$summary <- rbind(out$summary, transform(out$summary, residual_model = "free"))
  grid <- ravel_correct_grid(out, c(.3, .4))
  expect_true(all(grid$correction_status == "ok"))
  expect_true(all(is.infinite(grid$n_generations[grid$residual_model == "fixed"])))
  expect_true(all(is.na(grid$n_generations[grid$residual_model == "free"])))
  expect_equal(length(unique(grid$scenario_id)), 4L)
  s <- setNames(rep(0, 6), RAVEL:::.report_core_labels)
  expect_error(ravel_correct_grid(out, .4, pi_se = .02, core_pi_sensitivity = s), "multiple residual models")
  expect_true(all(ravel_correct_grid(out, .4, pi_se = .02,
    core_pi_sensitivity = list(fixed = s, free = s))$correction_status == "ok"))
})

test_that("multiple grid scenarios plot separately and render with the right inference caption", {
  estimated <- feature_grid_report("estimated_state")
  estimated$settings$outcome_type <- "binary_liability"
  estimated$packets$fixed$outcome$type <- "binary_liability"
  estimated$summary$outcome_scale <- "SEM_probit_liability"
  uncertain <- ravel_correct_grid(estimated, c(.3, .4, .5), pi_se = .02, ci_level = .9)
  p <- ravel_plot_report(uncertain)
  expect_equal(length(unique(p$data$model)), 3L)
  expect_true(p$binary)
  expect_match(p$caption, "90%")
  expect_match(p$caption, "external pi uncertainty")
  historical <- ravel_correct_grid(feature_grid_report(), c(.3, .4), c(1, 3, Inf))
  q <- ravel_plot_report(historical)
  expect_equal(length(unique(q$data$model)), 6L)
  bootstrap <- ravel_correct_grid(feature_grid_report("estimated_state", bootstrap = TRUE), .4)
  expect_match(ravel_plot_report(bootstrap)$caption, "bootstrap \\(percentile\\)")
  path <- tempfile(fileext = ".pdf")
  on.exit(unlink(path), add = TRUE)
  grDevices::pdf(path, width = 12, height = 9)
  device <- grDevices::dev.cur()
  tryCatch({
    expect_silent(print(p))
    expect_silent(print(q))
  }, finally = grDevices::dev.off(device))
  expect_gt(file.info(path)$size, 1000)
})

test_that("grid printing preserves scenario identifiers and incompatible-cell explanations", {
  out <- feature_grid_report("n_gen")
  grid <- suppressWarnings(ravel_correct_grid(out, c(.4, .5), c(3, 15)))
  for (printer in list(print_simple, print_ravel_results)) {
    text <- paste(capture.output(printer(grid)), collapse = "\n")
    expect_match(text, "scenario_id")
    expect_match(text, "hidden_fraction")
    expect_match(text, "n_generations")
    expect_match(text, "correction_status")
    expect_match(text, "failed")
    expect_match(text, "Refit|refit|match")
  }
})
