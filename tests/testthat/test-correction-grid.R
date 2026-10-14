grid_test_report <- function(bootstrap = FALSE) {
  make_packet <- function(model, shift) {
    theta <- setNames(c(.3 + shift, .1, .02, .1, .55, .55), RAVEL:::.report_core_labels)
    V <- matrix(.00002, 6, 6, dimnames = list(names(theta), names(theta)))
    diag(V) <- .0001
    packet <- list(schema = "ravel_visible_transfer_v1", diagnostics = list(admissible = TRUE),
      settings = list(uncertainty = if (bootstrap) "bootstrap" else "lavaan",
        residual_model = model, n_complete = 16000L, ci_level = .9),
      outcome = list(name = "trait"), theta = theta, vcov_core = V,
      software = list(transform_fingerprint = RAVEL:::.report_transform_fingerprint()))
    if (bootstrap) {
      draws <- sweep(outer(seq(-.015, .015, length.out = 20), c(1, .3, -.2, .1, .15, .2)),
        2, theta, "+")
      colnames(draws) <- names(theta)
      packet$bootstrap <- list(core_draws = draws, minimum_success = 16L)
    }
    packet
  }
  list(packets = list(fixed = make_packet("fixed", 0), free = make_packet("free", .02)),
    summary = data.frame(outcome_variable = "trait", regression = "RAVEL",
      residual_model = c("fixed", "free"), outcome_scale = "prepared_phenotype_SD", N = 16000L),
    settings = list(ci_level = .9))
}

test_that("scalar corrections and Cartesian grids match packet calls", {
  out <- grid_test_report()
  before <- out
  grid <- ravel_correct_grid(out, hidden_fraction = c(.4, .5), n_generations = c(1, 3, Inf))
  expect_identical(out, before)
  expect_equal(nrow(grid), 2 * 3 * 2 * 8)
  expect_equal(nrow(unique(grid[c("hidden_fraction", "n_generations")])), 6)
  expect_true(all(grid$correction_status == "ok"))
  expect_true(all(grid$ci_level == .9))
  expect_true(all(grid$N == 16000))
  columns <- c("parameter", "scale", "estimate", "se", "ci_low", "ci_high", "p_value", "inference")
  for (model in c("fixed", "free")) for (h in c(.4, .5)) for (g in c(1, 3, Inf)) {
    expected <- ravel_correct_packet(out$packets[[model]], h = h,
      am_correction = if (is.infinite(g)) "equilibrium" else "n_gen",
      am_n_generations = if (is.infinite(g)) 15L else g, ci_level = .9)
    got <- grid[grid$residual_model == model & grid$hidden_fraction == h & grid$n_generations == g, columns]
    rownames(got) <- NULL
    rownames(expected) <- NULL
    expect_equal(got, expected[columns])
  }
  scalar <- ravel_correct_grid(out, .4, 1, residual_models = "free", include_raw = FALSE)
  expect_equal(nrow(scalar), 4)
  expect_identical(scalar$parameter, c("beta_D", "beta_M", "beta_P", "rho_corrected"))
  deduplicated <- ravel_correct_grid(out, c(.4, .4, .5), c(Inf, 1, 1), include_raw = FALSE)
  expect_equal(nrow(deduplicated), 32)
  expect_equal(unique(deduplicated$n_generations), c(Inf, 1))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(out, path)
  expect_equal(ravel_correct_grid(path, .4, 1, residual_models = "free", include_raw = FALSE), scalar)
})

test_that("grid uncertainty preserves bootstrap draws and explicit model-based overrides", {
  out <- grid_test_report(bootstrap = TRUE)
  expected <- ravel_correct_packet(out$packets$fixed, .4, "equilibrium", ci_level = .95)
  got <- ravel_correct_grid(out, .4, Inf, residual_models = "fixed", ci_level = .95)
  expect_equal(got$estimate, expected$estimate)
  expect_equal(got$se, expected$se)
  expect_equal(got$ci_low, expected$ci_low)
  expect_equal(got$ci_high, expected$ci_high)
  expect_true(all(got$n_boot_success == 20))
  model <- ravel_correct_grid(out, .4, Inf, residual_models = "fixed", uncertainty = "lavaan", ci_level = .95)
  ref <- ravel_correct_packet(out$packets$fixed, .4, "equilibrium", uncertainty = "lavaan", ci_level = .95)
  expect_equal(model$estimate, got$estimate)
  expect_equal(model$se, ref$se)
  expect_true(all(model$inference == "lavaan_delta_wald"))
  out$packets$fixed$bootstrap$core_draws[1:5, ] <- NA
  expect_warning(insufficient <- ravel_correct_grid(out, .4, 1, residual_models = "fixed"), "unavailable")
  expect_true(all(is.na(insufficient$se)))
  expect_true(all(insufficient$correction_status == "insufficient_bootstrap"))
  expect_true(all(insufficient$n_boot_success == 15))
})

test_that("failed cells and absent fits are retained without replacing other cells", {
  out <- grid_test_report()
  expect_warning(grid <- ravel_correct_grid(out, c(.4, .9999), 1), "unavailable")
  expect_true(all(grid$correction_status[grid$hidden_fraction == .4] == "ok"))
  bad <- grid[grid$hidden_fraction == .9999, ]
  expect_equal(nrow(bad), 16)
  expect_true(all(is.na(bad$estimate)))
  expect_true(all(bad$correction_status == "failed"))
  expect_true(all(nzchar(bad$error)))
  expect_true(all(nzchar(bad$warning)))
  out$packets["fixed"] <- list(NULL)
  expect_warning(absent <- ravel_correct_grid(out, .4, 1, include_raw = FALSE), "unavailable")
  expect_equal(nrow(absent), 8)
  expect_true(all(absent$correction_status[absent$residual_model == "fixed"] == "unavailable_model"))
  expect_true(all(absent$correction_status[absent$residual_model == "free"] == "ok"))
  expect_true(all(absent$N == 16000))
})

test_that("unsafe inputs and missing prerequisites fail before evaluating the grid", {
  out <- grid_test_report()
  for (h in list(numeric(), c(.4, NA), -0.1, 1, Inf, "0.4"))
    expect_error(ravel_correct_grid(out, h, 1), "hidden_fraction")
  for (g in list(numeric(), c(1, NA), 0, -Inf, 1.5, .Machine$integer.max + 1, "Inf"))
    expect_error(ravel_correct_grid(out, .4, g), "n_generations")
  expect_error(ravel_correct_grid(out$summary, .4, 1), "full")
  expect_error(ravel_correct_grid(out, .4, 1, residual_models = "other"), "residual_models")
  expect_error(ravel_correct_grid(out, .4, 1, ci_level = 1), "ci_level")
  expect_error(ravel_correct_grid(out, .4, 1, include_raw = NA), "TRUE or FALSE")
  expect_error(ravel_correct_grid(out, .4, 1, uncertainty = "bootstrap"), "lacks joint")
  out$packets$fixed$settings$requires_covariate_mapping_review <- TRUE
  # The retired confirmation flag no longer gates the correction or changes it.
  expect_s3_class(ravel_correct_grid(out, .4, 1), "data.frame")
  expect_equal(ravel_correct_grid(out, .4, 1, covariate_mapping_confirmed = TRUE), ravel_correct_grid(out, .4, 1))
  out$packets$fixed$software$transform_fingerprint <- "different"
  expect_error(ravel_correct_grid(out, .4, 1, covariate_mapping_confirmed = TRUE), "align package versions")
})
