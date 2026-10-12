test_that("estimated-state reports retain their model and support aggregate correction", {
  data("ravel_example", package = "RAVEL", envir = environment())
  out <- ravel_dataframe_report(ravel_example, fit_constraints = "estimated_state",
    h = .4, am_n_generations = NA_integer_, make_plot = FALSE)
  expect_identical(out$model_status$status, "ok")
  expect_true(all(out$summary$scale[out$summary$regression == "RAVEL"] == "visible"))
  packet <- out$packet
  expect_true(packet$diagnostics$estimated_state_identified)
  expect_identical(packet$settings$constraint_spec$mode, "estimated_state")
  expect_false(packet$settings$constraint_spec$depends_on_pi)
  expect_true(is.na(packet$settings$constraint_spec$n_generations))
  expect_null(packet$constraint_refit)
  fixed <- ravel_correct_packet(packet, h = .4)
  expect_true(all(is.finite(fixed$estimate)))
  expect_true(all(is.finite(fixed$se)))
  expect_true(all(is.na(fixed$am_n_generations)))
  expect_type(attr(fixed, "estimated_state"), "list")
  uncertain <- ravel_correct_pi(packet, pi = .4, pi_se = .02)
  expect_true(all(uncertain$se >= fixed$se - 1e-8))
  expect_equal(uncertain$estimate, fixed$estimate, tolerance = 1e-8)
  expect_error(ravel_correct_packet(packet, h = .4, am_correction = "n_gen"), "regime")
  grid <- ravel_correct_grid(out, hidden_fraction = .4)
  expect_equal(grid$estimate, fixed$estimate)
  expect_equal(grid$se, fixed$se)
  expect_true(all(is.na(grid$n_generations)))
  unchecked <- packet
  unchecked$diagnostics$estimated_state_identified <- NULL
  expect_error(ravel_correct_packet(unchecked, h = .4), "identification check")
  expect_true(all(is.finite(ravel_correct_pi(packet, pi = .35, pi_se = .02)$se)))
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(packet, tmp)
  expect_equal(ravel_correct_pi(tmp, pi = .4, pi_se = .02), uncertain)
})

test_that("estimated-state report bootstrap resamples the new constraints", {
  data("ravel_example", package = "RAVEL", envir = environment())
  out <- ravel_dataframe_report(ravel_example, fit_constraints = "estimated_state",
    h = .4, n_boot = 10L, seed = 701L, make_plot = FALSE)
  expect_identical(out$model_status$status, "ok")
  draws <- out$packet$bootstrap$core_draws
  expect_equal(dim(draws), c(10L, 6L))
  valid <- apply(is.finite(draws), 1, all)
  expect_true(sum(valid) >= 8L)
  expect_equal(unname(draws[valid, "eta_m"]), unname(draws[valid, "eta_f"]), tolerance = 1e-5)
  growth <- ((1 + draws[valid, "rhoMF"]) / (2 * draws[valid, "eta_m"]))^2
  expect_true(stats::sd(growth) > 0)
  corrected <- ravel_correct_pi(out$packet, pi = .4, pi_se = .02)
  expect_true(all(is.finite(corrected$se)))
  expect_true(all(corrected$inference == "bootstrap_covariance_pi_delta_wald"))
})
