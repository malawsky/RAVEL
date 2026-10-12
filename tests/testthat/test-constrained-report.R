test_that("dataframe constraint options are validated before fitting", {
  expect_error(ravel_dataframe_report(data.frame(), fit_constraints = "equilibrium",
    residual_models = "free"), "requires residual_models")
  expect_error(ravel_dataframe_report(data.frame(), fit_constraints = "n_gen",
    h = 1), "h must")
  expect_error(ravel_dataframe_report(data.frame(), fit_constraints = "n_gen",
    am_n_generations = 1.5), "positive integer")
})

test_that("constrained reports preserve visible outputs and correction assumptions", {
  data("ravel_example", package = "RAVEL", envir = environment())
  for (mode in c("equilibrium", "n_gen")) {
    out <- ravel_dataframe_report(ravel_example, fit_constraints = mode,
      h = .4, am_n_generations = 3L, make_plot = FALSE)
    expect_identical(out$model_status$residual_model, "fixed")
    expect_identical(out$model_status$status, "ok")
    expect_true(all(out$summary$scale[out$summary$regression == "RAVEL"] == "visible"))
    packet <- out$packet
    expect_identical(packet$settings$constraint_spec$mode, mode)
    expect_identical(packet$settings$constraint_spec$depends_on_pi, mode == "n_gen")
    expect_null(packet$constraint_refit)
    corrected <- ravel_correct_packet(packet, h = .4,
      am_correction = mode, am_n_generations = 3L)
    expect_true(all(is.finite(corrected$se)))
    expect_equal(ravel_correct_packet(packet, h = .4), corrected)
    expect_error(ravel_correct_packet(packet, h = .4,
      am_correction = if (mode == "n_gen") "equilibrium" else "n_gen"), "regime")
    if (mode == "n_gen") {
      expect_error(ravel_correct_packet(packet, h = .5,
        am_correction = mode, am_n_generations = 3L), "refit")
      expect_error(ravel_correct_pi(packet, pi = .4, pi_ci = c(.3, .5)), "sensitivity")
    } else {
      with_pi <- ravel_correct_pi(packet, pi = .4, pi_ci = c(.3, .5))
      expect_true(all(is.finite(with_pi$se)))
      expect_true(all(with_pi$am_correction == "equilibrium"))
    }
  }
})

test_that("whole-trio report bootstrap retains the constrained model", {
  data("ravel_example", package = "RAVEL", envir = environment())
  out <- ravel_dataframe_report(ravel_example, fit_constraints = "equilibrium",
    h = .4, n_boot = 10L, seed = 733L, make_plot = FALSE)
  expect_identical(out$model_status$status, "ok")
  draws <- out$packet$bootstrap$core_draws
  expect_equal(dim(draws), c(10L, 6L))
  expect_true(all(is.finite(draws)))
  expect_equal(unname(draws[, "eta_m"]), unname((1 + draws[, "rhoMF"]) / 2),
    tolerance = 1e-5)
  expect_equal(unname(draws[, "eta_m"]), unname(draws[, "eta_f"]), tolerance = 1e-5)
  corrected <- ravel_correct_pi(out$packet, pi = .4, pi_ci = c(.3, .5))
  expect_true(all(is.finite(corrected$se)))
  expect_true(all(corrected$inference == "bootstrap_covariance_pi_delta_wald"))
})
