pi_test_packet <- function() {
  labels <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f")
  theta <- setNames(c(.3, .1, .02, .1, .55, .55), labels)
  V <- matrix(.00002, 6L, 6L, dimnames = list(labels, labels))
  diag(V) <- .0001
  list(schema = "ravel_visible_transfer_v1", diagnostics = list(admissible = TRUE),
    settings = list(uncertainty = "lavaan", residual_model = "fixed"),
    outcome = list(name = "trait"), theta = theta, vcov_core = V,
    software = list(transform_fingerprint = RAVEL:::.report_transform_fingerprint()))
}

test_that("direct raw-scale pi SE matches equivalent CI input without inventing CI metadata", {
  packet <- pi_test_packet()
  for (regime in c("n_gen", "equilibrium")) {
    direct <- ravel_correct_pi(packet, .4, pi_se = .02, am_correction = regime, am_n_generations = 1)
    interval <- ravel_correct_pi(packet, .4, .4 + c(-1, 1) * qnorm(.975) * .02,
                                 am_correction = regime, am_n_generations = 1)
    for (column in c("estimate", "se", "ci_low", "ci_high", "p_value", "pi_se"))
      expect_equal(direct[[column]], interval[[column]], tolerance = 1e-12)
    expect_equal(attr(direct, "covariance"), attr(interval, "covariance"), tolerance = 1e-12)
    expect_equal(attr(direct, "joint_input_covariance"), attr(interval, "joint_input_covariance"), tolerance = 1e-12)
    metadata <- attr(direct, "pi_uncertainty")
    expect_identical(metadata$source, "se")
    expect_identical(metadata$raw_scale_se, .02)
    expect_null(metadata$ci)
    expect_null(metadata$ci_level)
    expect_null(metadata$ci_scale)
    expect_true(all(is.na(direct$pi_ci_scale)))
    expect_true(all(direct$pi_uncertainty_source == "se"))
    expect_identical(attr(interval, "pi_uncertainty")$source, "ci")
  }
})

test_that("direct pi SE supports zero uncertainty and conditional overlap covariance", {
  packet <- pi_test_packet()
  for (pi in c(0, .4)) {
    got <- ravel_correct_pi(packet, pi, pi_se = 0, am_correction = "equilibrium")
    expected <- ravel_correct_packet(packet, pi, am_correction = "equilibrium")
    expect_equal(got$estimate, expected$estimate, tolerance = 1e-10)
    expect_equal(got$se, expected$se, tolerance = 1e-10)
    expect_equal(attr(got, "joint_input_covariance")["pi", "pi"], 0)
  }
  cross <- setNames(rep(0, 6), names(packet$theta))
  cross["bD"] <- -3e-5
  got <- ravel_correct_pi(packet, .4, pi_se = .02, core_pi_cov = cross, am_n_generations = 1)
  ref <- ravel_correct_pi(packet, .4, .4 + c(-1, 1) * qnorm(.975) * .02,
                          core_pi_cov = cross, am_n_generations = 1)
  expect_equal(attr(got, "covariance"), attr(ref, "covariance"), tolerance = 1e-12)
  expect_equal(attr(got, "joint_input_covariance")[1:6, "pi"], cross)
  expect_false(attr(got, "pi_uncertainty")$independent_of_core)
})

test_that("direct pi SE validates its scale, exclusivity, and boundaries", {
  packet <- pi_test_packet()
  expect_error(ravel_correct_pi(packet, .4), "exactly one")
  expect_error(ravel_correct_pi(packet, .4, c(.3, .5), pi_se = .02), "exactly one")
  for (bad in list(-.01, NA_real_, NaN, Inf, c(.01, .02), "0.02", .01 + 1i, TRUE))
    expect_error(ravel_correct_pi(packet, .4, pi_se = bad), "pi_se must")
  expect_error(ravel_correct_pi(packet, 0, pi_se = .01), "boundary-aware")
  expect_error(ravel_correct_pi(packet, 1, pi_se = .01), "pi must")
  expect_error(ravel_correct_pi(packet, .4, pi_se = 1e200), "finite sampling variance")
  expect_error(ravel_correct_pi(packet, .4, pi_se = .02, pi_ci_scale = "logit"), "only to pi_ci")
  expect_error(ravel_correct_pi(packet, .4, pi_se = .02, pi_ci_scale = "normal"), "only to pi_ci")
  expect_error(ravel_correct_pi(packet, .4, pi_se = .02, pi_ci_level = .95), "only to pi_ci")
})

test_that("a degenerate external CI reproduces fixed-pi correction uncertainty", {
  packet <- pi_test_packet()
  for (h in c(0, .4)) for (regime in c("n_gen", "equilibrium")) {
    expected <- ravel_correct_packet(packet, h, regime, am_n_generations = 1L)
    got <- ravel_correct_pi(packet, h, rep(h, 2), am_correction = regime, am_n_generations = 1L)
    for (column in c("estimate", "se", "ci_low", "ci_high", "p_value"))
      expect_equal(got[[column]], expected[[column]], tolerance = 1e-8)
    expect_equal(unname(attr(got, "jacobian")[, "pi"]), rep(0, 8))
    expect_equal(dim(attr(got, "covariance")), c(8L, 8L))
    expect_identical(packet, pi_test_packet())
  }
})

test_that("independent pi uncertainty adds the full rank-one covariance", {
  packet <- pi_test_packet()
  fixed <- ravel_correct_pi(packet, .4, c(.4, .4), am_correction = "equilibrium")
  got <- ravel_correct_pi(packet, .4, c(.3, .5), am_correction = "equilibrium")
  expected_se <- .1 / qnorm(.975)
  expect_equal(unique(got$pi_se), expected_se)
  J <- attr(got, "jacobian")
  increment <- tcrossprod(J[, "pi"]) * expected_se^2
  expect_equal(unname(attr(got, "covariance") - attr(fixed, "covariance")), unname(increment),
               tolerance = 1e-10)
  expect_true(all(got$se >= fixed$se))
  expect_gt(got$se[1], fixed$se[1])
  expect_equal(Reduce(`+`, attr(got, "variance_components")), attr(got, "covariance"),
               ignore_attr = TRUE, tolerance = 1e-12)
  expect_true(attr(got, "pi_uncertainty")$independent_of_core)
})

test_that("normal and logit-Wald CI conventions agree at first order", {
  packet <- pi_test_packet()
  h <- .4
  raw_se <- .025
  logit_se <- raw_se / (h * (1 - h))
  normal_ci <- h + c(-1, 1) * qnorm(.975) * raw_se
  logit_ci <- plogis(qlogis(h) + c(-1, 1) * qnorm(.975) * logit_se)
  a <- ravel_correct_pi(packet, h, normal_ci, am_n_generations = 1)
  b <- ravel_correct_pi(packet, h, logit_ci, pi_ci_scale = "logit", am_n_generations = 1)
  expect_equal(a$pi_se, b$pi_se, tolerance = 1e-12)
  expect_equal(a$se, b$se, tolerance = 1e-10)
  expect_false(attr(b, "pi_uncertainty")$boundary_truncation)
  expect_error(ravel_correct_pi(packet, .4, c(.3, .6)), "symmetric")
  expect_error(ravel_correct_pi(packet, .4, c(.3, .5), pi_ci_scale = "logit"), "symmetric")
  expect_error(ravel_correct_pi(packet, 0, c(0, .1)), "boundary-aware")
  expect_error(ravel_correct_pi(packet, .4, c(0, .8), pi_ci_scale = "logit"), "strictly")
  expect_error(ravel_correct_pi(packet, 1, c(.9, 1)), "pi must")
  expect_error(ravel_correct_pi(packet, .4, c(.5, .6)), "enclosing")
  expect_error(ravel_correct_pi(packet, .4, c(.3, .5), am_n_generations = 1.5), "integer")
  # Normal-Wald bounds are a sampling approximation, not a truncated prior.
  crossed <- ravel_correct_pi(packet, .05, c(-.03, .13), am_n_generations = 1)
  expect_equal(unique(crossed$pi_se), .08 / qnorm(.975))
})

test_that("external overlap covariance is retained and validated", {
  packet <- pi_test_packet()
  cross <- setNames(rep(0, 6), names(packet$theta))
  cross["bD"] <- -.3 * sqrt(packet$vcov_core[1, 1]) * .02
  ci <- .4 + c(-1, 1) * qnorm(.975) * .02
  got <- ravel_correct_pi(packet, .4, ci, core_pi_cov = cross, am_n_generations = 1)
  J <- attr(got, "jacobian")
  V <- attr(got, "joint_input_covariance")
  expect_equal(V[1:6, "pi"], cross)
  expect_equal(attr(got, "covariance"), J %*% V %*% t(J), tolerance = 1e-12)
  expect_false(attr(got, "pi_uncertainty")$independent_of_core)
  expect_error(ravel_correct_pi(packet, .4, ci, core_pi_cov = unname(cross)), "named")
  cross["bD"] <- 1
  expect_error(ravel_correct_pi(packet, .4, ci, core_pi_cov = cross), "positive semidefinite")
  packet$vcov_core[1, 2] <- .01
  expect_error(ravel_correct_pi(packet, .4, ci), "symmetric")
})

test_that("bootstrap covariance uses complete joint draws and honest labels", {
  packet <- pi_test_packet()
  set.seed(8903)
  packet$bootstrap <- list(core_draws = sweep(matrix(rnorm(1200), 200, 6) %*%
    chol(packet$vcov_core), 2, packet$theta, "+"), minimum_success = 160L)
  colnames(packet$bootstrap$core_draws) <- names(packet$theta)
  packet$settings$uncertainty <- "bootstrap"
  got <- ravel_correct_pi(packet, .4, c(.3, .5), am_n_generations = 1)
  expected <- packet
  expected$vcov_core <- cov(packet$bootstrap$core_draws)
  model <- ravel_correct_pi(expected, .4, c(.3, .5), am_n_generations = 1, uncertainty = "lavaan")
  expect_equal(got$se, model$se)
  expect_true(all(got$inference == "bootstrap_covariance_pi_delta_wald"))
  expect_true(all(got$n_boot_success == 200L))
  packet$bootstrap$core_draws[1:41, 1] <- NA_real_
  expect_error(ravel_correct_pi(packet, .4, c(.3, .5)), "Insufficient")
})

test_that("pi-dependent finite fits use the complete refitting derivative", {
  packet <- pi_test_packet()
  packet$settings$constraint_spec <- list(mode = "n_gen", h = .4, n_generations = 3L)
  sensitivity <- setNames(c(.01, -.02, .03, .04, .05, .05), names(packet$theta))
  expect_error(ravel_correct_pi(packet, .4, c(.3, .5), am_n_generations = 3), "core_pi_sensitivity")
  expect_error(ravel_correct_pi(packet, .3, c(.2, .4), am_n_generations = 3), "Refit")
  expect_error(ravel_correct_pi(packet, .4, c(.3, .5), am_correction = "equilibrium"), "mating history")
  got <- ravel_correct_pi(packet, .4, c(.3, .5), am_n_generations = 3,
                          core_pi_sensitivity = sensitivity)
  inherited <- ravel_correct_pi(packet, .4, c(.3, .5), core_pi_sensitivity = sensitivity)
  expect_equal(inherited, got)
  direct <- ravel_correct_pi(packet, .4, pi_se = .1 / qnorm(.975), core_pi_sensitivity = sensitivity)
  expect_equal(attr(direct, "covariance"), attr(got, "covariance"), tolerance = 1e-10)
  J <- attr(got, "jacobian")
  partial <- attr(got, "partial_jacobian")
  expect_equal(J[, "pi"], partial[, "pi"] + drop(partial[, 1:6] %*% sensitivity))
  evaluate <- function(delta) {
    adjusted <- pi_test_packet()
    adjusted$theta <- adjusted$theta + delta * sensitivity
    ravel_correct_pi(adjusted, .4 + delta, rep(.4 + delta, 2), am_n_generations = 3)$estimate
  }
  reference <- (evaluate(1e-5) - evaluate(-1e-5)) / 2e-5
  expect_equal(unname(J[, "pi"]), reference, tolerance = 2e-5)
})

test_that("joint delta SEs calibrate in a small independent moment simulation", {
  packet <- pi_test_packet()
  # This tests post-estimation propagation, with known moment-level covariance;
  # SEM fitting and its covariance calibration have separate integration tests.
  packet$vcov_core <- packet$vcov_core / 100
  h <- .4
  h_se <- .002
  got <- ravel_correct_pi(packet, h, h + c(-1, 1) * qnorm(.975) * h_se, am_n_generations = 1)
  set.seed(581033)
  B <- 1600L
  core <- sweep(matrix(rnorm(B * 6), B, 6) %*% chol(packet$vcov_core), 2, packet$theta, "+")
  hs <- rnorm(B, h, h_se)
  draws <- t(vapply(seq_len(B), function(i) {
    x <- core[i, ]
    value <- RAVEL:::.compute_scale_estimates(x[1:3], x[5], x[6], x[4], hs[i], "n_gen", 1L)
    c(value$beta_total_assumed_std, value$rho_total_assumed,
      value$beta_total_assumed, value$cov_total_assumed)
  }, numeric(8)))
  ratio <- apply(draws, 2, sd) / got$se
  expect_true(all(ratio > .92 & ratio < 1.08))
  covered <- colMeans(abs(sweep(draws, 2, got$estimate)) <=
                        matrix(qnorm(.975) * got$se, B, 8L, byrow = TRUE))
  expect_true(all(covered > .93 & covered < .97))
})

test_that("live legacy results retain all joint core uncertainty", {
  fixture <- cov_test_fixture("continuous", correlated = FALSE, with_covariates = FALSE, n = 5000L)
  fitted <- cov_test_fit(fixture, h = .4, regime = "equilibrium")
  expect_true(fitted$converged)
  corrected <- ravel_correct_pi(fitted, .4, c(.3, .5), am_correction = "equilibrium")
  expect_equal(ravel_correct_pi(fitted, .4, c(.3, .5)), corrected)
  direct <- ravel_correct_pi(fitted, .4, pi_se = .1 / qnorm(.975))
  expect_equal(direct$se, corrected$se, tolerance = 1e-10)
  expect_equal(corrected$estimate[1:3], unname(fitted$beta_total_assumed_std), tolerance = 1e-8)
  expect_true(all(corrected$se[1:3] > fitted$lavaan_inference$total_std$se[1:3]))
  bad <- fitted
  bad$converged <- FALSE
  expect_error(ravel_correct_pi(bad, .4, c(.3, .5)), "converged")
})

test_that("direct pi SE includes automatic finite-generation refitting sensitivity", {
  # This integration check compares uncertainty interfaces on the same fit;
  # population-model recovery is checked in the constrained-fitting tests.
  fixture <- cov_test_fixture("continuous", correlated = FALSE, with_covariates = FALSE, n = 5000L)
  fitted <- ravel_fit(fixture$dat, h = .4, fit_constraints = "n_gen", am_n_generations = 3L)
  expect_true(fitted$converged)
  direct <- ravel_correct_pi(fitted, .4, pi_se = .02)
  sensitivity <- attr(direct, "pi_uncertainty")$core_pi_sensitivity
  expect_gt(max(abs(sensitivity)), 1e-5)
  interval <- ravel_correct_pi(fitted, .4, .4 + c(-1, 1) * qnorm(.975) * .02,
                               core_pi_sensitivity = sensitivity)
  expect_equal(direct$estimate, interval$estimate, tolerance = 1e-12)
  expect_equal(direct$se, interval$se, tolerance = 1e-10)
  expect_equal(attr(direct, "covariance"), attr(interval, "covariance"), tolerance = 1e-10)
  expect_gt(max(abs(attr(direct, "jacobian")[, "pi"] -
                     attr(direct, "partial_jacobian")[, "pi"])), 1e-5)
})

test_that("live bootstrap core draws compose with pi correction without changing legacy summaries", {
  fixture <- cov_test_fixture("continuous", correlated = FALSE, with_covariates = FALSE, n = 5000L)
  fitted <- cov_test_fit(fixture, h = .4, regime = "equilibrium", n_boot = 12L)
  expect_true(fitted$converged)
  expect_equal(dim(fitted$bootstrap$core_draws), c(12L, 6L))
  expect_equal(fitted$bootstrap$minimum_success, 10L)
  draws <- fitted$bootstrap$core_draws
  complete <- apply(is.finite(draws), 1L, all)
  expect_gte(sum(complete), 10L)
  before <- ravel_to_dataframe(fitted, uncertainty = "bootstrap", include_raw = TRUE,
                               include_vanilla = FALSE)
  original_bootstrap <- fitted$bootstrap
  explicit <- ravel_correct_pi(fitted, .4, c(.3, .5), uncertainty = "bootstrap")
  automatic <- ravel_correct_pi(fitted, .4, c(.3, .5))
  expect_equal(automatic, explicit)
  direct <- ravel_correct_pi(fitted, .4, pi_se = .1 / qnorm(.975))
  expect_equal(direct$se, explicit$se, tolerance = 1e-10)
  expect_equal(attr(direct, "covariance"), attr(explicit, "covariance"), tolerance = 1e-10)
  expect_true(all(direct$inference == "bootstrap_covariance_pi_delta_wald"))
  expect_true(all(explicit$inference == "bootstrap_covariance_pi_delta_wald"))
  expect_true(all(explicit$n_boot_success == sum(complete)))
  expect_equal(attr(explicit, "joint_input_covariance")[1:6, 1:6],
               cov(draws[complete, , drop = FALSE]))
  J <- attr(explicit, "jacobian")
  expect_equal(attr(explicit, "covariance"),
               J %*% attr(explicit, "joint_input_covariance") %*% t(J), tolerance = 1e-12)
  expect_identical(fitted$bootstrap, original_bootstrap)
  expect_equal(ravel_to_dataframe(fitted, uncertainty = "bootstrap", include_raw = TRUE,
                                  include_vanilla = FALSE), before)

  # Removing the newly stored information gives the previous bootstrap object
  # schema. All pre-existing extraction and transformed-draw summaries agree.
  legacy <- fitted
  legacy$bootstrap$core_draws <- NULL
  legacy$bootstrap$minimum_success <- NULL
  expect_equal(ravel_to_dataframe(legacy, uncertainty = "bootstrap", include_raw = TRUE,
                                  include_vanilla = FALSE), before)
  transformed <- t(apply(draws[complete, , drop = FALSE], 1L, function(x) {
    value <- RAVEL:::.compute_scale_estimates(x[c("bD", "bm", "bf")],
      x["eta_m"], x["eta_f"], x["rhoMF"], .4, "equilibrium", 1L)
    c(value$beta_total_assumed_std, value$rho_total_assumed,
      value$beta_total_assumed, value$cov_total_assumed)
  }))
  expect_equal(unname(fitted$bootstrap$total_std$se),
               unname(apply(transformed[, 1:4, drop = FALSE], 2L, sd)), tolerance = 1e-10)
  expect_equal(unname(fitted$bootstrap$total_raw$se),
               unname(apply(transformed[, 5:8, drop = FALSE], 2L, sd)), tolerance = 1e-10)
  expected_ci <- apply(transformed[, 1:4, drop = FALSE], 2L, quantile,
                       probs = c(.025, .975), names = FALSE)
  expect_equal(unname(fitted$bootstrap$total_std$ci), unname(expected_ci), tolerance = 1e-10)
})
