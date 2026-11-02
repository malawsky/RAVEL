test_that("mixed audit exposes covariance estimates and conditional uncertainty", {
  spec <- list(mode = "estimated_state", pi = .4, generations = 1L)
  p <- c(state_rho = atanh(.3), state_variance = log(.8),
    loading1_mother = 0, loading1_father = 0, loading1_child = 0,
    loading2_mother = 0, loading2_father = 0, loading2_child = 0,
    noise1 = 0, noise2 = 0)
  V <- diag(.001, length(p)); dimnames(V) <- list(names(p), names(p))
  fit <- structure(list(parameters = p, settings = spec, vcov = V,
    means = numeric(), model_vcov = V, fit_statistics = c(logLik = -100, AIC = 220),
    N = 200, n_subjects = 200, n_families = 100, n_clusters = 90,
    converged = TRUE, optimizer_converged = TRUE, inference_status = "ok",
    table = data.frame(estimate = .2), variance_components = data.frame(estimate = .3)),
    class = "ravel_mixed_results")
  a <- ravel_mixed_audit(fit)
  covs <- a$genetic_covariances
  expect_equal(covs$estimate[covs$block == "total_founder_covariance" &
    covs$lhs == "mother" & covs$rhs == "father"], .3 * 1.3)
  expect_equal(covs$estimate[covs$block == "total_founder_covariance" &
    covs$lhs == "child1" & covs$rhs == "child2"], .5 * 1.3 * 1.3)
  expect_equal(covs$estimate[covs$block == "visible_founder_covariance" &
    covs$lhs == "child1" & covs$rhs == "child2"], .6 * .5 * 1.3 * 1.3 + .4 * .5)
  expect_equal(a$state$se[a$state$lhs == "rho"], (1 - .3^2) * sqrt(.001), tolerance = 1e-7)
  expect_true(all(is.finite(covs$se)))
  expect_true(all(grepl("conditional on supplied pi", covs$uncertainty_conditions)))
  expect_identical(a$coefficients, fit$table)
  expect_identical(a$parameter_vcov, V)
  expect_identical(dim(a$state_and_genetic_vcov), c(47L, 47L))
  expect_identical(a$reference_n_children, 2L)
  fit$vcov <- NULL; fit$inference_status <- "boundary"
  expect_true(all(is.na(ravel_mixed_audit(fit)$genetic_covariances$se)))
  expect_error(ravel_mixed_audit(fit, n_children = 0), "integer from 1 to 20")
  expect_error(ravel_mixed_audit(list()), "ravel_mixed_fit")
})
