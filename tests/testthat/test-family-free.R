test_that("least-restricted family measurement fits retain reduced estimands", {
  # Two-generation constant assortment has parental V=1+rho/2=1.15.
  # Exact population fixtures exercise free loadings and residual covariances,
  # rather than relying on constrained inheritance to recover measurement scales.
  p <- family_validation_population(pi = .4, rho = .3, parent_variance = 1.15)
  dat <- family_validation_data(p, N = 1200L, seed = 772L)
  configurations <- list(
    full = p$roles,
    trio = c("child", "mother", "father"),
    siblings = c("mother", "mat_sib"),
    adults = c("mother", "father", "mat_sib", "pat_sib"))
  for (roles in configurations) {
    target <- intersect(c("mother", "mat_sib"), roles)
    fit <- ravel_family_fit(dat, genetic_roles = roles, outcome_roles = target,
      fit_constraints = "none", crossperson = "free", pi = .4,
      am_correction = "n_gen", am_n_generations = 2L)
    expect_true(fit$converged)
    if (!isTRUE(fit$converged)) next
    oracle <- family_validation_oracle(p, target)
    expect_equal(fit$regression$coefficients$total_raw, oracle$beta_raw,
      tolerance = 2e-4)
    rows <- fit$regression$table$scale == "total_raw" &
      fit$regression$table$kind == "coefficient"
    expect_true(all(is.finite(fit$regression$table$se[rows])))
  }
})
