test_that("child-free inputs and fixed-state semantics are explicit", {
  p <- family_validation_population()
  d <- family_validation_data(p, N = 1200L)
  d <- d[, !grepl("child", names(d)), drop = FALSE]
  expect_error(ravel_family_fit(d, c("mother", "father"),
    fit_constraints = "estimated_state", pi = p$pi), "not identified")
  f <- ravel_family_fit(d, c("mother", "father"),
    fit_constraints = "estimated_state", pi = p$pi,
    parent_visible_variance = p$pi + (1 - p$pi) * p$parent_variance)
  expect_true(f$converged)
  expect_equal(f$regression$coefficients$total_raw,
    family_validation_oracle(p, c("mother", "father"))$beta_raw, tolerance = 2e-5)
  expect_equal(as.data.frame(f), f$regression$table)
  expect_equal(f$outcome_variable, "Y")
  expect_error(ravel_family_fit(d, c("mother", "mat_sib"),
    fit_constraints = "equilibrium", parent_visible_variance = 1.1), "cannot override")
  expect_error(ravel_family_fit(d, c("mother", "mat_sib"),
    crossperson = "equal_pc"), "not applicable")
  expect_error(ravel_family_fit(d, c("mother", "mat_sib"),
    fit_constraints = "estimated_state", fixed_crossperson_corr = c(e1_m_ms = .2)), "benchmark")
  expect_error(ravel_family_fit(d, c("mother", "mat_sib"), .4), "distinct, nonempty names")
  expect_error(ravel_family_fit(d, c("mother", "mat_sib"),
    pi = .4, pi = .4), "distinct, nonempty names")
})

test_that("adult covariate and binary routes retain appropriate moments", {
  p <- family_validation_population()
  d <- family_validation_data(p, N = 1600L)
  set.seed(44)
  d$age <- rnorm(nrow(d))
  d$Y <- d$Y + .2 * d$age
  f <- ravel_family_fit(d, c("mother", "mat_sib"), outcome_roles = "mother",
    fit_constraints = "estimated_state", pi = .4, fixed_effects = "age",
    covariate_mapping_confirmed = TRUE)
  expect_true(f$converged)
  expect_true(all(is.finite(f$regression$table$estimate)))
  expect_true(all(is.finite(f$regression$table$se)))
  d$Y <- as.numeric(d$Y > 0)
  b <- ravel_family_fit(d, c("mother", "mat_sib"),
    fit_constraints = "estimated_state", pi = .4, fixed_effects = "age",
    outcome = "binary_liability", covariate_mapping_confirmed = TRUE)
  expect_true(b$converged)
  expect_true(all(is.finite(b$regression$table$estimate)))
  expect_equal(b$family_descriptor$moments(b$family_descriptor$theta, .4)$outcome_variance, 1)
  expect_error(ravel_family_fit(d, c("mother", "mat_sib"),
    outcome = "binary_liability", cluster = "FID"), "continuous outcome")
})

test_that("adult clustering and resampling use declared family units", {
  p <- family_validation_population()
  d <- family_validation_data(p, N = 800L)
  d$family <- rep(seq_len(nrow(d) / 2L), each = 2L)
  expect_warning(f <- ravel_family_fit(d, c("mother", "mat_sib"),
    fit_constraints = "estimated_state", pi = .4, cluster = "family",
    n_boot = 3L, seed = 77L), "auto uses joint lavaan")
  expect_true(f$converged)
  expect_equal(f$n_clusters, 400L)
  expect_equal(nrow(f$family_descriptor$bootstrap_draws), 3L)
  expect_true(any(stats::complete.cases(f$family_descriptor$bootstrap_draws)))
  expect_error(ravel_family_fit(d, c("mother", "mat_sib"),
    cluster = "family", cluster_var = "FID"), "nest the cluster")
})

test_that("adult interaction models preserve role and moderator labels", {
  p <- family_validation_population()
  d <- family_validation_data(p, N = 2400L)
  set.seed(412)
  d$age <- rnorm(nrow(d), 30, 2)
  d$Y <- d$Y + .08 * (d$T_mother - .5 * d$T_mat_sib) * (d$age - 30)
  f <- ravel_family_fit(d, c("mother", "mat_sib"),
    fit_constraints = "estimated_state", pi = .4, interactions = "age",
    interaction_center = c(age = 30), covariate_mapping_confirmed = TRUE)
  expect_true(f$converged)
  expect_true(all(c("mother:age", "mat_sib:age", "delta_m:age") %in%
    f$regression$table$parameter))
  expect_true(all(is.finite(f$regression$table$estimate)))
  expect_error(ravel_family_fit(d, c("mother", "mat_sib"), pi = .4,
    interactions = "age", interaction_roles = "mother"), "every genetic role")
})

test_that("adult finite-history pi uncertainty refits the genetic model", {
  p <- family_validation_population(parent_variance = 1.2475)
  d <- family_validation_data(p, N = 1600L)
  f <- ravel_family_fit(d, c("mother", "father", "mat_sib"),
    fit_constraints = "n_gen", am_n_generations = 3L, pi = .4, pi_se = .02)
  expect_true(f$converged)
  expect_true(f$family_descriptor$pi_dependent)
  expect_length(f$regression$uncertainty_errors, 0L)
  expect_true(all(is.finite(f$regression$table$se)))
  expect_equal(f$regression$coefficients$total_raw,
    family_validation_oracle(p, c("mother", "father", "mat_sib"))$beta_raw,
    tolerance = 2e-5)
})

test_that("a null adult score-only reliability direction is not certified", {
  p <- family_validation_population(rho = 0, parent_variance = 1,
    beta = setNames(rep(0, 5), c("child", "mother", "father", "mat_sib", "pat_sib")),
    lineage = c(0, 0))
  d <- family_validation_data(p, N = 1200L)
  f <- ravel_family_fit(d, c("mother", "mat_sib"),
    fit_constraints = "estimated_state", pi = .4)
  expect_false(f$converged)
  expect_match(f$error, "locally identified")
})
