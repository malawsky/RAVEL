test_that("mixed coefficients and residual variances preserve input outcome units", {
  skip_on_cran()
  d <- mixed_validation_data(N_families = 500L, sib_sizes = 2L, visits = 1L,
    pi = .4, family_variance = .5, age_effect = 0, sex_effect = 0)
  base <- ravel_mixed_fit(d, pi = .4, diagnostic = FALSE)
  changed <- d
  score_names <- c("PGS1_mother", "PGS2_mother", "PGS1_father", "PGS2_father",
    "PGS1_child", "PGS2_child")
  factors <- c(1e-5, 1e3, .1, 17, .03, 12)
  for (j in seq_along(score_names)) changed[[score_names[j]]] <-
    d[[score_names[j]]] * factors[j] + 3 * factors[j]
  changed$Y <- 7 * d$Y + 13
  scaled <- ravel_mixed_fit(changed, pi = .4, diagnostic = FALSE)
  expect_true(base$converged && scaled$converged)
  expect_identical(base$inference_status, "ok")
  expect_identical(scaled$inference_status, "ok")
  expect_equal(scaled$table$estimate, base$table$estimate * 7, tolerance = 2e-4)
  expect_equal(scaled$table$se, base$table$se * 7, tolerance = 2e-4)
  expect_equal(scaled$variance_components$estimate,
    base$variance_components$estimate * 49, tolerance = 2e-4)
  expect_equal(scaled$variance_components$se,
    base$variance_components$se * 49, tolerance = 2e-4)
  expect_equal(scaled$state, base$state, tolerance = 2e-4)
  expect_equal(unname(scaled$fit_statistics), unname(base$fit_statistics), tolerance = 2e-5)
  a <- ravel_mixed_audit(base)
  b <- ravel_mixed_audit(scaled)
  multiplier <- c(PGS1_mother = 1e-5, PGS1_father = .1, PGS1_child1 = .03,
    PGS1_child2 = .03, PGS2_mother = 1e3, PGS2_father = 17,
    PGS2_child1 = 12, PGS2_child2 = 12, Y1 = 7, Y2 = 7)
  transform <- with(a$implied_covariances, unname(multiplier[lhs] * multiplier[rhs]))
  expect_equal(b$implied_covariances$estimate,
    a$implied_covariances$estimate * transform, tolerance = 3e-4)
  expect_equal(b$implied_covariances$se,
    a$implied_covariances$se * transform, tolerance = 3e-4)
  expect_identical(nrow(a$implied_covariances), 55L)
  expect_identical(a$implied_outcome_design$child, c("child1", "child2"))
  expect_identical(a$implied_status$status, "ok")
  expect_identical(ravel_mixed_audit(base, design_group = NULL)$implied_status$status, "not_extracted")
  scores_only <- base
  scores_only$prepared$groups <- scores_only$prepared$groups[1L]
  scores_only$prepared$groups[[1L]]$child <- integer()
  scores_only$prepared$groups[[1L]]$moderators <- matrix(numeric(), 0L, 1L)
  scores_only$prepared$groups[[1L]]$slope <- numeric()
  only <- ravel_mixed_audit(scores_only)
  expect_identical(nrow(only$implied_covariances), 36L)
  expect_identical(nrow(only$implied_outcome_design), 0L)
})
