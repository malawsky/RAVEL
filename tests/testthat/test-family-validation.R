test_that("independent family fixtures match their specified likelihood moments", {
  p <- family_validation_population()
  d <- family_validation_data(p, N = 800L)
  x <- as.matrix(d[, rownames(p$observed_map)])
  expect_equal(crossprod(x) / nrow(x), p$observed_cov, tolerance = 1e-10)
  full <- family_validation_oracle(p, p$roles)$beta_raw
  reduced <- family_validation_oracle(p, c("mother", "mat_sib"))$beta_raw
  # This fixture distinguishes marginal regression from slicing full-model
  # coefficients: omitted child and partner terms are deliberately nonzero.
  expect_gt(abs(full[["mother"]] - reduced[["mother"]]), .02)
})

test_that("trio-assisted adult projections recover all eight oracle regressions", {
  p <- family_validation_population()
  d <- family_validation_data(p)
  fit <- ravel_family_fit(d, genetic_roles = p$roles,
    outcome_roles = c("mother", "mat_sib"), pi = p$pi,
    fit_constraints = "estimated_state")
  expect_true(fit$converged)
  for (roles in family_validation_subsets()) {
    got <- ravel_family_regression(fit, outcome_roles = roles)
    truth <- family_validation_oracle(p, roles)
    expect_equal(got$coefficients$total_raw, truth$beta_raw,
                 tolerance = 2e-4, info = paste(roles, collapse = "+"))
    expect_equal(got$coefficients$total_std, truth$beta_std,
                 tolerance = 2e-4, info = paste(roles, collapse = "+"))
    expect_equal(unname(got$variance_explained$total_raw["genetic_variance"]),
                 truth$variance_explained, tolerance = 2e-4)
    expect_equal(unname(got$variance_explained$total_raw["genetic_R2"]),
                 truth$R2, tolerance = 2e-4)
    if (length(truth$contrasts_raw))
      expect_equal(got$contrasts$total_raw[names(truth$contrasts_raw)],
                   truth$contrasts_raw, tolerance = 2e-4)
  }
  roles <- c("mother", "father", "mat_sib")
  a <- ravel_family_regression(fit, outcome_roles = roles)
  b <- ravel_family_regression(fit, outcome_roles = rev(roles))
  expect_equal(a$coefficients$total_raw,
               b$coefficients$total_raw[names(a$coefficients$total_raw)],
               tolerance = 1e-9)
  V <- a$vcov
  expect_equal(V["total_raw:delta_m", "total_raw:delta_m"],
    V["total_raw:mother", "total_raw:mother"] +
      V["total_raw:mat_sib", "total_raw:mat_sib"] -
      2 * V["total_raw:mother", "total_raw:mat_sib"], tolerance = 1e-9)
  uncertain <- ravel_family_regression(fit, outcome_roles = roles, pi_se = .03)
  expect_equal(uncertain$coefficients$total_raw, a$coefficients$total_raw)
  visible <- grep("^visible:", rownames(V), value = TRUE)
  corrected <- grep("^total_raw:", rownames(V), value = TRUE)
  expect_equal(uncertain$vcov[visible, visible], V[visible, visible], tolerance = 1e-9)
  added <- diag(uncertain$vcov[corrected, corrected]) - diag(V[corrected, corrected])
  expect_true(all(added >= -1e-10))
  expect_gt(max(added), 1e-8)
})

test_that("trio-assisted projections recover zero-pi and changed-assortment states", {
  scenarios <- list(c(pi = 0, rho = .3, V = 1.15),
                    c(pi = .7, rho = -.2, V = .85),
                    c(pi = .2, rho = .1, V = 1.6))
  roles <- c("mother", "mat_sib", "father", "pat_sib")
  for (s in scenarios) {
    p <- family_validation_population(pi = s[["pi"]], rho = s[["rho"]],
                                      parent_variance = s[["V"]])
    d <- family_validation_data(p)
    fit <- ravel_family_fit(d, genetic_roles = p$roles, outcome_roles = roles,
                           pi = p$pi, fit_constraints = "estimated_state")
    expect_true(fit$converged)
    z <- ravel_family_regression(fit, outcome_roles = roles)
    expect_equal(z$coefficients$total_raw,
                 family_validation_oracle(p, roles)$beta_raw, tolerance = 3e-4)
  }
})

test_that("child-free sibling and partner configurations recover oracle projections", {
  p <- family_validation_population()
  d <- family_validation_data(p)
  # The child scores truly are unavailable, not fabricated or set to missing.
  d[c("PGS1_child", "PGS2_child")] <- NULL
  configurations <- list(c("mother", "mat_sib"),
                         c("mother", "father", "mat_sib"),
                         c("mother", "father", "pat_sib"),
                         c("mother", "father", "mat_sib", "pat_sib"))
  for (roles in configurations) {
    fit <- ravel_family_fit(d, genetic_roles = roles, outcome_roles = roles,
                           pi = p$pi, fit_constraints = "estimated_state")
    expect_true(fit$converged, info = paste(roles, collapse = "+"))
    z <- ravel_family_regression(fit, outcome_roles = roles)
    expect_equal(z$coefficients$total_raw,
                 family_validation_oracle(p, roles)$beta_raw, tolerance = 3e-4,
                 info = paste(roles, collapse = "+"))
  }
})

test_that("trio-only and one-sibling adapters recover independent adult projections", {
  for (sibling in list(character(0), "mat_sib", "pat_sib")) {
    genetic <- c("child", "mother", "father", sibling)
    beta <- c(child = .12, mother = .25, father = -.08,
              mat_sib = .03, pat_sib = -.02)
    beta[setdiff(c("mat_sib", "pat_sib"), sibling)] <- 0
    lineage <- c(maternal = if ("mat_sib" %in% sibling) .2 else 0,
                 paternal = if ("pat_sib" %in% sibling) -.1 else 0)
    p <- family_validation_population(beta = beta, lineage = lineage)
    d <- family_validation_data(p)
    unavailable <- setdiff(p$roles, genetic)
    for (j in 1:2) d[paste0("PGS", j, "_", unavailable)] <- NULL
    selected <- setdiff(genetic, "child")
    fit <- ravel_family_fit(d, genetic_roles = genetic, outcome_roles = selected,
                           pi = p$pi, fit_constraints = "estimated_state")
    expect_true(fit$converged)
    got <- ravel_family_regression(fit, selected)
    expect_equal(got$coefficients$total_raw,
                 family_validation_oracle(p, selected)$beta_raw,
                 tolerance = 3e-4, info = paste(genetic, collapse = "+"))
  }
})
