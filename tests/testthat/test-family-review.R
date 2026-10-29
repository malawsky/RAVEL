test_that("family bootstrap rejects a single or missing resampling unit", {
  p <- family_validation_population()
  d <- family_validation_data(p, N = 400L)
  d$FID <- 0L
  for (roles in list(c("mother", "mat_sib"), c("child", "mother", "father"))) {
    expect_error(ravel_family_fit(d, roles, fit_constraints = "estimated_state",
      pi = p$pi, n_boot = 10L), "at least two|at least 2|two.*units|2.*units")
  }
  d$FID <- seq_len(nrow(d))
  d$FID[1L] <- NA_integer_
  expect_error(ravel_family_fit(d, c("mother", "mat_sib"),
    fit_constraints = "estimated_state", pi = p$pi, n_boot = 10L),
    "missing|observed|complete")
  expect_error(ravel_family_fit(d, c("mother", "mat_sib"),
    fit_constraints = "estimated_state", pi = p$pi, n_boot = 10L,
    cluster_var = "absent_id"), "absent|not found|Missing|missing")
})

test_that("explicit row bootstrap bypasses a dummy FID without zero uncertainty", {
  p <- family_validation_population()
  d <- family_validation_data(p, N = 600L)
  d$FID <- 0L
  f <- ravel_family_fit(d, c("mother", "mat_sib"),
    fit_constraints = "estimated_state", pi = p$pi,
    cluster_var = NULL, n_boot = 10L, seed = 71L)
  expect_true(f$converged)
  target <- f$regression$table
  target <- target[target$kind == "coefficient", , drop = FALSE]
  expect_true(all(is.finite(target$se)))
  expect_true(all(target$se > 0))
  expect_true(all(target$ci_high > target$ci_low))
})

test_that("legacy fixed parent paths remain usable in complete-trio family fits", {
  p <- family_validation_population(pi = 0,
    beta = c(child = .3, mother = 0, father = 0, mat_sib = 0, pat_sib = 0),
    lineage = c(0, 0))
  d <- family_validation_data(p, N = 600L)
  f <- ravel_family_fit(d, c("child", "mother", "father"),
    fit_constraints = "estimated_state", pi = 0, fix_parent_paths = TRUE)
  expect_true(f$converged)
  expect_equal(unname(f$regression$coefficients$visible[c("mother", "father")]),
    c(0, 0), tolerance = 1e-10)
  expect_true(all(is.finite(f$regression$table$estimate)))
  expect_equal(f$regression$coefficients$total_raw[["child"]], .3,
    tolerance = 2e-5)
})

test_that("two corrected interaction blocks retain role order and moderator units", {
  roles <- c("mother", "father")
  pi <- .4
  G <- matrix(c(1, .3, .3, 1), 2, dimnames = list(roles, roles))
  T <- (1 - pi) * G + diag(pi, 2)
  means <- c(age = 1.3, cohort = -.7)
  W <- matrix(c(1, .25, .25, 1.5), 2, dimnames = list(names(means), names(means)))
  Z <- rbind(c(1, means), cbind(means, W + tcrossprod(means)))
  canonical <- c(roles, paste(roles, "age", sep = ":"), paste(roles, "cohort", sep = ":"))
  order <- c(1, 2, 4, 3, 6, 5)
  dn <- canonical[order]
  SV <- kronecker(Z, G)[order, order]
  ST <- kronecker(Z, T)[order, order]
  dimnames(SV) <- dimnames(ST) <- list(dn, dn)
  raw <- setNames(c(.12, -.05, .04, .08, -.02, .03), canonical)[dn]
  M <- diag(sqrt(1 - pi), 2)
  dimnames(M) <- list(roles, roles)
  m <- list(Sigma_GG = G, Sigma_TT = T, M = M, beta_total_raw = raw[roles],
    beta_visible = sqrt(1 - pi) * raw[roles],
    design_cov_visible = SV, design_beta_visible = sqrt(1 - pi) * raw,
    design_covariate_cov_visible = matrix(0, 6, 2),
    design_names = dn, design_roles = sub(":.*", "", dn),
    design_moderators = c("", "", "age", "age", "cohort", "cohort"),
    design_units = c(1, 1, 2, 2, 5, 5),
    moderator_means = means, moderator_covariance = W, Sigma_WW = W,
    outcome_variance = 2, outcome_covariate_cov = c(0, 0))
  out <- .ravel_family_projection_values(m, roles, "mother", pi, TRUE)$value
  idx <- c("mother", "mother:age", "mother:cohort")
  target <- drop(solve(ST[idx, idx], ST[idx, ] %*% raw)) / c(1, 2, 5)
  expect_equal(unname(out[paste0("total_raw:", idx)]), unname(target), tolerance = 1e-12)
  expect_equal(unname(out[paste0("total_std:", idx)]), unname(target), tolerance = 1e-12)
})
