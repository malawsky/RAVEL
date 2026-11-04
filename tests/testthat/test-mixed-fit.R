# Numerical likelihood checks independent of the fitted parameter reporting.
.mixed_central_gradient <- function(fun, x, step = 1e-5) {
  vapply(seq_along(x), function(j) {
    h <- step * max(1, abs(x[j])); a <- b <- x
    a[j] <- a[j] + h; b[j] <- b[j] - h
    (fun(a) - fun(b)) / (2 * h)
  }, numeric(1))
}

.mixed_preparation_fixture <- function(dat, outcome = "continuous", fixed = NULL) {
  .ravel_mixed_prepare(dat, outcome, "Y",
    c(child1 = "PGS1_child", child2 = "PGS2_child", mother1 = "PGS1_mother",
      mother2 = "PGS2_mother", father1 = "PGS1_father", father2 = "PGS2_father"),
    "M_ID", "IID", "FID", fixed, NULL, NULL, TRUE, FALSE, NULL,
    "estimated_state", 1L, .4, TRUE)
}

test_that("probit orthant derivatives agree with direct probability differences", {
  skip_if_not_installed("mvtnorm")
  for (m in 1:4) {
    mu <- c(-.4, .6, -.2, .3)[seq_len(m)]
    s <- c(.8, 1.2, .9, 1.1)[seq_len(m)]
    R <- matrix(.2, m, m); diag(R) <- 1
    S <- R * outer(s, s); y <- c(0, 1, 1, 0)[seq_len(m)]
    answer <- .ravel_mixed_probit_probability(mu, S, y)
    direct <- .mixed_central_gradient(function(x)
      .ravel_mixed_probit_probability(x, S, y, FALSE)$log_probability, mu)
    expect_equal(answer$mean_gradient, direct, tolerance = 2e-4)
    for (i in seq_len(m)) for (j in seq_len(i)) {
      H <- matrix(0, m, m); H[i, j] <- H[j, i] <- 1
      h <- 1e-5
      direct <- (.ravel_mixed_probit_probability(mu, S + h * H, y, FALSE)$log_probability -
        .ravel_mixed_probit_probability(mu, S - h * H, y, FALSE)$log_probability) / (2 * h)
      expect_equal(sum(answer$covariance_gradient * H), direct, tolerance = 3e-4)
    }
  }
})

test_that("binary family score retains its natural variance derivative at zero", {
  skip_if_not_installed("mvtnorm")
  dat <- mixed_validation_data(N_families = 60L, sib_sizes = 2L, visits = 1L,
    pi = .4, outcome = "binary_liability", residual_variance = 1,
    age_effect = 0, sex_effect = 0, prevalence_reference = .4)
  p <- .mixed_preparation_fixture(dat, "binary_liability")
  covariance <- setdiff(names(p$start), c("residual_sd", "family_sd"))
  par <- c(p$start[covariance], family_variance = .2,
    stats::setNames(rep(0, length(p$spec$mean_names)), p$spec$mean_names))
  covariance <- c(covariance, "family_variance")
  ev <- .ravel_mixed_binary_evaluator(p$groups, p$spec, covariance)
  actual <- ev(par)
  direct <- .mixed_central_gradient(function(x) ev(x, FALSE)$value, par)
  expect_true(actual$valid)
  expect_equal(unname(actual$gradient), direct, tolerance = 2e-4)
  expect_equal(-colSums(actual$family_scores), actual$gradient, tolerance = 1e-12)
  par["family_variance"] <- 0
  zero <- ev(par)
  plus <- par; plus["family_variance"] <- 1e-5
  one_sided <- (ev(plus, FALSE)$value - zero$value) / 1e-5
  expect_equal(unname(zero$gradient["family_variance"]), one_sided, tolerance = 2e-3)
  expect_equal(.ravel_mixed_binary_theta(par, covariance)[["family_sd"]], -Inf)
})

test_that("continuous sufficient statistics preserve differing fixed-effect designs", {
  dat <- mixed_validation_data(N_families = 60L, sib_sizes = 2L, visits = 1L, pi = .4)
  p <- .mixed_preparation_fixture(dat, fixed = c("age", "sex"))
  merged <- .ravel_mixed_compress_means(p$groups)
  a <- .ravel_mixed_evaluator(p$groups, p$spec)(p$start, TRUE)
  b <- .ravel_mixed_evaluator(merged, p$spec)(p$start, TRUE)
  expect_gt(length(p$groups), length(merged))
  expect_equal(a$value, b$value, tolerance = 1e-10)
  expect_equal(a$gradient, b$gradient, tolerance = 1e-9)
  expect_equal(a$means, b$means, tolerance = 1e-10)
  expect_equal(a$mean_information, b$mean_information, tolerance = 1e-10)
})

test_that("continuous variance scores include correlated individual slopes", {
  dat <- mixed_validation_data(N_families = 40L, sib_sizes = 2L, visits = 3L,
    pi = .4, age_pattern = "scheduled", individual_variance = .2,
    slope_variance = .003, intercept_slope_covariance = .012)
  p <- .ravel_mixed_prepare(dat, "continuous", "Y",
    c(child1 = "PGS1_child", child2 = "PGS2_child", mother1 = "PGS1_mother",
      mother2 = "PGS2_mother", father1 = "PGS1_father", father2 = "PGS2_father"),
    "M_ID", "IID", "FID", c("age", "sex"), NULL, NULL, TRUE, TRUE, "age",
    "estimated_state", 1L, .4, TRUE)
  p$groups <- .ravel_mixed_compress_means(p$groups)
  par <- p$start; par["intercept_slope_cor"] <- atanh(.3)
  ev <- .ravel_mixed_evaluator(p$groups, p$spec)
  analytic <- .ravel_mixed_variance_scores(par, p$groups, p$spec, ev(par, TRUE))
  for (name in names(analytic)) {
    variance <- exp(2 * par[[name]]); step <- variance * 1e-4
    a <- b <- par
    a[name] <- log(variance + step) / 2; b[name] <- log(variance - step) / 2
    direct <- (ev(a)$value - ev(b)$value) / (2 * step)
    expect_equal(unname(analytic[name]), direct, tolerance = 2e-5)
  }
})

test_that("mixed family identity checks apply after missing outcome exclusions", {
  dat <- mixed_validation_data(N_families = 60L, sib_sizes = 2L, visits = 2L, pi = .4)
  expect_error(ravel_mixed_fit(dat, fixed_effects = "age", diagnostic = FALSE,
    covariate_model = "independent"), "covariate_mapping_confirmed")
  broken <- dat; broken$FID[1] <- "another_extended_family"
  expect_error(ravel_mixed_fit(broken, diagnostic = FALSE), "one maternal family and one FID")
  broken <- dat; broken$PGS1_child[1] <- broken$PGS1_child[1] + .1
  expect_error(ravel_mixed_fit(broken, diagnostic = FALSE), "constant across visits")
  broken <- dat; broken$Y[duplicated(broken$IID)] <- NA_real_
  expect_error(ravel_mixed_fit(broken, individual_random = TRUE, diagnostic = FALSE), "repeated observed outcomes")
  expect_error(ravel_mixed_fit(dat, outcome = "binary_liability", diagnostic = FALSE), "binary outcome")
  broken <- dat
  families <- unique(broken$M_ID)
  one <- broken$M_ID == families[1L]; two <- broken$M_ID == families[2L]
  broken$F_ID[two] <- broken$F_ID[which(one)[1L]]
  expect_error(ravel_mixed_fit(broken, paternal_id = "F_ID", diagnostic = FALSE),
    "Each father must lie within one FID")
  broken$FID[two] <- broken$FID[which(one)[1L]]
  expect_error(ravel_mixed_fit(broken, paternal_id = "F_ID", diagnostic = FALSE),
    "Paternal PGSs disagree")
})

test_that("binary maternal-family boundary is distinguished from optimizer failure", {
  skip_on_cran()
  skip_if_not_installed("mvtnorm")
  dat <- mixed_validation_data(N_families = 100L, sib_sizes = 2L, visits = 1L,
    pi = .4, outcome = "binary_liability", residual_variance = 1,
    age_effect = 0, sex_effect = 0)
  fit <- ravel_mixed_fit(dat, outcome = "binary_liability", pi = .4, diagnostic = FALSE)
  expect_true(fit$optimizer_converged)
  expect_true(fit$converged)
  expect_gte(fit$parameters[["family_variance"]], 0)
  if (fit$parameters[["family_variance"]] == 0) {
    expect_identical(fit$inference_status, "boundary")
    expect_true("family_variance" %in% fit$boundary_parameters)
    expect_true(all(is.na(fit$table$se)))
    expect_gte(fit$gradient[["family_variance"]], -1e-6)
    expect_equal(fit$projected_gradient[["family_variance"]], 0)
  } else {
    expect_identical(fit$inference_status, "ok")
    expect_true(all(is.finite(fit$table$se)))
  }
  expect_equal(fit$variance_components$se[fit$variance_components$component == "residual"], 0)
  expect_true(all(fit$table$outcome_scale == "conditional probit residual SD = 1"))
})

test_that("continuous zero and near-zero family variance cannot falsely certify Wald inference", {
  skip_on_cran()
  for (variance in c(0, 1e-5)) {
    dat <- mixed_validation_data(N_families = 200L, sib_sizes = 2L, visits = 1L,
      pi = .4, family_variance = variance, age_effect = 0, sex_effect = 0, seed = 1830L)
    fit <- ravel_mixed_fit(dat, pi = .4, diagnostic = FALSE)
    expect_true(fit$optimizer_converged)
    expect_lt(max(abs(fit$projected_gradient)) / fit$n_families, 1e-5)
    expect_gt(fit$natural_variance_gradient[["family_sd"]] / fit$n_families, .01)
    expect_true(fit$inference_status %in% c("boundary", "not_converged"))
    expect_true(all(is.na(fit$table$se)))
    if (!length(fit$boundary_parameters)) expect_false(fit$converged)
  }
})

test_that("mixed-model sandwich uses scores aggregated at FID", {
  skip_on_cran()
  dat <- mixed_validation_data(N_families = 500L, sib_sizes = 2L, visits = 1L,
    pi = .4, family_variance = .5, age_effect = 0, sex_effect = 0)
  dat$FID <- paste0("outer_", (match(dat$M_ID, unique(dat$M_ID)) - 1L) %/% 3L)
  fit <- ravel_mixed_fit(dat, pi = .4, diagnostic = FALSE)
  expect_true(fit$converged)
  expect_identical(fit$inference_status, "ok", info = paste(fit$messages, collapse = "; "))
  u <- fit$uncertainty
  aggregated <- rowsum(u$family_scores, u$cluster_ids, reorder = FALSE)
  reference <- u$model_vcov %*% crossprod(aggregated) %*% u$model_vcov
  expect_equal(unname(u$vcov), unname(reference), tolerance = 1e-10)
  expect_equal(nrow(aggregated), fit$n_clusters)
  expect_equal(nrow(u$family_scores), fit$n_families)
})

test_that("binary uncertainty uses extended-family scores and the full joint Hessian", {
  skip_on_cran()
  skip_if_not_installed("mvtnorm")
  dat <- mixed_validation_data(N_families = 120L, sib_sizes = 2L, visits = 1L,
    pi = .4, outcome = "binary_liability", residual_variance = 1,
    age_effect = 0, sex_effect = 0)
  dat$FID <- paste0("outer_", (match(dat$M_ID, unique(dat$M_ID)) - 1L) %/% 3L)
  fit <- ravel_mixed_fit(dat, outcome = "binary_liability", pi = .4, diagnostic = FALSE)
  expect_true(fit$converged)
  expect_identical(fit$inference_status, "ok", info = paste(fit$messages, collapse = "; "))
  u <- fit$uncertainty
  aggregate <- rowsum(u$family_scores, u$cluster_ids, reorder = FALSE)
  reference <- u$model_vcov %*% crossprod(aggregate) %*% u$model_vcov
  expect_equal(unname(u$vcov), unname(reference), tolerance = 1e-10)
  expect_equal(nrow(aggregate), 40L)
  expect_equal(nrow(u$family_scores), 120L)
  expect_true(all(c(names(fit$parameters), names(fit$means)) %in% colnames(u$information)))
})
