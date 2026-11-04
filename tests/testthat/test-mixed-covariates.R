test_that("correlated child covariate fixtures preserve person identities and reference ages", {
  set.seed(831)
  previous <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  contract <- mixed_covariate_contract(sex_pc_association = .7)
  d <- mixed_covariate_data(N_families = 12L, visits = 3L,
    contract = contract, age_pattern = "random", seed = 8421L)
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), previous)
  expect_identical(d, mixed_covariate_data(N_families = 12L, visits = 3L,
    contract = contract, age_pattern = "random", seed = 8421L))
  for (x in split(d, d$IID)) {
    expect_equal(mean(x$age), x$baseline_age[1L], tolerance = 1e-12)
    expect_equal(nrow(unique(x[c("PC1", "PC2", "sex", "baseline_age")])), 1L)
    expect_gt(length(unique(x$age)), 1L)
  }
  expect_setequal(unique(d$sex), c(0, 1))
  expect_equal(attr(d, "truth")$population$pi, .4)
})

test_that("the independent Schur oracle retains marginal pi and AM state", {
  contract <- mixed_covariate_contract(sex_pc_association = .7)
  for (mode in c("estimated_state", "equilibrium", "n_gen")) {
    population <- mixed_validation_population(mode = mode)
    for (n in c(1L, 3L)) {
      x <- mixed_covariate_family(population, n, contract)
      expect_equal(unname(x$E), unname(x$projection %*% x$D), tolerance = 1e-11)
      expect_equal(x$conditional_TV, x$projection %*% x$conditional_V, tolerance = 1e-11)
      expect_equal(x$conditional_T,
        x$projection %*% x$conditional_V %*% t(x$projection) + x$hidden, tolerance = 1e-11)
      expect_equal(x$conditional_V + x$D %*% solve(x$Q, t(x$D)), x$VV, tolerance = 1e-11)
      expect_gt(min(eigen(x$conditional_V, symmetric = TRUE)$values), 0)
      expect_gt(min(eigen(x$hidden, symmetric = TRUE)$values), 0)
      expect_gt(max(abs(x$VV - x$conditional_V)), .02)
    }
  }
})

test_that("family covariate moments are projectively consistent across sibling counts", {
  contract <- mixed_covariate_contract(sex_pc_association = .5)
  two <- mixed_covariate_family(n_children = 2L, contract = contract)
  three <- mixed_covariate_family(n_children = 3L, contract = contract)
  q <- length(two$covariates)
  expect_equal(three$VV[1:4, 1:4], two$VV, tolerance = 1e-12)
  expect_equal(three$Q[seq_len(2L * q), seq_len(2L * q)], two$Q, tolerance = 1e-12)
  expect_equal(three$D[1:4, seq_len(2L * q)], two$D, tolerance = 1e-12)
  expect_equal(three$E[1:4, seq_len(2L * q)], two$E, tolerance = 1e-12)
  # Within this Mendelian generator, another sibling's covariates cannot
  # correlate with this child's independent segregation innovation.
  expect_equal(unname(two$D[3L, q + seq_len(q)]),
    unname(.5 * (two$D[1L, q + seq_len(q)] + two$D[2L, q + seq_len(q)])), tolerance = 1e-12)
})

test_that("sampled correlated PCs and binary sex match independent population moments", {
  contract <- mixed_covariate_contract(sex_pc_association = .7)
  d <- mixed_covariate_data(N_families = 20000L, contract = contract,
    fixed_pc = c(0, 0), age_effect = 0, sex_effect = 0, seed = 9131L)
  truth <- attr(d, "truth"); latent <- attr(d, "latent")
  family <- mixed_covariate_family(truth$population, 2L, contract)
  first <- seq.int(1L, nrow(latent), by = 2L)
  visible_sd <- stats::setNames(truth$population$visible_raw_sd, c("child", "mother", "father"))
  V <- cbind(latent$G_mother[first] * visible_sd["mother"],
    latent$G_father[first] * visible_sd["father"],
    latent$G_child[first] * visible_sd["child"],
    latent$G_child[first + 1L] * visible_sd["child"])
  T <- cbind(latent$T_mother[first], latent$T_father[first],
             latent$T_child[first], latent$T_child[first + 1L])
  C <- cbind(as.matrix(d[first, family$covariates]),
             as.matrix(d[first + 1L, family$covariates]))
  C <- sweep(C, 2L, rep(family$mean_C, 2L))
  expect_equal(unname(stats::cov(C)), unname(family$Q), tolerance = .04)
  expect_equal(unname(stats::cov(V, C)), unname(family$D), tolerance = .04)
  expect_equal(unname(stats::cov(T, C)), unname(family$E), tolerance = .04)
  mean_V <- C %*% solve(family$Q, t(family$D))
  mean_T <- C %*% solve(family$Q, t(family$E))
  expect_equal(unname(stats::cov(V - mean_V)), unname(family$conditional_V), tolerance = .04)
  expect_equal(unname(stats::cov(T - mean_T)), unname(family$conditional_T), tolerance = .04)
  expect_lt(max(abs(stats::cov(T - V %*% t(family$projection), C))), .04)
  expect_gt(abs(stats::cor(d$PGS1_child, d$PC1)), .1)
  expect_gt(abs(stats::cor(d$PGS1_child, d$sex)), .02)
  # PCs suffice for the genetic mean; sex and baseline age add no information
  # conditional on those PCs in this particular generator.
  coefficient <- family$D %*% solve(family$Q)
  redundant <- c(3L, 4L, 7L, 8L)
  expect_lt(max(abs(coefficient[, redundant])), 1e-12)
})

test_that("marginal pi alone does not identify a generic covariate-adjusted total map", {
  r <- sqrt(.6); d <- .4
  covariance <- function(e) matrix(c(1, r, d, r, 1, e, d, e, 1), 3L)
  projection_closure <- r * d
  trait_aligned_closure <- d / r
  expect_gt(min(eigen(covariance(projection_closure), symmetric = TRUE)$values), 0)
  expect_gt(min(eigen(covariance(trait_aligned_closure), symmetric = TRUE)$values), 0)
  multiplier <- function(e) (1 - d^2) / (r - d * e)
  expect_equal(multiplier(projection_closure), 1 / r, tolerance = 1e-12)
  expect_gt(multiplier(trait_aligned_closure) / multiplier(projection_closure), 1.14)
})

test_that("binary correlated-covariate data use the stated conditional-probit residual scale", {
  d <- mixed_covariate_data(N_families = 1000L, residual_variance = 1,
    outcome = "binary_liability", seed = 14132L)
  expect_setequal(unique(d$Y), c(0, 1))
  expect_equal(d$Y, as.numeric(attr(d, "liability") > 0))
  expect_equal(attr(d, "truth")$residual_variance, 1)
  expect_error(mixed_covariate_data(outcome = "binary_liability"), "residual variance")
})

test_that("fitted conditional covariance and mean match the independent primitive oracle", {
  contract <- mixed_covariate_contract(sex_pc_association = .7)
  values <- rbind(c(.3, -.4, 1, 10), c(-.2, .8, 0, 11))
  q <- ncol(values)
  triangle <- which(lower.tri(matrix(0, q, q), diag = TRUE), arr.ind = TRUE)
  for (mode in c("estimated_state", "equilibrium", "n_gen")) {
    population <- mixed_validation_population(mode = mode)
    family <- mixed_covariate_family(population, 2L, contract, c_values = values)
    parameter_names <- stats::setNames(lapply(c("mother", "father", "child_self", "child_sibling"),
      function(role) paste0("C_D_", role, "_", seq_len(q))),
      c("mother", "father", "child_self", "child_sibling"))
    n_second <- nrow(triangle)
    covariate <- list(q = q, triangle = triangle, has_siblings = TRUE,
      eta = c(family$mean_C, (family$Q_self + tcrossprod(family$mean_C))[triangle],
        (family$Q_sib + tcrossprod(family$mean_C))[triangle]),
      index = list(mean = seq_len(q), self = q + seq_len(n_second),
        sibling = q + n_second + seq_len(n_second)), parameter_names = parameter_names)
    names(covariate$eta) <- paste0("moment", seq_along(covariate$eta))
    beta <- c(child = .3, mother = .12, father = -.04)
    age_beta <- c(child = .025, mother = -.008, father = 0)
    sex_beta <- c(child = .12, mother = -.06, father = .02)
    beta_names <- unlist(lapply(c("main", "age", "sex"),
      function(z) paste0("genetic_", c("child", "mother", "father"), "_", z)))
    spec <- list(pi = .4, mode = mode, generations = 3L, beta_names = beta_names,
      family_random = TRUE, individual_random = TRUE, random_slope = "age", covariate = covariate)
    par <- c(state_rho = atanh(.3), state_variance = log(population$parent_variance - .5),
      noise1 = log(sqrt((1 - .65) / .65)), noise2 = log(sqrt((1 - .4) / .4)),
      residual_sd = log(sqrt(.75)), family_sd = log(sqrt(.25)),
      individual_sd = log(sqrt(.2)), slope_sd = log(sqrt(.003)),
      intercept_slope_cor = atanh(.012 / sqrt(.2 * .003)))
    scales <- matrix(c(.7, 1.2, 1.6, .9, 1.1, 1.8), 2L, 3L,
      dimnames = list(NULL, c("child", "mother", "father")))
    for (k in 1:2) for (role in c("child", "mother", "father"))
      par[paste0("loading", k, "_", role)] <- log(scales[k, role])
    par[beta_names] <- c(beta * population$total_sd, age_beta * population$total_sd,
                         sex_beta * population$total_sd)
    par[parameter_names$mother] <- family$D[1L, seq_len(q)]
    par[parameter_names$father] <- family$D[2L, seq_len(q)]
    par[parameter_names$child_self] <- family$D[3L, seq_len(q)]
    par[parameter_names$child_sibling] <- family$D[3L, q + seq_len(q)]
    child <- rep(1:2, each = 3L); age <- c(7, 10, 13, 8, 11, 14)
    sex <- rep(c(1, 0), each = 3L)
    group <- list(nchild = 2L, n = 1L, child = child,
      moderators = cbind(1, age - 10, sex), slope = age - 10,
      C = matrix(c(t(values)), 1L), observed = matrix(0, 1L, 14L))
    actual <- RAVEL:::.ravel_mixed_covariate_joint(par, group, spec)
    expected <- mixed_covariate_observed(family, child, age, sex, beta, age_beta, sex_beta,
      score_scales = scales, individual_variance = .2, slope_variance = .003,
      intercept_slope_covariance = .012)
    expect_equal(actual$S, expected$covariance, tolerance = 1e-10)
    expect_equal(drop(actual$offset), expected$mean, tolerance = 1e-10)
  }
})

mixed_covariate_cached_test_fit <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      d <- mixed_covariate_data(N_families = 300L,
        contract = mixed_covariate_contract(n_pc = 1L, sex_pc_association = .5), seed = 721L)
      value <<- ravel_mixed_fit(d, pi = .4, fixed_effects = c("PC1", "sex", "age"),
                               diagnostic = FALSE)
    }
    value
  }
})

test_that("full information, family scores and stacked covariate uncertainty agree with independent differences", {
  fit <- mixed_covariate_cached_test_fit()
  expect_true(fit$converged)
  expect_identical(fit$inference_status, "ok")
  audit <- mixed_covariate_engine_audit(fit)
  expect_lt(audit[["information_relative_error"]], 1e-4)
  expect_lt(audit[["family_score_max_error"]], 1e-5)
  expect_lt(audit[["moment_gradient_max_error"]], 1e-5)
  expect_lt(audit[["moment_influence_max_error"]], 1e-6)
  expect_lt(audit[["cross_derivative_relative_error"]], 1e-4)
  expect_lt(audit[["stacked_covariance_relative_error"]], 1e-3)
})

test_that("genetic estimates and stacked SEs are invariant to covariate affine units", {
  original <- mixed_covariate_cached_test_fit()
  d <- original$prepared$data
  d$PC1 <- 100 * d$PC1 + 7
  transformed <- ravel_mixed_fit(d, pi = .4, fixed_effects = c("PC1", "sex", "age"),
                                 diagnostic = FALSE)
  expect_true(transformed$converged)
  expect_identical(transformed$inference_status, "ok")
  expect_equal(original$table$estimate, transformed$table$estimate, tolerance = 2e-5)
  expect_equal(original$table$se, transformed$table$se, tolerance = 2e-5)
  expect_equal(unlist(original$state), unlist(transformed$state), tolerance = 2e-5)
})

test_that("no-covariate calls preserve the established mixed path exactly", {
  d <- mixed_validation_data(N_families = 200L, age_effect = 0, sex_effect = 0, seed = 172L)
  a <- ravel_mixed_fit(d, pi = .4, diagnostic = FALSE, se = FALSE, covariate_model = "joint")
  b <- ravel_mixed_fit(d, pi = .4, diagnostic = FALSE, se = FALSE, covariate_model = "independent")
  expect_identical(a$converged, b$converged)
  expect_equal(a$parameters, b$parameters, tolerance = 0)
  expect_equal(a$table$estimate, b$table$estimate, tolerance = 0)
  expect_equal(a$fit_statistics, b$fit_statistics, tolerance = 0)
})

test_that("joint optimizer polishes a failed interior certificate without relaxing it", {
  evaluate <- function(x, derivatives = FALSE) list(valid = TRUE,
    value = .5 * sum(x^2), gradient = if (derivatives) x else NULL)
  certify <- function(x, evaluated) list(gradient_ok = max(abs(evaluated$gradient)) < 1e-9,
                                        boundary_parameters = character())
  start <- c(small_variance_coordinate = 5e-7)
  control <- list(maxit = 100L, pgtol = 1e-6, factr = 100)
  ordinary <- RAVEL:::.ravel_mixed_covariate_optimize(start, evaluate, -1, 1, control, 1L)
  expect_false(certify(ordinary$par, evaluate(ordinary$par, TRUE))$gradient_ok)
  polished <- RAVEL:::.ravel_mixed_covariate_optimize(start, evaluate, -1, 1, control, 1L,
                                                     certify = certify)
  expect_true(certify(polished$par, evaluate(polished$par, TRUE))$gradient_ok)
  expect_true(any(vapply(polished$attempts, function(x)
    identical(x$method, "L-BFGS-B natural-variance polish"), logical(1))))
  boundary <- function(x, evaluated) list(gradient_ok = FALSE,
                                         boundary_parameters = "small_variance_coordinate")
  nonregular <- RAVEL:::.ravel_mixed_covariate_optimize(start, evaluate, -1, 1, control, 1L,
                                                      certify = boundary)
  expect_length(nonregular$attempts, 1L)
  expect_equal(nonregular$par, ordinary$par, tolerance = 0)
})
