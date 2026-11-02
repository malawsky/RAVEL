test_that("mixed validation preserves IDs, person-level scores and RNG state", {
  set.seed(912)
  original_seed <- .Random.seed
  d <- mixed_validation_data(N_families = 8L, sib_sizes = rep(c(1L, 3L), 4L),
    visits = 3L, families_per_FID = 2L, seed = 91)
  expect_identical(.Random.seed, original_seed)
  expect_equal(nrow(d), 48L)
  expect_equal(length(unique(d$IID)), 16L)
  expect_equal(length(unique(d$M_ID)), 8L)
  expect_equal(length(unique(d$FID)), 4L)
  expect_true(all(d$M_ID != d$FID))
  score <- grep("^PGS", names(d), value = TRUE)
  for (x in split(d, d$IID)) {
    expect_equal(nrow(unique(x[score])), 1L)
    expect_equal(length(unique(x$M_ID)), 1L)
    expect_equal(length(unique(x$FID)), 1L)
    expect_equal(length(unique(x$sex)), 1L)
  }
  for (x in split(d, d$M_ID)) {
    expect_equal(length(unique(x$FID)), 1L)
    expect_equal(nrow(unique(x[grep("mother|father", score, value = TRUE)])), 1L)
  }
  expect_identical(d, mixed_validation_data(N_families = 8L,
    sib_sizes = rep(c(1L, 3L), 4L), visits = 3L, families_per_FID = 2L, seed = 91))
})

test_that("independent component algebra reproduces genetic moments and Mendelian score errors", {
  d <- mixed_validation_data(N_families = 20000L, sib_sizes = 2L,
    family_variance = 0, age_effect = 0, sex_effect = 0, seed = 381)
  tr <- attr(d, "truth")
  latent <- attr(d, "latent")
  T <- as.matrix(latent[paste0("T_", c("child", "mother", "father"))])
  G <- as.matrix(latent[paste0("G_", c("child", "mother", "father"))])
  expect_equal(unname(stats::cov(T)), unname(tr$population$Sigma_TT), tolerance = .04)
  expect_equal(unname(stats::cov(G)), unname(tr$population$Sigma_GG), tolerance = .04)
  expect_equal(unname(stats::cov(G, T)), unname(tr$population$Sigma_GT), tolerance = .04)
  e <- lapply(1:2, function(k) {
    raw <- sweep(G, 2L, tr$population$visible_raw_sd, "*")
    (as.matrix(d[paste0("PGS", k, "_", c("child", "mother", "father"))]) - raw) /
      sqrt((1 - tr$reliability_founder[k]) / tr$reliability_founder[k])
  })
  target <- matrix(c(1, .5, .5, .5, 1, 0, .5, 0, 1), 3L)
  for (k in 1:2) expect_equal(unname(stats::cov(e[[k]])), target, tolerance = .04)
  expect_lt(max(abs(stats::cov(e[[1]], e[[2]]))), .035)
  first <- seq.int(1L, nrow(d), by = 2L)
  expected_sib_cov <- .5 * (1 + tr$population$rho_total) * tr$population$parent_variance
  expect_equal(stats::cov(T[first, 1L], T[first + 1L, 1L]), expected_sib_cov, tolerance = .04)
  expect_equal(stats::cov(e[[1]][first, 1L], e[[1]][first + 1L, 1L]), .5, tolerance = .04)
})

test_that("visible oracle is an independent covariance projection in each history mode", {
  for (mode in c("estimated_state", "n_gen", "equilibrium")) {
    d <- mixed_validation_data(N_families = 3L, mode = mode, pi = .45,
      age_beta = c(child = .02, mother = -.01, father = 0),
      sex_beta = c(child = .1, mother = 0, father = -.04))
    tr <- attr(d, "truth")
    p <- tr$population
    expect_equal(unname(p$Sigma_GG %*% p$total_to_visible), unname(p$Sigma_GT), tolerance = 1e-12)
    expect_equal(unname(tr$coefficients$corrected_std[1:3]),
      unname(tr$beta_raw * sqrt(diag(p$Sigma_TT))), tolerance = 1e-12)
    expect_equal(unname(tr$coefficients$visible[4:6]),
      unname(drop(p$total_to_visible %*% tr$age_beta_raw)), tolerance = 1e-12)
  }
  p0 <- mixed_validation_population(pi = 0)
  expect_equal(unname(p0$total_to_visible), diag(unname(p0$total_sd)), tolerance = 1e-12)
})

test_that("extended-family dependence preserves the declared maternal-family variance", {
  d <- mixed_validation_data(N_families = 20000L, sib_sizes = 1L,
    family_variance = .6, families_per_FID = 2L, extended_fraction = .5, seed = 192)
  u <- attr(d, "latent")$family_effect
  first <- seq.int(1L, length(u), by = 2L)
  expect_equal(stats::var(u), .6, tolerance = .025)
  expect_equal(stats::cov(u[first], u[first + 1L]), .3, tolerance = .025)
  expect_equal(attr(d, "truth")$family_variance, .6)
})

test_that("continuous outcome covariance separates genetic, family and individual terms", {
  d <- mixed_validation_data(N_families = 18000L, sib_sizes = 2L, visits = 2L,
    individual_variance = .2, family_variance = .3, residual_variance = .7,
    age_effect = 0, sex_effect = 0, seed = 82)
  tr <- attr(d, "truth")
  b <- tr$beta_raw
  c_m <- b["mother"] + .5 * b["child"]
  c_f <- b["father"] + .5 * b["child"]
  shared_genetic <- tr$population$parent_variance *
    (c_m^2 + c_f^2 + 2 * tr$population$rho_total * c_m * c_f)
  child_genetic <- .5 * b["child"]^2
  first <- seq.int(1L, nrow(d), by = 4L)
  expect_equal(stats::var(d$Y), unname(shared_genetic + child_genetic + .3 + .2 + .7), tolerance = .04)
  expect_equal(stats::cov(d$Y[first], d$Y[first + 1L]),
    unname(shared_genetic + child_genetic + .3 + .2), tolerance = .04)
  expect_equal(stats::cov(d$Y[first], d$Y[first + 2L]),
    unname(shared_genetic + .3), tolerance = .04)
})

test_that("binary data use a fixed conditional-probit scale and retain concordant sibships", {
  d <- mixed_validation_data(N_families = 20000L, sib_sizes = 2L,
    outcome = "binary_liability", residual_variance = 1,
    prevalence_reference = .15, age_effect = 0, sex_effect = 0, seed = 891)
  expect_setequal(unique(d$Y), c(0, 1))
  expect_equal(mean(d$Y), .15, tolerance = .015)
  first <- seq.int(1L, nrow(d), by = 2L)
  expect_true(any(d$Y[first] == d$Y[first + 1L]))
  expect_true(any(d$Y[first] != d$Y[first + 1L]))
  expect_equal(attr(d, "truth")$residual_variance, 1)
  expect_error(mixed_validation_data(outcome = "binary_liability"), "residual variance")
})

test_that("mixed likelihood covariance agrees with the independent primitive oracle", {
  skip_if_not(exists(".ravel_mixed_covariance", envir = asNamespace("RAVEL"), inherits = FALSE))
  for (mode in c("estimated_state", "n_gen", "equilibrium")) {
    p <- mixed_validation_population(mode = mode, pi = .4, n_generations = 3L)
    beta <- c(child = .3, mother = .12, father = -.04)
    age_beta <- c(child = .025, mother = -.008, father = 0)
    sex_beta <- c(child = .12, mother = -.06, father = .02)
    scales <- matrix(c(.7, 1.2, 1.6, .9, 1.1, 1.8), 2L, 3L,
                     dimnames = list(NULL, c("child", "mother", "father")))
    beta_names <- unlist(lapply(c("main", "age", "sex"),
      function(x) paste0("genetic_", c("child", "mother", "father"), "_", x)))
    spec <- list(pi = .4, mode = mode, generations = 3L,
      beta_names = beta_names, family_random = TRUE, individual_random = TRUE, random_slope = "age")
    par <- c(state_rho = atanh(p$rho_total), state_variance = log(p$parent_variance - .5),
             noise1 = log(sqrt((1 - .65) / .65)), noise2 = log(sqrt((1 - .4) / .4)),
             residual_sd = log(sqrt(.75)), family_sd = log(sqrt(.25)),
             individual_sd = log(sqrt(.2)), slope_sd = log(sqrt(.003)),
             intercept_slope_cor = atanh(.012 / sqrt(.2 * .003)))
    for (k in 1:2) for (role in c("child", "mother", "father"))
      par[paste0("loading", k, "_", role)] <- log(scales[k, role])
    par[beta_names] <- c(beta * p$total_sd, age_beta * p$total_sd, sex_beta * p$total_sd)
    child <- rep(1:3, each = 3L)
    age <- c(4, 10, 15, 6, 9, 16, 5, 11, 14)
    sex <- rep(c(0, 1, 1), each = 3L)
    group <- list(nchild = 3L, child = child,
      moderators = cbind(1, age - 10, sex), slope = age - 10)
    observed <- RAVEL:::.ravel_mixed_covariance(par, group, spec)
    expected <- mixed_validation_observed_covariance(p, n_children = 3L,
      child_index = child, age = age, sex = sex, score_scales = scales)
    expect_equal(observed, expected, tolerance = 1e-10)
  }
})

test_that("one-trio no-random-effects likelihood reproduces the original fitted effects", {
  skip_if_not(exists("ravel_mixed_fit", envir = asNamespace("RAVEL"), inherits = FALSE))
  d <- mixed_validation_data(N_families = 1200L, sib_sizes = 1L,
    family_variance = 0, individual_variance = 0, age_pattern = "scheduled",
    age_effect = 0, sex_effect = 0, seed = 183)
  original <- ravel_fit(d, pi = .4, fit_constraints = "estimated_state",
                        am_correction = "estimated_state")
  mixed <- ravel_mixed_fit(d, pi = .4, family_random = FALSE,
                          fit_constraints = "estimated_state", diagnostic = FALSE)
  expect_true(original$converged)
  expect_true(mixed$converged)
  a <- as.data.frame(original, include_raw = TRUE)
  wanted <- a$parameter %in% c("beta_D", "beta_M", "beta_F") &
    a$scale %in% c("visible", "corrected_std", "corrected_raw")
  a <- a[wanted, , drop = FALSE]
  idx <- match(paste(a$scale, a$parameter), paste(mixed$table$scale, mixed$table$parameter))
  expect_false(anyNA(idx))
  # The two paths use different optimizers; inferential methods also differ.
  # Equality of the three effect scales checks the statistical nesting itself.
  expect_equal(unname(a$estimate), unname(mixed$table$estimate[idx]), tolerance = 3e-5)
})
