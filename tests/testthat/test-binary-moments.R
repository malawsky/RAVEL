binary_test_keys <- function(labels) {
  vapply(labels, function(label) {
    if (!grepl("~~", label, fixed = TRUE)) return(label)
    paste(sort(strsplit(label, "~~", fixed = TRUE)[[1L]]), collapse = "~~")
  }, character(1), USE.NAMES = FALSE)
}

binary_test_multiplier <- function(labels, scales, reverse_y = FALSE) {
  vapply(labels, function(label) {
    if (label == "Y|t1") return(if (reverse_y) -1 else 1)
    if (endsWith(label, "~1")) return(unname(scales[sub("~1$", "", label)]))
    pair <- strsplit(label, "~~", fixed = TRUE)[[1L]]
    prod(scales[pair]) * if (reverse_y && "Y" %in% pair) -1 else 1
  }, numeric(1), USE.NAMES = FALSE)
}

test_that("binary moments match independent raw-probit and numerical-delta snapshots", {
  fixture <- cov_test_fixture(n = 2400L)
  fixture$dat$PC1 <- fixture$dat$PC1 + 2
  fixture$dat$age <- fixture$dat$age - 3
  group <- rep(c("A", "B", "C"), length.out = nrow(fixture$dat))
  fixture$dat$group_B <- as.numeric(group == "B")
  fixture$dat$group_C <- as.numeric(group == "C")
  fixture$covariates <- c(fixture$covariates, "group_B", "group_C")
  model <- cov_test_reference_model(fixture$covariates)
  before <- serialize(fixture$dat, NULL)
  expected <- cov_test_binary_reference(fixture$dat, model)
  got <- RAVEL:::.ravel_binary_moments(fixture$dat, model)
  expect_identical(names(got), names(expected$args))
  for (name in names(expected$args)) {
    expect_equal(got[[name]], expected$args[[name]], tolerance = 2e-6, info = name)
  }
  expect_identical(serialize(fixture$dat, NULL), before)
  expect_identical(rownames(got$NACOV), expected$labels)
  expect_identical(colnames(got$NACOV), expected$labels)
  expect_equal(got$sample.cov["Y", "Y"], 1)
  expect_equal(got$sample.mean["Y"], c(Y = 0))
  expect_true(all(is.finite(got$NACOV)))
  expect_true(all(is.finite(got$WLS.V)))
  expect_true(all(diag(got$WLS.V) > 0))
  expect_equal(got$NACOV, t(got$NACOV), tolerance = 1e-12)
  ev <- eigen(got$NACOV, symmetric = TRUE, only.values = TRUE)$values
  expect_gte(min(ev), -1e-9 * max(ev))
  expect_gte(sum(abs(ev) < 1e-9 * max(ev)), 4L)
  null <- setNames(numeric(length(expected$labels)), expected$labels)
  null["sex~1"] <- -(1 - 2 * mean(fixture$dat$sex))
  null["sex~~sex"] <- 1
  expect_lt(max(abs(got$NACOV %*% null)), 1e-8)
})

test_that("binary moment Jacobian and full Gamma agree by independent differentiation", {
  fixture <- cov_test_fixture(n = 2400L)
  model <- cov_test_reference_model(fixture$covariates)
  reference <- cov_test_binary_reference(fixture$dat, model)
  analytic <- cov_test_binary_analytic_jacobian(reference)
  fine <- cov_test_fd_jacobian(reference$map, reference$primitive, step = 2e-6)
  expect_equal(reference$jacobian, analytic, tolerance = 2e-7)
  expect_equal(fine, analytic, tolerance = 2e-7)
  gamma <- analytic %*% reference$primitive_gamma %*% t(analytic)
  expect_equal(reference$args$NACOV, gamma, tolerance = 2e-7)
  influence <- reference$primitive_if %*% t(analytic)
  expect_equal(reference$moment_if, influence, tolerance = 2e-7)
  expect_lt(max(abs(colMeans(influence))), 1e-10)
  got <- RAVEL:::.ravel_binary_moments(fixture$dat, model)
  expect_equal(got$NACOV, gamma, tolerance = 2e-6)

  # Dropping first-stage cross-covariances is not the same uncertainty calculation.
  diagonal_only <- analytic %*% diag(diag(reference$primitive_gamma)) %*% t(analytic)
  expect_gt(max(abs(diagonal_only - gamma)), 1e-4)
})

test_that("binary moment order follows model variables rather than dataframe order", {
  fixture <- cov_test_fixture(n = 2400L)
  model <- cov_test_reference_model(fixture$covariates)
  original <- RAVEL:::.ravel_binary_moments(fixture$dat, model)
  reordered_data <- fixture$dat[rev(seq_len(nrow(fixture$dat))), rev(names(fixture$dat))]
  got <- RAVEL:::.ravel_binary_moments(reordered_data, model)
  expect_equal(got, original, tolerance = 2e-7)

  statements <- strsplit(model, "\n", fixed = TRUE)[[1L]]
  reversed_model <- paste(rev(statements[nzchar(trimws(statements))]), collapse = "\n")
  reversed <- RAVEL:::.ravel_binary_moments(fixture$dat, reversed_model)
  independent <- cov_test_binary_reference(fixture$dat, reversed_model)
  expect_equal(reversed, independent$args, tolerance = 2e-6)
  keys <- binary_test_keys(rownames(original$NACOV))
  other_keys <- binary_test_keys(rownames(reversed$NACOV))
  expect_setequal(keys, other_keys)
  expect_equal(anyDuplicated(keys), 0L)
  expect_equal(anyDuplicated(other_keys), 0L)
  order <- match(keys, other_keys)
  vars <- rownames(original$sample.cov)
  expect_equal(reversed$sample.cov[vars, vars], original$sample.cov, tolerance = 2e-7)
  expect_equal(reversed$sample.mean[vars], original$sample.mean, tolerance = 2e-7)
  expect_equal(unname(reversed$NACOV[order, order]), unname(original$NACOV), tolerance = 2e-6)
  expect_equal(reversed$WLS.V[order, order], original$WLS.V, tolerance = 2e-7)
})

test_that("affine covariate recoding transforms all moments and Gamma consistently", {
  fixture <- cov_test_fixture(n = 2400L)
  model <- cov_test_reference_model(fixture$covariates)
  original <- RAVEL:::.ravel_binary_moments(fixture$dat, model)
  recoded <- fixture$dat
  recoded$PC1 <- recoded$PC1 - 3
  recoded$age <- .7 * recoded$age + 2
  recoded$sex <- 1 - recoded$sex
  got <- RAVEL:::.ravel_binary_moments(recoded, model)
  variables <- rownames(original$sample.cov)
  scales <- setNames(rep(1, length(variables)), variables)
  scales[c("age", "sex")] <- c(.7, -1)
  shifts <- setNames(rep(0, length(variables)), variables)
  shifts[c("PC1", "age", "sex")] <- c(-3, 2, 1)
  multiplier <- binary_test_multiplier(rownames(original$NACOV), scales)
  expect_equal(got$sample.cov, original$sample.cov * outer(scales, scales), tolerance = 2e-7)
  expect_equal(got$sample.mean, original$sample.mean * scales + shifts, tolerance = 2e-7)
  expect_equal(unname(got$sample.th["Y"]), unname(original$sample.th["Y"]), tolerance = 2e-7)
  expect_equal(got$NACOV, original$NACOV * outer(multiplier, multiplier), tolerance = 2e-6)
  expect_equal(got$WLS.V, original$WLS.V / outer(multiplier, multiplier), tolerance = 2e-7)
  expect_equal(got, cov_test_binary_reference(recoded, model)$args, tolerance = 2e-6)
})

test_that("binary reversal changes liability signs and preserves covariance uncertainty", {
  fixture <- cov_test_fixture(n = 2400L)
  model <- cov_test_reference_model(fixture$covariates)
  original <- RAVEL:::.ravel_binary_moments(fixture$dat, model)
  reversed <- fixture$dat
  reversed$Y <- 1 - reversed$Y
  got <- RAVEL:::.ravel_binary_moments(reversed, model)
  scales <- setNames(rep(1, nrow(original$sample.cov)), rownames(original$sample.cov))
  multiplier <- binary_test_multiplier(rownames(original$NACOV), scales, reverse_y = TRUE)
  signs <- scales
  signs["Y"] <- -1
  expect_equal(got$sample.cov, original$sample.cov * outer(signs, signs), tolerance = 2e-7)
  expect_equal(got$sample.mean, original$sample.mean, tolerance = 2e-7)
  expect_equal(unname(got$sample.th["Y"]), -unname(original$sample.th["Y"]), tolerance = 2e-7)
  expect_equal(got$NACOV, original$NACOV * outer(multiplier, multiplier), tolerance = 2e-6)
  expect_equal(got$WLS.V, original$WLS.V, tolerance = 2e-7)
})

test_that("binary moments use complete model rows without dropping unused missing columns", {
  fixture <- cov_test_fixture(n = 2400L)
  fixture$dat$unused <- NA_real_
  model <- cov_test_reference_model(fixture$covariates)
  expected <- RAVEL:::.ravel_binary_moments(fixture$dat, model)
  expect_equal(expected$sample.nobs, nrow(fixture$dat))
  fixture$dat$sex[1:3] <- NA_real_
  fixture$dat$Y[4:5] <- NA_real_
  fixture$dat$PGS2_mother[6:7] <- NA_real_
  got <- RAVEL:::.ravel_binary_moments(fixture$dat, model)
  expect_equal(got$sample.nobs, nrow(fixture$dat) - 7L)
  expect_equal(got, cov_test_binary_reference(fixture$dat, model)$args, tolerance = 2e-6)
  expect_equal(got, RAVEL:::.ravel_binary_moments(fixture$dat[-(1:7), ], model), tolerance = 2e-7)
})

test_that("the free binary SEM uses the supplied matrices and independent sandwich covariance", {
  fixture <- cov_test_fixture()
  model <- cov_test_reference_model(fixture$covariates, fixed_corr = FALSE)
  reference <- cov_test_binary_reference(fixture$dat, model)
  expected <- cov_test_fit(fixture, "reference", fixed_corr = FALSE)
  result <- cov_test_fit(fixture, fixed_corr = FALSE)
  cov_test_expect_reference(result, expected)
  if (!isTRUE(result$converged) || !isTRUE(expected$converged)) {
    stop("Both free SEM fits must converge for the independent sandwich check.")
  }
  fit <- result$fit
  observed <- lavaan::lavInspect(fit, "wls.obs")
  D <- lavaan::lavTech(fit, "delta")[[1L]]
  Gamma <- lavaan::lavTech(fit, "gamma")[[1L]]
  W <- lavaan::lavTech(fit, "wls.v")[[1L]]
  if (is.null(dim(W))) W <- diag(W)
  expect_identical(names(observed), reference$labels)
  expect_equal(unname(unclass(observed)), unname(reference$moments), tolerance = 2e-6)
  expect_equal(unname(unclass(Gamma)), unname(reference$args$NACOV), tolerance = 2e-6)
  expect_equal(unname(unclass(W)), unname(reference$args$WLS.V), tolerance = 2e-7)
  expect_equal(nrow(D), length(reference$labels))
  expect_equal(ncol(D), length(lavaan::coef(fit)))
  expect_equal(qr(D)$rank, ncol(D))
  n <- reference$args$sample.nobs
  A <- solve(crossprod(D, reference$args$WLS.V %*% D))
  WD <- reference$args$WLS.V %*% D
  sandwich <- A %*% crossprod(WD, reference$args$NACOV %*% WD) %*% A
  options <- lavaan::lavInspect(fit, "options")
  # lavaan 0.6-21 divides non-ML covariance by N-1, with an optional extra factor.
  denominator <- n - 1L
  if (isTRUE(options$gamma.vcov.mplus)) sandwich <- sandwich * n / (n - 1L)
  sandwich <- sandwich / denominator
  expect_equal(unname(unclass(lavaan::vcov(fit))), unname(sandwich), tolerance = 2e-5)
  coefficients <- lavaan::coef(fit)
  rownames(sandwich) <- colnames(sandwich) <- names(coefficients)
  core <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f")
  expect_true(all(diag(sandwich)[core] > 0))
  expect_equal(unname(unclass(lavaan::vcov(fit)[core, core])), unname(sandwich[core, core]),
               tolerance = 2e-5)

  # Dummy moment redundancy must not remove any tested residual direction.
  d <- nrow(D) - ncol(D)
  expect_gt(d, 0L)
  Z <- qr.Q(qr(D), complete = TRUE)[, seq.int(ncol(D) + 1L, nrow(D)), drop = FALSE]
  residual_variance <- crossprod(Z, reference$args$NACOV %*% Z)
  eigenvalues <- eigen(residual_variance, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(length(eigenvalues), d)
  expect_gt(min(eigenvalues), 1e-9 * max(eigenvalues))
  expect_equal(as.numeric(lavaan::fitMeasures(fit, "df")), as.numeric(d))
})
