test_that("exactly shared covariates retain their mean-space conditioning", {
  W <- diag(c(0, .4))
  B <- matrix(c(.8, .2, .2, .5), 2L)
  moments <- list(mean = c(0, 0), Qself = B + W, Qsib = B, shared = 1L)
  for (n in c(1L, 2L, 4L)) {
    Q <- kronecker(diag(n), W) + kronecker(matrix(1, n, n), B)
    actual <- .ravel_mixed_covariate_inverse(moments, n)$Qinv
    eigen_Q <- eigen(Q, symmetric = TRUE)
    keep <- eigen_Q$values > 1e-10
    U <- eigen_Q$vectors[, keep, drop = FALSE]
    oracle <- sweep(U, 2L, eigen_Q$values[keep], "/") %*% t(U)
    expect_equal(actual, oracle, tolerance = 1e-11)
    expect_equal(Q %*% actual %*% Q, Q, tolerance = 1e-11)
    # Repeating a shared value does not create extra independent information.
    shared <- rep(c(1, 0), n)
    expect_true(is.finite(drop(crossprod(shared, actual %*% shared))))
    if (n > 1L) {
      contrast <- c(1, 0, -1, 0, rep(0, 2L * (n - 2L)))
      expect_equal(drop(actual %*% contrast), rep(0, 2L * n), tolerance = 1e-12)
    }
  }
})

test_that("a single entirely family-shared feature has estimable moments and influence", {
  set.seed(1824)
  value <- rnorm(60L, mean = 2)
  sizes <- rep(c(1L, 1L, 2L), 20L)
  family <- lapply(seq_along(sizes), function(i) matrix(rep(value[i], sizes[i]), ncol = 1L))
  fitted <- .ravel_mixed_covariate_moments(family)
  expect_identical(fitted$shared, 1L)
  expected_mean <- sum(value * sizes) / sum(sizes)
  expected_variance <- mean((value - expected_mean)^2)
  expect_equal(fitted$mean, expected_mean, tolerance = 1e-12)
  expect_equal(drop(fitted$Qself), expected_variance, tolerance = 1e-8)
  expect_equal(fitted$Qsib, fitted$Qself, tolerance = 1e-12)
  plus <- minus <- rep(1, length(family))
  plus[3L] <- plus[3L] + 1e-4; minus[3L] <- minus[3L] - 1e-4
  ep <- .ravel_mixed_covariate_moments(family, weights = plus, influence = FALSE)$eta
  em <- .ravel_mixed_covariate_moments(family, weights = minus, influence = FALSE)$eta
  expect_equal(unname(fitted$moment_scores[3L, ] / length(family)),
    unname((ep - em) / 2e-4), tolerance = 1e-6)
})

test_that("shared-covariate genetic paths are tied and moment derivatives respect support", {
  dat <- mixed_covariate_data(N_families = 100L, sib_sizes = 2L, seed = 1823)
  dat$PC1 <- ave(dat$PC1, dat$M_ID, FUN = mean)
  p <- .ravel_mixed_prepare(dat, "continuous", "Y",
    c(child1 = "PGS1_child", child2 = "PGS2_child", mother1 = "PGS1_mother",
      mother2 = "PGS2_mother", father1 = "PGS1_father", father2 = "PGS2_father"),
    "M_ID", "IID", "FID", c("PC1", "PC2", "sex", "age"), NULL, NULL,
    FALSE, FALSE, NULL, "estimated_state", 1L, .4, FALSE,
    covariate_model = "joint")
  p <- .ravel_mixed_covariate_prepare(p)
  cv <- p$spec$covariate
  shared <- match("PC1", cv$columns)
  expect_true(shared %in% cv$shared)
  expect_true("PC1" %in% p$spec$covariates)
  expect_identical(cv$parameter_names$child_self[shared],
    cv$parameter_names$child_sibling[shared])
  expect_false(anyDuplicated(names(p$start)) > 0L)
  values <- .ravel_mixed_covariate_unpack(cv)
  expect_equal((values$Qself - values$Qsib)[shared, ], rep(0, cv$q), tolerance = 1e-12)
  par <- p$start
  names_D <- unique(unlist(cv$parameter_names, use.names = FALSE))
  par[names_D] <- .025 * sin(seq_along(names_D))
  evaluate <- .ravel_mixed_covariate_evaluator(p, p$spec)
  base <- evaluate(par, moments = TRUE)
  expect_true(base$valid)
  # A family reweighting supplies a tangent that preserves the exact support.
  direction <- cv$moment_scores[3L, ] / max(1, max(abs(cv$moment_scores[3L, ])))
  step <- 1e-5
  plus <- minus <- p$spec
  plus$covariate$eta <- cv$eta + step * direction
  minus$covariate$eta <- cv$eta - step * direction
  ep <- .ravel_mixed_covariate_evaluator(p, plus)(par)
  em <- .ravel_mixed_covariate_evaluator(p, minus)(par)
  expect_true(ep$valid && em$valid)
  expect_equal(sum(base$moment_gradient * direction),
    (ep$value - em$value) / (2 * step), tolerance = 1e-6)
})
