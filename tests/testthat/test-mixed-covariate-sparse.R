.sparse_covariate_fixture <- function(sizes, q = 3L, seed = 173L, shared = integer()) {
  set.seed(seed)
  family <- rep(seq_along(sizes), sizes)
  between <- matrix(rnorm(length(sizes) * q), length(sizes), q)
  C <- sqrt(.8) * between[family, , drop = FALSE] +
    sqrt(.2) * matrix(rnorm(length(family) * q), length(family), q)
  if (length(shared)) C[, shared] <- between[family, shared, drop = FALSE]
  C <- scale(C)
  lapply(split(seq_len(nrow(C)), family), function(i) C[i, , drop = FALSE])
}

test_that("sparse sibling covariance estimation gives a coherent pair distribution", {
  x <- .sparse_covariate_fixture(c(rep(1L, 380L), rep(2L, 20L)))
  m <- .ravel_mixed_covariate_moments(x)
  expect_gt(min(eigen(m$Qself - m$Qsib, symmetric = TRUE)$values), 0)
  expect_gt(min(eigen(m$Qself + m$Qsib, symmetric = TRUE)$values), 0)
  expect_true(all(is.finite(m$moment_scores)))
  expect_lt(max(abs(colSums(m$moment_scores))), 1e-4)
  expect_equal(m$mean, colMeans(do.call(rbind, x)), tolerance = 1e-10)
})

test_that("covariance influence includes weighted singleton and sibling contributions", {
  x <- .sparse_covariate_fixture(c(rep(1L, 380L), rep(2L, 20L)), q = 2L)
  m <- .ravel_mixed_covariate_moments(x)
  n <- length(x); sizes <- vapply(x, nrow, integer(1))
  for (i in c(which(sizes == 1L)[3L], which(sizes == 2L)[3L])) {
    h <- 1e-3; plus <- minus <- rep(1, n)
    plus[i] <- plus[i] + h; minus[i] <- minus[i] - h
    numerical <- (.ravel_mixed_covariate_moments(x, weights = plus, influence = FALSE)$eta -
      .ravel_mixed_covariate_moments(x, weights = minus, influence = FALSE)$eta) / (2 * h)
    expect_equal(unname(m$moment_scores[i, ] / n), unname(numerical), tolerance = 3e-5)
  }
})

test_that("balanced sibling covariances agree with their direct empirical solution", {
  x <- .sparse_covariate_fixture(rep(2L, 100L), q = 2L, seed = 175L)
  m <- .ravel_mixed_covariate_moments(x)
  C <- do.call(rbind, x); mean <- colMeans(C)
  self <- crossprod(C) / nrow(C) - tcrossprod(mean)
  sibling <- Reduce(`+`, lapply(x, function(z)
    tcrossprod(colSums(z)) - crossprod(z))) / nrow(C) - tcrossprod(mean)
  expect_equal(unname(m$Qself), unname(self), tolerance = 1e-5)
  expect_equal(unname(m$Qsib), unname(sibling), tolerance = 1e-5)
})

test_that("a truly family-shared covariate retains its zero within-family support", {
  x <- .sparse_covariate_fixture(c(rep(1L, 380L), rep(2L, 20L)), shared = 2L)
  m <- .ravel_mixed_covariate_moments(x)
  expect_identical(as.integer(m$shared), 2L)
  expect_lt(max(abs((m$Qself - m$Qsib)[2L, ])), 1e-10)
  expect_true(all(is.finite(m$moment_scores)))
  expect_gt(m$Qself[2L, 2L], 0)
})
