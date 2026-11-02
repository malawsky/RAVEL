# Known-measurement positive control for the binary auxiliary diagnostic.
# Run from the repository root; regenerates the same 500 independent samples
# as dev/run-binary-score-calibration.R, but conditions on the true A matrix.
C <- matrix(c(1, .6, .6, .6, 1, .2, .6, .2, 1), 3L)
K <- matrix(c(1, .5, .5, .5, 1, 0, .5, 0, 1), 3L)
C_chol <- chol(C); K_chol <- chol(K)
score_covariance <- rbind(cbind(C + .6 * K, C), cbind(C, C + 1.3 * K))
A <- solve(score_covariance, rbind(C, C))
restriction <- cbind(-A[4:6, , drop = FALSE] %*% solve(A[1:3, , drop = FALSE]), diag(3L))
n <- 1500L
p_values <- numeric(500L)
for (replicate in seq_len(500L)) {
  set.seed(620000L + replicate)
  genetic <- matrix(stats::rnorm(n * 3L), n, 3L) %*% C_chol
  error1 <- matrix(stats::rnorm(n * 3L), n, 3L) %*% K_chol
  error2 <- matrix(stats::rnorm(n * 3L), n, 3L) %*% K_chol
  y <- as.numeric(drop(genetic %*% c(.4, .1, .05)) + stats::rnorm(n) > .6)
  scores <- cbind(genetic + sqrt(.6) * error1, genetic + sqrt(1.3) * error2)
  X <- cbind(intercept = 1, scores)
  fit <- stats::glm.fit(X, y, family = stats::binomial("probit"))
  stopifnot(fit$converged, fit$rank == ncol(X), all(is.finite(fit$coefficients)))
  eta <- drop(X %*% fit$coefficients)
  probability <- stats::pnorm(eta); density <- stats::dnorm(eta)
  variance <- probability * (1 - probability); residual <- y - probability
  score <- X * as.vector(density * residual / variance)
  weights <- density^2 / variance + eta * density * residual / variance +
    density^2 * residual * (1 - 2 * probability) / variance^2
  bread <- solve(crossprod(X, X * as.vector(weights)))
  covariance <- bread %*% crossprod(score) %*% bread
  discrepancy <- drop(restriction %*% fit$coefficients[-1L])
  tested_covariance <- restriction %*% covariance[-1L, -1L] %*% t(restriction)
  statistic <- drop(crossprod(discrepancy, solve(tested_covariance, discrepancy)))
  p_values[replicate] <- stats::pchisq(statistic, 3L, lower.tail = FALSE)
}
output <- "dev/mixed-validation/binary-diagnostic-oracle"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(replicate = seq_len(500L), p = p_values),
  file.path(output, "attempts.csv"), row.names = FALSE)
interval <- stats::binom.test(sum(p_values < .05), 500L)$conf.int
summary <- data.frame(attempts = 500L, rejected = sum(p_values < .05), rate = mean(p_values < .05),
  ci_low = interval[1L], ci_high = interval[2L])
utils::write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
print(summary)
