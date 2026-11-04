# Independent review check: unbalanced families with negative sibling covariance.
# Run from the package root: Rscript dev/check-negative-sibling-covariates.R
source("R/ravel_mixed_covariate_moments.R")
destination <- "dev/mixed-validation/sparse-sibling-covariates/negative-sibling-review"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
set.seed(661)
records <- list(); covariances <- list()
for (q in c(1L, 2L)) {
  sizes <- c(rep(1L, 475L), rep(2L, 15L), rep(3L, 10L))
  family_values <- lapply(sizes, function(n) {
    Q <- kronecker(diag(n), 1.3 * diag(q)) +
      kronecker(matrix(1, n, n), -.3 * diag(q))
    matrix(drop(t(chol(Q)) %*% rnorm(n * q)), n, q, byrow = TRUE)
  })
  fitted <- .ravel_mixed_covariate_moments(family_values)
  for (i in c(1L, 476L, 491L)) {
    step <- .001
    plus <- minus <- rep(1, length(family_values))
    plus[i] <- 1 + step; minus[i] <- 1 - step
    ep <- .ravel_mixed_covariate_moments(family_values, plus, FALSE)$eta
    em <- .ravel_mixed_covariate_moments(family_values, minus, FALSE)$eta
    numerical <- (ep - em) / (2 * step)
    analytical <- fitted$moment_scores[i, ] / length(family_values)
    error <- max(abs(numerical - analytical))
    records[[length(records) + 1L]] <- data.frame(
      q = q, family_size = sizes[i], family_index = i,
      max_abs_influence_error = error, tolerance = 1e-6,
      status = if (is.finite(error) && error < 1e-6) "pass" else "FAIL")
  }
  ev <- min(eigen(fitted$Qself + 2 * fitted$Qsib, symmetric = TRUE)$values)
  covariances[[length(covariances) + 1L]] <- data.frame(
    q = q, covariate = seq_len(q), true_sibling_covariance = -.3,
    estimated_sibling_covariance = diag(fitted$Qsib),
    maximum_family_mean_min_eigenvalue = ev,
    status = if (all(diag(fitted$Qsib) < 0) && ev > 0) "pass" else "FAIL")
}
checks <- do.call(rbind, records)
covariances <- do.call(rbind, covariances)
write.csv(checks, file.path(destination, "influence-checks.csv"), row.names = FALSE)
write.csv(covariances, file.path(destination, "covariances.csv"), row.names = FALSE)
hashes <- tools::md5sum(c("R/ravel_mixed_covariate_moments.R",
  "dev/check-negative-sibling-covariates.R"))
writeLines(sub("[[:space:]]+$", "", c(capture.output(sessionInfo()), capture.output(hashes))),
  file.path(destination, "session.txt"))
print(checks, row.names = FALSE)
print(covariances, row.names = FALSE)
stopifnot(all(checks$status == "pass"), all(covariances$status == "pass"))
