# Consolidate completed runs without dropping boundary or unsuccessful attempts.
# Rscript dev/summarize-mixed-validation.R
root <- file.path("dev", "mixed-validation")
jobs <- list(
  list(path = "mc-continuous", expected = c(cross_state = 500L, extended_family = 500L)),
  list(path = "pilot-interactions", expected = c(long_age_sex_random_slope = 50L)),
  list(path = "pilot-binary", expected = c(binary_siblings = 20L)))
status <- recovery <- list()
for (job in jobs) {
  path <- file.path(root, job$path)
  a <- utils::read.csv(file.path(path, "attempts.csv"))
  e <- utils::read.csv(file.path(path, "estimates.csv"))
  for (scenario in names(job$expected)) {
    aa <- a[a$scenario == scenario, ]
    if (nrow(aa) != job$expected[[scenario]])
      stop(scenario, " is incomplete: ", nrow(aa), "/", job$expected[[scenario]], " attempts.")
    status[[length(status) + 1L]] <- data.frame(scenario = scenario,
      evidence_directory = job$path, attempted = nrow(aa), errors = sum(aa$error),
      optimizer_converged = sum(aa$optimizer_converged), converged = sum(aa$converged),
      finite_inference = sum(aa$inference_status == "ok", na.rm = TRUE),
      boundary = sum(aa$inference_status == "boundary", na.rm = TRUE),
      other_unavailable = sum(is.na(aa$inference_status) | !aa$inference_status %in% c("ok", "boundary")),
      total_fit_seconds = sum(aa$elapsed_seconds), median_fit_seconds = stats::median(aa$elapsed_seconds))
    ee <- e[e$scenario == scenario & e$scale %in% c("corrected_std", "variance"), ]
    for (key in unique(paste(ee$scale, ee$parameter))) {
      x <- ee[paste(ee$scale, ee$parameter) == key, ]
      point <- is.finite(x$estimate)
      valid <- x$available
      npoint <- sum(point); nvalid <- sum(valid)
      stopifnot(nrow(x) == job$expected[[scenario]])
      error <- x$estimate[point] - x$truth[point]
      coverage <- if (nvalid) mean(x$covered[valid]) else NA_real_
      ci <- if (nvalid) stats::binom.test(sum(x$covered[valid]), nvalid)$conf.int else c(NA, NA)
      null <- all(x$truth == 0)
      rejected <- if (null && nvalid) sum(x$reject_zero[valid]) else NA_integer_
      null_ci <- if (null && nvalid) stats::binom.test(rejected, nvalid)$conf.int else c(NA, NA)
      recovery[[length(recovery) + 1L]] <- data.frame(scenario = scenario,
        scale = x$scale[1L], parameter = x$parameter[1L], truth = x$truth[1L],
        attempted = nrow(x), finite_point_estimates = npoint, finite_SE = nvalid,
        bias_all_finite_points = if (npoint) mean(error) else NA_real_,
        bias_MCSE_all_finite_points = if (npoint > 1L) stats::sd(error) / sqrt(npoint) else NA_real_,
        empirical_SD_all_finite_points = if (npoint > 1L) stats::sd(x$estimate[point]) else NA_real_,
        mean_SE = if (nvalid) mean(x$se[valid]) else NA_real_,
        coverage_among_available = coverage,
        coverage_MCSE = if (nvalid) sqrt(coverage * (1 - coverage) / nvalid) else NA_real_,
        coverage_exact_low = ci[1L], coverage_exact_high = ci[2L],
        coverage_all_attempts = mean(x$covered),
        null_rejections = rejected,
        null_rejection_among_available = if (null && nvalid) rejected / nvalid else NA_real_,
        null_rejection_exact_low = null_ci[1L], null_rejection_exact_high = null_ci[2L])
    }
  }
}
destination <- file.path(root, "consolidated")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
status <- do.call(rbind, status)
recovery <- do.call(rbind, recovery)
utils::write.csv(status, file.path(destination, "fit-status.csv"), row.names = FALSE)
utils::write.csv(recovery, file.path(destination, "recovery-coverage.csv"), row.names = FALSE)
print(status, row.names = FALSE)
cat("Covariance guards differ by recorded source version; consult each run's session.txt\n",
    "and final-certification replay tables before extrapolating inference availability.\n", sep = "")
