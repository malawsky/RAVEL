# Deterministic independent checks for the optional mixed likelihood.
# Run from the package root: Rscript dev/check-mixed-nesting-pi.R
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-mixed-validation.R")
destination <- file.path("dev", "mixed-validation", "nesting-pi-check")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
source_hashes <- tools::md5sum(c("dev/check-mixed-nesting-pi.R",
  "tests/testthat/helper-mixed-validation.R", list.files("R", pattern = "mixed", full.names = TRUE)))

checks <- list()
modes <- c("estimated_state", "equilibrium", "n_gen")
seeds <- c(328001L, 358001L, 368001L)
for (i in seq_along(modes)) {
  d <- mixed_validation_data(N_families = 500L, seed = seeds[i], mode = modes[i],
    n_generations = 3L, age_pattern = "scheduled", age_effect = 0, sex_effect = 0)
  set.seed(seeds[i] + 1L)
  supplied_pi <- stats::rnorm(1L, .4, .03)
  f <- ravel_mixed_fit(d, pi = supplied_pi, pi_se = .03, fit_constraints = modes[i],
    am_n_generations = 3L, diagnostic = FALSE)
  stopifnot(f$converged, identical(f$inference_status, "ok"))
  # A wider step and tighter independent optimization check sensitivity to
  # numerical optimizer noise in the default local pi derivative.
  sides <- lapply(c(-.005, .005), function(delta)
    ravel_mixed_fit(d, pi = supplied_pi + delta, fit_constraints = modes[i],
      am_n_generations = 3L, diagnostic = FALSE, se = FALSE,
      control = list(factr = 10), start = f$parameters))
  stopifnot(all(vapply(sides, function(x) x$converged, logical(1))))
  derivative <- (sides[[2L]]$table$estimate - sides[[1L]]$table$estimate) / .01
  default <- unlist(f$pi_sensitivity, use.names = FALSE)
  checks[[i]] <- data.frame(mode = modes[i], parameter = f$table$parameter,
    scale = f$table$scale, default_derivative = default,
    independent_derivative = derivative, difference = default - derivative,
    SE_discrepancy_ratio = .03 * abs(default - derivative) / f$table$se)
}
checks <- do.call(rbind, checks)
utils::write.csv(checks, file.path(destination, "pi-sensitivity.csv"), row.names = FALSE)
stopifnot(max(checks$SE_discrepancy_ratio) < .01)

d <- mixed_validation_data(N_families = 1200L, sib_sizes = 1L,
  family_variance = 0, individual_variance = 0, age_pattern = "scheduled",
  age_effect = 0, sex_effect = 0, seed = 183L)
original <- ravel_fit(d, pi = .4, fit_constraints = "estimated_state",
                      am_correction = "estimated_state")
mixed <- ravel_mixed_fit(d, pi = .4, family_random = FALSE,
                        fit_constraints = "estimated_state", diagnostic = FALSE)
stopifnot(original$converged, mixed$converged)
a <- as.data.frame(original, include_raw = TRUE)
a <- a[a$parameter %in% c("beta_D", "beta_M", "beta_F") &
  a$scale %in% c("visible", "corrected_std", "corrected_raw"), ]
idx <- match(paste(a$scale, a$parameter), paste(mixed$table$scale, mixed$table$parameter))
comparison <- data.frame(scale = a$scale, parameter = a$parameter,
  original_estimate = a$estimate, mixed_estimate = mixed$table$estimate[idx],
  difference = mixed$table$estimate[idx] - a$estimate)
utils::write.csv(comparison, file.path(destination, "no-random-original-comparison.csv"), row.names = FALSE)
stopifnot(max(abs(comparison$difference)) < 3e-5)
writeLines(c(capture.output(sessionInfo()), "", "Source hashes:", capture.output(source_hashes)),
           file.path(destination, "session.txt"))
cat("Maximum pi-derivative discrepancy / reported SE:", max(checks$SE_discrepancy_ratio), "\n")
cat("Maximum old/new no-random effect discrepancy:", max(abs(comparison$difference)), "\n")
