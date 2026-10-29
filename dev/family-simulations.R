# Independent family-model Monte Carlo validation.
# Run from the package root, for example:
# Rscript dev/family-simulations.R --configuration=trio --reps=100 --n=2000
# Every replicate, including failures, is retained in the output CSV.
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", tail(hit, 1))
}
configuration <- option("configuration", "trio")
scenario <- option("scenario", "baseline")
reps <- as.integer(option("reps", "10"))
N <- as.integer(option("n", "2000"))
seed0 <- as.integer(option("seed", "20261029"))
pi_se <- as.numeric(option("pi-se", "0"))
outdir <- option("output", "dev/family-simulation-results")
stopifnot(reps >= 1, N >= 100, is.finite(pi_se), pi_se >= 0)

pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-family-validation.R")
population_args <- switch(scenario,
  baseline = list(pi = .4, rho = .3, parent_variance = 1.15),
  changed = list(pi = .4, rho = -.2, parent_variance = 1.5),
  random_mating = list(pi = .4, rho = 0, parent_variance = 1),
  high_pi = list(pi = .7, rho = .3, parent_variance = 1.15),
  stop("Unknown scenario: ", scenario))
p <- do.call(family_validation_population, population_args)
genetic_roles <- switch(configuration,
  trio = p$roles,
  siblings = c("mother", "mat_sib"),
  adult_one = c("mother", "father", "mat_sib"),
  adult_both = c("mother", "father", "mat_sib", "pat_sib"),
  stop("Unknown configuration: ", configuration))
targets <- list(sibling_pair = c("mother", "mat_sib"))
if (length(genetic_roles) > 2L)
  targets$all_adults <- setdiff(genetic_roles, "child")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
stem <- sprintf("%s-%s-N%d-pise%s-seed%d", configuration, scenario, N,
                format(pi_se, scientific = FALSE, trim = TRUE), seed0)
replicate_path <- file.path(outdir, paste0(stem, "-replicates.csv"))
summary_path <- file.path(outdir, paste0(stem, "-summary.csv"))
metadata_path <- file.path(outdir, paste0(stem, "-metadata.txt"))
writeLines(c(sprintf("configuration: %s", configuration),
             sprintf("scenario: %s", scenario), sprintf("reps: %d", reps),
             sprintf("N: %d", N), sprintf("seed0: %d", seed0),
             sprintf("pi: %.8g; pi_se: %.8g", p$pi, pi_se),
             sprintf("rho_total: %.8g; parent_variance: %.8g", p$rho,
                     p$parent_variance),
             "Generator: independent Gaussian pedigree innovations; exact=FALSE",
             "R source hashes:",
             capture.output(tools::md5sum(list.files("R", pattern = "\\.R$",
                                                      full.names = TRUE))),
             capture.output(sessionInfo())), metadata_path)

records <- vector("list", reps)
for (rep in seq_len(reps)) {
  seed <- seed0 + rep - 1L
  dat <- family_validation_data(p, N = N, seed = seed, exact = FALSE)
  # Simulate the independent calibration study when pi uncertainty is used.
  set.seed(seed0 + 100000L + rep)
  fitted_pi <- if (pi_se > 0) stats::rnorm(1, p$pi, pi_se) else p$pi
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  fitted <- if (fitted_pi < 0 || fitted_pi >= 1) {
    simpleError("Independent pi estimate was outside [0,1); replicate retained")
  } else tryCatch(withCallingHandlers(
    ravel_family_fit(dat, genetic_roles = genetic_roles,
      outcome_roles = targets[[1]], pi = fitted_pi,
      pi_se = if (pi_se > 0) pi_se else NULL,
      fit_constraints = "estimated_state"),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = identity)
  fit_time <- proc.time()[["elapsed"]] - start
  rows <- list()
  for (target in names(targets)) {
    roles <- targets[[target]]
    oracle <- family_validation_oracle(p, roles)
    truth <- c(oracle$beta_raw, oracle$contrasts_raw)
    status <- if (inherits(fitted, "error")) "fit_error" else
      if (!isTRUE(fitted$converged)) "nonconverged" else "ok"
    reason <- if (inherits(fitted, "error")) conditionMessage(fitted) else
      if (!isTRUE(fitted$converged)) paste(fitted$error, collapse = "; ") else ""
    z <- if (status == "ok") tryCatch(withCallingHandlers(
      ravel_family_regression(fitted, outcome_roles = roles),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = identity) else NULL
    if (inherits(z, "error")) {
      status <- "projection_error"
      reason <- conditionMessage(z)
    }
    for (parameter in names(truth)) {
      estimate <- se <- NA_real_
      parameter_status <- status
      if (status == "ok") {
        values <- c(z$coefficients$total_raw, z$contrasts$total_raw)
        estimate <- unname(values[parameter])
        key <- paste0("total_raw:", parameter)
        if (!is.null(z$vcov) && key %in% rownames(z$vcov))
          se <- sqrt(z$vcov[key, key])
        if (!is.finite(estimate)) parameter_status <- "correction_failed" else
          if (!is.finite(se) || se <= 0) parameter_status <- "inference_failed"
      }
      rows[[length(rows) + 1L]] <- data.frame(
        configuration = configuration, scenario = scenario,
        outcome_model = target, N = N, pi = p$pi, pi_se = pi_se,
        replicate = rep, seed = seed, fitted_pi = fitted_pi,
        parameter = parameter, truth = unname(truth[parameter]),
        estimate = estimate, se = se, status = parameter_status,
        reason = reason, warnings = paste(unique(warnings), collapse = " | "),
        fit_seconds = fit_time, stringsAsFactors = FALSE)
    }
  }
  records[[rep]] <- do.call(rbind, rows)
  utils::write.csv(do.call(rbind, records[seq_len(rep)]), replicate_path,
                   row.names = FALSE)
  cat(sprintf("%s / %s: %d/%d, fit %.1fs, statuses %s\n", configuration,
              scenario, rep, reps, fit_time,
              paste(unique(records[[rep]]$status), collapse = ",")))
  flush.console()
}

draws <- do.call(rbind, records)
groups <- split(draws, interaction(draws$outcome_model, draws$parameter,
                                   drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(d) {
  keep <- d$status == "ok" & is.finite(d$estimate) & is.finite(d$se)
  good <- d[keep, , drop = FALSE]
  n <- nrow(good)
  errors <- good$estimate - good$truth
  covered <- abs(errors) <= stats::qnorm(.975) * good$se
  coverage <- if (n) mean(covered) else NA_real_
  coverage_ci <- if (n) stats::binom.test(sum(covered), n)$conf.int else c(NA, NA)
  data.frame(configuration = configuration, scenario = scenario,
    outcome_model = d$outcome_model[1], parameter = d$parameter[1],
    N = N, pi = p$pi, pi_se = pi_se, attempted = nrow(d), valid = n,
    failed = nrow(d) - n, failure_rate = 1 - n / nrow(d),
    truth = d$truth[1], bias = if (n) mean(errors) else NA_real_,
    bias_mcse = if (n > 1L) stats::sd(errors) / sqrt(n) else NA_real_,
    rmse = if (n) sqrt(mean(errors^2)) else NA_real_,
    empirical_sd = if (n > 1L) stats::sd(good$estimate) else NA_real_,
    mean_se = if (n) mean(good$se) else NA_real_, coverage = coverage,
    coverage_mcse = if (n) sqrt(coverage * (1 - coverage) / n) else NA_real_,
    coverage_ci_low = coverage_ci[1], coverage_ci_high = coverage_ci[2],
    successful_coverage_all_attempts = sum(covered) / nrow(d),
    median_fit_seconds = stats::median(d$fit_seconds), row.names = NULL)
}))
utils::write.csv(summary, summary_path, row.names = FALSE)
print(summary, row.names = FALSE)
cat("Replicates:", replicate_path, "\nSummary:", summary_path, "\n")
