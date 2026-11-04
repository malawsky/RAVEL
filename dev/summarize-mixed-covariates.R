# Summarize every scheduled simulation attempt, including failed inference.
# The source directories retain original results; explicit retries are reported
# separately and never overwrite the original Monte Carlo denominator.
directories <- commandArgs(trailingOnly = TRUE)
if (!length(directories)) directories <- file.path("dev", "mixed-covariate-validation",
  c("mc-continuous", "longitudinal", "binary-pilot"))
read_results <- function(directory, table) {
  path <- file.path(directory, paste0(table, ".csv"))
  if (!file.exists(path)) return(NULL)
  out <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (table == "attempts" && !"source_run_id" %in% names(out))
    out$source_run_id <- "historical_pre_natural_polish"
  out$run <- basename(directory)
  out
}
attempts <- do.call(rbind, lapply(directories, read_results, table = "attempts"))
estimates <- do.call(rbind, lapply(directories, read_results, table = "estimates"))
states <- do.call(rbind, lapply(directories, read_results, table = "states"))
stopifnot(!is.null(attempts), !is.null(estimates), !is.null(states))
interval <- function(success, n) {
  if (!n) return(c(estimate = NA_real_, mcse = NA_real_, lower = NA_real_, upper = NA_real_))
  p <- success / n
  c(estimate = p, mcse = sqrt(p * (1 - p) / n),
    lower = if (!success) 0 else stats::qbeta(.025, success, n - success + 1),
    upper = if (success == n) 1 else stats::qbeta(.975, success + 1, n - success))
}
attempt_summary <- do.call(rbind, lapply(split(attempts,
  interaction(attempts$run, attempts$scenario, drop = TRUE)), function(x) {
  data.frame(run = x$run[1L], scenario = x$scenario[1L], attempted = nrow(x),
    errors = sum(x$error), converged = sum(x$converged),
    regular_inference = sum(x$inference_status == "ok"),
    statuses = paste(names(table(x$inference_status)), as.integer(table(x$inference_status)),
      sep = ":", collapse = "; "),
    total_fit_seconds = sum(x$elapsed_seconds), median_fit_seconds = median(x$elapsed_seconds))
}))
source_summary <- do.call(rbind, lapply(split(attempts,
  interaction(attempts$run, attempts$scenario, attempts$source_run_id, drop = TRUE)), function(x) {
  data.frame(run = x$run[1L], scenario = x$scenario[1L], source_run_id = x$source_run_id[1L],
    attempted = nrow(x), errors = sum(x$error), converged = sum(x$converged),
    regular_inference = sum(x$inference_status == "ok"),
    statuses = paste(names(table(x$inference_status)), as.integer(table(x$inference_status)),
      sep = ":", collapse = "; "),
    first_seed = min(x$seed), last_seed = max(x$seed))
}))
coefficient_summary <- do.call(rbind, lapply(split(estimates,
  interaction(estimates$run, estimates$scenario, estimates$scale, estimates$parameter,
    drop = TRUE)), function(x) {
  finite <- x$finite_point & is.finite(x$estimate)
  inference <- x$finite_inference & finite
  error <- x$estimate[finite] - x$truth[finite]
  coverage_all <- interval(sum(x$covered), nrow(x))
  coverage_regular <- interval(sum(x$covered), sum(inference))
  row <- data.frame(run = x$run[1L], scenario = x$scenario[1L], scale = x$scale[1L],
    parameter = x$parameter[1L], truth = x$truth[1L], attempted = nrow(x),
    finite_point = sum(finite), finite_inference = sum(inference), covered = sum(x$covered),
    bias = if (length(error)) mean(error) else NA_real_,
    bias_mcse = if (length(error) > 1L) sd(error) / sqrt(length(error)) else NA_real_,
    empirical_sd = if (sum(finite) > 1L) sd(x$estimate[finite]) else NA_real_,
    mean_se = if (any(inference)) mean(x$se[inference]) else NA_real_,
    rmse = if (length(error)) sqrt(mean(error^2)) else NA_real_)
  row[paste0("all_attempt_coverage_", names(coverage_all))] <- as.list(coverage_all)
  row[paste0("conditional_coverage_", names(coverage_regular))] <- as.list(coverage_regular)
  row
}))
states <- merge(states, attempts[c("run", "scenario", "replicate", "converged")],
  by = c("run", "scenario", "replicate"), all.x = TRUE)
state_summary <- do.call(rbind, lapply(split(states,
  interaction(states$run, states$scenario, states$parameter, drop = TRUE)), function(x) {
  good <- x$converged & is.finite(x$estimate)
  error <- x$estimate[good] - x$truth[good]
  data.frame(run = x$run[1L], scenario = x$scenario[1L], parameter = x$parameter[1L],
    truth = x$truth[1L], attempted = nrow(x), finite_point = sum(good),
    bias = if (length(error)) mean(error) else NA_real_,
    bias_mcse = if (length(error) > 1L) sd(error) / sqrt(length(error)) else NA_real_,
    rmse = if (length(error)) sqrt(mean(error^2)) else NA_real_)
}))
destination <- file.path("dev", "mixed-covariate-validation", "summary")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
for (name in c("attempt_summary", "source_summary", "coefficient_summary", "state_summary"))
  utils::write.csv(get(name), file.path(destination, paste0(name, ".csv")), row.names = FALSE)
writeLines(c(paste("Inputs:", paste(directories, collapse = ", ")),
  "All-attempt coverage counts failed/nonregular fits as uncovered.",
  "Conditional coverage uses only finite, regular Wald intervals.",
  "Coverage intervals are exact binomial 95% intervals; MCSE uses observed binomial variance.",
  "Bias and empirical SD use converged finite point estimates and retain their denominator.",
  "Historical guards and final-source retries must be identified using attempts source_run_id.",
  capture.output(sessionInfo())), file.path(destination, "session.txt"))
print(attempt_summary, row.names = FALSE)
print(coefficient_summary[coefficient_summary$scale == "corrected_std", ], row.names = FALSE)
