# Replay selected historical Monte Carlo seeds on the current inference guards.
# The original attempt tables and source hashes are never overwritten.
# Rscript dev/recheck-mixed-certification.R
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-mixed-validation.R")
arguments <- commandArgs(trailingOnly = TRUE)
output_argument <- arguments[startsWith(arguments, "--output=")]
destination <- if (length(output_argument)) sub("^--output=", "", output_argument[1L]) else
  file.path("dev", "mixed-validation", "final-certification-replay")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
source_hashes <- tools::md5sum(c("dev/recheck-mixed-certification.R",
  "tests/testthat/helper-mixed-validation.R", list.files("R", pattern = "mixed", full.names = TRUE)))
scenarios <- list(
  cross_state = list(path = "mc-continuous", data = list(), fit = list()),
  extended_family = list(path = "mc-continuous",
    data = list(families_per_FID = 3L, extended_fraction = .6), fit = list()),
  long_age_sex_random_slope = list(path = "pilot-interactions",
    data = list(visits = 4L, individual_variance = .2, slope_variance = .003,
      intercept_slope_covariance = .012,
      age_beta = c(child = .025, mother = -.008, father = 0),
      sex_beta = c(child = .12, mother = -.06, father = .02)),
    fit = list(individual_random = TRUE, random_slope = "age", interactions = c("age", "sex"),
      interaction_center = c(age = 10, sex = 0))))
chosen <- arguments[startsWith(arguments, "--scenarios=")]
if (length(chosen)) scenarios <- scenarios[strsplit(sub("^--scenarios=", "", chosen[1L]), ",", fixed = TRUE)[[1L]]]
stopifnot(length(scenarios) > 0L, all(vapply(scenarios, is.list, logical(1))))
results <- parameters <- list()
for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  path <- file.path("dev", "mixed-validation", scenario$path)
  attempts <- utils::read.csv(file.path(path, "attempts.csv"))
  attempts <- attempts[attempts$scenario == scenario_name, ]
  if (!nrow(attempts)) next
  historical <- utils::read.csv(file.path(path, "estimates.csv"))
  historical <- historical[historical$scenario == scenario_name, ]
  selected <- unique(c(1L, nrow(attempts),
    which(is.na(attempts$inference_status) | attempts$inference_status != "ok"),
    head(order(attempts$projected_gradient_max, decreasing = TRUE), 2L)))
  variances <- historical[historical$scale == "variance" & historical$parameter !=
    "individual_intercept_slope_covariance" & is.finite(historical$estimate), ]
  for (component in unique(variances$parameter)) {
    v <- variances[variances$parameter == component, ]
    selected <- union(selected, match(v$replicate[which.min(v$estimate)], attempts$replicate))
  }
  for (row in selected) {
    a <- attempts[row, ]
    d <- do.call(mixed_validation_data, modifyList(list(N_families = a$N_families,
      seed = a$seed, age_pattern = "scheduled", age_effect = 0, sex_effect = 0), scenario$data))
    args <- modifyList(list(dat = d, pi = a$pi_supplied,
      covariate_mapping_confirmed = TRUE, diagnostic = FALSE), scenario$fit)
    f <- tryCatch(do.call(ravel_mixed_fit, args), error = identity)
    error <- inherits(f, "error")
    difference <- NA_real_
    if (!error) {
      old <- historical[historical$replicate == a$replicate & historical$scale != "variance", ]
      idx <- match(paste(old$scale, old$parameter), paste(f$table$scale, f$table$parameter))
      differences <- abs(old$estimate - f$table$estimate[idx])
      if (any(is.finite(differences))) difference <- max(differences, na.rm = TRUE)
      parameters[[length(parameters) + 1L]] <- data.frame(scenario = scenario_name,
        replicate = a$replicate, parameter = names(f$parameters), estimate = unname(f$parameters))
    }
    results[[length(results) + 1L]] <- data.frame(scenario = scenario_name,
      replicate = a$replicate, seed = a$seed, historical_status = a$inference_status,
      current_status = if (error) "error" else f$inference_status,
      current_converged = !error && f$converged,
      natural_variance_score_max = if (error || !length(f$projected_natural_variance_gradient))
        NA_real_ else max(abs(f$projected_natural_variance_gradient)),
      boundary_parameters = if (error) "" else paste(f$boundary_parameters, collapse = ";"),
      max_effect_difference = difference,
      note = if (error) conditionMessage(f) else paste(f$messages, collapse = ";"))
    utils::write.csv(do.call(rbind, results), file.path(destination, "attempts.csv"), row.names = FALSE)
    if (length(parameters)) utils::write.csv(do.call(rbind, parameters),
      file.path(destination, "parameters.csv"), row.names = FALSE)
    cat(scenario_name, a$replicate, "historical", a$inference_status,
        "current", if (error) "error" else f$inference_status, "\n")
  }
}
writeLines(c(capture.output(sessionInfo()), "", "Source hashes:", capture.output(source_hashes)),
           file.path(destination, "session.txt"))
