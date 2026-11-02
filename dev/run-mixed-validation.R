# Independent, reproducible validation of the optional family mixed likelihood.
# Examples:
# Rscript dev/run-mixed-validation.R --stage=smoke
# Rscript dev/run-mixed-validation.R --stage=pilot --scenarios=cross_state --replicates=50
# Rscript dev/run-mixed-validation.R --stage=mc --replicates=500 --families=1000
# Every attempted fit is retained; failures never disappear from denominators.

pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-mixed-validation.R")
source_hashes <- tools::md5sum(c("tests/testthat/helper-mixed-validation.R",
  "dev/run-mixed-validation.R", list.files("R", pattern = "mixed", full.names = TRUE)))
arguments <- commandArgs(trailingOnly = TRUE)
option <- function(key, default) {
  value <- arguments[startsWith(arguments, paste0("--", key, "="))]
  if (!length(value)) default else sub(paste0("^--", key, "="), "", value[[length(value)]])
}
stage <- option("stage", "smoke")
stopifnot(stage %in% c("smoke", "pilot", "mc"))
replicates <- as.integer(option("replicates", switch(stage, smoke = "1", pilot = "50", mc = "500")))
families <- as.integer(option("families", switch(stage, smoke = "500", pilot = "500", mc = "1000")))
output <- option("output", file.path("dev", "mixed-validation", stage))
stopifnot(is.finite(replicates), replicates >= 1L, is.finite(families), families >= 50L)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
scenarios <- list(
  cross_state = list(data = list(), fit = list()),
  cross_equilibrium = list(data = list(mode = "equilibrium"), fit = list(fit_constraints = "equilibrium")),
  cross_n_gen = list(data = list(mode = "n_gen", n_generations = 3L),
    fit = list(fit_constraints = "n_gen", am_n_generations = 3L)),
  cross_sex = list(data = list(sex_beta = c(child = .12, mother = -.06, father = .02)),
    fit = list(interactions = "sex", interaction_center = c(sex = 0))),
  long_intercept = list(data = list(visits = 3L, individual_variance = .2),
    fit = list(individual_random = TRUE)),
  long_irregular = list(data = list(visits = 3L, individual_variance = .2,
    age_pattern = "random", age_beta = c(child = .025, mother = -.008, father = 0)),
    fit = list(individual_random = TRUE, interactions = "age", interaction_center = c(age = 10))),
  long_age = list(data = list(visits = 3L, individual_variance = .2,
    age_beta = c(child = .025, mother = -.008, father = 0)),
    fit = list(individual_random = TRUE, interactions = "age", interaction_center = c(age = 10))),
  long_age_sex = list(data = list(visits = 3L, individual_variance = .2,
    age_beta = c(child = .025, mother = -.008, father = 0),
    sex_beta = c(child = .12, mother = -.06, father = .02)),
    fit = list(individual_random = TRUE, interactions = c("age", "sex"),
      interaction_center = c(age = 10, sex = 0))),
  long_random_slope = list(data = list(visits = 4L, individual_variance = .2,
    slope_variance = .003, intercept_slope_covariance = .012,
    age_beta = c(child = .025, mother = -.008, father = 0)),
    fit = list(individual_random = TRUE, random_slope = "age", interactions = "age",
      interaction_center = c(age = 10))),
  extended_family = list(data = list(families_per_FID = 3L, extended_fraction = .6), fit = list()),
  pi_uncertainty = list(data = list(), fit = list(pi_se = .03), external_pi_sd = .03),
  binary_siblings = list(data = list(outcome = "binary_liability", residual_variance = 1,
    prevalence_reference = .15), fit = list(outcome = "binary_liability")),
  binary_sex = list(data = list(outcome = "binary_liability", residual_variance = 1,
    prevalence_reference = .15, sex_beta = c(child = .12, mother = -.06, father = .02)),
    fit = list(outcome = "binary_liability", interactions = "sex", interaction_center = c(sex = 0))),
  pi_uncertainty_equilibrium = list(data = list(mode = "equilibrium"),
    fit = list(fit_constraints = "equilibrium", pi_se = .03), external_pi_sd = .03),
  pi_uncertainty_n_gen = list(data = list(mode = "n_gen", n_generations = 3L),
    fit = list(fit_constraints = "n_gen", am_n_generations = 3L, pi_se = .03), external_pi_sd = .03)
)
selected <- strsplit(option("scenarios", if (stage == "smoke")
  "cross_state,cross_equilibrium,cross_n_gen,cross_sex,long_age_sex,extended_family" else "cross_state"),
  ",", fixed = TRUE)[[1L]]
stopifnot(length(selected) > 0L, all(selected %in% names(scenarios)))
if (!exists("ravel_mixed_fit", mode = "function"))
  stop("ravel_mixed_fit is not yet available; the independent generator tests can run separately.")

coefficient_table <- function(fit) {
  if (is.data.frame(fit$table)) return(fit$table)
  if (is.data.frame(fit$coefficients)) return(fit$coefficients)
  tryCatch(as.data.frame(fit, include_raw = TRUE), error = function(e) data.frame())
}
expected_coefficients <- function(truth, interactions) {
  wanted <- c("beta_D", "beta_M", "beta_F")
  if (length(interactions)) wanted <- c(wanted,
    unlist(lapply(interactions, function(x) paste0(c("beta_D", "beta_M", "beta_F"), ":", x))))
  do.call(rbind, lapply(names(truth$coefficients), function(scale) {
    x <- truth$coefficients[[scale]][wanted]
    data.frame(scale = scale, parameter = names(x), truth = as.numeric(x), stringsAsFactors = FALSE)
  }))
}
attempts <- list(); estimates <- list(); variances <- list()
write_progress <- function() {
  if (length(attempts)) utils::write.csv(do.call(rbind, attempts), file.path(output, "attempts.csv"), row.names = FALSE)
  if (length(estimates)) utils::write.csv(do.call(rbind, estimates), file.path(output, "estimates.csv"), row.names = FALSE)
  if (length(variances)) utils::write.csv(do.call(rbind, variances), file.path(output, "variance-components.csv"), row.names = FALSE)
}
for (scenario_name in selected) for (replicate in seq_len(replicates)) {
  scenario <- scenarios[[scenario_name]]
  seed <- 218000L + match(scenario_name, names(scenarios)) * 10000L + replicate
  data_arguments <- modifyList(list(N_families = families, seed = seed,
    age_pattern = "scheduled", age_effect = 0, sex_effect = 0), scenario$data)
  d <- do.call(mixed_validation_data, data_arguments)
  truth <- attr(d, "truth")
  fit_arguments <- modifyList(list(dat = d, family_id = "M_ID", subject_id = "IID",
    cluster = "FID", fit_constraints = "estimated_state", family_random = TRUE,
    fixed_effects = NULL, pi = truth$population$pi,
    covariate_mapping_confirmed = TRUE, diagnostic = FALSE), scenario$fit)
  if (!is.null(scenario$external_pi_sd)) {
    set.seed(seed + 1L)
    supplied_pi <- stats::rnorm(1L, truth$population$pi, scenario$external_pi_sd)
    while (supplied_pi <= 0 || supplied_pi >= 1)
      supplied_pi <- stats::rnorm(1L, truth$population$pi, scenario$external_pi_sd)
    fit_arguments$pi <- supplied_pi
  }
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(do.call(ravel_mixed_fit, fit_arguments),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity)
  elapsed <- proc.time()[["elapsed"]] - start
  failed <- inherits(fit, "error")
  converged <- !failed && isTRUE(fit$converged)
  note <- paste(unique(c(if (failed) conditionMessage(fit) else c(fit$error, fit$messages), warnings)), collapse = "; ")
  attempts[[length(attempts) + 1L]] <- data.frame(scenario = scenario_name,
    replicate = replicate, seed = seed, N_families = families, N_rows = nrow(d),
    N_FIDs = length(unique(d$FID)), pi_supplied = fit_arguments$pi,
    error = failed, converged = converged,
    optimizer_converged = if (failed) FALSE else isTRUE(fit$optimizer_converged),
    inference_status = if (failed || is.null(fit$inference_status)) NA_character_ else fit$inference_status,
    projected_gradient_max = if (failed || !length(fit$projected_gradient)) NA_real_ else max(abs(fit$projected_gradient)),
    elapsed_seconds = elapsed, note = note,
    stringsAsFactors = FALSE)
  expected <- expected_coefficients(truth, fit_arguments$interactions)
  tab <- if (converged) coefficient_table(fit) else data.frame()
  expected$estimate <- expected$se <- NA_real_
  if (nrow(tab) && all(c("scale", "parameter", "estimate", "se") %in% names(tab))) {
    key <- paste(expected$scale, expected$parameter)
    idx <- match(key, paste(tab$scale, tab$parameter))
    expected$estimate <- tab$estimate[idx]
    expected$se <- tab$se[idx]
  }
  expected$available <- converged & is.finite(expected$estimate) &
    is.finite(expected$se) & expected$se > 0
  expected$covered <- expected$available &
    abs(expected$estimate - expected$truth) <= stats::qnorm(.975) * expected$se
  expected$reject_zero <- expected$available & abs(expected$estimate / expected$se) > stats::qnorm(.975)
  estimates[[length(estimates) + 1L]] <- cbind(scenario = scenario_name,
    replicate = replicate, expected, stringsAsFactors = FALSE)
  if (converged && is.data.frame(fit$variance_components)) {
    v <- fit$variance_components
    v$scenario <- scenario_name; v$replicate <- replicate
    target_variance <- c(residual = truth$residual_variance,
      family_intercept = truth$family_variance, individual_intercept = truth$individual_variance,
      individual_slope = truth$slope_variance,
      individual_intercept_slope_covariance = truth$intercept_slope_covariance)
    v$truth <- unname(target_variance[v$component])
    v$is_fixed <- truth$outcome == "binary_liability" & v$component == "residual"
    variances[[length(variances) + 1L]] <- v
    estimated_v <- v[!v$is_fixed, , drop = FALSE]
    ev <- data.frame(scale = "variance", parameter = estimated_v$component, truth = estimated_v$truth,
      estimate = estimated_v$estimate, se = estimated_v$se)
    ev$available <- is.finite(ev$estimate) & is.finite(ev$se) & ev$se > 0
    ev$covered <- ev$available & abs(ev$estimate - ev$truth) <= stats::qnorm(.975) * ev$se
    ev$reject_zero <- ev$available & abs(ev$estimate / ev$se) > stats::qnorm(.975)
    estimates[[length(estimates) + 1L]] <- cbind(scenario = scenario_name,
      replicate = replicate, ev, stringsAsFactors = FALSE)
  } else {
    names_variance <- c(if (truth$outcome != "binary_liability") "residual",
      if (isTRUE(fit_arguments$family_random)) "family_intercept",
      if (isTRUE(fit_arguments$individual_random)) "individual_intercept",
      if (!is.null(fit_arguments$random_slope)) c("individual_slope",
        if (isTRUE(fit_arguments$individual_random)) "individual_intercept_slope_covariance"))
    targets <- c(residual = truth$residual_variance, family_intercept = truth$family_variance,
      individual_intercept = truth$individual_variance, individual_slope = truth$slope_variance,
      individual_intercept_slope_covariance = truth$intercept_slope_covariance)
    ev <- data.frame(scale = "variance", parameter = names_variance,
      truth = unname(targets[names_variance]), estimate = NA_real_, se = NA_real_,
      available = FALSE, covered = FALSE, reject_zero = FALSE)
    estimates[[length(estimates) + 1L]] <- cbind(scenario = scenario_name,
      replicate = replicate, ev, stringsAsFactors = FALSE)
  }
  cat(sprintf("%s %d/%d converged=%s %.2fs %s\n", scenario_name,
              replicate, replicates, converged, elapsed, note))
  if (stage == "smoke") saveRDS(fit, file.path(output, paste0(scenario_name, "-fit.rds")))
  write_progress()
}

all_estimates <- do.call(rbind, estimates)
summarize <- function(x) {
  good <- x$available
  n <- nrow(x); ng <- sum(good)
  error <- x$estimate[good] - x$truth[good]
  empirical_sd <- if (ng > 1L) stats::sd(x$estimate[good]) else NA_real_
  coverage <- if (ng) mean(x$covered[good]) else NA_real_
  binomial <- if (ng) stats::binom.test(sum(x$covered[good]), ng)$conf.int else c(NA_real_, NA_real_)
  rejection <- if (ng && all(x$truth == 0)) mean(x$reject_zero[good]) else NA_real_
  data.frame(scenario = x$scenario[1L], scale = x$scale[1L], parameter = x$parameter[1L],
    truth = x$truth[1L], attempted = n, available = ng, unavailable = n - ng,
    bias = if (ng) mean(error) else NA_real_,
    bias_MCSE = if (ng > 1L) stats::sd(error) / sqrt(ng) else NA_real_,
    empirical_sd = empirical_sd, mean_SE = if (ng) mean(x$se[good]) else NA_real_,
    SE_over_empirical_SD = if (ng > 1L) mean(x$se[good]) / empirical_sd else NA_real_,
    coverage_among_available = coverage,
    coverage_MCSE = if (ng) sqrt(coverage * (1 - coverage) / ng) else NA_real_,
    coverage_binomial_low = binomial[1L], coverage_binomial_high = binomial[2L],
    coverage_all_attempts = mean(x$covered),
    null_rejection_among_available = rejection,
    null_rejection_MCSE = if (ng && is.finite(rejection)) sqrt(rejection * (1 - rejection) / ng) else NA_real_)
}
summary <- do.call(rbind, lapply(split(all_estimates,
  interaction(all_estimates$scenario, all_estimates$scale, all_estimates$parameter, drop = TRUE)), summarize))
rownames(summary) <- NULL
utils::write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
writeLines(c(paste("stage:", stage), paste("replicates:", replicates),
  paste("families:", families), paste("scenarios:", paste(selected, collapse = ", ")),
  "", capture.output(sessionInfo()), "", "Source hashes:",
  capture.output(source_hashes)),
  file.path(output, "session.txt"))
print(summary[summary$scale == "corrected_std", ], row.names = FALSE)
cat("All attempted fits and estimation failures retained in", normalizePath(output), "\n")
