# Recover longitudinal variance point estimates from saved parameters. This
# regenerates each original data set only to recover its outcome/age scaling;
# it does not refit, alter failure status, or claim variance-interval coverage.
source("tests/testthat/helper-mixed-validation.R")
source("tests/testthat/helper-mixed-covariates.R")
directory <- file.path("dev", "mixed-covariate-validation", "longitudinal")
attempts <- utils::read.csv(file.path(directory, "attempts.csv"), stringsAsFactors = FALSE)
parameters <- utils::read.csv(file.path(directory, "parameters.csv"), stringsAsFactors = FALSE)
contract <- mixed_covariate_contract(sex_pc_association = .5,
  baseline_pc_association = 0, baseline_between = 0, baseline_within = 0)
results <- lapply(seq_len(nrow(attempts)), function(i) {
  a <- attempts[i, ]
  d <- mixed_covariate_data(N_families = a$N_families, seed = a$seed,
    contract = contract, visits = 4L, individual_variance = .2,
    slope_variance = .003, intercept_slope_covariance = .012,
    age_beta = c(child = .025, mother = -.008, father = 0),
    sex_beta = c(child = .12, mother = -.06, father = .02))
  z <- parameters[parameters$scenario == a$scenario & parameters$replicate == a$replicate, ]
  value <- stats::setNames(z$estimate, z$parameter)
  sy <- stats::sd(d$Y); sx <- stats::sd(d$age)
  estimate <- c(residual = exp(2 * value["residual_sd"]),
    family_intercept = exp(2 * value["family_sd"]),
    individual_intercept = exp(2 * value["individual_sd"]),
    individual_slope = exp(2 * value["slope_sd"]) / sx^2,
    individual_intercept_slope_covariance =
      exp(value["individual_sd"] + value["slope_sd"]) *
      tanh(value["intercept_slope_cor"]) / sx) * sy^2
  truth <- attr(d, "truth")
  expected <- c(residual = truth$residual_variance, family_intercept = truth$family_variance,
    individual_intercept = truth$individual_variance, individual_slope = truth$slope_variance,
    individual_intercept_slope_covariance = truth$intercept_slope_covariance)
  data.frame(scenario = a$scenario, replicate = a$replicate, seed = a$seed,
    inference_status = a$inference_status, converged = a$converged,
    component = names(expected), truth = unname(expected), estimate = unname(estimate))
})
results <- do.call(rbind, results)
summary <- do.call(rbind, lapply(split(results, results$component), function(x) {
  good <- x$converged & is.finite(x$estimate)
  error <- x$estimate[good] - x$truth[good]
  data.frame(component = x$component[1L], truth = x$truth[1L], attempted = nrow(x),
    finite_point = sum(good), bias = mean(error), bias_mcse = stats::sd(error) / sqrt(sum(good)),
    rmse = sqrt(mean(error^2)), empirical_sd = stats::sd(x$estimate[good]))
}))
utils::write.csv(results, file.path(directory, "variance-point-estimates.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(directory, "variance-point-summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
