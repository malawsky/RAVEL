# Independent step-size check of external-pi propagation for the joint child-
# covariate model. Run from the repository root; one R process, modes sequential.
output <- "dev/mixed-covariate-validation/pi-sensitivity"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(output, "attempts.csv")))
  stop("Existing attempted fits are retained; use a fresh output directory for a rerun.")
files <- c(list.files("R", pattern = "\\.R$", full.names = TRUE),
  "tests/testthat/helper-mixed-validation.R", "tests/testthat/helper-mixed-covariates.R",
  file.path(output, "check.R"))
utils::write.csv(data.frame(file = files, md5 = unname(tools::md5sum(files))),
  file.path(output, "source-md5.csv"), row.names = FALSE)
for (file in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(file)
source("tests/testthat/helper-mixed-validation.R")
source("tests/testthat/helper-mixed-covariates.R")
modes <- c("estimated_state", "equilibrium", "n_gen")
pi_value <- .4; pi_se <- .03; wide_step <- .003
attempts <- comparisons <- uncertainty_checks <- list()
write_outputs <- function() {
  if (length(attempts)) utils::write.csv(do.call(rbind, attempts),
    file.path(output, "attempts.csv"), row.names = FALSE)
  if (length(comparisons)) utils::write.csv(do.call(rbind, comparisons),
    file.path(output, "derivatives.csv"), row.names = FALSE)
  if (length(uncertainty_checks)) utils::write.csv(do.call(rbind, uncertainty_checks),
    file.path(output, "uncertainty-addition.csv"), row.names = FALSE)
}
capture_fit <- function(arguments, mode, stage, seed) {
  warnings <- character(); started <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(do.call(ravel_mixed_fit, arguments),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = identity)
  error <- inherits(f, "error")
  valid <- !error && isTRUE(f$converged) && !length(f$boundary_parameters)
  attempts[[length(attempts) + 1L]] <<- data.frame(mode = mode, stage = stage, seed = seed,
    N_families = 400L, pi = arguments$pi, pi_se = arguments$pi_se %||% 0,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    converged = !error && isTRUE(f$converged), interior_valid = valid,
    inference_status = if (error) "error" else f$inference_status,
    max_gradient_per_family = if (error) NA_real_ else max(abs(f$gradient)) / f$n_families,
    max_natural_gradient_per_family = if (error) NA_real_ else
      max(abs(f$projected_natural_variance_gradient)) / f$n_families,
    boundaries = if (error) "" else paste(f$boundary_parameters, collapse = "; "),
    note = paste(unique(c(if (error) conditionMessage(f) else f$messages, warnings)), collapse = "; "))
  write_outputs()
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), mode, stage,
    if (error) conditionMessage(f) else f$inference_status, "\n"); flush.console()
  list(fit = if (error) NULL else f, valid = valid)
}
for (i in seq_along(modes)) {
  mode <- modes[i]; seed <- 741400L + i
  d <- mixed_covariate_data(N_families = 400L, sib_sizes = 2L, visits = 1L,
    seed = seed, mode = mode, n_generations = 3L,
    contract = mixed_covariate_contract(n_pc = 2L, sex_pc_association = .5))
  arguments <- list(dat = d, pi = pi_value, pi_se = pi_se,
    fit_constraints = mode, am_n_generations = 3L,
    fixed_effects = c("PC1", "PC2", "sex", "age"),
    covariate_model = "joint", diagnostic = FALSE, se = TRUE, keep_data = FALSE)
  main <- capture_fit(arguments, mode, "main_with_pi_se", seed)
  arguments$pi_se <- NULL; arguments$se <- FALSE
  if (!is.null(main$fit)) arguments$start <- main$fit$parameters
  arguments$pi <- pi_value - wide_step
  minus <- capture_fit(arguments, mode, "wider_minus", seed)
  arguments$pi <- pi_value + wide_step
  plus <- capture_fit(arguments, mode, "wider_plus", seed)
  available <- main$valid && identical(main$fit$inference_status, "ok") &&
    !is.null(main$fit$pi_sensitivity) && minus$valid && plus$valid
  if (!available) {
    comparisons[[length(comparisons) + 1L]] <- data.frame(mode = mode, kind = "unavailable",
      scale = NA_character_, parameter = NA_character_, internal_step = .001,
      wider_step = wide_step, stored_derivative = NA_real_, wider_derivative = NA_real_,
      absolute_difference = NA_real_, conditional_se = NA_real_,
      pi_se_weighted_difference = NA_real_, fraction_conditional_se = NA_real_,
      material_agreement = FALSE)
    write_outputs(); next
  }
  f <- main$fit
  append_check <- function(kind, scale, stored, wider, conditional_se) {
    difference <- abs(stored - wider)
    comparisons[[length(comparisons) + 1L]] <<- data.frame(mode = mode,
      kind = kind, scale = scale, parameter = names(stored), internal_step = .001,
      wider_step = wide_step, stored_derivative = unname(stored),
      wider_derivative = unname(wider[names(stored)]), absolute_difference = unname(difference),
      conditional_se = unname(conditional_se), pi_se_weighted_difference = pi_se * unname(difference),
      fraction_conditional_se = pi_se * unname(difference) / pmax(unname(conditional_se), 1e-12),
      material_agreement = pi_se * unname(difference) <= .05 * unname(conditional_se) + 1e-6)
  }
  for (scale in names(f$pi_sensitivity)) {
    estimates <- function(x) {
      table <- x$table[x$table$scale == scale, ]
      stats::setNames(table$estimate, table$parameter)
    }
    wider <- (estimates(plus$fit) - estimates(minus$fit)) / (2 * wide_step)
    map <- .ravel_mixed_delta(function(x) .ravel_mixed_effects(x, f$settings, scale), f$parameters)
    conditional <- map$jacobian %*% f$vcov[names(f$parameters), names(f$parameters), drop = FALSE] %*%
      t(map$jacobian)
    stored <- f$pi_sensitivity[[scale]]
    expected <- conditional + pi_se^2 * tcrossprod(stored)
    uncertainty_checks[[length(uncertainty_checks) + 1L]] <- data.frame(mode = mode,
      scale = scale, maximum_addition_error = max(abs(expected - f$coefficient_vcov[[scale]])))
    append_check("coefficient", scale, stored, wider, sqrt(pmax(diag(conditional), 0)))
  }
  variance_values <- function(x) stats::setNames(x$variance_components$estimate,
    x$variance_components$component)
  wider <- (variance_values(plus$fit) - variance_values(minus$fit)) / (2 * wide_step)
  vm <- .ravel_mixed_delta(function(x) .ravel_mixed_variances(x, f$settings), f$parameters)
  conditional <- vm$jacobian %*% f$vcov[names(f$parameters), names(f$parameters), drop = FALSE] %*%
    t(vm$jacobian)
  append_check("residual_variance", "input_outcome_variance", f$variance_pi_sensitivity,
    wider, sqrt(pmax(diag(conditional), 0)))
  write_outputs()
}
results <- do.call(rbind, comparisons)
summary <- do.call(rbind, lapply(modes, function(mode) {
  d <- results[results$mode == mode, ]
  data.frame(mode = mode, comparisons = nrow(d), all_available = all(d$kind != "unavailable"),
    all_material_agreement = all(d$material_agreement),
    maximum_absolute_derivative_difference = max(d$absolute_difference),
    maximum_pi_se_weighted_difference = max(d$pi_se_weighted_difference),
    maximum_fraction_conditional_se = max(d$fraction_conditional_se))
}))
utils::write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())), file.path(output, "session.txt"))
stopifnot(all(summary$all_available), all(summary$all_material_agreement),
  max(do.call(rbind, uncertainty_checks)$maximum_addition_error) < 1e-10)
cat("All three external-pi sensitivity checks passed.\n")
