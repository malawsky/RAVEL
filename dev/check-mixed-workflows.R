# Reproducible workflow integration checks; these are not calibration simulations.
# From the package root: Rscript dev/check-mixed-workflows.R [--binary | --binary-only]
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-mixed-validation.R")
binary_only <- "--binary-only" %in% commandArgs(trailingOnly = TRUE)
destination <- file.path("dev", "mixed-validation", if (binary_only) "workflow-binary-check" else "workflow-check")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
workflow_path <- file.path("analysis", "parental_deviation", "ravel_parental_deviation.R")
source_hashes <- tools::md5sum(c(workflow_path, "dev/check-mixed-workflows.R",
  "tests/testthat/helper-mixed-validation.R", list.files("R", pattern = "mixed|score_diagnostic", full.names = TRUE)))
previous_option <- getOption("ravel.workflow_functions_only")
options(ravel.workflow_functions_only = TRUE)
workflow <- new.env(parent = globalenv())
sys.source(workflow_path, envir = workflow)
options(ravel.workflow_functions_only = previous_option)

# Load the ACTUAL collector function from the script's guarded RUN block,
# without executing its data-import or output loop and without reproducing it.
assignments <- list()
for (expression in parse(workflow_path)) {
  if (!is.call(expression) || !identical(expression[[1L]], as.name("if")) ||
      !any(grepl("ravel.workflow_functions_only", deparse(expression[[2L]]), fixed = TRUE))) next
  block <- expression[[3L]]
  if (!is.call(block) || !identical(block[[1L]], as.name("{"))) next
  for (statement in as.list(block)[-1L])
    if (is.call(statement) && identical(statement[[1L]], as.name("<-")) &&
        identical(statement[[2L]], as.name("run_outcome")))
      assignments[[length(assignments) + 1L]] <- statement
}
stopifnot(length(assignments) == 1L)
eval(assignments[[1L]], envir = workflow)

checks <- list()
check <- function(condition, name, detail = "") {
  ok <- isTRUE(condition)
  checks[[length(checks) + 1L]] <<- data.frame(check = name,
    status = if (ok) "pass" else "FAIL", detail = detail)
  utils::write.csv(do.call(rbind, checks), file.path(destination, "checks.csv"), row.names = FALSE)
  if (!ok) stop("Workflow check failed: ", name, ": ", detail)
  invisible(TRUE)
}
pgs <- c(child1 = "PGS1_child", child2 = "PGS2_child", mother1 = "PGS1_mother",
  mother2 = "PGS2_mother", father1 = "PGS1_father", father2 = "PGS2_father")
configure <- function(dat, binary = FALSE) {
  settings <- list(dat = dat, binary_outcome = binary, sib_configs = "none",
    models = list("Estimated state" = list(fit_constraints = "estimated_state"),
      "Equilibrium" = list(fit_constraints = "equilibrium")),
    design = "cross_sectional", pgs_columns = pgs, covariates = NULL,
    family_id = "FID", subject_id = "IID", diagnostic_cluster = "auto",
    maternal_id = "M_ID", random_family = FALSE, random_individual = FALSE,
    random_slope = NULL, genetic_interactions = NULL, covariate_mapping_confirmed = FALSE,
    age = "age", age_center = NULL, pi = .4, pi_se = NULL, pi_grid = .45,
    grid_pi_se = NULL, include_raw = TRUE, drop_missing_outcome = TRUE,
    sib_columns = NULL, mixed_settings = list(maternal_id = "M_ID",
      subject_id = "IID", paternal_id = "F_ID"))
  list2env(settings, envir = workflow)
  invisible(settings)
}

w <- c(1, -1)
if (!binary_only) {
# NULL dispatch must execute the legacy function body. Strip only the resolver
# and mixed dispatch before its original design statement; no alternative
# preprocessing is written and changes to prefix length do not change the oracle.
reference <- workflow$ravel_parental_deviation
statements <- as.list(body(reference))[-1L]
legacy_start <- which(vapply(statements, function(statement)
  is.call(statement) && identical(statement[[1L]], as.name("<-")) &&
    identical(statement[[2L]], as.name("design")) &&
    identical(statement[[3L]], quote(match.arg(design))), logical(1)))
stopifnot(length(legacy_start) == 1L)
body(reference) <- as.call(c(list(as.name("{")), statements[seq.int(legacy_start, length(statements))]))
formals(reference)[["mixed_settings"]] <- NULL
workflow$main_calls <- list()
workflow$ravel_fit <- function(...) {
  workflow$main_calls[[length(workflow$main_calls) + 1L]] <- list(...)
  stop("Controlled primary failure for dispatch comparison")
}
workflow$ravel_score_diagnostic <- function(...) list(status = "routing_fixture",
  reason = "Diagnostic fit omitted for dispatch comparison", N = NA_integer_)
d <- mixed_validation_data(N_families = 250L, sib_sizes = 2L, visits = 1L,
  pi = .4, family_variance = .5, age_effect = 0, sex_effect = 0)
legacy_args <- list(dat = d, pgs_columns = pgs, family_id = "FID", pi = .4,
  pi_se = NULL, pi_grid = NULL, fit_constraints = "estimated_state", include_comparison = FALSE)
old <- suppressMessages(do.call(reference, legacy_args))
old_calls <- workflow$main_calls; workflow$main_calls <- list()
omitted <- suppressMessages(do.call(workflow$ravel_parental_deviation, legacy_args))
omitted_calls <- workflow$main_calls; workflow$main_calls <- list()
explicit <- suppressMessages(do.call(workflow$ravel_parental_deviation,
  c(legacy_args, list(mixed_settings = NULL))))
check(identical(old_calls, omitted_calls) && identical(old_calls, workflow$main_calls),
  "default_NULL_preserves_legacy_primary_calls")
check(identical(old, omitted) && identical(old, explicit), "default_NULL_preserves_legacy_return")
configure(d); workflow$mixed_settings <- NULL
workflow$models <- list("Estimated state" = list(fit_constraints = "estimated_state"))
workflow$pi_grid <- NULL; workflow$main_calls <- list()
legacy_collected <- suppressMessages(workflow$run_outcome("Y"))
check(identical(old_calls, workflow$main_calls) &&
  !length(attr(legacy_collected, "mixed_fit_records", exact = TRUE)) &&
  !length(attr(legacy_collected, "mixed_audits", exact = TRUE)), "default_NULL_collector_stays_on_legacy_estimator")
rm("ravel_fit", "ravel_score_diagnostic", envir = workflow)

# Real collector: two model labels, each with a main and a full pi-grid refit.
configure(d)
continuous <- suppressMessages(workflow$run_outcome("Y"))
records <- attr(continuous, "mixed_fit_records", exact = TRUE)
audits <- attr(continuous, "mixed_audits", exact = TRUE)
expected_keys <- c("none|Estimated state|main", "none|Estimated state|grid_1",
  "none|Equilibrium|main", "none|Equilibrium|grid_1")
check(identical(names(records), expected_keys) && identical(names(audits), expected_keys),
  "collector_retains_every_model_main_and_grid")
check(all(vapply(records, function(record) inherits(record$fit, "ravel_mixed_results"), logical(1))),
  "collector_preserves_fitted_classes")
check(all(vapply(records, function(record) isTRUE(record$fit$converged), logical(1))),
  "continuous_collected_fits_converged")
check(all(vapply(records, function(record) if (record$kind == "main")
  !is.null(record$fit$one_factor) else is.null(record$fit$one_factor), logical(1))),
  "one_diagnostic_per_main_no_grid_diagnostics")
fit <- records[[1L]]$fit
corrected <- continuous[continuous$model == "Estimated state" &
  continuous$analysis == "Corrected for pi" & continuous$parameter == "delta_M_mid", ]
w <- c(1, -1); V <- fit$coefficient_vcov$corrected_std[c("beta_M", "beta_F"), c("beta_M", "beta_F")]
check(nrow(corrected) == 1L && abs(corrected$se - sqrt(drop(t(w) %*% V %*% w))) < 1e-10,
  "parent_difference_uses_joint_covariance")
raw <- continuous[continuous$model == "Estimated state" &
  continuous$analysis == "Corrected for pi (founder scale)" & continuous$parameter == "cov_corrected", ]
check(nrow(raw) == 1L && abs(raw$estimate - fit$state$rho * fit$state$parent) < 1e-12,
  "founder_spouse_row_is_covariance_not_correlation")
check(any(grepl("condition on supplied pi", raw$note)), "spouse_uncertainty_scope_is_explicit")
check(!any(grepl("Plain PGS regressions", continuous$model)), "no_unimplemented_plain_PGS_mixed_comparator")
check(!any(vapply(records, function(record) is.null(record$fit$prepared), logical(1))),
  "main_and_grid_keep_prepared_designs")

saveRDS(continuous, file.path(destination, "continuous-result.rds"))
restored <- readRDS(file.path(destination, "continuous-result.rds"))
check(identical(names(attr(restored, "mixed_fit_records", exact = TRUE)), expected_keys),
  "RDS_retains_every_fit_record")
unlink(file.path(destination, "continuous-result.rds"))
paths <- workflow$write_ravel_mixed_workflow_audits(audits, file.path(destination, "continuous"))
check(length(paths) >= 19L && all(file.exists(paths)) &&
  all(c("implied_covariances", "implied_covariance_vcov", "state_and_genetic_vcov",
    "parameter_vcov", "coefficient_vcov") %in% names(paths)), "all_detailed_audit_CSVs_written",
  paste(length(paths), "tables; empty covariate-scaling table is omitted when no covariates are fitted"))
check(all(vapply(paths, function(path) !anyDuplicated(names(utils::read.csv(path, check.names = FALSE))), logical(1))),
  "audit_CSV_columns_remain_unique")
state_rows <- utils::read.csv(paths[["state"]], check.names = FALSE)
check(setequal(unique(state_rows$source_fit), expected_keys), "CSV_source_metadata_identifies_all_fits")

# Test missingness routing with the actual public selector, but no redundant
# estimator or auxiliary SEM refits. Selection completes before the controlled
# auxiliary-fit failure, so its selected/used row indices remain inspectable.
workflow$models <- list("Estimated state" = list(fit_constraints = "estimated_state"))
workflow$pi_grid <- NULL
auxiliary_stop <- RAVEL::ravel_fit
body(auxiliary_stop) <- quote(stop("Controlled auxiliary-fit stop after selection"))
selection_run <- function() testthat::with_mocked_bindings(
  suppressMessages(workflow$run_outcome("Y")),
  ravel_mixed_fit = function(...) fit,
  ravel_fit = auxiliary_stop, .package = "RAVEL")
before <- selection_run()
original_selection <- attr(before, "diagnostics", exact = TRUE)$none[[1L]]$selection
chosen <- original_selection$rows_selected
workflow$dat$Y[chosen[1:8]] <- NA_real_
workflow$dat[workflow$dat$IID == workflow$dat$IID[chosen[9L]], unname(pgs)] <- NA_real_
after <- selection_run()
selection <- attr(after, "diagnostics", exact = TRUE)$none[[1L]]$selection
check(identical(selection$rows_selected, chosen) && selection$N_incomplete_selected == 9L,
  "missing_selected_Y_and_PGS_do_not_trigger_replacement")
check(selection$N_input == nrow(d), "diagnostic_sees_original_pre_exclusion_rows")

# Actual collector also retains failed main/grid records and their audits.
configure(d); workflow$models <- list("Estimated state" = list(fit_constraints = "estimated_state"))
failed <- testthat::with_mocked_bindings(suppressMessages(workflow$run_outcome("Y")),
  ravel_mixed_fit = function(...) stop("Controlled primary fit failure"),
  ravel_score_diagnostic = function(...) list(status = "routing_fixture", reason = "Captured failed-primary diagnostic", N = nrow(d)),
  .package = "RAVEL")
failed_records <- attr(failed, "mixed_fit_records", exact = TRUE)
failed_audits <- attr(failed, "mixed_audits", exact = TRUE)
check(length(failed_records) == 2L && all(vapply(failed_records, function(x) is.null(x$fit), logical(1))),
  "failed_main_and_grid_records_retained")
check(length(failed_audits) == 2L && all(vapply(failed_audits,
  function(x) identical(x$status$inference_status, "fit_error"), logical(1))), "failed_audit_statuses_retained")
check(is.null(attr(failed, "fits", exact = TRUE)$none[[1L]]), "failed_main_fit_uses_exact_attribute_lookup")
failed_paths <- workflow$write_ravel_mixed_workflow_audits(failed_audits, file.path(destination, "failed"))
check(all(file.exists(failed_paths)), "failed_fit_audit_CSVs_written")

for (setting in list(list(parental_sibs = "maternal"), list(sib_columns = pgs), list(fit_constraints = "none"))) {
  error <- tryCatch(do.call(workflow$ravel_parental_deviation,
    c(list(dat = d, pgs_columns = pgs, mixed_settings = workflow$mixed_settings), setting)), error = identity)
  check(inherits(error, "error"), paste0("unsupported_request_rejected_", names(setting)))
}

}

if (any(c("--binary", "--binary-only") %in% commandArgs(trailingOnly = TRUE))) {
  binary_data <- mixed_validation_data(N_families = 300L, sib_sizes = 2L, visits = 1L,
    pi = .4, family_variance = .5, age_effect = 0, sex_effect = .12,
    sex_beta = c(child = .12, mother = -.06, father = .02), outcome = "binary_liability",
    residual_variance = 1, prevalence_reference = .3)
  configure(binary_data, TRUE)
  workflow$models <- list("Estimated state" = list(fit_constraints = "estimated_state"))
  workflow$pi_grid <- NULL
  workflow$mixed_settings$interactions <- "sex"
  workflow$mixed_settings$covariate_mapping_confirmed <- TRUE
  binary <- suppressMessages(workflow$run_outcome("Y"))
  binary_fit <- attr(binary, "fits", exact = TRUE)$none[[1L]]
  check(isTRUE(binary_fit$converged) && identical(binary_fit$inference_status, "ok"),
    "binary_sex_interaction_workflow_inference", "300 maternal families; 600 outcomes; smoke only")
  check(all(c("beta_D:sex", "beta_M:sex", "beta_F:sex", "beta_mid:sex", "delta_M_mid:sex") %in% binary$parameter),
    "binary_all_genetic_sex_slopes_and_parent_contrasts")
  r <- binary[binary$analysis == "Corrected for pi" & binary$parameter == "delta_M_mid:sex", ]
  V <- binary_fit$coefficient_vcov$corrected_std[c("beta_M:sex", "beta_F:sex"), c("beta_M:sex", "beta_F:sex")]
  check(abs(r$se - sqrt(drop(t(w) %*% V %*% w))) < 1e-10, "binary_interaction_contrast_uses_joint_covariance")
  printed <- capture.output(workflow$print_parental_deviation(binary))
  check(any(grepl("conditional probit scale", printed, fixed = TRUE)), "binary_console_scale_is_explicit")
  check(!any(grepl("plain regressions keep ordinary SEs", printed, ignore.case = TRUE)),
    "binary_console_does_not_invent_comparator")
  binary_diagnostic <- attr(binary, "diagnostics", exact = TRUE)$none[[1L]]
  diagnostic_rows <- binary[binary$parameter == "one_factor_test", , drop = FALSE]
  check(identical(binary_diagnostic$status, "ok") &&
    identical(binary_diagnostic$calibration_status, "uncalibrated") &&
    is.finite(binary_diagnostic$p_value) &&
    all(diagnostic_rows$p_value == binary_diagnostic$p_value),
    "binary_nominal_diagnostic_p_returned_with_uncalibrated_status")
  check(is.finite(binary_diagnostic$statistic) &&
    identical(binary_diagnostic$diagnostics$exploratory$calibrated, FALSE) &&
    identical(binary_diagnostic$p_value, binary_diagnostic$diagnostics$exploratory$nominal_p_value),
    "binary_diagnostic_preserves_explicitly_exploratory_details")
  check(any(grepl("8.0%", binary_diagnostic$auxiliary$warnings, fixed = TRUE)) &&
    any(grepl("nominal 5% threshold", binary_diagnostic$workflow_warnings, fixed = TRUE)) &&
    all(grepl("8.0%", diagnostic_rows$note, fixed = TRUE)),
    "binary_calibration_warning_retained_in_auxiliary_workflow_and_table")
  csv_file <- file.path(destination, "binary-sex-result.csv")
  utils::write.csv(binary, csv_file, row.names = FALSE)
  csv_result <- utils::read.csv(csv_file, stringsAsFactors = FALSE)
  csv_diagnostic <- csv_result[csv_result$parameter == "one_factor_test", , drop = FALSE]
  check(nrow(csv_diagnostic) == 1L && is.finite(csv_diagnostic$p_value) &&
    abs(csv_diagnostic$p_value - binary_diagnostic$p_value) < 1e-14 &&
    grepl("8.0%", csv_diagnostic$note, fixed = TRUE) &&
    grepl("nominal 5% threshold", csv_diagnostic$note, fixed = TRUE),
    "CSV_reports_binary_nominal_p_and_actual_calibration_warning")
  saveRDS(binary, file.path(destination, "binary-sex-result.rds"))
  restored_binary <- readRDS(file.path(destination, "binary-sex-result.rds"))
  restored_diagnostic <- attr(restored_binary, "diagnostics", exact = TRUE)$none[[1L]]
  check(identical(restored_diagnostic$status, "ok") &&
    identical(restored_diagnostic$calibration_status, "uncalibrated") &&
    identical(restored_diagnostic$p_value, binary_diagnostic$p_value) &&
    any(grepl("8.0%", restored_diagnostic$auxiliary$warnings, fixed = TRUE)),
    "RDS_preserves_binary_nominal_p_and_uncalibrated_warning")
  unlink(file.path(destination, "binary-sex-result.rds"))
  binary_paths <- workflow$write_ravel_mixed_workflow_audits(
    attr(binary, "mixed_audits", exact = TRUE), file.path(destination, "binary-sex"))
  check(length(binary_paths) >= 20L, "binary_detailed_audit_CSVs_written")
}
writeLines(c(capture.output(sessionInfo()), "", "Source hashes:", capture.output(source_hashes)),
  file.path(destination, "session.txt"))
cat("Workflow integration checks passed:", length(checks), "\n")
