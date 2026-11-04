# Actual child-workflow collector checks for correlated child covariates.
# These exercise integration and reporting, not Monte Carlo calibration.
# Run from the repository root: Rscript dev/check-mixed-covariate-workflow.R
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-mixed-validation.R")
source("tests/testthat/helper-mixed-covariates.R")
destination <- file.path("dev", "mixed-covariate-validation", "workflow")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
workflow_path <- file.path("analysis", "parental_deviation", "ravel_parental_deviation.R")
source_files <- c(workflow_path, "dev/check-mixed-covariate-workflow.R",
  "tests/testthat/helper-mixed-validation.R", "tests/testthat/helper-mixed-covariates.R",
  list.files("R", pattern = "mixed|score_diagnostic", full.names = TRUE))
source_hashes <- tools::md5sum(source_files)
utils::write.csv(data.frame(file = names(source_hashes), md5 = unname(source_hashes)),
  file.path(destination, "source-md5.csv"), row.names = FALSE)
previous <- getOption("ravel.workflow_functions_only")
options(ravel.workflow_functions_only = TRUE)
workflow <- new.env(parent = globalenv())
sys.source(workflow_path, envir = workflow)
options(ravel.workflow_functions_only = previous)

# Extract the actual run_outcome definition from the guarded RUN block. Do not
# duplicate the collector, run imports, or evaluate the script's example data.
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
checks <- list(); fit_status <- list()
check <- function(value, label, detail = "") {
  ok <- isTRUE(value)
  checks[[length(checks) + 1L]] <<- data.frame(check = label,
    status = if (ok) "pass" else "FAIL", detail = detail)
  utils::write.csv(do.call(rbind, checks), file.path(destination, "checks.csv"), row.names = FALSE)
  if (!ok) stop("Correlated-covariate workflow check failed: ", label, ": ", detail)
}
pgs <- c(child1 = "PGS1_child", child2 = "PGS2_child", mother1 = "PGS1_mother",
  mother2 = "PGS2_mother", father1 = "PGS1_father", father2 = "PGS2_father")
rename_inputs <- function(d) {
  names(d)[match(c("M_ID", "IID", "FID", "age"), names(d))] <-
    c("Nuclear", "Person", "Extended", "VisitAge")
  d
}
configure <- function(d, longitudinal = FALSE) {
  invisible(list2env(list(dat = d, binary_outcome = FALSE, sib_configs = "none",
    models = list("Estimated state" = list(fit_constraints = "estimated_state")),
    design = if (longitudinal) "longitudinal" else "cross_sectional",
    pgs_columns = pgs, covariates = c("PC1", "sex", if (!longitudinal) "VisitAge"),
    family_id = "Extended", subject_id = "Person", diagnostic_cluster = "auto",
    maternal_id = "Nuclear", random_family = TRUE, random_individual = longitudinal,
    random_slope = NULL, genetic_interactions = if (longitudinal) c("VisitAge", "sex") else NULL,
    covariate_mapping_confirmed = FALSE, mixed_settings = NULL,
    age = "VisitAge", age_center = 10, pi = .4, pi_se = NULL,
    pi_grid = if (longitudinal) NULL else .45, grid_pi_se = NULL,
    include_raw = TRUE, drop_missing_outcome = TRUE, sib_columns = NULL), envir = workflow))
}
inspect <- function(result, scenario, expected_keys, elapsed, require_diagnostic_p = TRUE) {
  records <- attr(result, "mixed_fit_records", exact = TRUE)
  audits <- attr(result, "mixed_audits", exact = TRUE)
  check(identical(names(records), expected_keys) && identical(names(audits), expected_keys),
    paste0(scenario, "_actual_collector_retains_main_and_grid"))
  for (key in names(records)) {
    fit <- records[[key]]$fit
    fit_status[[length(fit_status) + 1L]] <<- data.frame(scenario = scenario, source_fit = key,
      pi = records[[key]]$pi, converged = isTRUE(fit$converged),
      inference_status = if (is.null(fit)) "fit_error" else fit$inference_status,
      N = if (is.null(fit)) NA_integer_ else fit$N,
      families = if (is.null(fit)) NA_integer_ else fit$n_families,
      clusters = if (is.null(fit)) NA_integer_ else fit$n_clusters,
      scenario_elapsed_seconds = elapsed,
      note = paste(records[[key]]$note, collapse = "; "))
  }
  utils::write.csv(do.call(rbind, fit_status), file.path(destination, "fit-status.csv"), row.names = FALSE)
  check(all(vapply(records, function(r) inherits(r$fit, "ravel_mixed_results") &&
    isTRUE(r$fit$converged) && identical(r$fit$inference_status, "ok"), logical(1))),
    paste0(scenario, "_fits_have_regular_inference"))
  check(all(vapply(records, function(r) identical(r$fit$settings$covariate_model, "joint") &&
    all(c("PC1", "sex") %in% r$fit$settings$covariate$columns), logical(1))),
    paste0(scenario, "_default_joint_accepts_PC_and_sex_without_confirmation"))
  check(all(vapply(records, function(r) identical(r$fit$family_id, "Nuclear") &&
    identical(r$fit$subject_id, "Person") && identical(r$fit$cluster, "Extended"), logical(1))),
    paste0(scenario, "_custom_nuclear_person_and_FID_ids_retained"))
  check(all(vapply(records, function(r) if (r$kind == "main")
    !is.null(r$fit$one_factor) else is.null(r$fit$one_factor), logical(1))),
    paste0(scenario, "_diagnostic_only_for_main_fit"))
  diagnostic <- records[[1L]]$fit$one_factor
  diagnostic_row <- result[result$parameter == "one_factor_test", , drop = FALSE]
  warning_text <- c(diagnostic$workflow_warnings, diagnostic$auxiliary$warnings,
    diagnostic_row$note)
  utils::write.csv(data.frame(status = diagnostic$status,
    p_value = if (is.null(diagnostic$p_value)) NA_real_ else diagnostic$p_value,
    reason = paste(diagnostic$reason, collapse = "; "),
    warnings = paste(warning_text, collapse = "; ")),
    file.path(destination, paste0(scenario, "-diagnostic.csv")), row.names = FALSE)
  successful_diagnostic <- identical(diagnostic$status, "ok") && is.finite(diagnostic$p_value) &&
    nrow(diagnostic_row) == 1L && diagnostic_row$p_value == diagnostic$p_value
  check(successful_diagnostic || (!require_diagnostic_p &&
    !identical(diagnostic$status, "ok") && nrow(diagnostic_row) == 1L &&
    is.na(diagnostic_row$p_value) && nzchar(diagnostic$reason)),
    paste0(scenario, if (require_diagnostic_p) "_one_factor_p_retained_with_covariates" else
      "_one_factor_result_or_explicit_failure_retained"),
    paste(diagnostic$status, diagnostic$reason, collapse = "; "))
  check(any(grepl("Covariates were not used", warning_text, fixed = TRUE)) &&
    grepl("Covariates were not used", diagnostic_row$note, fixed = TRUE),
    paste0(scenario, "_unadjusted_diagnostic_warning_in_saved_table"))
  check(identical(diagnostic$workflow$original_cluster, "Extended"),
    paste0(scenario, "_diagnostic_uses_custom_outer_cluster"))
  check(all(vapply(audits, function(a) identical(a$covariate_audit_status$status, "ok") &&
    all(is.finite(a$covariate_moments$se)) &&
    all(is.finite(a$genetic_covariate_covariances$se)) &&
    all(is.finite(a$conditional_implied_covariances$se)), logical(1))),
    paste0(scenario, "_new_audits_have_finite_joint_uncertainty"))
  check(all(vapply(records, function(r) is.null(r$fit$settings$covariate$moment_scores) &&
    is.null(r$fit$prepared$spec$covariate$moment_scores), logical(1))),
    paste0(scenario, "_audit_works_after_large_moment_scores_removed"))
  tmp <- tempfile(paste0("ravel-", scenario, "-")); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  saved <- file.path(tmp, "fit.rds"); saveRDS(result, saved)
  restored <- readRDS(saved)
  check(identical(names(attr(restored, "mixed_fit_records", exact = TRUE)), expected_keys) &&
    identical(attr(restored, "mixed_audits", exact = TRUE)[[1L]]$covariate_audit_vcov,
      audits[[1L]]$covariate_audit_vcov), paste0(scenario, "_RDS_retains_nested_audit_matrices"))
  paths <- workflow$write_ravel_mixed_workflow_audits(audits, file.path(tmp, scenario))
  required <- c("covariate_moments", "genetic_covariate_covariances",
    "conditional_genetic_covariances", "conditional_implied_covariances",
    "marginal_implied_covariances", "covariate_audit_vcov", "covariate_features",
    "covariate_audit_status")
  check(all(required %in% names(paths)) && all(file.exists(paths)),
    paste0(scenario, "_all_new_audit_CSVs_written"))
  check(all(vapply(paths, function(path) !anyDuplicated(names(utils::read.csv(path,
    check.names = FALSE))), logical(1))), paste0(scenario, "_audit_CSV_columns_unique"))
  for (kind in required[1:4]) {
    frame <- utils::read.csv(paths[[kind]], check.names = FALSE)
    check(setequal(unique(frame$source_fit), expected_keys) && all(is.finite(frame$se)),
      paste0(scenario, "_", kind, "_CSV_sources_and_SEs"))
  }
  matrix_rows <- utils::read.csv(paths[["covariate_audit_vcov"]], check.names = FALSE)
  check(setequal(unique(matrix_rows$covariance_scale), c("moments", "cross", "genetic", "implied")) &&
    setequal(unique(matrix_rows$source_fit), expected_keys) && all(is.finite(matrix_rows$covariance)),
    paste0(scenario, "_list_of_covariance_matrices_exported"))
  kept <- paths[c(required, "fit_statistics", "status")]
  for (path in kept) file.copy(path, file.path(destination, basename(path)), overwrite = TRUE)
  utils::write.csv(result, file.path(destination, paste0(scenario, "-effects.csv")), row.names = FALSE)
  records
}

cross <- rename_inputs(mixed_covariate_data(N_families = 1200L, sib_sizes = 2L,
  contract = mixed_covariate_contract(sex_pc_association = .6), fixed_pc = c(.12, 0),
  family_variance = .45, families_per_FID = 2L, extended_fraction = .4, seed = 724201L))
configure(cross)
started <- proc.time()[["elapsed"]]
cross_result <- suppressMessages(workflow$run_outcome("Y"))
cross_records <- inspect(cross_result, "cross", c("none|Estimated state|main",
  "none|Estimated state|grid_1"), proc.time()[["elapsed"]] - started)
check(identical(vapply(cross_records, function(r) r$fit$settings$pi, numeric(1)),
  stats::setNames(c(.4, .45), names(cross_records))), "cross_pi_grid_full_refit_settings")
check(cross_records[[1L]]$fit$n_families == 2L * cross_records[[1L]]$fit$n_clusters,
  "cross_nuclear_groups_differ_from_outer_FID_groups")
cat("Cross-sectional joint-covariate collector passed.\n")

long <- rename_inputs(mixed_covariate_data(N_families = 200L, sib_sizes = 2L, visits = 3L,
  contract = mixed_covariate_contract(sex_pc_association = .6), fixed_pc = c(.12, 0),
  age_beta = c(child = .025, mother = -.008, father = 0),
  sex_beta = c(child = .12, mother = -.06, father = .02),
  family_variance = .4, individual_variance = .25, seed = 724202L))
configure(long, longitudinal = TRUE)
started <- proc.time()[["elapsed"]]
long_result <- suppressMessages(workflow$run_outcome("Y"))
long_records <- inspect(long_result, "longitudinal", "none|Estimated state|main",
  proc.time()[["elapsed"]] - started, require_diagnostic_p = FALSE)
fit <- long_records[[1L]]$fit
check(identical(fit$settings$person_mean_covariates, "VisitAge") &&
  fit$covariate_features$aggregation[fit$covariate_features$column == "VisitAge"] ==
    "person_mean_of_supplied_visits", "longitudinal_custom_age_person_mean_mapped_automatically")
check(all(unlist(lapply(c("VisitAge", "sex"), function(x)
  paste0(c("beta_D", "beta_M", "beta_F"), ":", x))) %in% fit$table$parameter),
  "longitudinal_age_and_sex_genetic_slopes_retained")
check(all(c("family_intercept", "individual_intercept") %in% fit$variance_components$component) &&
  is.null(fit$settings$random_slope), "longitudinal_both_intercepts_without_unrequested_random_slope")
check(fit$N == nrow(long) && fit$n_subjects == length(unique(long$Person)),
  "longitudinal_repeated_outcomes_and_unique_people_counted")
writeLines(sub("[[:space:]]+$", "", c(capture.output(sessionInfo()), "", "Source hashes:",
  capture.output(source_hashes))), file.path(destination, "session.txt"))
cat("Correlated-covariate child-workflow checks passed:", length(checks), "\n")
