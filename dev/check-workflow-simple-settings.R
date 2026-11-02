# Routing and one real collector check for the simplified child-workflow controls.
# Run from the package root: Rscript dev/check-workflow-simple-settings.R
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-mixed-validation.R")
path <- file.path("analysis", "parental_deviation", "ravel_parental_deviation.R")
destination <- file.path("dev", "mixed-validation", "simple-settings-check")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
hashes <- tools::md5sum(c(path, "dev/check-workflow-simple-settings.R",
  "tests/testthat/helper-mixed-validation.R"))
previous <- getOption("ravel.workflow_functions_only")
options(ravel.workflow_functions_only = TRUE)
workflow <- new.env(parent = globalenv()); sys.source(path, envir = workflow)
options(ravel.workflow_functions_only = previous)
checks <- list()
check <- function(value, name) {
  ok <- isTRUE(value)
  checks[[length(checks) + 1L]] <<- data.frame(check = name, status = if (ok) "pass" else "FAIL")
  utils::write.csv(do.call(rbind, checks), file.path(destination, "checks.csv"), row.names = FALSE)
  if (!ok) stop("Simple-setting check failed: ", name)
}
pgs <- c(child1 = "PGS1_child", child2 = "PGS2_child", mother1 = "PGS1_mother",
  mother2 = "PGS2_mother", father1 = "PGS1_father", father2 = "PGS2_father")
dat <- mixed_validation_data(N_families = 250L, sib_sizes = 2L, visits = 1L,
  family_variance = .5, age_effect = 0, sex_effect = 0)
dat$Nuclear <- dat$M_ID; dat$Person <- dat$IID; dat$Extended <- dat$FID
dat$BroadDiagnostic <- paste0("B", (match(dat$Extended, unique(dat$Extended)) - 1L) %/% 2L)
mixed_function <- workflow$ravel_parental_deviation_mixed
workflow$ravel_parental_deviation_mixed <- function(...) list(...)
call <- function(...) workflow$ravel_parental_deviation(dat, pgs_columns = pgs, ...)

family <- call(random_family = TRUE)
check(family$mixed_settings$family_random && !family$mixed_settings$individual_random &&
  family$mixed_settings$maternal_id == "M_ID" && !anyDuplicated(family$dat$IID),
  "cross_sectional_one_row_per_child_family_intercept")
individual <- call(design = "longitudinal", random_individual = TRUE)
check(!individual$mixed_settings$family_random && individual$mixed_settings$individual_random &&
  individual$mixed_settings$maternal_id == "M_ID" && individual$mixed_settings$subject_id == "IID",
  "individual_only_keeps_nuclear_genetic_group_without_family_intercept")
check(is.null(individual$mixed_settings$interactions), "NULL_preserves_helper_automatic_longitudinal_age")
both <- call(design = "longitudinal", random_family = TRUE, random_individual = TRUE)
check(both$mixed_settings$family_random && both$mixed_settings$individual_random,
  "both_random_intercepts_are_requested_explicitly")
slope <- call(design = "longitudinal", random_individual = TRUE, random_slope = "age",
  genetic_interactions = c("age", "sex"), covariate_mapping_confirmed = TRUE)
check(identical(slope$mixed_settings$random_slope, "age") &&
  identical(slope$mixed_settings$interactions, c("age", "sex")) &&
  slope$mixed_settings$covariate_mapping_confirmed, "age_random_slope_and_age_sex_genetic_slopes_are_distinct")
custom <- call(random_family = TRUE, maternal_id = "Nuclear", subject_id = "Person",
  family_id = "Extended", diagnostic_cluster = "BroadDiagnostic")
check(custom$mixed_settings$maternal_id == "Nuclear" && custom$mixed_settings$subject_id == "Person" &&
  custom$family_id == "Extended" && custom$diagnostic_cluster == "BroadDiagnostic", "custom_IDs_forward_without_substitution")
interactions <- call(genetic_interactions = "sex", covariate_mapping_confirmed = TRUE)
check(!interactions$mixed_settings$family_random && !interactions$mixed_settings$individual_random &&
  identical(interactions$mixed_settings$interactions, "sex"), "interaction_only_uses_existing_no_random_mixed_model")
no_interactions <- call(genetic_interactions = character(0))
check(identical(no_interactions$mixed_settings$interactions, character(0)) &&
  !no_interactions$mixed_settings$family_random, "empty_interactions_explicitly_selects_no_interaction_mixed_model")
advanced <- list(maternal_id = "Nuclear", subject_id = "Person", family_random = FALSE,
  individual_random = TRUE, covariate_mapping_confirmed = TRUE)
check(identical(call(design = "longitudinal", mixed_settings = advanced)$mixed_settings, advanced),
  "advanced_settings_remain_unchanged_when_simple_options_are_inactive")
for (arguments in list(list(random_family = TRUE, mixed_settings = advanced),
    list(maternal_id = "OtherMother", mixed_settings = advanced),
    list(subject_id = "OtherPerson", mixed_settings = advanced),
    list(random_individual = TRUE), list(random_slope = "age"),
    list(random_family = TRUE, maternal_id = "missing_mother"))) {
  error <- tryCatch(do.call(call, arguments), error = identity)
  check(inherits(error, "error"), paste0("explicit_error_", paste(names(arguments), collapse = "_")))
}
workflow$ravel_parental_deviation_mixed <- mixed_function

# An all-off call must use the exact existing legacy body and main arguments.
reference <- workflow$ravel_parental_deviation
statements <- as.list(body(reference))[-1L]
start <- which(vapply(statements, function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
  identical(x[[2L]], as.name("design")), logical(1)))[1L]
body(reference) <- as.call(c(list(as.name("{")), statements[start:length(statements)]))
workflow$captured <- list()
workflow$ravel_fit <- function(...) {
  workflow$captured[[length(workflow$captured) + 1L]] <- list(...)
  stop("Controlled primary failure for unchanged legacy routing")
}
workflow$ravel_score_diagnostic <- function(...) list(status = "routing_fixture", reason = "Routing fixture", N = NA_integer_)
base <- list(dat = dat, pgs_columns = pgs, pi_se = NULL, include_comparison = FALSE)
old <- suppressMessages(do.call(reference, base)); old_calls <- workflow$captured
workflow$captured <- list()
now <- suppressMessages(do.call(workflow$ravel_parental_deviation, base))
check(identical(old_calls, workflow$captured) && identical(old, now), "all_off_NULL_preserves_existing_legacy_route_and_return")
rm("ravel_fit", "ravel_score_diagnostic", envir = workflow)

# Use the full helper to check primary/auxiliary routing while controlled fit
# failures avoid unnecessary estimator reruns. Diagnostic grouping must never
# split a known primary extended family into independent rows.
captured <- list(); captured_diagnostic <- list()
inspect <- function(...) testthat::with_mocked_bindings(
  suppressMessages(workflow$ravel_parental_deviation(dat, pgs_columns = pgs, pi_se = NULL, ...)),
  ravel_mixed_fit = function(...) { captured[[length(captured) + 1L]] <<- list(...); stop("Controlled primary stop") },
  ravel_score_diagnostic = function(...) { captured_diagnostic[[length(captured_diagnostic) + 1L]] <<- list(...)
    list(status = "routing_fixture", reason = "Routing fixture", N = NA_integer_) }, .package = "RAVEL")
z <- inspect(random_family = TRUE, maternal_id = "Nuclear", subject_id = "Person",
  family_id = "Extended", diagnostic_cluster = "BroadDiagnostic")
check(captured[[1L]]$family_id == "Nuclear" && captured[[1L]]$subject_id == "Person" &&
  captured[[1L]]$cluster == "Extended" && captured_diagnostic[[1L]]$cluster == "BroadDiagnostic",
  "primary_sandwich_and_diagnostic_clusters_are_independent_controls")
check(attr(z, "one_factor", exact = TRUE)$workflow$original_cluster == "BroadDiagnostic",
  "diagnostic_records_actual_requested_cluster")
z <- inspect(random_family = TRUE, family_id = "BroadDiagnostic", diagnostic_cluster = "auto")
check(identical(attr(z, "one_factor", exact = TRUE)$status, "invalid_design") &&
  grepl("split known extended families", attr(z, "one_factor", exact = TRUE)$reason),
  "auto_original_FID_cannot_split_a_coarser_primary_FID")
z <- inspect(random_family = TRUE, diagnostic_cluster = "IID")
check(identical(attr(z, "one_factor", exact = TRUE)$status, "invalid_design"), "custom_diagnostic_cluster_cannot_split_primary_FID")
z <- inspect(random_family = TRUE, diagnostic_cluster = NULL)
check(identical(attr(z, "one_factor", exact = TRUE)$status, "invalid_design"), "NULL_diagnostic_cluster_rejects_related_subjects")
error <- tryCatch(inspect(random_family = TRUE, diagnostic_cluster = "missing_cluster"), error = identity)
check(inherits(error, "error") && grepl("diagnostic_cluster", conditionMessage(error)), "missing_custom_diagnostic_column_rejected")

# Extract the actual RUN-block collector and perform one real simple-settings
# family-intercept fit with custom nuclear/person/extended-family columns.
assignment <- list()
for (expression in parse(path)) {
  if (!is.call(expression) || !identical(expression[[1L]], as.name("if")) ||
      !any(grepl("ravel.workflow_functions_only", deparse(expression[[2L]]), fixed = TRUE))) next
  block <- expression[[3L]]
  if (!is.call(block) || !identical(block[[1L]], as.name("{"))) next
  for (statement in as.list(block)[-1L]) if (is.call(statement) &&
    identical(statement[[1L]], as.name("<-")) && identical(statement[[2L]], as.name("run_outcome")))
    assignment[[length(assignment) + 1L]] <- statement
}
stopifnot(length(assignment) == 1L); eval(assignment[[1L]], workflow)
invisible(list2env(list(dat = dat, binary_outcome = FALSE, sib_configs = "none", sib_columns = NULL,
  models = list("Estimated state" = list(fit_constraints = "estimated_state")),
  design = "cross_sectional", pgs_columns = pgs, covariates = NULL, family_id = "Extended",
  subject_id = "Person", diagnostic_cluster = "BroadDiagnostic", maternal_id = "Nuclear",
  random_family = TRUE, random_individual = FALSE, random_slope = NULL,
  genetic_interactions = NULL, covariate_mapping_confirmed = FALSE, mixed_settings = NULL,
  age = "age", age_center = NULL, pi = .4, pi_se = NULL, pi_grid = NULL, grid_pi_se = NULL,
  include_raw = TRUE, drop_missing_outcome = TRUE), workflow))
result <- suppressMessages(workflow$run_outcome("Y"))
fit <- attr(result, "fits", exact = TRUE)$none[[1L]]
check(inherits(fit, "ravel_mixed_results") && fit$converged && fit$inference_status == "ok",
  "actual_simple_settings_collector_fit_is_regular")
check(fit$family_id == "Nuclear" && fit$subject_id == "Person" && fit$cluster == "Extended" &&
  fit$settings$family_random && !fit$settings$individual_random, "actual_simple_settings_fit_uses_requested_IDs_and_variance_components")
diagnostic <- attr(result, "diagnostics", exact = TRUE)$none[[1L]]
check(diagnostic$workflow$original_cluster == "BroadDiagnostic", "actual_auxiliary_fit_uses_custom_coarser_cluster")
collector_status <- data.frame(scenario = "cross_sectional_family_only",
  N = fit$N, maternal_families = fit$n_families,
  extended_families = fit$n_clusters, converged = fit$converged,
  inference_status = fit$inference_status, diagnostic_status = diagnostic$status,
  maternal_id = fit$family_id, subject_id = fit$subject_id, sandwich_cluster = fit$cluster,
  diagnostic_cluster = diagnostic$workflow$original_cluster)

# The requested individual-only longitudinal combination retains the same
# maternal genetic block even though its residual family intercept is disabled.
long_data <- mixed_validation_data(N_families = 250L, sib_sizes = 2L, visits = 3L,
  family_variance = 0, individual_variance = .3, age_pattern = "scheduled",
  age_effect = .03, sex_effect = 0)
long_data$Nuclear <- long_data$M_ID; long_data$Person <- long_data$IID
long_data$Extended <- long_data$FID
long_data$BroadDiagnostic <- paste0("B", (match(long_data$Extended, unique(long_data$Extended)) - 1L) %/% 2L)
workflow$dat <- long_data; workflow$design <- "longitudinal"
workflow$random_family <- FALSE; workflow$random_individual <- TRUE
workflow$age_center <- 10; workflow$covariate_mapping_confirmed <- TRUE
long_result <- suppressMessages(workflow$run_outcome("Y"))
long_fit <- attr(long_result, "fits", exact = TRUE)$none[[1L]]
check(inherits(long_fit, "ravel_mixed_results") && long_fit$converged &&
  long_fit$inference_status == "ok", "actual_individual_only_longitudinal_collector_fit_is_regular")
check(!long_fit$settings$family_random && long_fit$settings$individual_random &&
  long_fit$family_id == "Nuclear" && long_fit$subject_id == "Person" && long_fit$cluster == "Extended" &&
  "individual_intercept" %in% long_fit$variance_components$component &&
  !"family_intercept" %in% long_fit$variance_components$component,
  "actual_individual_only_fit_retains_genetic_family_but_omits_residual_family_variance")
check(identical(long_fit$settings$interactions, "age") &&
  all(c("beta_D:age", "beta_M:age", "beta_F:age") %in% long_result$parameter),
  "actual_simple_longitudinal_fit_preserves_automatic_genetic_age_slopes")
long_diagnostic <- attr(long_result, "diagnostics", exact = TRUE)$none[[1L]]
collector_status <- rbind(collector_status, data.frame(scenario = "longitudinal_individual_only",
  N = long_fit$N, maternal_families = long_fit$n_families, extended_families = long_fit$n_clusters,
  converged = long_fit$converged, inference_status = long_fit$inference_status,
  diagnostic_status = long_diagnostic$status, maternal_id = long_fit$family_id,
  subject_id = long_fit$subject_id, sandwich_cluster = long_fit$cluster,
  diagnostic_cluster = long_diagnostic$workflow$original_cluster))
utils::write.csv(collector_status, file.path(destination, "collector-status.csv"), row.names = FALSE)
writeLines(sub("[[:blank:]]+$", "", c(capture.output(sessionInfo()), "", "Source hashes:", capture.output(hashes))),
  file.path(destination, "session.txt"))
cat("Simple workflow setting checks passed:", length(checks), "\n")
