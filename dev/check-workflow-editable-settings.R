# Focused checks for the short SETTINGS block; no changes to the estimators.
# Run from the package root: Rscript dev/check-workflow-editable-settings.R
pkgload::load_all(".", quiet = TRUE)
RAVEL:::.ensure_lavaan_runtime()
source("tests/testthat/helper-mixed-validation.R")
source("tests/testthat/helper-mixed-covariates.R")
path <- "analysis/parental_deviation/ravel_parental_deviation.R"
destination <- "dev/mixed-validation/editable-settings-check"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
previous <- getOption("ravel.workflow_functions_only")
options(ravel.workflow_functions_only = TRUE)
w <- new.env(parent = globalenv()); sys.source(path, w)
options(ravel.workflow_functions_only = previous)
checks <- list()
check <- function(ok, name) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name,
    status = if (isTRUE(ok)) "pass" else "FAIL")
  write.csv(do.call(rbind, checks), file.path(destination, "checks.csv"), row.names = FALSE)
  if (!isTRUE(ok)) stop(name)
}
controls <- function(...) do.call(w$.ravel_script_controls, modifyList(list(
  design = "cross_sectional", se_cluster = "Outer", iid = "Person",
  family_random_effect = NULL, individual_random_effect = FALSE,
  interactions = NULL, nuclear_family_id = "Nuclear", age = "VisitAge",
  random_slope = NULL, sib_columns = NULL), list(...), keep.null = TRUE))
x <- controls()
check(is.null(x$genetic_interactions) && !x$random_family && !x$random_individual,
  "cross_all_off_keeps_legacy_route")
check(x$family_id == "Outer" && x$subject_id == "Person" && x$diagnostic_cluster == "Outer",
  "one_cluster_setting_and_one_person_setting_reach_both_uses")
x <- controls(family_random_effect = "Mother", nuclear_family_id = "unused")
check(x$random_family && x$maternal_id == "Mother" && identical(x$genetic_interactions, character()),
  "family_column_enables_intercept_and_identifies_shared_parents")
x <- controls(design = "longitudinal", individual_random_effect = TRUE)
check(!x$random_family && x$random_individual && x$maternal_id == "Nuclear",
  "individual_only_uses_explicit_genetic_family_fallback")
check(identical(x$genetic_interactions, character()), "NULL_never_adds_an_implicit_genetic_age_slope")
check(identical(controls(design = "long")$design, "longitudinal"), "design_abbreviations_are_canonicalized")
x <- controls(design = "longitudinal", interactions = c("VisitAge", "sex"))
check(identical(x$genetic_interactions, c("VisitAge", "sex")), "requested_interactions_are_explicit")
x <- controls(design = "longitudinal", individual_random_effect = TRUE, random_slope = "VisitAge")
check(x$random_individual && identical(x$genetic_interactions, character()),
  "individual_random_slope_does_not_enable_genetic_interactions")
x <- controls(design = "longitudinal", interactions = "VisitAge", sib_columns = c(a = "sib"))
check(is.null(x$genetic_interactions) && is.null(x$mixed_settings), "existing_parental_sibling_age_workflow_retained")
for (args in list(list(family_random_effect = TRUE), list(individual_random_effect = "Person"),
    list(individual_random_effect = TRUE), list(interactions = c("sex", "sex")),
    list(design = "longitudinal", random_slope = "VisitAge"),
    list(design = "longitudinal", se_cluster = NULL),
    list(design = "longitudinal", sib_columns = c(a = "sib")))) {
  z <- tryCatch(do.call(controls, args), error = identity)
  check(inherits(z, "error"), paste0("invalid_setting_rejected_", length(checks) + 1L))
}

# Execute the real RUN-block bridge and collector, not copies of their logic.
blocks <- lapply(parse(path), function(x) {
  if (is.call(x) && identical(x[[1L]], as.name("if")) &&
      any(grepl("ravel.workflow_functions_only", deparse(x[[2L]]), fixed = TRUE)))
    as.list(x[[3L]])[-1L] else list()
})
statements <- unlist(blocks, recursive = FALSE)
bridge <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("invisible")) &&
  any(grepl(".ravel_script_controls", deparse(x), fixed = TRUE)), statements)
age_fixed <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("if")) &&
  identical(x[[2L]], quote(design == "longitudinal")), statements)
collector <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
  identical(x[[2L]], as.name("run_outcome")), statements)
stopifnot(length(bridge) == 1L, length(age_fixed) == 1L, length(collector) == 1L)
eval(collector[[1L]], w)
invisible(list2env(list(design = "long", se_cluster = "Outer", iid = "Person",
  family_random_effect = NULL, individual_random_effect = TRUE, interactions = NULL,
  nuclear_family_id = "Nuclear", age = "VisitAge", random_slope = NULL,
  sib_columns = NULL, covariates = "PC1"), w))
eval(bridge[[1L]], w); eval(age_fixed[[1L]], w)
check(w$design == "longitudinal" && "VisitAge" %in% w$covariates,
  "actual_bridge_adds_fixed_age_for_abbreviated_long_design")
w$individual_random_effect <- FALSE; w$interactions <- "VisitAge"
w$sib_columns <- c(a = "sib"); w$covariates <- c("PC1", "VisitAge")
eval(bridge[[1L]], w); eval(age_fixed[[1L]], w)
check(identical(w$covariates, "PC1") && is.null(w$genetic_interactions),
  "legacy_sibling_interaction_keeps_age_out_of_inplace_covariate_scaling")
run <- function(long) {
  dat <- mixed_covariate_data(N_families = 300L, visits = if (long) 3L else 1L,
    contract = mixed_covariate_contract(n_pc = 1L),
    family_variance = if (long) 0 else .5, individual_variance = if (long) .3 else 0,
    sex_beta = if (long) c(child = .12, mother = -.06, father = .02) else c(child = 0, mother = 0, father = 0),
    seed = if (long) 1963L else 721L)
  dat$Nuclear <- dat$M_ID; dat$Person <- dat$IID
  dat$Outer <- dat$FID; dat$VisitAge <- dat$age
  settings <- list(dat = dat, design = if (long) "longitudinal" else "cross_sectional",
    se_cluster = "Outer", iid = "Person", family_random_effect = if (long) NULL else "Nuclear",
    nuclear_family_id = "Nuclear", individual_random_effect = long,
    interactions = if (long) "sex" else NULL, age = "VisitAge", age_center = 10,
    random_slope = NULL, covariates = c("sex", "PC1"), sib_columns = NULL,
    sib_configs = "none", binary_outcome = FALSE, pi = .4, pi_se = NULL,
    pi_grid = NULL, grid_pi_se = NULL, include_raw = TRUE, drop_missing_outcome = TRUE,
    models = list("Estimated state" = list(fit_constraints = "estimated_state")),
    pgs_columns = c(child1 = "PGS1_child", child2 = "PGS2_child", mother1 = "PGS1_mother",
      mother2 = "PGS2_mother", father1 = "PGS1_father", father2 = "PGS2_father"))
  list2env(settings, w); eval(bridge[[1L]], w); eval(age_fixed[[1L]], w)
  result <- suppressMessages(w$run_outcome("Y"))
  fit <- attr(result, "fits", exact = TRUE)$none[[1L]]
  label <- if (long) "longitudinal" else "cross_sectional"
  check(isTRUE(fit$converged) && fit$inference_status == "ok", paste0(label, "_real_fit_regular"))
  check(fit$cluster == "Outer" && fit$family_id == "Nuclear" && fit$subject_id == "Person",
    paste0(label, "_custom_IDs_preserved"))
  check(identical(fit$settings$family_random, !long) && identical(fit$settings$individual_random, long),
    paste0(label, "_only_requested_random_intercepts"))
  check(!any(grepl(":VisitAge$", result$parameter)) &&
    identical(any(grepl(":sex$", result$parameter)), long), paste0(label, "_only_requested_genetic_slopes_reported"))
  diagnostic <- attr(result, "diagnostics", exact = TRUE)$none[[1L]]
  check(diagnostic$workflow$subject_id == "Person" && diagnostic$workflow$original_cluster == "Outer",
    paste0(label, "_one_factor_uses_requested_IDs"))
  if (long) check("VisitAge" %in% fit$settings$covariates,
    "measurement_age_is_a_fixed_covariate_without_a_genetic_age_slope")
  data.frame(design = label, N = fit$N, converged = fit$converged,
    inference_status = fit$inference_status, diagnostic_status = diagnostic$status)
}
status <- rbind(run(FALSE), run(TRUE))
write.csv(status, file.path(destination, "fits.csv"), row.names = FALSE)
hashes <- tools::md5sum(c(path, "dev/check-workflow-editable-settings.R"))
writeLines(sub("[[:space:]]+$", "", c(capture.output(sessionInfo()), capture.output(hashes))),
  file.path(destination, "session.txt"))
print(do.call(rbind, checks), row.names = FALSE)
