# Synthetic end-to-end reporting checks; no statistical estimator changes.
# Run from the package root: Rscript dev/check-paper-reporting.R
pkgload::load_all(".", quiet = TRUE)
RAVEL:::.ensure_lavaan_runtime()
source("tests/testthat/helper-family-validation.R")
source("tests/testthat/helper-mixed-validation.R")
source("tests/testthat/helper-mixed-covariates.R")
destination <- "dev/paper-reporting-validation"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
reuse <- "--reuse-fits" %in% args
cache_directory <- "/tmp/ravel-paper-reporting-cache"
dir.create(cache_directory, recursive = TRUE, showWarnings = FALSE)
export_root <- tempfile("ravel-paper-reporting-exports-", tmpdir = "/tmp")
dir.create(export_root)
child_path <- "analysis/parental_deviation/ravel_parental_deviation.R"
parent_path <- "analysis/parental_deviation/ravel_parental_phenotypes.R"
source_hashes <- tools::md5sum(c(child_path, parent_path, "dev/check-paper-reporting.R"))
previous <- options(ravel.workflow_functions_only = TRUE)
child <- new.env(parent = globalenv()); sys.source(child_path, child)
parent <- new.env(parent = globalenv()); source(parent_path, local = parent)
options(previous)
statements <- unlist(lapply(parse(child_path), function(x) {
  if (is.call(x) && identical(x[[1L]], as.name("if")) &&
      any(grepl("ravel.workflow_functions_only", deparse(x[[2L]]), fixed = TRUE)))
    as.list(x[[3L]])[-1L] else list()
}), recursive = FALSE)
collector <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
  identical(x[[2L]], as.name("run_outcome")), statements)
stopifnot(length(collector) == 1L)
eval(collector[[1L]], child)
checks <- list()
check <- function(ok, label, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(check = label,
    status = if (isTRUE(ok)) "pass" else "FAIL", detail = detail)
  write.csv(do.call(rbind, checks), file.path(destination, "checks.csv"), row.names = FALSE)
  if (!isTRUE(ok)) stop(label, ": ", detail)
}
equal <- function(x, y, tolerance = 1e-8) isTRUE(all.equal(x, y,
  tolerance = tolerance, check.attributes = FALSE))
pgs <- c(child1 = "PGS1_child", child2 = "PGS2_child", mother1 = "PGS1_mother",
  mother2 = "PGS2_mother", father1 = "PGS1_father", father2 = "PGS2_father")
sibs <- c(mat_sib1 = "PGS1_mat_sib", mat_sib2 = "PGS2_mat_sib",
  pat_sib1 = "PGS1_pat_sib", pat_sib2 = "PGS2_pat_sib")
model <- list("Estimated state" = list(fit_constraints = "estimated_state"))
configure <- function(dat, mixed = FALSE, covariates = NULL, sib_columns = NULL,
    pi_grid = NULL) {
  invisible(list2env(list(dat = dat, binary_outcome = FALSE,
    sib_configs = if (is.null(sib_columns)) "none" else c("none", "both"),
    models = model, design = "cross_sectional", pgs_columns = pgs,
    covariates = covariates, family_id = "FID", subject_id = "IID",
    diagnostic_cluster = "FID", maternal_id = "M_ID", random_family = mixed,
    random_individual = FALSE, random_slope = NULL, genetic_interactions = NULL,
    covariate_mapping_confirmed = FALSE, mixed_settings = NULL,
    age = "age", age_center = 10, pi = .4, pi_se = NULL, pi_grid = pi_grid,
    grid_pi_se = NULL, include_raw = TRUE, drop_missing_outcome = TRUE,
    sib_columns = sib_columns), child))
}
cache <- function(name, expression) {
  path <- file.path(cache_directory, paste0(name, ".rds"))
  if (reuse && file.exists(path)) return(readRDS(path))
  answer <- force(expression); saveRDS(answer, path); answer
}
family <- family_validation_data(family_validation_population(), N = 500L, seed = 420L)
family$IID <- family$FID
legacy <- cache("child-sem", {
  configure(family, sib_columns = sibs, pi_grid = .45)
  suppressMessages(child$run_outcome("Y"))
})
adult <- cache("parent-sem", suppressMessages(parent$ravel_parental_phenotypes(
  family, outcome = "Y", pgs_columns = pgs, sib_columns = sibs,
  models = model, pi = .4, pi_se = NULL, pi_grid = .45, grid_pi_se = NULL,
  subject_id = "IID", include_raw = TRUE, include_sem_audit = TRUE)))
mixed_data <- mixed_covariate_data(N_families = 300L,
  contract = mixed_covariate_contract(n_pc = 1L), family_variance = .5, seed = 721L)
mixed <- cache("child-mixed", {
  configure(mixed_data, mixed = TRUE, covariates = c("PC1", "sex"), pi_grid = .45)
  suppressMessages(child$run_outcome("Y"))
})
plain_data <- mixed_validation_data(N_families = 300L, family_variance = .5,
  sex_effect = 0, age_effect = 0, seed = 715L)
plain <- cache("child-mixed-no-covariates", {
  configure(plain_data, mixed = TRUE)
  suppressMessages(child$run_outcome("Y"))
})
objects <- list(child_sem = legacy, parent_sem = adult, child_mixed = mixed,
  child_mixed_no_covariates = plain)
if ("--prepare-only" %in% args) {
  cat("Synthetic model caches complete.\n")
  quit(status = 0L)
}
reports <- lapply(names(objects), function(name) child$ravel_reporting_tables(
  objects[[name]], workflow = if (name == "parent_sem") "parental" else "child"))
names(reports) <- names(objects)
required <- c("fit_statistics", "parameters", "genetic_structure", "diagnostics", "status")
for (name in names(reports)) {
  report <- reports[[name]]
  check(all(required %in% names(report)), paste0(name, "_common_report_tables"))
  check(all(vapply(report[required], is.data.frame, logical(1))), paste0(name, "_tables_are_data_frames"))
  check(!any(grepl("vcov|implied|sampling", names(report))), paste0(name, "_compact_exports_omit_large_matrices"))
  check(nrow(report$fit_statistics) > 0 && nrow(report$parameters) > 0,
    paste0(name, "_fit_and_parameter_results_present"))
}
# Additional schema-specific assertions are kept below; failures must remain
# visible rather than dropping a source fit or replacing an unavailable metric.
all_records <- function(result, workflow) {
  if (workflow == "parental") attr(result, "fits", exact = TRUE) else
    attr(result, "report_fit_records", exact = TRUE)
}
for (name in names(objects)) {
  result <- objects[[name]]; report <- reports[[name]]
  records <- all_records(result, if (name == "parent_sem") "parental" else "child")
  check(length(records) > 0, paste0(name, "_retains_source_records"))
  check(!anyDuplicated(names(records)) && all(nzchar(names(records))),
    paste0(name, "_source_keys_are_unique_and_nonempty"))
  check(all(vapply(records, function(z) !is.null(z$fit) && length(z$model) == 1L &&
    is.finite(z$pi), logical(1))), paste0(name, "_every_main_and_grid_record_has_fit_model_and_pi"))
  check(setequal(unique(report$fit_statistics$source_fit), names(records)),
    paste0(name, "_all_main_and_grid_fits_reported"))
  for (key in names(records)) {
    fit <- records[[key]]$fit
    check(!is.null(fit) && isTRUE(fit$converged), paste0(name, "_fixture_converged_", key))
    stats <- report$fit_statistics[report$fit_statistics$source_fit == key, , drop = FALSE]
    value <- function(stat) stats$value[tolower(stats$statistic) == tolower(stat)]
    if (inherits(fit$fit, "lavaan")) {
      native <- lavaan::fitMeasures(fit$fit)
      for (stat in c("aic", "bic", "chisq", "df", "cfi", "rmsea", "srmr"))
        check(length(value(stat)) == 1L && equal(value(stat), unname(native[stat])),
          paste0(name, "_native_", stat, "_", key))
    } else {
      for (stat in names(fit$fit_statistics))
        check(length(value(stat)) == 1L && equal(value(stat), unname(fit$fit_statistics[stat])),
          paste0(name, "_native_", stat, "_", key))
      if (!is.null(fit$settings$covariate)) {
        aic <- stats[tolower(stats$statistic) == "aic", , drop = FALSE]
        descriptive <- stats[stats$statistic == "AIC_conditional_moments_fixed", , drop = FALSE]
        check(is.na(aic$value) && nzchar(aic$note), paste0(name, "_ordinary_AIC_has_reason_", key))
        check(is.finite(descriptive$value) && grepl("descriptive", descriptive$status, ignore.case = TRUE),
          paste0(name, "_fixed_moments_AIC_is_explicitly_descriptive_", key))
      } else check(is.finite(value("AIC")), paste0(name, "_ordinary_AIC_retained_", key))
    }
  }
}
check(length(attr(legacy, "report_fit_records")) == 4L &&
  sum(vapply(attr(legacy, "report_fit_records"), function(z) z$kind == "grid", logical(1))) == 2L,
  "child_SEM_retains_unnamed_grid_fits_for_both_genetic_blocks")
check(length(attr(adult, "outcome_models")$both) == 15L,
  "parental_workflow_keeps_all_15_adult_subsets")
check(length(unique(reports$parent_sem$fit_statistics$source_fit)) == 2L,
  "parental_fit_statistics_describe_two_source_fits_not_thirty_subsets")

# Every implied liability-pair correlation is present on both visible and
# total scales, with uncertainty conditional on the supplied pi.
population <- family_validation_population()
for (name in c("child_sem", "parent_sem")) {
  records <- all_records(objects[[name]], if (name == "parent_sem") "parental" else "child")
  selected <- which(vapply(records, function(r) identical(r$kind, "main") &&
    identical(r$parental_sibs, "both"), logical(1)))
  check(length(selected) == 1L, paste0(name, "_one_full_five_role_source"))
  key <- names(records)[selected]; record <- records[[selected]]
  descriptor <- record$fit$family_descriptor
  if (is.null(descriptor)) descriptor <- RAVEL:::.ravel_family_legacy_descriptor(record$fit)
  for (block in c("visible_correlation", "total_correlation")) {
    rows <- reports[[name]]$genetic_structure
    rows <- rows[rows$source_fit == key & rows$block == block & rows$lhs != rows$rhs, , drop = FALSE]
    check(nrow(rows) == 10L && !anyDuplicated(paste(rows$lhs, rows$rhs)),
      paste0(name, "_ten_unique_", block, "_pairs"))
    expected <- stats::cov2cor(if (block == "visible_correlation") population$Sigma_GG else population$Sigma_TT)
    oracle <- expected[cbind(match(rows$lhs, rownames(expected)), match(rows$rhs, colnames(expected)))]
    check(equal(rows$estimate, oracle, tolerance = 3e-4), paste0(name, "_", block, "_independent_population_oracle"))
    check(all(is.finite(rows$se)) && all(is.finite(rows$ci_low)) && all(is.finite(rows$ci_high)),
      paste0(name, "_", block, "_all_pair_SEs_and_CIs"))
    for (pair in list(c("mother", "mat_sib"), c("father", "pat_sib"))) {
      row <- rows[(rows$lhs == pair[1L] & rows$rhs == pair[2L]) |
        (rows$rhs == pair[1L] & rows$lhs == pair[2L]), , drop = FALSE]
      check(nrow(row) == 1L, paste0(name, "_", block, "_", paste(pair, collapse = "_"), "_included"))
      fun <- function(theta) {
        m <- descriptor$moments(theta, record$pi)
        stats::cov2cor(if (block == "visible_correlation") m$Sigma_GG else m$Sigma_TT)[pair[1L], pair[2L]]
      }
      gradient <- numDeriv::grad(fun, descriptor$theta)
      se <- sqrt(drop(gradient %*% descriptor$vcov %*% gradient))
      check(equal(row$se, se, tolerance = 1e-6),
        paste0(name, "_", block, "_", paste(pair, collapse = "_"), "_joint_delta_SE"))
    }
  }
}
for (name in c("child_mixed", "child_mixed_no_covariates")) {
  records <- all_records(objects[[name]], "child")
  key <- names(records)[1L]; fit <- records[[1L]]$fit
  parameters <- reports[[name]]$parameters
  parameters <- parameters[parameters$source_fit == key, , drop = FALSE]
  check(!any(parameters$kind == "internal_parameter"), paste0(name, "_no_optimizer_coordinate_dump"))
  roles <- c("mother", "father", "child")
  load_names <- unlist(lapply(1:2, function(k) paste0("loading", k, "_", roles)))
  theta <- c(fit$parameters, fit$means[setdiff(names(fit$means), names(fit$parameters))])
  state <- RAVEL:::.ravel_mixed_state(theta, fit$settings)
  load <- exp(theta[load_names]) * rep(sqrt(c(state$visible_parent, state$visible_parent, state$visible_child)), 2L)
  residual <- exp(2 * theta[load_names] + 2 * rep(theta[c("noise1", "noise2")], each = 3L))
  rows <- parameters[match(load_names, parameters$parameter), , drop = FALSE]
  check(equal(rows$estimate, load) && all(is.finite(rows$se)), paste0(name, "_natural_measurement_loadings"))
  residual_names <- sub("loading", "score_residual_variance", load_names, fixed = TRUE)
  rows <- parameters[match(residual_names, parameters$parameter), , drop = FALSE]
  check(equal(rows$estimate, residual) && all(is.finite(rows$se)), paste0(name, "_natural_measurement_residual_variances"))
  if (length(fit$settings$covariates)) {
    covariates <- fit$settings$covariates; fixed_names <- paste0("fixed_", covariates)
    expected <- theta[fixed_names] * fit$settings$outcome_scale / fit$settings$scales[covariates]
    rows <- parameters[match(fixed_names, parameters$parameter), , drop = FALSE]
    check(equal(rows$estimate, expected) && all(is.finite(rows$se)), paste0(name, "_fixed_slopes_in_original_covariate_units"))
  }
  row <- parameters[parameters$parameter == "loading1_mother", , drop = FALSE]
  gradient <- numDeriv::grad(function(x) exp(x["loading1_mother"]) *
    sqrt(RAVEL:::.ravel_mixed_state(x, fit$settings)$visible_parent), theta)
  se <- sqrt(drop(gradient %*% fit$vcov[names(theta), names(theta)] %*% gradient))
  check(equal(row$se, se, tolerance = 1e-6), paste0(name, "_loading_SE_uses_joint_state_measurement_covariance"))
}

# Unavailable source fits and diagnostic p-values remain explicit.
failed <- legacy
record <- attr(legacy, "report_fit_records", exact = TRUE)[[1L]]
record$fit <- NULL; record$note <- "Controlled failed source fit"
attr(failed, "report_fit_records") <- list(controlled_failure = record)
failed_report <- child$ravel_reporting_tables(failed)
check(any(failed_report$status$source_fit == "controlled_failure") &&
  any(grepl("Controlled failed source fit", failed_report$status$note, fixed = TRUE)),
  "failed_source_fit_retained_with_reason")
not_converged <- legacy
records <- attr(not_converged, "report_fit_records", exact = TRUE)
records[[1L]]$fit$converged <- FALSE
attr(not_converged, "report_fit_records") <- records
bad_report <- child$ravel_reporting_tables(not_converged)
check(any(bad_report$status$source_fit == names(records)[1L] &
  bad_report$status$status == "not_converged"), "nonconverged_source_fit_not_presented_as_success")
for (name in names(objects)) {
  native_rows <- objects[[name]][which(objects[[name]]$parameter == "one_factor_test"), , drop = FALSE]
  diagnostic <- reports[[name]]$diagnostics
  check(nrow(diagnostic) == nrow(native_rows), paste0(name, "_diagnostic_count_matches_main_results"))
  check(equal(sort(diagnostic$p_value, na.last = TRUE), sort(native_rows$p_value, na.last = TRUE)),
    paste0(name, "_diagnostic_p_values_preserved"))
  check(all(grepl("unadjusted", diagnostic$note, ignore.case = TRUE)),
    paste0(name, "_diagnostic_scope_is_explicit"))
}

# Use the exact writer selected by both scripts, and verify decimal conventions.
for (name in names(objects)) {
  for (decimal in c(".", ",")) {
    code <- if (decimal == ".") "point" else "comma"
    directory <- file.path(export_root, name, code)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    result <- objects[[name]]; before <- result
    prefix <- file.path(directory, "result")
    paths <- child$write_ravel_workflow_outputs(result, prefix,
      workflow = if (name == "parent_sem") "parental" else "child", decimal = decimal)
    files <- list.files(directory, full.names = TRUE)
    csv <- files[grepl("[.]csv$", files)]
    check(length(csv) == 6L && all(grepl("_excel[.]csv$", csv)),
      paste0(name, "_", code, "_only_six_Excel_CSVs"))
    sep <- if (decimal == ",") ";" else ","
    for (path in csv) check(identical(readLines(path, n = 1L, warn = FALSE), paste0("sep=", sep)),
      paste0(name, "_", code, "_separator_header_", basename(path)))
    main <- read.table(paste0(prefix, "_excel.csv"), header = TRUE, sep = sep,
      dec = decimal, skip = 1L, quote = "\"", comment.char = "", check.names = FALSE,
      na.strings = "", stringsAsFactors = FALSE)
    check(identical(names(main), names(result)) && nrow(main) == nrow(result),
      paste0(name, "_", code, "_main_schema_unchanged"))
    check(equal(main$estimate, result$estimate) && equal(main$se, result$se),
      paste0(name, "_", code, "_coefficient_numbers_roundtrip"))
    saved <- readRDS(paste0(prefix, "_fits.rds"))
    check(length(all_records(saved, if (name == "parent_sem") "parental" else "child")) > 0,
      paste0(name, "_", code, "_RDS_retains_full_fit_objects"))
    check(identical(result, before), paste0(name, "_", code, "_writer_does_not_mutate_primary_results"))
  }
}
text_fixture <- data.frame(label = c('Comma, quote " and semicolon;', "Two\nlines", "ordinary"),
  value = c(1.25, -2.75, NA_real_), stringsAsFactors = FALSE)
for (decimal in c(".", ",")) {
  path <- file.path(export_root, paste0("quoted-", if (decimal == ".") "point" else "comma", "_excel.csv"))
  child$write_excel_csv(text_fixture, path, decimal = decimal)
  sep <- if (decimal == ",") ";" else ","
  roundtrip <- read.table(path, header = TRUE, sep = sep, dec = decimal, skip = 1L,
    quote = "\"", comment.char = "", check.names = FALSE, na.strings = "", stringsAsFactors = FALSE)
  check(equal(roundtrip, text_fixture), paste0("Excel_quoting_and_multiline_roundtrip_", decimal))
}
writeLines(sub("[[:space:]]+$", "", c(capture.output(sessionInfo()),
  capture.output(source_hashes))),
  file.path(destination, "session.txt"))
writeLines(export_root, file.path(destination, "export-location.txt"))
print(do.call(rbind, checks), row.names = FALSE)
