# Run from the repository root: Rscript dev/one-factor-workflow-check.R
# Only bundled or synthetic data are used. Sourcing the workflow definitions
# does not run their SETTINGS/RUN blocks or write user analysis outputs.
pkgload::load_all(".", quiet = TRUE)
local({
  previous <- options(ravel.workflow_functions_only = TRUE)
  on.exit(options(previous), add = TRUE)
  source("analysis/parental_deviation/ravel_parental_deviation.R", local = TRUE)
  source("tests/testthat/helper-family-validation.R", local = TRUE)
  dir.create("dev/one-factor-check", showWarnings = FALSE, recursive = TRUE)
  cols <- c(child1 = "PGS1_child", child2 = "PGS2_child",
    mother1 = "PGS1_mother", mother2 = "PGS2_mother",
    father1 = "PGS1_father", father2 = "PGS2_father")
  sib_cols <- c(mat_sib1 = "PGS1_mat_sib", mat_sib2 = "PGS2_mat_sib",
    pat_sib1 = "PGS1_pat_sib", pat_sib2 = "PGS2_pat_sib")
  answers <- diagnostics <- list()
  run_case <- function(label, dat, expected, ...) {
    message("Checking one-factor workflow: ", label)
    answer <- ravel_parental_deviation(dat, pgs_columns = cols,
      pi = .528, pi_se = NULL, pi_grid = NULL, include_comparison = FALSE, ...)
    diagnostic <- attr(answer, "one_factor")
    stopifnot(is.list(diagnostic), identical(diagnostic$status, expected))
    row <- answer[answer$parameter == "one_factor_test", , drop = FALSE]
    stopifnot(nrow(row) == 1L, identical(row$diagnostic_status, expected),
      isTRUE(all.equal(row$diagnostic_df, diagnostic$df %||% NA_real_)))
    if (expected == "ok") {
      stopifnot(is.finite(row$estimate), is.finite(row$p_value),
        identical(row$estimate, diagnostic$statistic),
        identical(row$p_value, diagnostic$p_value), length(diagnostic$fits) == 3L)
    } else stopifnot(is.na(row$estimate), is.na(row$p_value))
    diagnostics[[label]] <<- diagnostic
    answers[[label]] <<- data.frame(case = label, expected_status = expected,
      status = diagnostic$status, reason = diagnostic$reason,
      statistic = row$estimate, df = row$diagnostic_df, p_value = row$p_value,
      N = row$diagnostic_N, stringsAsFactors = FALSE)
    answer
  }
  `%||%` <- function(x, y) if (is.null(x)) y else x
  data("ravel_example", package = "RAVEL", envir = environment())
  example <- ravel_example
  names(example) <- sub("^pgs", "PGS", names(example))
  for (mode in c("estimated_state", "n_gen", "equilibrium"))
    run_case(paste0("bundled_", mode), example, "ok",
      fit_constraints = mode, am_n_generations = 3L)
  run_case("bundled_free_measurement", example, "unsupported_free_measurement",
    fit_constraints = "none")

  # Explicitly identify a known dummy FID for the auxiliary diagnostic. Auto
  # mode must preserve a genuine extended-family FID even when main family_id
  # uses a different column; group size cannot distinguish these situations.
  alias <- example
  alias$independent_family <- alias$FID
  alias$FID <- 0L
  aliased <- run_case("named_ID_with_dummy_FID", alias, "ok",
    fit_constraints = "estimated_state", family_id = "independent_family",
    diagnostic_cluster = "independent_family")
  stopifnot(!isTRUE(attr(aliased, "fit")$score_signal_input$dependent_families),
    !"FID" %in% names(attr(aliased, "analysis_data")))

  # No true child effect: the phenotype belongs to the mother. Scores for the
  # child and every adult still belong to the joint diagnostic source model.
  population <- family_validation_population(pi = .528, rho = .35,
    parent_variance = 1.25,
    beta = c(child = 0, mother = .3, father = .04, mat_sib = .02, pat_sib = 0),
    lineage = c(maternal = .1, paternal = 0))
  parent <- family_validation_data(population, N = 2000L, seed = 4317L, exact = FALSE)
  parent$mom_height <- parent$Y
  run_case("maternal_phenotype_complete", parent, "ok", outcome = "mom_height",
    fit_constraints = "estimated_state", parental_sibs = "both", sib_columns = sib_cols)
  missing_sib <- parent
  missing_sib$PGS1_mat_sib[seq_len(40L)] <- NA_real_
  missing_sib_result <- run_case("maternal_phenotype_missing_sibling", missing_sib, "ok",
    outcome = "mom_height", fit_constraints = "estimated_state",
    parental_sibs = "both", sib_columns = sib_cols)
  stopifnot(attr(missing_sib_result, "one_factor")$selection$N_incomplete_selected == 40L,
    attr(missing_sib_result, "one_factor")$selection$N_used == nrow(parent) - 40L)
  missing_y <- example
  missing_y$Y[seq_len(20L)] <- NA_real_
  missing_y_result <- run_case("retained_missing_outcome", missing_y, "ok",
    fit_constraints = "estimated_state", drop_missing_outcome = FALSE)
  stopifnot(attr(missing_y_result, "one_factor")$selection$N_incomplete_selected == 20L,
    nrow(attr(missing_y_result, "analysis_data")) == nrow(example))
  repeated <- example
  repeated$FID <- rep(seq_len(ceiling(nrow(repeated) / 2)), each = 2L)[seq_len(nrow(repeated))]
  repeated$IID <- sprintf("I%07d", seq_len(nrow(repeated)))
  repeated_result <- run_case("repeated_family", repeated, "ok",
    fit_constraints = "estimated_state", subject_id = "IID")
  repeated_diagnostic <- attr(repeated_result, "one_factor")
  stopifnot(repeated_diagnostic$selection$N_selected == length(unique(repeated$FID)),
    !anyDuplicated(repeated_diagnostic$selection$clusters_used))
  adjusted <- example
  set.seed(8201L)
  adjusted$sex <- stats::rbinom(nrow(adjusted), 1L, .5)
  adjusted_result <- run_case("covariates_omitted", adjusted, "ok",
    fit_constraints = "estimated_state", covariates = "sex")
  diagnostic_row <- adjusted_result[adjusted_result$parameter == "one_factor_test", , drop = FALSE]
  stopifnot(grepl("Covariates were not used to compute the one-factor test", diagnostic_row$note, fixed = TRUE),
    identical(attr(adjusted_result, "one_factor")$auxiliary$covariates_omitted, "sex"))

  summary <- do.call(rbind, answers)
  rownames(summary) <- NULL
  utils::write.csv(summary, "dev/one-factor-check/workflow-summary.csv", row.names = FALSE)
  writeLines(c(capture.output(sessionInfo()), "", "Source hashes:",
    capture.output(tools::md5sum(c("R/ravel_score_signal.R", "R/ravel_score_diagnostic.R",
      "analysis/parental_deviation/ravel_parental_deviation.R",
      "dev/one-factor-workflow-check.R")))),
    "dev/one-factor-check/workflow-session.txt")
  print(summary)
  cat("One-factor workflow checks passed: ", nrow(summary), " cases.\n", sep = "")
})
