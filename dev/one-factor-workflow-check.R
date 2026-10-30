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

  # An explicitly selected independent-family ID takes precedence over a
  # dummy PLINK FID. The workflow must remove the unused FID before fitting.
  alias <- example
  alias$independent_family <- alias$FID
  alias$FID <- 0L
  aliased <- run_case("named_ID_with_dummy_FID", alias, "ok",
    fit_constraints = "estimated_state", family_id = "independent_family")
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
  run_case("maternal_phenotype_missing_sibling", missing_sib, "unsupported_missingness",
    outcome = "mom_height", fit_constraints = "estimated_state",
    parental_sibs = "both", sib_columns = sib_cols)
  missing_y <- example
  missing_y$Y[seq_len(20L)] <- NA_real_
  run_case("retained_missing_outcome", missing_y, "unsupported_missingness",
    fit_constraints = "estimated_state", drop_missing_outcome = FALSE)
  repeated <- example
  repeated$FID <- rep(seq_len(ceiling(nrow(repeated) / 2)), each = 2L)[seq_len(nrow(repeated))]
  run_case("repeated_family", repeated, "unsupported_design",
    fit_constraints = "estimated_state")

  summary <- do.call(rbind, answers)
  rownames(summary) <- NULL
  utils::write.csv(summary, "dev/one-factor-check/workflow-summary.csv", row.names = FALSE)
  saveRDS(diagnostics, "dev/one-factor-check/workflow-diagnostics.rds")
  writeLines(c(capture.output(sessionInfo()), "", "Source hashes:",
    capture.output(tools::md5sum(c("R/ravel_score_signal.R",
      "analysis/parental_deviation/ravel_parental_deviation.R",
      "dev/one-factor-workflow-check.R")))),
    "dev/one-factor-check/workflow-session.txt")
  print(summary)
  cat("One-factor workflow checks passed: ", nrow(summary), " cases.\n", sep = "")
})
