# Independent checks of the parental-phenotype analysis workflow.
# Run from the repository root: Rscript dev/run-parental-phenotype-checks.R
# Uses synthetic data only; no settings-driven analysis is run while sourcing.
pkgload::load_all(".", quiet = TRUE)
RAVEL:::.ensure_lavaan_runtime()
source("tests/testthat/helper-family-validation.R")

workflow_env <- new.env(parent = globalenv())
previous_options <- options(ravel.workflow_functions_only = TRUE)
source("analysis/parental_deviation/ravel_parental_phenotypes.R", local = workflow_env)
options(previous_options)
stopifnot(is.function(workflow_env$ravel_parental_phenotypes),
          !exists("results", envir = workflow_env, inherits = FALSE))
workflow <- workflow_env$ravel_parental_phenotypes

checks_passed <- 0L
check <- function(label, expression) {
  force(expression)
  checks_passed <<- checks_passed + 1L
  cat(sprintf("PASS %02d: %s\n", checks_passed, label))
}
equal <- function(observed, expected, tolerance = 3e-4) {
  result <- all.equal(observed, expected, tolerance = tolerance)
  if (!isTRUE(result)) stop(paste(result, collapse = "\n"), call. = FALSE)
}
rejects <- function(expression, pattern) {
  answer <- tryCatch(force(expression), error = identity)
  stopifnot(inherits(answer, "error"), grepl(pattern, conditionMessage(answer), ignore.case = TRUE))
}
roles <- c("child", "mother", "father", "mat_sib", "pat_sib")
pgs <- stats::setNames(unlist(lapply(roles[1:3], function(r) paste0("PGS", 1:2, "_", r))),
                       unlist(lapply(roles[1:3], function(r) paste0(r, 1:2))))
sibs <- stats::setNames(unlist(lapply(roles[4:5], function(r) paste0("PGS", 1:2, "_", r))),
                        unlist(lapply(roles[4:5], function(r) paste0(r, 1:2))))
population <- family_validation_population()
dat <- family_validation_data(population, N = 1000L, seed = 420L)
state_model <- list("Estimated state" = list(fit_constraints = "estimated_state"))
run <- function(data = dat, ...) {
  arguments <- list(dat = data, outcome = "Y", pgs_columns = pgs, sib_columns = sibs,
                    models = state_model, pi = .4, subject_id = "FID", include_raw = TRUE)
  more <- list(...)
  for (name in names(more)) arguments[name] <- more[name]
  do.call(workflow, arguments)
}
parameter <- c(mother = "beta_M", father = "beta_F", mat_sib = "beta_MS", pat_sib = "beta_FS")
coefficient_rows <- function(x) x[x$parameter %in% unname(parameter), , drop = FALSE]
estimates <- function(x, analysis, model, scale_roles, pi = NULL) {
  keep <- x$analysis == analysis & x$outcome_model == model & x$parameter %in% parameter[scale_roles]
  if (!is.null(pi)) keep <- keep & !is.na(x$pi) & abs(x$pi - pi) < 1e-12
  y <- x[which(keep), , drop = FALSE]
  stopifnot(nrow(y) == length(scale_roles), !anyDuplicated(y$parameter))
  stats::setNames(y$estimate[match(parameter[scale_roles], y$parameter)], scale_roles)
}
sets <- function(table, config = "both") attr(table, "outcome_models")[[config]]

main <- run(pi_se = .03, pi_grid = c(.4, .5), grid_pi_se = NULL)
check("legacy table schema and fitted-object attributes survive the new workflow", {
  stopifnot(is.data.frame(main), all(c("siblings", "model", "analysis", "pi", "pi_se",
    "effect", "parameter", "estimate", "se", "ci_low", "ci_high", "p_value", "note",
    "outcome_model", "outcome_predictors", "genetic_members", "status", "inference_status", "N") %in% names(main)),
    all(vapply(c("fits", "settings", "diagnostics", "projections"),
      function(nm) is.list(attr(main, nm)), logical(1))))
})
check("default maternal analysis contains precisely eight focal-preserving subsets", {
  stopifnot(length(sets(main)) == 8L,
    all(vapply(sets(main), function(x) "mother" %in% x && !"child" %in% x, logical(1))),
    setequal(vapply(sets(main), function(x) paste(sort(x), collapse = "+"), ""),
      vapply(family_validation_subsets("mother"), function(x) paste(sort(x), collapse = "+"), "")))
})
check("every corrected maternal regression matches the independent population oracle", {
  outcome_sd <- stats::sd(dat$Y)
  for (model in names(sets(main))) {
    selected <- sets(main)[[model]]
    oracle <- family_validation_oracle(population, selected)
    equal(estimates(main, "Corrected for pi (founder scale)", model, selected),
          oracle$beta_raw / outcome_sd)
    equal(estimates(main, "Corrected for pi", model, selected), oracle$beta_std / outcome_sd)
  }
})
check("child scores and nonzero nuisance association remain, but child coefficients are not reported", {
  records <- attr(main, "fits")
  source <- Filter(function(x) identical(x$kind, "main"), records)[[1L]]$fit
  stopifnot(source$converged, all(c("PGS1_child", "PGS2_child") %in% lavaan::lavNames(source$fit, "ov")),
    abs(source$beta_total_assumed[["bD"]]) > .08,
    !any(grepl("^(beta_D|beta_child|child)($|:)", main$parameter)),
    all(grepl("child", coefficient_rows(main)$genetic_members, fixed = TRUE)))
  full_name <- names(sets(main))[vapply(sets(main), length, integer(1)) == 4L]
  adult_mother <- estimates(main, "Corrected for pi (founder scale)", full_name, roles[-1])[["mother"]]
  stopifnot(abs(adult_mother - source$beta_total_assumed[["bm"]]) > .05)
})
check("plain PGS comparators refit the selected adults with their own contrast covariance", {
  for (model in names(sets(main))) for (version in 1:2) {
    selected <- sets(main)[[model]]
    x <- as.data.frame(scale(dat[paste0("PGS", version, "_", selected)]))
    names(x) <- selected
    x$Y <- as.numeric(scale(dat$Y))
    ordinary <- stats::lm(stats::reformulate(selected, "Y"), data = x)
    rows <- main[main$analysis == paste0("Plain PGS", version, " regression") &
                   main$outcome_model == model, , drop = FALSE]
    equal(stats::setNames(rows$estimate[match(parameter[selected], rows$parameter)], selected),
          stats::coef(ordinary)[selected], tolerance = 1e-8)
    equal(stats::setNames(rows$se[match(parameter[selected], rows$parameter)], selected),
          sqrt(diag(stats::vcov(ordinary)))[selected], tolerance = 1e-8)
    if (all(c("mother", "mat_sib") %in% selected)) {
      V <- stats::vcov(ordinary)
      delta <- rows[rows$parameter == "delta_M", ]
      equal(delta$se, sqrt(V["mother", "mother"] + V["mat_sib", "mat_sib"] -
                            2 * V["mother", "mat_sib"]), tolerance = 1e-8)
    }
  }
})
check("main pi uncertainty and exact grid values remain distinct at the same pi", {
  rows <- coefficient_rows(main)
  point <- rows[rows$analysis == "Corrected for pi" & rows$pi == .4, ]
  grid <- rows[rows$analysis == "Grid" & rows$pi == .4, ]
  key <- function(x) paste(x$outcome_model, x$parameter)
  grid <- grid[match(key(point), key(grid)), ]
  equal(point$estimate, grid$estimate, tolerance = 1e-7)
  stopifnot(all(point$pi_se == .03), all(is.na(grid$pi_se) | grid$pi_se == 0),
    all(point$se >= grid$se - 1e-7), max(point$se - grid$se) > 1e-5)
  records <- attr(main, "fits")
  stopifnot(sum(vapply(records, function(x) x$kind == "main", logical(1))) == 1L,
    setequal(vapply(Filter(function(x) x$kind == "grid", records), `[[`, numeric(1), "pi"), c(.4, .5)))
})
check("genetic diagnostics are retained once per shared source model", {
  diagnostic <- main[main$parameter == "one_factor_test", , drop = FALSE]
  stopifnot(nrow(diagnostic) == 1L, diagnostic$outcome_model == "[genetic block]",
    all(c("diagnostic_status", "diagnostic_df", "diagnostic_N") %in% names(main)),
    length(attr(main, "diagnostics")) == 1L)
})

father <- run(focal = "father", include_comparison = FALSE)
check("paternal analysis keeps father in all eight subsets and recovers its oracle targets", {
  stopifnot(length(sets(father)) == 8L,
    all(vapply(sets(father), function(x) "father" %in% x, logical(1))))
  for (model in names(sets(father))) {
    selected <- sets(father)[[model]]
    equal(estimates(father, "Corrected for pi (founder scale)", model, selected),
          family_validation_oracle(population, selected)$beta_raw / stats::sd(dat$Y))
  }
})

nested <- run(genetic_configurations = "nested", include_comparison = FALSE)
check("nested genetic configurations have two, four, four and eight maternal subsets", {
  combinations <- attr(nested, "outcome_models")
  equal(vapply(combinations[c("none", "maternal", "paternal", "both")], length, integer(1)),
        c(none = 2L, maternal = 4L, paternal = 4L, both = 8L), tolerance = 0)
  stopifnot(all(vapply(unlist(combinations, recursive = FALSE), function(x)
    "mother" %in% x && !"child" %in% x, logical(1))),
    all(vapply(attr(nested, "fits"), function(x) isTRUE(x$fit$converged), logical(1))))
})

check("invalid focal choices, omitted focal predictors and duplicate subjects are rejected", {
  rejects(run(focal = NULL), "focal")
  rejects(run(focal = "mat_sib"), "focal")
  rejects(run(outcome_models = list(omit_mother = "father")), "focal|mother")
  rejects(run(outcome_models = list(child = c("mother", "child"))), "adult|child")
  rejects(run(outcome_models = list(a = c("mother", "father"), b = c("father", "mother"))), "same subset|repeat")
  rejects(run(rbind(dat, dat[1L, ])), "subject_id|unique")
  bad_models <- list(list(fit_constraints = "estimated_state"), list(fit_constraints = "equilibrium"))
  names(bad_models) <- c("same", "same")
  rejects(run(models = bad_models), "distinct")
})

renamed <- dat
new_pgs <- pgs; new_sibs <- sibs
for (name in c(names(pgs), names(sibs))) {
  old <- c(pgs, sibs)[[name]]
  new <- paste0("custom_", name)
  names(renamed)[names(renamed) == old] <- new
  if (name %in% names(new_pgs)) new_pgs[[name]] <- new else new_sibs[[name]] <- new
}
names(renamed)[names(renamed) == "Y"] <- "maternal_height"
custom <- run(renamed, outcome = "maternal_height", pgs_columns = new_pgs,
  sib_columns = new_sibs, outcome_models = list(own = "mother", pair = c("mother", "mat_sib")),
  include_comparison = FALSE)
check("custom score and outcome names preserve preprocessing and regression units", {
  stopifnot(setequal(names(sets(custom)), c("own", "pair")))
  for (model in names(sets(custom))) {
    selected <- sets(custom)[[model]]
    equal(estimates(custom, "Corrected for pi (founder scale)", model, selected),
          family_validation_oracle(population, selected)$beta_raw / stats::sd(dat$Y))
  }
})

missing <- dat
missing$PGS1_child[1:12] <- NA_real_
missing$Y[13:17] <- NA_real_
missing[18:47, c("PGS1_mat_sib", "PGS2_mat_sib")] <- NA_real_
incomplete <- run(missing, outcome_models = list(own = "mother"), include_comparison = FALSE)
check("original complete-trio filtering and sibling FIML sample accounting are preserved", {
  settings <- attr(incomplete, "settings")$source_settings[[1L]]
  fit <- attr(incomplete, "fits")[[1L]]$fit
  stopifnot(settings$N == nrow(dat) - 17L, settings$N_dropped == 17L,
    settings$dropped[["incomplete_trio"]] == 12L, settings$dropped[["missing_outcome"]] == 5L,
    fit$N == nrow(dat) - 17L, fit$parental_sibs$n_with_sib[["mat_sib"]] == nrow(dat) - 47L,
    fit$parental_sibs$n_with_sib[["pat_sib"]] == nrow(dat) - 17L)
})

finite <- run(models = list("Two generations" = list(fit_constraints = "n_gen", am_n_generations = 2L)),
  outcome_models = list(own = "mother"), pi_grid = c(.35, .45), include_comparison = FALSE)
check("finite-generation grid fits actually use the corresponding pi", {
  records <- attr(finite, "fits")
  grid <- Filter(function(x) x$kind == "grid", records)
  stopifnot(length(grid) == 2L)
  for (record in grid) {
    stopifnot(record$fit$converged, isTRUE(record$fit$constraint_spec$depends_on_pi))
    equal(record$fit$h, record$pi, tolerance = 0)
  }
  direct_dat <- dat
  direct_dat$Y <- as.numeric(scale(dat$Y))
  direct <- RAVEL::ravel_family_fit(direct_dat, genetic_roles = roles, outcome_roles = "mother",
    pi = .35, fit_constraints = "n_gen", am_n_generations = 2L)
  equal(estimates(finite, "Grid (founder scale)", "own", "mother", pi = .35),
        direct$regression$coefficients$total_raw, tolerance = 4e-4)
})

cat(sprintf("\nParental-phenotype workflow: all %d independent check groups passed.\n", checks_passed))
