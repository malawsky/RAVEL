## =============================================================================
## RAVEL: parental phenotypes, with the child retained for genetic identification
## =============================================================================
## Edit SETTINGS and source this file, just as for ravel_parental_deviation.R.
## Keep these scripts, ravel_sem_audit.R and the RAVEL package bundle together. The original script
## supplies the same preprocessing, source fits, pi grid and one-factor check.
##
## The child's scores stay in every genetic model. Reported outcome regressions
## contain adults only and always include the focal parent: eight subsets with
## four adults, four with three adults, or two with the two parents. Set
## focal = "mother" for maternal outcomes or "father" for paternal outcomes.
## These are reduced regressions of the corrected joint moments. The child's
## source outcome association is a nuisance parameter, not fixed to zero.
##
## The one-factor diagnostic concerns the FULL source model, once per main fit.
## It does not test whether the child's outcome coefficient is zero, or validate
## every reduced regression. Full diagnostics and fits are saved in the RDS file.
## Separate SEM audit CSVs contain fit statistics, all parameter estimates,
## variances/covariances with SEs, and the parameter uncertainty covariance.
## To load functions without running SETTINGS: set
## options(ravel.workflow_functions_only = TRUE) before sourcing, then reset it.
## =============================================================================

.ravel_parental_script_dir <- local({
  files <- vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else x$ofile, "")
  files <- files[nzchar(files)]
  cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  path <- if (length(files)) tail(files, 1L) else if (length(cli)) sub("^--file=", "", cli[1L]) else "."
  if (identical(path, ".")) getwd() else dirname(normalizePath(path, mustWork = TRUE))
})
.ravel_parental_engine <- local({
  engine <- NULL
  function() {
    if (is.null(engine)) {
      path <- file.path(.ravel_parental_script_dir, "ravel_parental_deviation.R")
      if (!file.exists(path)) stop("Keep ravel_parental_deviation.R beside this script.")
      old <- options(ravel.workflow_functions_only = TRUE)
      on.exit(options(old), add = TRUE)
      cwd <- getwd(); on.exit(setwd(cwd), add = TRUE)
      setwd(.ravel_parental_script_dir)
      candidate <- new.env(parent = globalenv())
      sys.source(path, envir = candidate)
      engine <<- candidate
    }
    engine
  }
})

## ---- SETTINGS: edit these ---------------------------------------------------
if (!isTRUE(getOption("ravel.workflow_functions_only", FALSE))) {
  data_file <- NULL                 # .csv, .tsv/.txt or .rds; NULL = synthetic adult example
  design <- "cross_sectional"      # or "longitudinal" for repeated adult measurements
  outcome <- "mom_height"          # one or more parental outcome columns
  binary_outcome <- FALSE          # TRUE for all outcomes, or names of binary outcomes
  drop_missing_outcome <- TRUE

  pgs_columns <- c(
    child1 = "PGS1_child", child2 = "PGS2_child",
    mother1 = "PGS1_mother", mother2 = "PGS2_mother",
    father1 = "PGS1_father", father2 = "PGS2_father")
  sib_columns <- c(mat_sib1 = "PGS1_mat_sib", mat_sib2 = "PGS2_mat_sib",
                   pat_sib1 = "PGS1_pat_sib", pat_sib2 = "PGS2_pat_sib")
  # Set sib_columns = NULL without parental siblings, or supply only one pair.

  genetic_configurations <- "available"  # all supplied scores inform every subset
  # "nested" also runs trio-only, + maternal sibling, + paternal sibling, + both,
  # for the configurations your score columns support, as in the original script.
  focal <- "mother"                # required in every subset; use "father" for paternal outcomes
  outcome_models <- NULL           # optional named list, e.g.
  # outcome_models <- list(own = "mother", sib = c("mother", "mat_sib"),
  #                       partners = c("mother", "father"))

  covariates <- NULL               # numerical covariates, e.g. c("sex", "PC1", "PC2")
  covariate_mapping_confirmed <- FALSE
  # Corrected reduced regressions with covariates/interactions additionally assume
  # Cov(T - M'G, W) = 0. Set TRUE only when using this assumption (family-models.md).
  # Otherwise visible estimates remain available and corrected rows explain why not.
  family_id <- "FID"              # true family ID; NULL if rows are independent
  subject_id <- "FID"             # cross-sectional: unique focal-adult ID, or NULL
  # Choose one row per adult: clustering repeated copies once per child does not
  # remove duplicate-adult weighting. Longitudinal measurements may repeat adults.
  age <- "age"                   # longitudinal adult age column
  age_center <- NULL              # NULL uses its mean

  pi <- 0.528
  pi_se <- 0.021
  pi_grid <- seq(0.40, 0.60, by = 0.05)
  grid_pi_se <- NULL              # each sensitivity value is treated as fixed
  models <- list(
    "Estimated state" = list(fit_constraints = "estimated_state"),
    "3 generations" = list(fit_constraints = "n_gen", am_n_generations = 3L),
    "Equilibrium" = list(fit_constraints = "equilibrium"))
  include_raw <- FALSE
  include_sem_audit <- TRUE        # extra detailed SEM tables; primary results stay unchanged
  output_csv <- "ravel_parental_phenotypes.csv"
}

## ---- ANALYSIS FUNCTIONS (no need to edit) -----------------------------------

.ravel_adult_outcome_models <- function(adults, focal = "mother", outcome_models = NULL) {
  if (is.null(outcome_models)) {
    ans <- list()
    for (mask in seq_len(2^length(adults) - 1L)) {
      roles <- adults[as.logical(intToBits(mask)[seq_along(adults)])]
      if (focal %in% roles) ans[[paste(roles, collapse = "+")]] <- roles
    }
    return(ans)
  }
  Filter(function(x) all(x %in% adults) && focal %in% x, outcome_models)
}

.ravel_adult_parameter <- function(x) {
  map <- c(mother = "beta_M", father = "beta_F", mat_sib = "beta_MS", pat_sib = "beta_FS",
           delta_m = "delta_M", delta_f = "delta_F")
  base <- sub(":.*$", "", x)
  hit <- base %in% names(map)
  x[hit] <- paste0(unname(map[base[hit]]), substring(x[hit], nchar(base[hit]) + 1L))
  x
}
.ravel_adult_effect <- function(x) {
  labels <- c(beta_M = "Maternal (beta_M)", beta_F = "Paternal (beta_F)",
    beta_MS = "Maternal sibling (beta_MS)", beta_FS = "Paternal sibling (beta_FS)",
    beta_mid = "Parental midpoint (beta_M + beta_F)",
    delta_M_mid = "Mother - father difference (beta_M - beta_F)",
    delta_P_mid = "Father - mother difference (beta_F - beta_M)",
    delta_M = "Mother - maternal sibling (delta_M)",
    delta_F = "Father - paternal sibling (delta_F)",
    rho_visible = "Spouse correlation", rho_corrected = "Spouse correlation",
    genetic_variance = "Variance of the genetic prediction", genetic_R2 = "Genetic R-squared",
    conditional_genetic_variance = "Conditional variance of the genetic prediction",
    partial_R2 = "Partial genetic R-squared")
  base <- sub(":.*$", "", x)
  ans <- unname(labels[base]); ans[is.na(ans)] <- base[is.na(ans)]
  suffix <- substring(x, nchar(base) + 1L)
  paste0(ans, ifelse(nzchar(suffix), paste0(" x ", substring(suffix, 2L)), ""))
}

# Add midpoint/deviation contrasts using the FULL covariance, including external
# pi uncertainty. Sibling contrasts are already provided by the family API.
.ravel_adult_projection_table <- function(projection) {
  tab <- projection$table
  extra <- list()
  for (sc in unique(tab$scale)) {
    params <- tab$parameter[tab$scale == sc & tab$kind == "coefficient"]
    mods <- unique(ifelse(grepl(":", params, fixed = TRUE), sub("^[^:]+", "", params), ""))
    for (mod in mods) {
      pair <- paste0(c("mother", "father"), mod)
      if (!all(pair %in% params)) next
      rows <- match(paste0(sc, ":", pair), paste0(tab$scale, ":", tab$parameter))
      b <- tab$estimate[rows]
      V <- projection$vcov[paste0(sc, ":", pair), paste0(sc, ":", pair), drop = FALSE]
      W <- rbind(beta_mid = c(1, 1), delta_M_mid = c(1, -1), delta_P_mid = c(-1, 1))
      est <- drop(W %*% b)
      se <- sqrt(pmax(diag(W %*% V %*% t(W)), 0))
      z <- stats::qnorm(1 - (1 - projection$ci_level) / 2)
      out <- tab[rep(rows[1L], 3L), , drop = FALSE]
      out$parameter <- paste0(rownames(W), mod); out$estimate <- est; out$se <- se
      out$ci_low <- est - z * se; out$ci_high <- est + z * se
      out$p_value <- ifelse(is.finite(se) & se > 0, 2 * stats::pnorm(-abs(est / se)), NA_real_)
      out$kind <- "contrast"
      out$inference_method <- ifelse(is.finite(se), "joint_covariance_wald", "unavailable")
      out$inference_status <- ifelse(is.finite(se), "ok", "unavailable")
      extra[[length(extra) + 1L]] <- out
    }
  }
  if (length(extra)) tab <- do.call(rbind, c(list(tab), extra))
  tab$parameter <- .ravel_adult_parameter(tab$parameter)
  rownames(tab) <- NULL
  tab
}

# Refit each plain PGS comparator on precisely the requested adult predictors.
# Selecting rows from the full trio regression would target a different estimand.
.ravel_adult_plain <- function(dat, outcome, roles, pgs_columns, covariates,
                               binary = FALSE, version = 1L) {
  d <- data.frame(.Y = if (binary) dat[[outcome]] else as.numeric(scale(dat[[outcome]])))
  predictors <- paste0(".PGS_", roles)
  for (i in seq_along(roles))
    d[[predictors[i]]] <- as.numeric(scale(dat[[pgs_columns[[paste0(roles[i], version)]]]]))
  for (nm in covariates) d[[nm]] <- dat[[nm]]
  formula <- stats::reformulate(c(predictors, covariates), response = ".Y")
  fit <- if (binary) stats::glm(formula, data = d, family = stats::binomial("probit")) else
    stats::lm(formula, data = d)
  b <- stats::coef(fit); V <- stats::vcov(fit)
  W <- matrix(0, length(roles), length(b), dimnames = list(roles, names(b)))
  for (i in seq_along(roles)) W[i, predictors[i]] <- 1
  pairs <- list(beta_mid = c("mother", "father"), delta_M_mid = c("mother", "father"),
    delta_P_mid = c("father", "mother"), delta_m = c("mother", "mat_sib"),
    delta_f = c("father", "pat_sib"))
  for (nm in names(pairs)) if (all(pairs[[nm]] %in% roles)) {
    w <- stats::setNames(rep(0, length(b)), names(b))
    w[paste0(".PGS_", pairs[[nm]])] <- if (nm == "beta_mid") c(1, 1) else c(1, -1)
    W <- rbind(W, w); rownames(W)[nrow(W)] <- nm
  }
  est <- vapply(seq_len(nrow(W)), function(i) {
    use <- which(W[i, ] != 0)
    sum(W[i, use] * b[use])
  }, numeric(1))
  se <- vapply(seq_len(nrow(W)), function(i) {
    use <- which(W[i, ] != 0)
    sqrt(max(drop(crossprod(W[i, use], V[use, use, drop = FALSE] %*% W[i, use])), 0))
  }, numeric(1))
  statistic <- est / se
  p <- if (binary) 2 * stats::pnorm(-abs(statistic)) else 2 * stats::pt(-abs(statistic), fit$df.residual)
  critical <- if (binary) stats::qnorm(.975) else stats::qt(.975, fit$df.residual)
  converged <- !binary || isTRUE(fit$converged)
  status <- ifelse(is.finite(est) & converged, "ok", if (converged) "aliased_coefficient" else "not_converged")
  data.frame(parameter = .ravel_adult_parameter(rownames(W)), estimate = est, se = se,
    ci_low = est - critical * se, ci_high = est + critical * se,
    p_value = p, kind = c(rep("coefficient", length(roles)), rep("contrast", nrow(W) - length(roles))),
    scale = "vanilla", status = status, corrected_status = "not_applicable",
    inference_method = if (binary) "glm_probit" else "ordinary_lm",
    inference_status = ifelse(is.finite(se) & converged, "ok", "unavailable"),
    N = stats::nobs(fit), stringsAsFactors = FALSE, row.names = NULL)
}

ravel_parental_phenotypes <- function(
    dat, outcome = "Y", pgs_columns, sib_columns = NULL,
    genetic_configurations = c("available", "nested"), focal = "mother", outcome_models = NULL,
    models = list("Estimated state" = list(fit_constraints = "estimated_state"),
      "3 generations" = list(fit_constraints = "n_gen", am_n_generations = 3L),
      "Equilibrium" = list(fit_constraints = "equilibrium")),
    pi = .4, pi_se = NULL, pi_grid = NULL, grid_pi_se = NULL,
    covariates = NULL, covariate_mapping_confirmed = FALSE, family_id = "FID", subject_id = NULL,
    design = c("cross_sectional", "longitudinal"), age = "age", age_center = NULL,
    binary_outcome = FALSE, include_raw = FALSE, include_comparison = TRUE,
    drop_missing_outcome = TRUE, include_sem_audit = TRUE) {
  design <- match.arg(design); genetic_configurations <- match.arg(genetic_configurations)
  if (!is.data.frame(dat)) stop("dat must be a data.frame.")
  if (!is.character(outcome) || length(outcome) != 1L || is.na(outcome))
    stop("Supply one outcome column per function call; the script loops over outcome columns.")
  if (!is.logical(binary_outcome) || length(binary_outcome) != 1L || is.na(binary_outcome))
    stop("binary_outcome must be TRUE or FALSE for this outcome.")
  trio_keys <- c("child1", "child2", "mother1", "mother2", "father1", "father2")
  if (is.null(names(pgs_columns)) || anyDuplicated(names(pgs_columns)) ||
      !all(trio_keys %in% names(pgs_columns))) stop("pgs_columns must map all six trio scores.")
  pgs_columns <- pgs_columns[trio_keys]
  sib_roles <- character(0)
  if (length(sib_columns)) {
    if (is.null(names(sib_columns)) || anyDuplicated(names(sib_columns)) ||
        any(!names(sib_columns) %in% c("mat_sib1", "mat_sib2", "pat_sib1", "pat_sib2")))
      stop("sib_columns must map mat_sib1/mat_sib2 and/or pat_sib1/pat_sib2.")
    for (r in c("mat_sib", "pat_sib")) {
      present <- paste0(r, 1:2) %in% names(sib_columns)
      if (any(present) && !all(present)) stop("Supply both scores for ", r, ".")
      if (all(present)) sib_roles <- c(sib_roles, r)
    }
  }
  all_pgs <- c(pgs_columns, sib_columns)
  if (!is.character(all_pgs) || anyNA(all_pgs) || any(!nzchar(all_pgs)) || anyDuplicated(all_pgs))
    stop("PGS mappings must name distinct, nonmissing columns.")
  absent <- setdiff(c(unname(all_pgs), outcome, covariates, subject_id), names(dat))
  if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))
  numeric_cols <- unique(c(unname(all_pgs), outcome, covariates))
  if (any(!vapply(dat[numeric_cols], is.numeric, logical(1))))
    stop("Scores, outcome and covariates must be numeric.")
  if (any(vapply(dat[numeric_cols], function(x) any(!is.na(x) & !is.finite(x)), logical(1))))
    stop("Scores, outcome or covariates contain infinite values.")
  if (!is.null(subject_id)) {
    if (!is.character(subject_id) || length(subject_id) != 1L) stop("subject_id must name one column.")
    if (anyNA(dat[[subject_id]])) stop("subject_id must be observed.")
    if (design == "cross_sectional" && anyDuplicated(dat[[subject_id]]))
      stop("subject_id must be unique: choose one record per focal adult before a cross-sectional analysis.")
  }
  adults <- c("mother", "father", sib_roles)
  if (!is.character(focal) || length(focal) != 1L || is.na(focal) || !focal %in% c("mother", "father"))
    stop("focal must be mother or father; every outcome regression includes that parent.")
  if (!is.null(outcome_models)) {
    if (!is.list(outcome_models) || !length(outcome_models) ||
        any(vapply(outcome_models, function(x) !is.character(x) || !length(x) || anyNA(x) ||
          anyDuplicated(x) || !all(x %in% adults) || !focal %in% x, logical(1))))
      stop("outcome_models must contain adult-role subsets including the focal parent; child is never an outcome predictor.")
    if (is.null(names(outcome_models)))
      names(outcome_models) <- vapply(outcome_models, paste, collapse = "+", character(1))
    if (anyNA(names(outcome_models)) || any(!nzchar(names(outcome_models))) || anyDuplicated(names(outcome_models)))
      stop("outcome_models needs distinct, nonempty names.")
    canonical <- vapply(outcome_models, function(x) paste(sort(x), collapse = "+"), "")
    if (anyDuplicated(canonical)) stop("outcome_models repeats the same subset.")
  }
  if (!is.list(models) || !length(models) || is.null(names(models)) || anyNA(names(models)) ||
      any(!nzchar(names(models))) || anyDuplicated(names(models)))
    stop("models must have distinct, nonempty names.")
  allowed_model_args <- c("fit_constraints", "am_n_generations")
  if (any(vapply(models, function(x) !is.list(x) || is.null(names(x)) || anyNA(names(x)) ||
    anyDuplicated(names(x)) || any(!names(x) %in% allowed_model_args), logical(1))))
    stop("Each model supplies fit_constraints and optionally am_n_generations.")
  engine <- .ravel_parental_engine()
  available <- if (length(sib_roles) == 2L) "both" else if (identical(sib_roles, "mat_sib"))
    "maternal" else if (identical(sib_roles, "pat_sib")) "paternal" else "none"
  configs <- if (genetic_configurations == "available") available else
    c("none", if ("mat_sib" %in% sib_roles) "maternal", if ("pat_sib" %in% sib_roles) "paternal",
      if (length(sib_roles) == 2L) "both")
  rows <- fits <- diagnostics <- projections <- source_tables <- settings <- subset_records <- list()
  add <- function(tab, ps, model_name, analysis, p = NA_real_, p_se = NA_real_,
                  subset = "[genetic block]", selected = character(), genetic, N, note = "") {
    if (is.null(tab) || !nrow(tab)) return(invisible(NULL))
    value <- function(nm, default) if (nm %in% names(tab)) tab[[nm]] else rep_len(default, nrow(tab))
    row_note <- value("note", "")
    row_note[is.na(row_note)] <- ""
    combined <- vapply(row_note, function(x) paste(c(note[nzchar(note)], x[nzchar(x)]), collapse = "; "), "")
    status <- value("status", ifelse(is.finite(tab$estimate), "ok", "unavailable"))
    if ("diagnostic_status" %in% names(tab)) {
      diagnostic <- !is.na(tab$diagnostic_status) & nzchar(tab$diagnostic_status)
      status[diagnostic] <- tab$diagnostic_status[diagnostic]
    }
    tab <- data.frame(outcome = outcome, siblings = c(none = "Trio only", maternal = "+ maternal sibling",
      paternal = "+ paternal sibling", both = "+ both siblings")[[ps]], model = model_name,
      outcome_model = subset, outcome_predictors = paste(selected, collapse = "+"),
      genetic_members = paste(genetic, collapse = "+"), analysis = analysis,
      pi = p, pi_se = if (is.null(p_se)) NA_real_ else p_se,
      effect = value("effect", .ravel_adult_effect(tab$parameter)), parameter = tab$parameter,
      estimate = tab$estimate, se = value("se", NA_real_), ci_low = value("ci_low", NA_real_),
      ci_high = value("ci_high", NA_real_), p_value = value("p_value", NA_real_), note = combined,
      scale = value("scale", NA_character_), kind = value("kind", "diagnostic"),
      status = status, corrected_status = value("corrected_status", "not_applicable"),
      inference_method = value("inference_method", ""), inference_status = value("inference_status", ""),
      N = value("N", N), diagnostic_status = value("diagnostic_status", NA_character_),
      diagnostic_df = value("diagnostic_df", NA_real_), diagnostic_N = value("diagnostic_N", NA_real_),
      stringsAsFactors = FALSE, row.names = NULL)
    rows[[length(rows) + 1L]] <<- tab
  }
  failure <- function(message) data.frame(parameter = NA_character_, estimate = NA_real_,
    status = "failed", note = message, stringsAsFactors = FALSE)
  for (ps in configs) {
    siblings <- switch(ps, none = character(), maternal = "mat_sib", paternal = "pat_sib",
                       both = c("mat_sib", "pat_sib"))
    genetic <- c("child", "mother", "father", siblings)
    subsets <- .ravel_adult_outcome_models(setdiff(genetic, "child"), focal, outcome_models)
    subset_records[[ps]] <- subsets
    if (!length(subsets)) {
      message("No requested adult subset is supported by genetic configuration '", ps, "'; skipped.")
      next
    }
    comparator_done <- FALSE
    for (model_name in names(models)) {
      source <- tryCatch(do.call(engine$ravel_parental_deviation,
        c(list(dat = dat, design = design, outcome = outcome, pgs_columns = pgs_columns,
          covariates = covariates, binary_outcome = binary_outcome, family_id = family_id,
          age = age, age_center = age_center, pi = pi, pi_se = pi_se, pi_grid = pi_grid,
          grid_pi_se = grid_pi_se, model = model_name, include_comparison = FALSE,
          include_raw = include_raw, drop_missing_outcome = drop_missing_outcome,
          parental_sibs = ps, sib_columns = sib_columns), models[[model_name]])), error = identity)
      source_key <- paste(ps, model_name, sep = "|")
      if (inherits(source, "error")) {
        for (nm in names(subsets)) add(failure(conditionMessage(source)), ps, model_name,
          "Source fit failed", pi, pi_se, nm, subsets[[nm]], genetic, NA_integer_)
        next
      }
      source_tables[[source_key]] <- source
      st <- attr(source, "settings"); settings[[source_key]] <- st
      diagnostics[[source_key]] <- attr(source, "one_factor")
      records <- attr(source, "fit_records")
      if (is.null(records)) stop("The companion ravel_parental_deviation.R is outdated: fit_records are missing.")
      # Source-level correlations and one-factor check are shown only once, not
      # repeated as though they independently validate each adult regression.
      keep <- source$parameter %in% c("rho_visible", "rho_corrected", "one_factor_test")
      for (i in which(keep)) {
        row <- source[i, , drop = FALSE]
        scope_note <- if (row$parameter == "one_factor_test")
          "Full joint source model; not a test of excluding child liability" else "Shared genetic model"
        add(row, ps, model_name, row$analysis, row$pi, row$pi_se, genetic = genetic, N = st$N,
          note = scope_note)
      }
      for (j in seq_along(records)) {
        record <- records[[j]]
        key <- paste(source_key, record$kind, format(record$pi, digits = 15), j, sep = "|")
        fit <- record$fit
        # Cache the source adapter once: its pi-sensitivity evaluator is shared
        # across subsets, so finite-history uncertainty does not repeat its refits.
        if (!is.null(fit) && isTRUE(fit$converged) && is.null(fit$family_descriptor)) {
          descriptor <- tryCatch(getFromNamespace(".ravel_family_legacy_descriptor", "RAVEL")(fit),
            error = identity)
          if (!inherits(descriptor, "error")) fit$family_descriptor <- descriptor
        }
        record$fit <- fit
        record$model <- model_name
        record$parental_sibs <- ps
        record$genetic_roles <- genetic
        record$outcome <- outcome
        fits[[key]] <- record
        projections[[key]] <- list()
        for (nm in names(subsets)) {
          selected <- subsets[[nm]]
          if (is.null(fit) || !isTRUE(fit$converged)) {
            add(failure(record$note), ps, model_name, if (record$kind == "main") "Corrected for pi" else "Grid",
                record$pi, record$pi_se, nm, selected, genetic, st$N)
            next
          }
          projection_warnings <- character()
          projected <- tryCatch(withCallingHandlers(RAVEL::ravel_family_regression(fit, selected,
            pi = record$pi, pi_se = record$pi_se,
            covariate_mapping_confirmed = covariate_mapping_confirmed), warning = function(w) {
              projection_warnings <<- c(projection_warnings, conditionMessage(w))
              invokeRestart("muffleWarning")
            }), error = identity)
          if (inherits(projected, "error")) {
            add(failure(conditionMessage(projected)), ps, model_name,
              if (record$kind == "main") "Corrected for pi" else "Grid", record$pi, record$pi_se,
              nm, selected, genetic, st$N)
            next
          }
          projections[[key]][[nm]] <- projected
          tab <- tryCatch(.ravel_adult_projection_table(projected), error = identity)
          if (inherits(tab, "error")) {
            add(failure(conditionMessage(tab)), ps, model_name, "Projection table failed",
              record$pi, record$pi_se, nm, selected, genetic, st$N)
            next
          }
          scales <- c(if (record$kind == "main") "visible", "total_std", if (include_raw) "total_raw")
          for (sc in scales) {
            selected_rows <- tab[tab$scale == sc, , drop = FALSE]
            label <- if (sc == "visible") "Visible (measured by the PGSs)" else
              paste0(if (record$kind == "main") "Corrected for pi" else "Grid",
                if (sc == "total_raw") " (founder scale)")
            note <- c(record$note, projection_warnings, projected$uncertainty_note, projected$uncertainty_errors,
              if (binary_outcome) "Reduced regression on the fitted probit liability; not a refitted marginal probit model",
              if (sc != "visible" && !identical(projected$correction_status, "ok")) projected$correction_status)
            note <- unique(note[!is.na(note) & nzchar(note)])
            add(selected_rows, ps, model_name, label,
              if (sc == "visible") NA_real_ else record$pi,
              if (sc == "visible") NA_real_ else record$pi_se, nm, selected, genetic, st$N,
              paste(note, collapse = "; "))
          }
        }
      }
      if (include_comparison && design == "cross_sectional" && !comparator_done) {
        prepared <- attr(source, "analysis_data")
        if (is.null(prepared)) stop("The companion ravel_parental_deviation.R is outdated: analysis_data are missing.")
        for (nm in names(subsets)) for (version in 1:2) {
          msgs <- character()
          plain <- tryCatch(withCallingHandlers(.ravel_adult_plain(prepared, outcome, subsets[[nm]],
            all_pgs, st$covariates, binary_outcome, version), warning = function(w) {
              msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
            }), error = identity)
          if (inherits(plain, "error")) plain <- failure(conditionMessage(plain))
          note <- c(if (!is.null(st$cluster_by)) "Ordinary comparison SEs treat rows as independent; RAVEL uses clustered SEs",
                    if (binary_outcome) "Conditional probit coefficients; SEM reduced-liability coefficients use a different normalization", msgs)
          add(plain, ps, "Plain PGS regressions (comparison)", paste0("Plain PGS", version, " regression"),
            subset = nm, selected = subsets[[nm]], genetic = genetic, N = st$N, note = paste(note, collapse = "; "))
        }
        comparator_done <- TRUE
      }
    }
  }
  if (!length(rows)) stop("No requested adult regression could be analysed.")
  table <- do.call(rbind, rows); rownames(table) <- NULL
  attr(table, "fits") <- fits
  attr(table, "diagnostics") <- diagnostics
  attr(table, "projections") <- projections
  attr(table, "settings") <- list(outcome = outcome, design = design,
    genetic_configurations = genetic_configurations, focal = focal, source_settings = settings,
    covariate_mapping_confirmed = covariate_mapping_confirmed, subject_id = subject_id,
    note = "Child scores inform the joint genetic model; all reported outcome predictors are adults.")
  attr(table, "outcome_models") <- subset_records
  attr(table, "source_tables") <- source_tables
  if (isTRUE(include_sem_audit) && length(fits)) {
    audit_path <- file.path(.ravel_parental_script_dir, "ravel_sem_audit.R")
    if (!file.exists(audit_path)) stop("Keep ravel_sem_audit.R beside this script to produce SEM audit tables.")
    audit_env <- new.env(parent = baseenv())
    sys.source(audit_path, envir = audit_env)
    attr(table, "sem_audit") <- audit_env$ravel_sem_audit(table)
  }
  table
}

print_parental_phenotypes <- function(table, digits = 3L) {
  f <- function(x) if (!length(x) || is.na(x)) "NA" else formatC(x, digits = digits, format = "f")
  fp <- function(x) if (!length(x) || is.na(x)) "NA" else if (x < .001) "<0.001" else
    formatC(x, digits = 3, format = "f")
  same <- function(x, y) if (is.na(y)) is.na(x) else !is.na(x) & abs(x - y) < 1e-12
  cat("\nRAVEL parental phenotype: ", attr(table, "settings")$outcome,
      "\nChild scores remain in the genetic block. Cells are estimate (SE); contrasts also show p.\n", sep = "")
  diagnostics <- table[table$parameter %in% "one_factor_test", , drop = FALSE]
  if (nrow(diagnostics)) {
    cat("\nOne-factor checks of the shared source models:\n")
    for (i in seq_len(nrow(diagnostics))) {
      d <- diagnostics[i, ]
      cat("  ", d$siblings, ", ", d$model, ": ",
        if (is.finite(d$estimate)) paste0("chi-square = ", f(d$estimate), ", df = ", d$diagnostic_df,
          ", p = ", fp(d$p_value)) else d$status,
        if (nzchar(d$note)) paste0("; ", d$note), "\n", sep = "")
    }
  }
  keys <- unique(table[table$outcome_model != "[genetic block]", c("siblings", "outcome_model")])
  short <- c(beta_M = "Mother", beta_F = "Father", beta_MS = "Mat sib", beta_FS = "Pat sib",
    delta_M_mid = "M - F [p]", delta_M = "M - mat sib [p]", delta_F = "F - pat sib [p]")
  for (i in seq_len(nrow(keys))) {
    sub <- table[table$siblings == keys$siblings[i] & table$outcome_model == keys$outcome_model[i], ]
    cat("\n", keys$siblings[i], ": Y ~ ", sub$outcome_predictors[1L],
      if (keys$outcome_model[i] != sub$outcome_predictors[1L]) paste0(" [", keys$outcome_model[i], "]"),
      "\n", sep = "")
    pars <- intersect(names(short), sub$parameter)
    width <- ifelse(grepl("^delta", pars), 27L, 17L)
    header <- paste0(sprintf("  %-29s %-29s %6s %6s %6s", "Model", "Analysis", "pi", "pi SE", "N"),
      paste(vapply(seq_along(pars), function(j) sprintf(paste0(" %", width[j], "s"), short[pars[j]]), ""), collapse = ""))
    cat(header, "\n", sep = "")
    groups <- unique(sub[, c("model", "analysis", "pi", "pi_se")])
    for (k in seq_len(nrow(groups))) {
      key <- groups[k, ]
      r <- sub[sub$model == key$model & sub$analysis == key$analysis & same(sub$pi, key$pi) &
        same(sub$pi_se, key$pi_se), , drop = FALSE]
      cell <- function(par) {
        row <- r[r$parameter %in% par, , drop = FALSE]
        if (!nrow(row)) return("")
        paste0(f(row$estimate[1]), " (", f(row$se[1]), ")",
          if (grepl("^delta", par)) paste0(" [", fp(row$p_value[1]), "]"))
      }
      n <- unique(r$N[!is.na(r$N)])
      cat(sprintf("  %-29s %-29s %6s %6s %6s", substr(key$model, 1L, 29L), substr(key$analysis, 1L, 29L),
        if (is.na(key$pi)) "" else f(key$pi), if (is.na(key$pi_se)) "" else f(key$pi_se),
        if (length(n) == 1L) as.character(n) else ""),
        paste(vapply(seq_along(pars), function(j) sprintf(paste0(" %", width[j], "s"), cell(pars[j])), ""), collapse = ""),
        "\n", sep = "")
      slopes <- r$parameter[grepl(":", r$parameter, fixed = TRUE)]
      if (length(slopes)) cat("      slopes: ", paste0(slopes, " ", vapply(slopes, cell, ""), collapse = ", "), "\n", sep = "")
      if (all(is.na(r$estimate))) cat("      Unavailable: ", paste(unique(r$status), collapse = "; "), " (see notes).\n", sep = "")
    }
    notes <- unique(sub$note[!is.na(sub$note) & nzchar(sub$note)])
    if (length(notes)) cat(paste0("  Note: ", notes, "\n"), sep = "")
  }
  cat("\nCSV tables also contain parental midpoint, reverse contrasts, variance summaries and source spouse correlations.\n")
  invisible(table)
}

## ---- RUN -------------------------------------------------------------------
if (!isTRUE(getOption("ravel.workflow_functions_only", FALSE))) {
  engine <- .ravel_parental_engine()
  if (is.null(data_file)) {
    example <- file.path(.ravel_parental_script_dir, "synthetic_family_example.rds")
    if (!file.exists(example)) stop("Synthetic example not found: ", example)
    dat <- readRDS(example)
    if (!"mom_height" %in% names(dat) && "Y" %in% names(dat)) dat$mom_height <- dat$Y
    message("No data_file set: using synthetic adult outcomes, not real family observations.")
    design <- "cross_sectional"
  } else dat <- engine$read_trio_data(data_file)
  if (anyDuplicated(outcome)) stop("Each outcome should be listed once.")
  results <- list()
  for (y in outcome) {
    is_binary <- if (is.logical(binary_outcome)) isTRUE(binary_outcome) else y %in% binary_outcome
    res <- tryCatch(ravel_parental_phenotypes(dat, y, pgs_columns, sib_columns,
      genetic_configurations = genetic_configurations, focal = focal, outcome_models = outcome_models,
      models = models, pi = pi, pi_se = pi_se, pi_grid = pi_grid, grid_pi_se = grid_pi_se,
      covariates = covariates, covariate_mapping_confirmed = covariate_mapping_confirmed,
      family_id = family_id, subject_id = subject_id, design = design, age = age, age_center = age_center,
      binary_outcome = is_binary, include_raw = include_raw, drop_missing_outcome = drop_missing_outcome,
      include_sem_audit = include_sem_audit), error = identity)
    if (inherits(res, "error")) {
      message("Outcome ", y, " could not be analysed: ", conditionMessage(res))
      next
    }
    results[[y]] <- res
    print_parental_phenotypes(res)
    base <- sub("[.]csv$", "", output_csv)
    path <- paste0(base, "_", gsub("[^A-Za-z0-9_.-]+", "_", y), ".csv")
    utils::write.csv(res, path, row.names = FALSE)
    engine$write_excel_csv(res, sub("[.]csv$", "_excel.csv", path))
    audit <- attr(res, "sem_audit")
    if (length(audit)) for (nm in names(audit))
      utils::write.csv(audit[[nm]], sub("[.]csv$", paste0("_", nm, ".csv"), path), row.names = FALSE)
    saveRDS(res, sub("[.]csv$", "_fits.rds", path))
    message("Saved table, Excel copy, ", length(audit), " SEM audit tables and fitted objects/diagnostics for ",
      y, " at ", normalizePath(path))
  }
  message("Done. Analysed: ", if (length(results)) paste(names(results), collapse = ", ") else "none")
}
