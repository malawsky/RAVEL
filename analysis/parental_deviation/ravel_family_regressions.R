# Source this file to define the analysis function; sourcing never runs an analysis.
# Requires RAVEL >= 0.2.0.9001. See README.md in this directory for an example.

ravel_family_regressions <- function(
    dat, outcomes, genetic_roles, focal = "mother", pgs_columns = NULL,
    outcome_models = NULL, covariates = NULL, cluster = NULL,
    subject_id = NULL, pi = .4, pi_se = NULL, pi_grid = NULL,
    models = list(
      estimated_state = list(fit_constraints = "estimated_state"),
      n_gen_3 = list(fit_constraints = "n_gen", am_n_generations = 3L),
      equilibrium = list(fit_constraints = "equilibrium")),
    binary_outcomes = character(0), standardize_outcome = FALSE,
    complete_cases = FALSE, covariate_mapping_confirmed = FALSE,
    output_prefix = NULL, ...) {
  if (!requireNamespace("RAVEL", quietly = TRUE) ||
      !all(c("ravel_family_fit", "ravel_family_regression") %in%
           getNamespaceExports("RAVEL")))
    stop("Install the bundled RAVEL 0.2.0.9001 or newer before running this function.")
  if (!is.data.frame(dat)) stop("dat must be a data.frame.")
  allowed <- c("child", "mother", "father", "mat_sib", "pat_sib")
  if (!length(genetic_roles) || anyDuplicated(genetic_roles) ||
      any(!genetic_roles %in% allowed)) stop("Invalid or duplicated genetic_roles.")
  if (length(focal) != 1L || !focal %in% genetic_roles || focal == "child")
    stop("focal must name one adult in genetic_roles.")
  if (!is.character(outcomes) || !length(outcomes) || anyNA(outcomes) || anyDuplicated(outcomes))
    stop("outcomes must contain unique outcome column names.")
  if (!is.list(models) || !length(models) || is.null(names(models)) ||
      anyNA(names(models)) || anyDuplicated(names(models)) || any(!nzchar(names(models))))
    stop("models must be a list with distinct, nonempty names for its fitting settings.")
  if (is.null(pgs_columns)) {
    keys <- unlist(lapply(genetic_roles, function(r) paste0(r, 1:2)))
    pgs_columns <- stats::setNames(unlist(lapply(genetic_roles,
      function(r) paste0("PGS", 1:2, "_", r))), keys)
  }
  keys <- unlist(lapply(genetic_roles, function(r) paste0(r, 1:2)))
  if (is.null(names(pgs_columns)) || anyDuplicated(names(pgs_columns)) ||
      any(!keys %in% names(pgs_columns)))
    stop("pgs_columns must map each selected role's two scores, e.g. mother1 and mother2.")
  score_cols <- unname(pgs_columns[keys])
  if (anyNA(score_cols) || anyDuplicated(score_cols)) stop("PGS columns must be distinct.")
  needed <- unique(c(score_cols, outcomes, covariates, cluster, subject_id))
  absent <- setdiff(needed, names(dat))
  if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))
  numeric_cols <- unique(c(score_cols, outcomes, covariates))
  if (any(!vapply(dat[numeric_cols], is.numeric, logical(1))))
    stop("Scores, outcomes and covariates must be numeric.")
  if (any(vapply(dat[numeric_cols], function(x) any(!is.na(x) & !is.finite(x)), logical(1))))
    stop("Numeric analysis columns contain infinite values.")
  if (!is.null(subject_id) && (anyNA(dat[[subject_id]]) || anyDuplicated(dat[[subject_id]])))
    stop("subject_id must be observed and unique. Choose one record per focal adult before ",
         "this cross-sectional analysis; clustering alone does not remove duplicate-adult weighting.")
  if (!is.null(cluster) && anyNA(dat[[cluster]])) stop("Cluster IDs must be observed.")
  if (is.null(outcome_models)) {
    optional <- setdiff(genetic_roles, c("child", focal))
    outcome_models <- list()
    for (mask in 0:(2^length(optional) - 1L)) {
      take <- if (length(optional)) as.logical(intToBits(mask)[seq_along(optional)]) else logical(0)
      roles <- c(focal, optional[take])
      outcome_models[[paste(roles, collapse = "+")]] <- roles
    }
  }
  if (!is.list(outcome_models) || !length(outcome_models) ||
      any(vapply(outcome_models, function(r) !length(r) || anyDuplicated(r) ||
        !all(r %in% genetic_roles), logical(1))))
    stop("outcome_models must be a list of unique subsets of genetic_roles.")
  if (is.null(names(outcome_models)))
    names(outcome_models) <- vapply(outcome_models, paste, collapse = "+", character(1))
  if (anyNA(names(outcome_models)) || anyDuplicated(names(outcome_models)) || any(!nzchar(names(outcome_models))))
    stop("Each outcome model needs a distinct, nonempty name.")
  col_args <- list()
  for (r in genetic_roles) for (j in 1:2)
    col_args[[paste0("col_pgs", j, "_", r)]] <- pgs_columns[[paste0(r, j)]]
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || anyNA(names(dots)) ||
      anyDuplicated(names(dots)) || any(!nzchar(names(dots)))))
    stop("Additional fitting arguments must have distinct, nonempty names.")
  protected <- c("dat", "genetic_roles", "outcome_roles", "col_outcome", "outcome",
                 "pi", "pi_se", "fixed_effects", "cluster", "covariate_mapping_confirmed", names(col_args))
  if (any(names(dots) %in% protected)) stop("A fitting argument duplicates an analysis setting.")
  for (settings in models) if (!is.list(settings) || any(names(settings) %in% protected) ||
      (length(settings) && (is.null(names(settings)) || anyNA(names(settings)) ||
        anyDuplicated(names(settings)) || any(!nzchar(names(settings))))))
    stop("Model settings need distinct names and must not override data, roles, columns or pi.")
  if (!is.numeric(pi) || length(pi) != 1L ||
      (!is.null(pi_se) && (!is.numeric(pi_se) || length(pi_se) != 1L)))
    stop("Supply one main pi value and at most one pi_se; use pi_grid for sensitivity values.")
  assumptions <- data.frame(kind = "main", pi = pi,
    pi_se = if (is.null(pi_se)) 0 else pi_se, stringsAsFactors = FALSE)
  if (length(pi_grid)) assumptions <- rbind(assumptions,
    data.frame(kind = "grid", pi = unique(pi_grid), pi_se = 0))
  if (any(!is.finite(assumptions$pi) | assumptions$pi < 0 | assumptions$pi >= 1) ||
      any(!is.finite(assumptions$pi_se) | assumptions$pi_se < 0))
    stop("pi must lie in [0,1) and pi_se must be finite and nonnegative.")
  rows <- fits <- list()
  failures <- data.frame(outcome = character(), model = character(), assumption = character(),
    pi = numeric(), regression = character(), message = character())
  samples <- list()
  add_failure <- function(y, m, a, p, r, msg) {
    failures <<- rbind(failures, data.frame(outcome = y, model = m,
      assumption = a, pi = p, regression = r, message = msg))
  }
  for (y in outcomes) {
    d <- dat
    if (complete_cases) d <- d[stats::complete.cases(d[c(score_cols, y, covariates)]), , drop = FALSE]
    if (standardize_outcome && !y %in% binary_outcomes) {
      sy <- stats::sd(d[[y]], na.rm = TRUE)
      if (!is.finite(sy) || sy <= 0) stop("Outcome has no observed variation: ", y)
      d[[y]] <- (d[[y]] - mean(d[[y]], na.rm = TRUE)) / sy
    }
    samples[[y]] <- data.frame(outcome = y, role = genetic_roles, rows = nrow(d),
      n_with_both_scores = vapply(genetic_roles, function(r)
        sum(stats::complete.cases(d[unname(pgs_columns[paste0(r, 1:2)])])), integer(1)),
      n_with_outcome = sum(!is.na(d[[y]])))
    for (m in names(models)) for (a in seq_len(nrow(assumptions))) {
      setting <- assumptions[a, ]
      label <- paste(y, m, setting$kind, format(setting$pi, digits = 15), sep = "|")
      message("Fitting ", label)
      options <- utils::modifyList(dots, models[[m]])
      args <- c(list(dat = d, genetic_roles = genetic_roles, outcome_roles = genetic_roles,
        col_outcome = y, outcome = if (y %in% binary_outcomes) "binary_liability" else "continuous",
        fixed_effects = covariates, cluster = cluster, pi = setting$pi,
        pi_se = if (setting$pi_se > 0) setting$pi_se else NULL,
        covariate_mapping_confirmed = covariate_mapping_confirmed), col_args, options)
      # Always refit a new pi value. This also handles history constraints whose
      # measurement fit depends on pi, rather than treating pi as post-processing.
      fit <- tryCatch(do.call(RAVEL::ravel_family_fit, args), error = identity)
      if (inherits(fit, "error") || !isTRUE(fit$converged)) {
        msg <- if (inherits(fit, "error")) conditionMessage(fit) else fit$error
        add_failure(y, m, setting$kind, setting$pi, "all", msg)
        next
      }
      fits[[label]] <- fit
      for (r in names(outcome_models)) {
        projected <- tryCatch(RAVEL::ravel_family_regression(fit,
          outcome_roles = outcome_models[[r]], pi = setting$pi,
          pi_se = if (setting$pi_se > 0) setting$pi_se else NULL,
          covariate_mapping_confirmed = covariate_mapping_confirmed), error = identity)
        if (inherits(projected, "error")) {
          add_failure(y, m, setting$kind, setting$pi, r, conditionMessage(projected))
          next
        }
        tab <- projected$table
        tab$analysis_outcome <- y
        tab$analysis_model <- m
        tab$assumption <- setting$kind
        tab$pi <- setting$pi
        tab$pi_se <- setting$pi_se
        tab$regression <- r
        tab$genetic_members <- paste(genetic_roles, collapse = "+")
        tab$outcome_predictors <- paste(outcome_models[[r]], collapse = "+")
        tab$outcome_units <- if (y %in% binary_outcomes) "probit_liability" else
          if (standardize_outcome) "phenotype_SD" else "input_units"
        rows[[paste(label, r, sep = "|")]] <- tab
      }
    }
  }
  table <- if (length(rows)) do.call(rbind, rows) else data.frame()
  rownames(table) <- NULL
  sample_counts <- do.call(rbind, samples)
  if (!is.null(output_prefix)) {
    utils::write.csv(table, paste0(output_prefix, ".csv"), row.names = FALSE)
    utils::write.csv(failures, paste0(output_prefix, "_failures.csv"), row.names = FALSE)
    utils::write.csv(sample_counts, paste0(output_prefix, "_sample_counts.csv"), row.names = FALSE)
    saveRDS(fits, paste0(output_prefix, "_fits.rds"))
  }
  if (nrow(failures)) warning(nrow(failures), " fitting/projection failures; inspect $failures.")
  list(table = table, fits = fits, failures = failures, sample_counts = sample_counts,
    outcome_models = outcome_models, genetic_roles = genetic_roles,
    note = "Regressions share the fitted genetic model and sample within each outcome/assumption.")
}
