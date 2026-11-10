# Auxiliary, unadjusted score diagnostic. The scientific fit is never modified.

#' Automatically attempt an unadjusted shared-factor diagnostic
#'
#' Refit an auxiliary model without covariates or genetic interactions and run
#' [ravel_test_score_signal()]. Related/repeated records are reduced to one
#' reproducibly selected individual and one occasion per independent cluster.
#' This is a marginal diagnostic, not a test of the adjusted or longitudinal
#' model. The original diagnostic's identification restrictions still apply.
#'
#' @param dat Data frame for the auxiliary analysis. When `primary_fit` supplies
#'   fixed residual covariances, its PGS units must match these data.
#' @param ... Named arguments to [ravel_fit()], including column mappings and
#'   genetic restrictions. These override settings recovered from `primary_fit`.
#' @param primary_fit Optional result of [ravel_fit()] supplying its column maps
#'   and genetic settings. It is not modified or required to have converged.
#' @param family_id Maternal-family ID column, used to check nesting when present.
#' @param subject_id Individual ID column. Required when `cluster` repeats.
#' @param cluster Extended-family ID column. One individual is selected per
#'   cluster. Set `NULL` only when the supplied individuals are independent.
#' @param covariates Adjustment variables in the primary model, omitted here.
#' @param interactions Genetic moderators in the primary model, omitted here.
#' @param age Optional numeric age column. With repeated observations, select
#'   the closest age to `target_age`, or the earliest age if no target is given.
#' @param target_age Optional finite scalar specifying the diagnostic age.
#' @param occasion Optional visit column, used to select a prespecified visit or
#'   the earliest ordered visit when no age is supplied, and to break age ties.
#' @param occasion_value Optional exact visit value. Selection does not replace
#'   a selected individual who lacks this visit with another family member.
#' @param selection_seed Integer seed for selecting individuals, uniformly over
#'   sorted unique IDs within each sorted cluster. The caller's RNG is preserved.
#' @return A `ravel_score_signal_test` with `auxiliary`, `auxiliary_fit`, and
#'   `selection` metadata. Unsupported or failed analyses retain missing inference
#'   and an explicit status. Binary outcomes use a separate exploratory
#'   conditional-probit calculation. Its nominal p-values are returned with
#'   an explicit warning because it has not passed null calibration; binary
#'   outcomes are never treated as continuous outcomes.
#' @details Selection occurs before removal of incomplete score/outcome records
#'   and never uses outcome values. Missing selected records are not replaced.
#'   Both original and selected sample counts, omissions, warnings, and row IDs
#'   are retained. A test on selected independent records does not assess random
#'   effects or genetic interaction slopes. Removing adjustment can itself alter
#'   the tested marginal associations. External pi is held fixed in this test.
#' @export
ravel_score_diagnostic <- function(dat, ..., primary_fit = NULL,
                                   family_id = "M_ID", subject_id = "IID",
                                   cluster = "FID", covariates = NULL,
                                   interactions = NULL, age = NULL,
                                   target_age = NULL, occasion = NULL,
                                   occasion_value = NULL, selection_seed = 1L) {
  out <- structure(list(status = "not_run", reason = "", statistic = NA_real_,
    df = NA_real_, p_value = NA_real_, raw_statistic = NA_real_, test = NA_character_,
    N = if (is.data.frame(dat)) nrow(dat) else NA_integer_,
    diagnostics = list(), fits = list(), auxiliary_fit = NULL,
    auxiliary = list(type = "unadjusted_independent_cross_sectional",
      primary_model_test = FALSE, covariates_omitted = character(),
      interactions_omitted = character(), warnings = character()),
    selection = list()), class = c("ravel_score_diagnostic", "ravel_score_signal_test"))
  finish <- function(status, reason) { out$status <- status; out$reason <- reason; out }
  dots <- list(...)
  if (!is.data.frame(dat) || !nrow(dat)) return(finish("invalid_data", "dat must be a nonempty data frame."))
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) || anyDuplicated(names(dots))))
    return(finish("invalid_arguments", "Every fitting argument in ... must have a unique name."))
  if ("dat" %in% names(dots)) return(finish("invalid_arguments", "Supply dat only once."))
  unknown <- setdiff(names(dots), names(formals(ravel_fit)))
  if (length(unknown)) return(finish("invalid_arguments", paste("Unknown ravel_fit argument(s):", paste(unknown, collapse = ", "))))
  args <- list()
  if (!is.null(primary_fit)) {
    if (!inherits(primary_fit, "ravel_results"))
      return(finish("invalid_arguments", "primary_fit must be a ravel_fit result; pass explicit settings for other fit classes."))
    args <- list(outcome = primary_fit$outcome, measurement = primary_fit$measurement,
      crossperson = primary_fit$crossperson_structure,
      fixed_crossperson = primary_fit[["fixed_crossperson"]],
      fixed_crossperson_corr = primary_fit[["fixed_crossperson_corr"]],
      fit_constraints = primary_fit$constraint_spec$mode,
      am_correction = primary_fit$am_correction,
      pi = primary_fit$pi %||% primary_fit$h,
      parental_sibs = primary_fit$parental_sibs$mode %||% "none")
    if (identical(args$fit_constraints, "n_gen"))
      args$am_n_generations <- primary_fit$constraint_spec$n_generations
    map <- primary_fit$score_signal_input$column_map
    for (canonical in names(map)) {
      argument <- if (canonical == "Y") "col_outcome" else
        sub("^PGS([12])_", "col_pgs\\1_", canonical)
      if (argument %in% names(formals(ravel_fit))) args[[argument]] <- unname(map[[canonical]])
    }
    if (inherits(primary_fit$fit, "lavaan")) {
      pt <- lavaan::parameterTable(primary_fit$fit)
      paths <- pt[pt$lhs == "Y" & pt$op == "~" & pt$rhs %in% c("Gm", "Gf"), , drop = FALSE]
      args$fix_parent_paths <- nrow(paths) == 2L && all(paths$free == 0L & paths$est == 0)
    }
    if (missing(cluster) && !is.null(primary_fit$cluster)) cluster <- primary_fit$cluster
  }
  args <- utils::modifyList(args, dots, keep.null = TRUE)
  if ("h" %in% names(dots) && !"pi" %in% names(dots)) args$pi <- NULL
  covariates <- unique(c(covariates, primary_fit$fixed_effects, args$fixed_effects))
  interactions <- unique(c(interactions, primary_fit$interactions$moderators, args$interactions))
  out$auxiliary$covariates_omitted <- covariates
  out$auxiliary$interactions_omitted <- interactions
  announce <- function(message) {
    out$auxiliary$warnings <<- c(out$auxiliary$warnings, message)
    warning(message, call. = FALSE)
  }
  if (length(covariates)) announce("Covariates were not used to compute the one-factor test.")
  if (length(interactions)) announce(paste("Genetic interactions were not used to compute the one-factor test;",
    "this auxiliary marginal test does not test interaction slopes."))
  args$outcome <- args$outcome %||% "continuous"
  args$fit_constraints <- args$fit_constraints %||% "none"
  out$auxiliary$fit_constraints <- args$fit_constraints
  out$auxiliary$pi <- args$pi %||% args$h %||% 0
  if (!args$outcome %in% c("continuous", "binary_liability"))
    return(finish("unsupported_outcome", "Only continuous or binary_liability outcomes are supported."))
  column <- function(name, optional = FALSE) {
    is.null(name) || (is.character(name) && length(name) == 1L && !is.na(name) && nzchar(name) &&
      (optional || name %in% names(dat)))
  }
  if (!column(cluster) || !column(subject_id, TRUE) || !column(family_id, TRUE) ||
      !column(age) || !column(occasion))
    return(finish("invalid_arguments", "ID and occasion arguments must be single column names; the requested cluster/age/occasion columns must exist."))
  if (!is.numeric(selection_seed) || length(selection_seed) != 1L || !is.finite(selection_seed) ||
      selection_seed < 0 || selection_seed > .Machine$integer.max || selection_seed != as.integer(selection_seed))
    return(finish("invalid_arguments", "selection_seed must be one nonnegative integer."))
  if (!is.null(target_age) && (is.null(age) || !is.numeric(target_age) || length(target_age) != 1L || !is.finite(target_age)))
    return(finish("invalid_arguments", "target_age requires an age column and one finite numeric value."))
  if (!is.null(age) && (!is.numeric(dat[[age]]) || any(!is.na(dat[[age]]) & !is.finite(dat[[age]]))))
    return(finish("invalid_data", "The age column must be numeric and finite or missing."))
  if (!is.null(occasion_value) && (is.null(occasion) || length(occasion_value) != 1L || is.na(occasion_value)))
    return(finish("invalid_arguments", "occasion_value requires an occasion column and one nonmissing value."))

  family_present <- !is.null(family_id) && family_id %in% names(dat)
  subject_present <- !is.null(subject_id) && subject_id %in% names(dat)
  id_columns <- unique(c(cluster, if (family_present) family_id, if (subject_present) subject_id))
  if (any(vapply(dat[id_columns], function(z) anyNA(z) || any(!nzchar(as.character(z))), logical(1))))
    return(finish("invalid_design", "Identification columns must not contain missing or empty IDs."))
  if (is.null(cluster)) {
    if (family_present && anyDuplicated(dat[[family_id]]) && !subject_present)
      return(finish("needs_subject_id", "Repeated families require individual IDs and an independent-cluster specification."))
    group <- if (subject_present) as.character(dat[[subject_id]]) else as.character(seq_len(nrow(dat)))
    if (family_present) {
      units <- unique(data.frame(f = as.character(dat[[family_id]]), s = group))
      if (anyDuplicated(units$f)) return(finish("needs_cluster", "Multiple individuals share a maternal family; supply an extended-family cluster column."))
    }
  } else group <- as.character(dat[[cluster]])
  if (anyDuplicated(group) && !subject_present)
    return(finish("needs_subject_id", "Repeated clusters require subject_id so individuals can be selected without favouring those with more visits."))
  subject <- if (subject_present) as.character(dat[[subject_id]]) else as.character(seq_len(nrow(dat)))
  nested <- function(lower, upper) !anyDuplicated(unique(data.frame(lower = lower, upper = upper))$lower)
  if (!nested(subject, group) || (family_present && (!nested(subject, as.character(dat[[family_id]])) ||
      !is.null(cluster) && !nested(as.character(dat[[family_id]]), group))))
    return(finish("invalid_design", "Each individual must belong to one maternal family and cluster, and each maternal family must lie within one FID."))

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(as.integer(selection_seed))
  groups <- sort(unique(group))
  selected <- integer()
  absent_visit <- character()
  for (g in groups) {
    ids <- sort(unique(subject[group == g]))
    chosen <- ids[if (length(ids) == 1L) 1L else sample.int(length(ids), 1L)]
    rows <- which(group == g & subject == chosen)
    if (!is.null(occasion_value)) rows <- rows[!is.na(dat[[occasion]][rows]) &
      as.character(dat[[occasion]][rows]) == as.character(occasion_value)]
    if (!is.null(age)) rows <- rows[!is.na(dat[[age]][rows])]
    if (!is.null(occasion) && is.null(age) && is.null(occasion_value)) rows <- rows[!is.na(dat[[occasion]][rows])]
    if (!length(rows)) { absent_visit <- c(absent_visit, g); next }
    if (length(rows) > 1L) {
      if (is.null(age) && is.null(occasion))
        return(finish("needs_occasion", "Repeated observations require age or occasion to define the auxiliary cross-sectional diagnostic."))
      if (!is.null(age)) {
        a <- dat[[age]][rows]
        distance <- if (is.null(target_age)) a else abs(a - target_age)
        rows <- rows[distance == min(distance)]
        rows <- rows[dat[[age]][rows] == min(dat[[age]][rows])]
      }
      if (length(rows) > 1L && !is.null(occasion)) {
        visit <- dat[[occasion]][rows]
        if (anyNA(visit)) return(finish("ambiguous_occasion", "Tied diagnostic ages have missing occasion IDs."))
        rows <- rows[visit == sort(unique(visit))[1L]]
      }
      if (length(rows) > 1L)
        return(finish("ambiguous_occasion", "More than one row remains at the selected individual/occasion; duplicate observations must be resolved explicitly."))
    }
    selected <- c(selected, rows)
  }
  out$selection <- list(cluster = cluster, family_id = if (family_present) family_id else NULL,
    subject_id = if (subject_present) subject_id else NULL, seed = as.integer(selection_seed),
    age = age, target_age = target_age, occasion = occasion, occasion_value = occasion_value,
    N_input = nrow(dat), N_subjects = length(unique(subject)), N_clusters = length(groups),
    N_selected = length(selected), clusters_without_selected_visit = absent_visit,
    rule = paste("One uniformly sampled individual per sorted independent cluster, selected before complete-case filtering;",
      if (!is.null(occasion_value)) "prespecified occasion;" else "",
      if (!is.null(age)) if (is.null(target_age)) "earliest observed age." else "nearest target age, earlier age breaks ties." else
        if (!is.null(occasion)) "earliest ordered occasion." else "one observation per individual."),
    rows_selected = selected)

  sibs <- args$parental_sibs %||% "none"
  sib_roles <- switch(sibs, none = character(), maternal = "mat_sib", paternal = "pat_sib", both = c("mat_sib", "pat_sib"), NULL)
  if (is.null(sib_roles)) return(finish("invalid_arguments", "Unknown parental_sibs setting."))
  roles <- c("child", "mother", "father", sib_roles)
  canonical <- c("Y", unlist(lapply(roles, function(r) paste0("PGS", 1:2, "_", r))))
  arguments <- c("col_outcome", sub("^PGS([12])_", "col_pgs\\1_", canonical[-1L]))
  defaults <- c("Y", canonical[-1L])
  mapping <- stats::setNames(vapply(seq_along(arguments), function(j) args[[arguments[j]]] %||% defaults[j], ""), canonical)
  if (anyDuplicated(mapping) || any(!mapping %in% names(dat)))
    return(finish("invalid_data", "Every outcome/PGS mapping must name a distinct existing data column."))
  d <- dat[selected, unname(mapping), drop = FALSE]
  names(d) <- canonical
  if (any(!vapply(d, is.numeric, logical(1))))
    return(finish("invalid_data", "The diagnostic outcome and all PGSs must be numeric."))
  complete <- stats::complete.cases(d) & apply(as.matrix(d), 1L, function(z) all(is.finite(z)))
  out$selection$N_incomplete_selected <- sum(!complete)
  out$selection$rows_used <- selected[complete]
  out$selection$rows_incomplete <- selected[!complete]
  out$selection$N_used <- sum(complete)
  out$selection$individuals_used <- subject[selected[complete]]
  out$selection$clusters_used <- group[selected[complete]]
  out$auxiliary$column_map <- mapping
  out$N <- sum(complete)
  if (out$N != nrow(dat)) announce(sprintf(paste(
    "The one-factor test uses an auxiliary independent cross-sectional sample (%d complete observations from %d input rows).",
    "It does not test random effects or genetic interaction slopes."), out$N, nrow(dat)))
  if (sum(complete) < 50L) return(finish("insufficient_data", "Fewer than 50 complete independently selected observations remain."))
  d <- d[complete, , drop = FALSE]
  sample_variances <- vapply(d, stats::var, numeric(1))
  invalid_variances <- names(sample_variances)[!is.finite(sample_variances) | sample_variances <= 0]
  if (length(invalid_variances)) return(finish("invalid_data", paste(
    "Selected diagnostic variables must have finite positive variance:",
    paste(invalid_variances, collapse = ", "))))
  for (j in seq_along(arguments)) args[[arguments[j]]] <- canonical[j]
  args$fixed_effects <- NULL
  args$interactions <- NULL
  args$interaction_center <- NULL
  args$interaction_roles <- NULL
  args$cluster <- NULL
  args$cluster_var <- NULL
  args$n_boot <- 0
  args$pi_se <- NULL
  args$core_pi_cov <- NULL
  args$verbose <- FALSE
  capture <- .ravel_score_signal_capture(do.call(ravel_fit, c(list(dat = d), args)))
  out$auxiliary$fit_warnings <- capture$warnings
  out$auxiliary_fit <- capture$value
  out$auxiliary$refit_converged <- isTRUE(capture$value$converged)
  if (is.null(capture$value)) return(finish("auxiliary_fit_error", capture$error))
  if (identical(args$outcome, "binary_liability")) {
    if (!exists(".ravel_test_binary_score_signal", mode = "function", inherits = TRUE))
      return(finish("unsupported_outcome", "A conditional-probit shared-factor diagnostic is not available; binary outcomes are not treated as continuous."))
    # The binary direction test estimates the measurement model and an
    # unrestricted probit model. An inadmissible outcome null is not needed to
    # construct those fits, but its failure remains recorded above.
    binary_input <- capture$value
    binary_input$outcome <- "binary_liability"
    binary_input$constraint_spec <- binary_input$constraint_spec %||% list(
      mode = args$fit_constraints, h = args$pi %||% args$h %||% 0,
      n_generations = args$am_n_generations %||% 1L,
      depends_on_pi = identical(args$fit_constraints, "n_gen") && (args$am_n_generations %||% 1L) > 1L)
    binary_input$fixed_crossperson <- binary_input[["fixed_crossperson"]] %||% args$fixed_crossperson
    binary_input$fixed_crossperson_corr <- binary_input[["fixed_crossperson_corr"]] %||% args$fixed_crossperson_corr
    binary_input$crossperson_structure <- args$crossperson %||% "free"
    binary_input$parental_sibs <- binary_input$parental_sibs %||% list(mode = sibs,
      roles = .parental_sib_roles(sibs))
    binary_input$score_signal_input <- binary_input$score_signal_input %||% list(
      column_map = stats::setNames(canonical, canonical), dependent_families = FALSE)
    test <- .ravel_score_signal_capture(.ravel_test_binary_score_signal(binary_input, d))
  } else test <- .ravel_score_signal_capture(ravel_test_score_signal(capture$value))
  if (is.null(test$value)) return(finish("diagnostic_error", test$error))
  diagnostic <- test$value
  if (identical(args$outcome, "binary_liability") && identical(diagnostic$status, "ok")) {
    # Report the requested nominal result, distinguishing successful numerical
    # calculation from calibration of its reference distribution.
    diagnostic$diagnostics$exploratory <- list(
      calibrated = FALSE,
      label = "Nominal asymptotic p-value; null calibration showed excess rejections",
      nominal_p_value = diagnostic$p_value,
      statistic = diagnostic$statistic, df = diagnostic$df,
      null_calibration = list(attempts = 500L, valid = 499L,
        rejected = 40L, nominal_level = .05, rejection_rate = 40 / 499))
    diagnostic$calibration_status <- "uncalibrated"
    diagnostic$reason <- paste(
      "Binary one-factor diagnostic: the nominal p-value is reported with a calibration warning.",
      "Null simulations rejected 40/499 (8.0%) at a nominal 5% threshold; interpret this p-value cautiously.",
      "This warning concerns the auxiliary diagnostic; primary model inference is separate.")
    announce(diagnostic$reason)
    diagnostic$test <- paste("Exploratory, uncalibrated:", diagnostic$test)
  }
  diagnostic$auxiliary <- out$auxiliary
  diagnostic$auxiliary$diagnostic_warnings <- test$warnings
  diagnostic$auxiliary_fit <- capture$value
  diagnostic$selection <- out$selection
  class(diagnostic) <- c("ravel_score_diagnostic", class(diagnostic))
  diagnostic
}

#' Print an auxiliary unadjusted score diagnostic
#' @param x An auxiliary score diagnostic returned by [ravel_score_diagnostic()].
#' @param ... Passed to the underlying diagnostic print method.
#' @export
print.ravel_score_diagnostic <- function(x, ...) {
  print.ravel_score_signal_test(x, ...)
  cat("Scope: auxiliary unadjusted independent cross-sectional diagnostic.\n")
  if (length(x$selection$N_input)) cat("Input rows:", x$selection$N_input,
    "; independent complete observations:", x$selection$N_used %||% NA_integer_, "\n")
  for (message in x$auxiliary$warnings) cat("Note:", message, "\n")
  invisible(x)
}
