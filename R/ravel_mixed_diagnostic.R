# Automatic auxiliary diagnostics for the new optional fitter. The primary
# family likelihood and its fitted random effects are never modified here.
.ravel_mixed_run_diagnostic <- function(prepared, diagnostic_control = list()) {
  allowed <- c("age", "target_age", "occasion", "occasion_value", "selection_seed")
  if (!is.list(diagnostic_control) || (length(diagnostic_control) &&
      (is.null(names(diagnostic_control)) || anyDuplicated(names(diagnostic_control)) ||
       any(!names(diagnostic_control) %in% allowed))))
    stop("diagnostic_control must be a named list of age, target_age, occasion, occasion_value or selection_seed.")
  spec <- prepared$spec
  # A separate auxiliary fit starts in well-conditioned score units, just as
  # the mixed likelihood does. No parameters or constraints of a primary SEM
  # are altered. This matters for raw PGSs with extremely small variances.
  dat <- prepared$data
  for (name in unname(prepared$pgs_columns)) dat[[name]] <- prepared$adjusted_data[[name]]
  age <- spec$random_slope
  if (is.null(age) && "age" %in% names(dat) && is.numeric(dat$age)) age <- "age"
  args <- list(dat = dat, family_id = prepared$family_id, subject_id = prepared$subject_id,
    cluster = prepared$cluster, covariates = spec$covariates, interactions = spec$interactions,
    age = age, outcome = spec$outcome, col_outcome = prepared$col_outcome,
    fit_constraints = spec$mode, am_n_generations = spec$generations, pi = spec$pi,
    pi_se = NULL)
  if (!is.null(age) && age %in% names(spec$centers)) args$target_age <- unname(spec$centers[[age]])
  for (role in c("child", "mother", "father")) for (k in 1:2)
    args[[paste0("col_pgs", k, "_", role)]] <- prepared$pgs_columns[[paste0(role, k)]]
  args <- utils::modifyList(args, diagnostic_control, keep.null = TRUE)
  result <- do.call(ravel_score_diagnostic, args)
  result$auxiliary$score_units <- "Primary mixed fitter's role-specific standardized PGS units"
  result$auxiliary$random_effects_omitted <- c(if (spec$family_random) "maternal-family intercept",
    if (spec$individual_random) "individual intercept", if (!is.null(spec$random_slope)) "individual slope")
  result
}
