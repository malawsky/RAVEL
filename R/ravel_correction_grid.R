#' @export
ravel_correct_grid <- function(out, hidden_fraction, n_generations = 15L,
                               residual_models = NULL,
                               uncertainty = c("auto", "lavaan", "bootstrap"),
                               ci_level = NULL, covariate_mapping_confirmed = FALSE,
                               include_raw = TRUE) {
  if (is.character(out) && length(out) == 1L && !is.na(out)) out <- readRDS(out)
  if (!is.list(out) || is.data.frame(out) || !is.list(out$packets) || !length(out$packets) ||
      !is.data.frame(out$summary) || !all(c("outcome_variable", "regression", "residual_model",
        "outcome_scale", "N") %in% names(out$summary)))
    stop("Supply the full ravel_dataframe_report() output, not out$summary; joint parameter uncertainty is required")
  uncertainty <- match.arg(uncertainty)
  if (!is.numeric(hidden_fraction) || !length(hidden_fraction) || any(!is.finite(hidden_fraction)) ||
      any(hidden_fraction < 0 | hidden_fraction >= 1))
    stop("hidden_fraction must contain finite values in [0, 1)")
  if (!is.numeric(n_generations) || !length(n_generations) || anyNA(n_generations) ||
      any(n_generations < 1) || any(is.finite(n_generations) &
        (n_generations > .Machine$integer.max | n_generations != floor(n_generations))))
    stop("n_generations must contain positive integers, or Inf for equilibrium")
  for (value in list(include_raw, covariate_mapping_confirmed))
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop("include_raw and covariate_mapping_confirmed must be TRUE or FALSE")
  available <- names(out$packets)
  if (is.null(available) || anyNA(available) || anyDuplicated(available) ||
      !all(available %in% c("fixed", "free"))) stop("Invalid residual-model names in out$packets")
  if (is.null(residual_models)) residual_models <- available
  if (!is.character(residual_models) || !length(residual_models) || anyNA(residual_models) ||
      anyDuplicated(residual_models) || !all(residual_models %in% available))
    stop("residual_models must select names present in out$packets without duplicates")
  packets <- out$packets[residual_models]
  usable <- packets[!vapply(packets, is.null, logical(1))]
  if (is.null(ci_level)) {
    ci_level <- out$settings$ci_level
    if (is.null(ci_level) && length(usable)) ci_level <- usable[[1L]]$settings$ci_level
    if (is.null(ci_level)) ci_level <- .95
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L || !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1)
    stop("ci_level must lie between zero and one")
  fingerprint <- .report_transform_fingerprint()
  for (packet in usable) {
    if (!identical(packet$schema, "ravel_visible_transfer_v1") || !isTRUE(packet$diagnostics$admissible))
      stop("A supplied correction packet is not an admissible visible-scale packet")
    if (!identical(packet$software$transform_fingerprint, fingerprint))
      stop("Correction code differs from the report; align package versions before proceeding")
    if (isTRUE(packet$settings$requires_covariate_mapping_review) && !covariate_mapping_confirmed)
      stop("Confirm that the visible-to-total mapping applies conditional on these covariates before correcting")
    bootstrap <- uncertainty == "bootstrap" ||
      (uncertainty == "auto" && identical(packet$settings$uncertainty, "bootstrap"))
    if (bootstrap && is.null(packet$bootstrap$core_draws))
      stop("Bootstrap requested but a selected packet lacks joint core-parameter draws")
  }
  grid <- expand.grid(hidden_fraction = unique(unname(hidden_fraction)),
    n_generations = unique(unname(n_generations)), KEEP.OUT.ATTRS = FALSE)
  grid$scenario_id <- seq_len(nrow(grid))
  scales <- if (include_raw) c("corrected_std", "corrected_raw") else "corrected_std"
  values <- c("estimate", "se", "ci_low", "ci_high", "p_value")
  results <- lapply(residual_models, function(model) {
    packet <- packets[[model]]
    original <- out$summary[out$summary$regression == "RAVEL" & out$summary$residual_model == model, , drop = FALSE]
    outcome <- if (!is.null(packet)) packet$outcome$name else
      if (nrow(original)) original$outcome_variable[1L] else NA_character_
    outcome_scale <- if (nrow(original)) original$outcome_scale[1L] else NA_character_
    n <- if (!is.null(packet)) packet$settings$n_complete else
      if (nrow(original)) original$N[1L] else NA_integer_
    if (is.null(n)) n <- NA_integer_
    do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
      h <- grid$hidden_fraction[i]
      generations <- grid$n_generations[i]
      regime <- if (is.infinite(generations)) "equilibrium" else "n_gen"
      messages <- character()
      error <- ""
      result <- if (is.null(packet)) NULL else tryCatch(withCallingHandlers(
        ravel_correct_packet(packet, h = h, am_correction = regime,
          am_n_generations = if (regime == "n_gen") generations else 15L,
          ci_level = ci_level, covariate_mapping_confirmed = covariate_mapping_confirmed,
          uncertainty = uncertainty), warning = function(w) {
            messages <<- c(messages, conditionMessage(w))
            invokeRestart("muffleWarning")
          }), error = function(e) {
            error <<- conditionMessage(e)
            NULL
          })
      if (is.null(result)) {
        if (is.null(packet)) error <- "No usable original SEM fit for this residual model"
        result <- data.frame(scale = rep(scales, each = 4L),
          parameter = unlist(lapply(scales, function(scale)
            c("beta_D", "beta_M", "beta_P", if (scale == "corrected_std") "rho_corrected" else "cov_corrected"))),
          stringsAsFactors = FALSE)
        for (column in values) result[[column]] <- NA_real_
        result$inference <- "unavailable"
        result$n_boot_success <- NA_integer_
        result$uncertainty_status <- "unavailable"
        result$correction_status <- if (is.null(packet)) "unavailable_model" else "failed"
      } else {
        result <- result[result$scale %in% scales, , drop = FALSE]
        result$correction_status <- ifelse(result$uncertainty_status == "ok", "ok", result$uncertainty_status)
      }
      result$outcome_variable <- outcome
      result$regression <- "RAVEL"
      result$residual_model <- model
      result$outcome_scale <- outcome_scale
      result$N <- n
      result$scenario_id <- grid$scenario_id[i]
      result$hidden_fraction <- h
      result$n_generations <- generations
      result$am_correction <- regime
      result$ci_level <- ci_level
      result$fit_status <- if (is.null(packet)) "failed" else "ok"
      result$error <- error
      result$warning <- paste(unique(messages), collapse = "; ")
      result[c("outcome_variable", "regression", "residual_model", "scale", "outcome_scale",
        "parameter", values, "N", "inference", "scenario_id", "hidden_fraction", "n_generations",
        "am_correction", "ci_level", "n_boot_success", "uncertainty_status", "fit_status",
        "correction_status", "error", "warning")]
    }))
  })
  result <- do.call(rbind, results)
  rownames(result) <- NULL
  if (any(result$correction_status != "ok"))
    warning("Some corrections or intervals are unavailable; inspect correction_status, error and n_boot_success")
  result
}
