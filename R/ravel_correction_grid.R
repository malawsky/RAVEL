#' @export
ravel_correct_grid <- function(out, hidden_fraction, n_generations = 15L,
                               residual_models = NULL,
                               uncertainty = c("auto", "lavaan", "bootstrap"),
                               ci_level = NULL, covariate_mapping_confirmed = FALSE,
                               include_raw = TRUE, pi_se = NULL,
                               core_pi_sensitivity = NULL) {
  generations_supplied <- !missing(n_generations)
  if (is.character(out) && length(out) == 1L && !is.na(out)) out <- readRDS(out)
  if (!is.list(out) || is.data.frame(out) || !is.list(out$packets) || !length(out$packets) ||
      !is.data.frame(out$summary) || !all(c("outcome_variable", "regression", "residual_model",
        "outcome_scale", "N") %in% names(out$summary)))
    stop("Supply the full ravel_dataframe_report() output, not out$summary; joint parameter uncertainty is required")
  uncertainty <- match.arg(uncertainty)
  if (!is.numeric(hidden_fraction) || is.complex(hidden_fraction) ||
      !length(hidden_fraction) || any(!is.finite(hidden_fraction)) ||
      any(hidden_fraction < 0 | hidden_fraction >= 1))
    stop("hidden_fraction must contain finite values in [0, 1)")
  valid_generations <- function(x) is.numeric(x) && !is.complex(x) && length(x) > 0L &&
    !anyNA(x) && all(x >= 1) && all(!is.finite(x) |
      (x <= .Machine$integer.max & x == floor(x)))
  if (!valid_generations(n_generations))
    stop("n_generations must contain positive integers, or Inf for equilibrium")
  if (!is.null(pi_se) && (!is.numeric(pi_se) || is.complex(pi_se) ||
      length(pi_se) != 1L || !is.finite(pi_se) || pi_se < 0 || !is.finite(pi_se^2)))
    stop("pi_se must be a finite nonnegative scalar on the raw probability scale")
  if (is.null(pi_se) && !is.null(core_pi_sensitivity))
    stop("core_pi_sensitivity requires pi_se")
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
  sensitivities <- setNames(vector("list", length(residual_models)), residual_models)
  if (!is.null(core_pi_sensitivity)) {
    if (is.list(core_pi_sensitivity)) {
      if (is.null(names(core_pi_sensitivity)) || anyNA(names(core_pi_sensitivity)) ||
          anyDuplicated(names(core_pi_sensitivity)) ||
          !all(names(core_pi_sensitivity) %in% residual_models))
        stop("core_pi_sensitivity lists must be named by selected residual models")
      sensitivities[names(core_pi_sensitivity)] <- core_pi_sensitivity
    } else {
      if (length(residual_models) != 1L)
        stop("For multiple residual models, supply core_pi_sensitivity as a named list")
      sensitivities[[1L]] <- core_pi_sensitivity
    }
    for (model in names(sensitivities)) {
      s <- sensitivities[[model]]
      labels <- .report_packet_labels(packets[[model]])
      if (!is.null(s) && (!is.numeric(s) || is.complex(s) || any(!is.finite(s)) ||
          is.null(names(s)) || anyDuplicated(names(s)) || !setequal(names(s), labels)))
        stop(if (length(labels) == 6L)
          "Each core_pi_sensitivity must be a finite vector named for all six core parameters" else
          "Each core_pi_sensitivity must be a finite vector named for all core parameters in its packet")
    }
  }
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
    if (identical(packet$settings$constraint_spec$mode, "estimated_state") &&
        !identical(packet$diagnostics$estimated_state_identified, TRUE))
      stop("Estimated-state packets require a validated identification check; regenerate the report")
    if (!identical(packet$software$transform_fingerprint, fingerprint))
      stop("Correction code differs from the report; align package versions before proceeding")
    if (isTRUE(packet$settings$requires_covariate_mapping_review) && !covariate_mapping_confirmed)
      stop("Confirm that the visible-to-total mapping applies conditional on these covariates before correcting")
    bootstrap <- uncertainty == "bootstrap" ||
      (uncertainty == "auto" && identical(packet$settings$uncertainty, "bootstrap"))
    if (bootstrap && is.null(packet$bootstrap$core_draws))
      stop("Bootstrap requested but a selected packet lacks joint core-parameter draws")
  }
  grids <- lapply(residual_models, function(model) {
    history <- packets[[model]]$settings$constraint_spec
    mode <- history$mode %||% out$settings$fit_constraints %||% "none"
    if (identical(mode, "estimated_state")) {
      if (generations_supplied)
        stop("Do not supply n_generations for estimated-state packets; the inferred state has no imposed generation count")
      generations <- NA_real_
    } else if (!generations_supplied && identical(mode, "equilibrium")) {
      generations <- Inf
    } else if (!generations_supplied && identical(mode, "n_gen")) {
      generations <- history$n_generations %||% out$settings$am_n_generations
      if (!valid_generations(generations) || length(generations) != 1L || !is.finite(generations))
        stop("The constrained packet must contain its fitted positive integer generation count")
    } else generations <- n_generations
    grid <- expand.grid(hidden_fraction = unique(unname(hidden_fraction)),
      n_generations = unique(unname(generations)), KEEP.OUT.ATTRS = FALSE)
    grid$am_correction <- if (identical(mode, "estimated_state")) "estimated_state" else
      ifelse(is.infinite(grid$n_generations), "equilibrium", "n_gen")
    grid
  })
  names(grids) <- residual_models
  # Equal scenarios share an identifier even if selected packets inherit
  # different histories. Existing unconstrained Cartesian ordering is retained.
  scenario_key <- function(grid) paste(sprintf("%.17g", grid$hidden_fraction),
    grid$n_generations, grid$am_correction, sep = ":")
  scenario_keys <- unique(unlist(lapply(grids, scenario_key), use.names = FALSE))
  grids <- lapply(grids, function(grid) {
    grid$scenario_id <- match(scenario_key(grid), scenario_keys)
    grid
  })
  scales <- if (include_raw) c("corrected_std", "corrected_raw") else "corrected_std"
  values <- c("estimate", "se", "ci_low", "ci_high", "p_value")
  results <- lapply(residual_models, function(model) {
    packet <- packets[[model]]
    grid <- grids[[model]]
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
      regime <- grid$am_correction[i]
      messages <- character()
      error <- ""
      result <- if (is.null(packet)) NULL else tryCatch(withCallingHandlers({
        args <- list(am_correction = regime,
          ci_level = ci_level, covariate_mapping_confirmed = covariate_mapping_confirmed,
          uncertainty = uncertainty)
        if (regime == "n_gen") args$am_n_generations <- generations
        if (is.null(pi_se)) do.call(ravel_correct_packet, c(list(packet = packet, h = h), args)) else
          do.call(ravel_correct_pi, c(list(object = packet, pi = h, pi_se = pi_se,
            core_pi_sensitivity = sensitivities[[model]]), args))
      }, warning = function(w) {
            messages <<- c(messages, conditionMessage(w))
            invokeRestart("muffleWarning")
          }), error = function(e) {
            error <<- conditionMessage(e)
            NULL
          })
      if (is.null(result)) {
        if (is.null(packet)) error <- "No usable original SEM fit for this residual model"
        roles <- .parental_sib_roles(packet$settings$parental_sibs %||% out$settings$parental_sibs %||% "none")
        parameters <- lapply(scales, function(scale) .report_parameters(roles, scale))
        result <- data.frame(scale = rep(scales, lengths(parameters)),
          parameter = unlist(parameters, use.names = FALSE), stringsAsFactors = FALSE)
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
      if (!is.na(outcome_scale) || is.null(result$outcome_scale))
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
      if (!is.null(pi_se)) {
        result$pi_se <- pi_se
        result$pi_uncertainty_source <- "se"
      }
      result[c("outcome_variable", "regression", "residual_model", "scale", "outcome_scale",
        "parameter", values, "N", "inference", "scenario_id", "hidden_fraction", "n_generations",
        "am_correction", "ci_level", "n_boot_success", "uncertainty_status", "fit_status",
        "correction_status", "error", "warning",
        if (!is.null(pi_se)) c("pi_se", "pi_uncertainty_source"))]
    }))
  })
  result <- do.call(rbind, results)
  # A scalar correction carries covariance/state attributes for one cell.
  # They cannot describe the combined grid and must not leak from its first row.
  attributes(result) <- attributes(result)[c("names", "row.names", "class")]
  rownames(result) <- NULL
  attr(result, "ci_level") <- ci_level
  settings <- out$settings %||% list()
  settings$ci_level <- ci_level
  inference <- unique(result$inference[result$inference != "unavailable"])
  selected <- unique(ifelse(grepl("bootstrap", inference), "bootstrap", "lavaan"))
  settings$uncertainty <- if (length(selected) == 1L) selected else if (length(selected) > 1L)
    "mixed" else if (uncertainty != "auto") uncertainty else settings$uncertainty %||% "unavailable"
  settings$pi_se <- pi_se
  attr(result, "settings") <- settings
  if (any(result$correction_status != "ok"))
    warning("Some corrections or intervals are unavailable; inspect correction_status, error and n_boot_success")
  result
}
