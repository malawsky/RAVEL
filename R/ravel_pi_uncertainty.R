# Utilities used only by the external-pi uncertainty extension. Existing
# fixed-h correction entry points deliberately retain their existing behavior.
.ravel_pi_ci_variance <- function(pi, pi_ci, pi_ci_scale, pi_ci_level) {
  scalar <- function(x) is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x)
  if (!scalar(pi) || pi < 0 || pi >= 1) stop("pi must be a finite scalar in [0, 1)")
  if (!scalar(pi_ci_level) || pi_ci_level <= 0 || pi_ci_level >= 1)
    stop("pi_ci_level must lie strictly between zero and one")
  if (!is.numeric(pi_ci) || is.complex(pi_ci) || length(pi_ci) != 2L ||
      any(!is.finite(pi_ci)) || pi_ci[1L] > pi || pi_ci[2L] < pi)
    stop("pi_ci must contain finite lower and upper bounds enclosing pi")
  if (all(pi_ci == pi)) return(list(variance = 0, working_se = 0))
  if (pi == 0) stop("A nondegenerate Wald interval at pi = 0 is not supported; supply a boundary-aware analysis")
  working <- if (pi_ci_scale == "logit") {
    if (any(pi_ci <= 0 | pi_ci >= 1))
      stop("Logit-Wald pi_ci bounds must lie strictly between zero and one")
    stats::qlogis(c(pi_ci[1L], pi, pi_ci[2L]))
  } else c(pi_ci[1L], pi, pi_ci[2L])
  widths <- c(working[2L] - working[1L], working[3L] - working[2L])
  if (min(widths) <= 0 || abs(diff(widths)) > max(1e-8, 1e-4 * max(widths)))
    stop("pi_ci must be symmetric about pi on the selected scale; arbitrary percentile or asymmetric intervals do not identify a Wald variance")
  working_se <- mean(widths) / stats::qnorm((1 - pi_ci_level) / 2, lower.tail = FALSE)
  multiplier <- if (pi_ci_scale == "logit") pi * (1 - pi) else 1
  list(variance = (working_se * multiplier)^2, working_se = working_se)
}

.ravel_pi_se_variance <- function(pi, pi_se) {
  scalar <- function(x) is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x)
  if (!scalar(pi) || pi < 0 || pi >= 1) stop("pi must be a finite scalar in [0, 1)")
  if (!scalar(pi_se) || pi_se < 0)
    stop("pi_se must be a finite nonnegative scalar on the raw probability scale")
  if (pi == 0 && pi_se > 0)
    stop("A positive pi_se at pi = 0 is not supported; supply a boundary-aware analysis")
  if (!is.finite(pi_se^2)) stop("pi_se is too large to produce a finite sampling variance")
  list(variance = as.numeric(pi_se)^2, working_se = as.numeric(pi_se), raw_se = as.numeric(pi_se))
}

.ravel_pi_validate_covariance <- function(V, labels, what) {
  if (!is.matrix(V) || !is.numeric(V) || is.complex(V) ||
      !identical(dim(V), rep(length(labels), 2L)) || any(!is.finite(V)) ||
      is.null(rownames(V)) || is.null(colnames(V)) ||
      anyDuplicated(rownames(V)) || anyDuplicated(colnames(V)) ||
      !setequal(rownames(V), labels) || !setequal(colnames(V), labels))
    stop(what, " must be a finite numeric matrix with unique core-parameter names")
  V <- V[labels, labels, drop = FALSE]
  tolerance <- 1e-8 * max(abs(V), .Machine$double.eps)
  if (max(abs(V - t(V))) > tolerance) stop(what, " must be symmetric")
  V <- (V + t(V)) / 2
  if (min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) < -tolerance)
    stop(what, " must be positive semidefinite")
  V
}

.ravel_pi_extract <- function(object, uncertainty, covariate_mapping_confirmed) {
  if (is.character(object) && length(object) == 1L) object <- readRDS(object)
  labels <- .report_core_labels
  sib_roles <- character(0)
  if (inherits(object, "ravel_results") && !is.null(object$parental_sibs)) {
    sib_roles <- object$parental_sibs$roles
    labels <- c(labels, .sib_core_labels(sib_roles))
  }
  if (inherits(object, "ravel_results")) {
    if (!isTRUE(object$converged) || is.null(object$fit) ||
        !suppressWarnings(isTRUE(lavaan::lavInspect(object$fit, "post.check"))))
      stop("ravel_results must contain a converged, admissible SEM fit")
    pt <- lavaan::parTable(object$fit)
    estimated_state_fit <- identical(object$constraint_spec$mode, "estimated_state") ||
      identical(object$fit_constraints, "estimated_state") ||
      identical(object$am_correction, "estimated_state") ||
      "am_inverse_vG_parent" %in% pt$label
    # Older live objects may predate the identification flag. Recheck the full
    # fitted model: a valid six-core submatrix cannot certify the measurement
    # parameters, and bootstrap covariance must not bypass identification.
    if (estimated_state_fit && !.ravel_estimated_state_fit_valid(object$fit))
      stop("Estimated-state correction requires a locally identified full SEM fit with valid covariance")
    theta <- c(object$beta_vis[c("bD", "bm", "bf")],
               object$latent_cov_vis[c("rhoMF", "eta_m", "eta_f")],
               .sib_theta(object$parental_sibs, sib_roles))
    V <- object$vcov_core
    if (is.null(V)) {
      full <- lavaan::vcov(object$fit)
      V <- matrix(0, length(labels), length(labels), dimnames = list(labels, labels))
      shared <- intersect(labels, rownames(full))
      V[shared, shared] <- full[shared, shared, drop = FALSE]
      # Missing labels may be fixed parameters, but must never silently become
      # zero uncertainty when they are estimated/defined by a constrained model.
      missing <- setdiff(labels, shared)
      if (any(pt$label %in% missing & (pt$free > 0L | pt$op == ":=")))
        stop("This fit needs its complete derived core covariance (vcov_core)")
    }
    packet <- list(theta = theta, vcov_core = V,
      diagnostics = list(estimated_state_identified = if (estimated_state_fit) TRUE else NULL),
      settings = list(uncertainty = if (is.null(object$bootstrap)) "lavaan" else "bootstrap",
        requires_covariate_mapping_review = length(object$fixed_effects) > 0L,
        residual_model = object$crossperson_structure,
        parental_sibs = object$parental_sibs$mode %||% "none",
        constraint_spec = object$constraint_spec),
      outcome = list(name = object$outcome_variable, type = object$outcome,
        scale = if (identical(object$outcome, "binary_liability")) "SEM_probit_liability" else "input_phenotype_units"),
      bootstrap = object$bootstrap)
  } else {
    if (!is.list(object) || !identical(object$schema, "ravel_visible_transfer_v1") ||
        !isTRUE(object$diagnostics$admissible))
      stop("object must be an admissible ravel_results or visible-transfer packet")
    if (!identical(object$software$transform_fingerprint, .report_transform_fingerprint()))
      stop("The correction code differs from the packet's; align package versions before proceeding")
    if (identical(object$settings$constraint_spec$mode, "estimated_state") &&
        !identical(object$diagnostics$estimated_state_identified, TRUE))
      stop("Estimated-state packets require a successful full-fit identification check; regenerate the packet from a validated fit")
    packet <- object
    sib_roles <- .parental_sib_roles(packet$settings$parental_sibs %||% "none")
    labels <- c(.report_core_labels, .sib_core_labels(sib_roles))
  }
  if (!is.numeric(packet$theta) || is.complex(packet$theta) ||
      is.null(names(packet$theta)) || anyDuplicated(names(packet$theta)) ||
      !all(labels %in% names(packet$theta)) || any(!is.finite(packet$theta[labels])))
    stop(if (length(labels) == 6L) "The object must contain all six finite, uniquely named core estimates" else
      "The object must contain all finite, uniquely named core estimates")
  packet$theta <- packet$theta[labels]
  if (uncertainty == "auto") uncertainty <- if (identical(packet$settings$uncertainty, "bootstrap"))
    "bootstrap" else "lavaan"
  n_success <- NA_integer_
  if (uncertainty == "bootstrap") {
    draws <- packet$bootstrap$core_draws
    if (!is.matrix(draws) || !is.numeric(draws) || is.complex(draws) ||
        is.null(colnames(draws)) || anyDuplicated(colnames(draws)) ||
        !setequal(colnames(draws), labels))
      stop("Bootstrap uncertainty requires stored joint ", if (length(labels) == 6L) "six-parameter" else
        "core-parameter", " core_draws; use uncertainty = 'lavaan' explicitly for model-based inference")
    draws <- draws[, labels, drop = FALSE]
    complete <- apply(is.finite(draws), 1L, all)
    n_success <- sum(complete)
    stored_minimum <- packet$bootstrap$minimum_success %||% 0L
    if (!is.numeric(stored_minimum) || length(stored_minimum) != 1L ||
        !is.finite(stored_minimum) || stored_minimum < 0)
      stop("The stored bootstrap minimum_success must be finite and nonnegative")
    minimum <- max(10L, ceiling(.8 * nrow(draws)), stored_minimum)
    if (n_success < minimum) stop("Insufficient complete joint core bootstrap draws")
    packet$vcov_core <- stats::cov(draws[complete, , drop = FALSE])
  }
  packet$vcov_core <- .ravel_pi_validate_covariance(packet$vcov_core, labels, "Core covariance")
  list(packet = packet, uncertainty = uncertainty, n_success = n_success,
       live_result = if (inherits(object, "ravel_results")) object else NULL,
       labels = labels, sib_roles = sib_roles)
}

#' Correct RAVEL Estimates with Uncertainty in an External Hidden Fraction
#'
#' @export
ravel_correct_pi <- function(object, pi, pi_ci = NULL,
                             pi_ci_scale = c("normal", "logit"), pi_ci_level = .95,
                             am_correction = c("n_gen", "equilibrium", "estimated_state"),
                             am_n_generations = 15L, ci_level = .95,
                             uncertainty = c("auto", "lavaan", "bootstrap"),
                             core_pi_cov = NULL, core_pi_sensitivity = NULL,
                             covariate_mapping_confirmed = FALSE, pi_se = NULL) {
  if (is.null(pi_ci) == is.null(pi_se))
    stop("Supply exactly one of pi_se or pi_ci")
  external_source <- if (is.null(pi_se)) "ci" else "se"
  if (external_source == "se" && (!missing(pi_ci_scale) || !missing(pi_ci_level)))
    stop("pi_ci_scale and pi_ci_level apply only to pi_ci; pi_se is always on the raw probability scale")
  if (is.character(object) && length(object) == 1L) object <- readRDS(object)
  if (is.list(object)) {
    history <- if (inherits(object, "ravel_results"))
      list(mode = object$am_correction, n_generations = object$am_n_generations) else
      object$settings$constraint_spec
    if (missing(am_correction) && !is.null(history$mode) && history$mode %in% c("n_gen", "equilibrium", "estimated_state"))
      am_correction <- history$mode
    if (missing(am_n_generations) && identical(history$mode, "n_gen") && !is.null(history$n_generations))
      am_n_generations <- history$n_generations
  }
  pi_ci_scale <- match.arg(pi_ci_scale)
  am_correction <- match.arg(am_correction)
  uncertainty <- match.arg(uncertainty)
  external <- if (external_source == "se") .ravel_pi_se_variance(pi, pi_se) else
    .ravel_pi_ci_variance(pi, pi_ci, pi_ci_scale, pi_ci_level)
  if (external_source == "ci") external$raw_se <- sqrt(external$variance)
  if (!is.numeric(ci_level) || is.complex(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1)
    stop("ci_level must lie strictly between zero and one")
  if (am_correction == "n_gen" && (!is.numeric(am_n_generations) ||
      is.complex(am_n_generations) || length(am_n_generations) != 1L ||
      !is.finite(am_n_generations) || am_n_generations < 1 ||
      am_n_generations > .Machine$integer.max || am_n_generations != as.integer(am_n_generations)))
    stop("am_n_generations must be an integer of at least one")
  input <- .ravel_pi_extract(object, uncertainty, covariate_mapping_confirmed)
  packet <- input$packet
  labels <- input$labels
  sib_roles <- input$sib_roles
  n_core <- length(labels)
  # Finite-generation fits may use pi in their measurement model. The total
  # derivative includes the fitted core's response to changing that input.
  spec <- packet$settings$constraint_spec
  constrained <- !is.null(spec$mode) && !identical(spec$mode, "none")
  depends <- constrained && identical(spec$mode, "n_gen") && spec$n_generations > 1L
  if (am_correction == "estimated_state" &&
      (!constrained || !identical(spec$mode, "estimated_state")))
    stop("Estimated-state correction requires a matching fit_constraints = 'estimated_state' fit or packet")
  if (constrained && (!identical(am_correction, spec$mode) ||
      (am_correction == "n_gen" && am_n_generations != spec$n_generations)))
    stop("The correction's mating history must match the constrained SEM fit")
  if (depends && (is.null(spec$h) || abs(pi - spec$h) > 1e-12))
    stop("Refit the finite-generation constrained SEM at the supplied pi before correcting")
  if (am_correction == "estimated_state" &&
      abs(packet$theta["eta_m"] - packet$theta["eta_f"]) > 1e-5)
    stop("Estimated-state correction requires equal maternal and paternal latent parent-child correlations")
  sensitivity <- setNames(numeric(n_core), labels)
  if (depends && external$variance > 0 && is.null(core_pi_sensitivity)) {
    if (is.null(input$live_result))
      stop("This pi-dependent packet requires core_pi_sensitivity from refitting the SEM at perturbed pi")
    core_pi_sensitivity <- .ravel_constraint_core_h_sensitivity(input$live_result)
  }
  if (!is.null(core_pi_sensitivity)) {
    if (!is.numeric(core_pi_sensitivity) || is.complex(core_pi_sensitivity) ||
        length(core_pi_sensitivity) != n_core || any(!is.finite(core_pi_sensitivity)) ||
        is.null(names(core_pi_sensitivity)) || anyDuplicated(names(core_pi_sensitivity)) ||
        !setequal(names(core_pi_sensitivity), labels))
      stop(if (n_core == 6L) "core_pi_sensitivity must be a finite vector named for all six core parameters" else
        paste0("core_pi_sensitivity must be a finite vector named for all core parameters: ",
               paste(labels, collapse = ", ")))
    if (!depends && any(core_pi_sensitivity != 0))
      stop("Nonzero core_pi_sensitivity applies only to a pi-dependent constrained fit")
    sensitivity <- core_pi_sensitivity[labels]
  }
  cross <- setNames(numeric(n_core), labels)
  independent <- is.null(core_pi_cov)
  .ravel_pi_validate_cross(core_pi_cov, external$variance, labels)
  if (!independent) cross <- core_pi_cov[labels]
  joint_labels <- c(labels, "pi")
  V <- matrix(0, n_core + 1L, n_core + 1L, dimnames = list(joint_labels, joint_labels))
  V[labels, labels] <- packet$vcov_core
  V["pi", "pi"] <- external$variance
  V[labels, "pi"] <- V["pi", labels] <- cross
  V <- .ravel_pi_validate_covariance(V, joint_labels, "Joint core/pi covariance")
  x <- c(packet$theta, pi = pi)
  sib_map <- .parental_sib_coef_map(sib_roles)
  sib_map_names <- names(sib_map)
  trio_std <- c("beta_D", "beta_M", "beta_P", "rho_corrected")
  trio_raw <- c("beta_D", "beta_M", "beta_P", "cov_corrected")
  row_scale <- c(rep("corrected_std", 4L), rep("corrected_raw", 4L),
                 rep("corrected_std", length(sib_map)), rep("corrected_raw", length(sib_map)))
  row_parameter <- c(trio_std, trio_raw, unname(sib_map), unname(sib_map))
  row_names <- paste0(row_scale, ":", row_parameter)
  evaluate <- function(value) {
    names(value) <- joint_labels
    S <- matrix(c(1, value["eta_m"], value["eta_f"],
                  value["eta_m"], 1, value["rhoMF"],
                  value["eta_f"], value["rhoMF"], 1), 3L, 3L)
    if (min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) <= 0)
      stop("The fitted latent trio covariance must be positive definite")
    result <- .compute_scale_estimates(value[c("bD", "bm", "bf")],
      value["eta_m"], value["eta_f"], value["rhoMF"], value["pi"],
      am_correction, am_n_generations, sib = .sib_input_from_theta(value, sib_roles))
    answer <- c(result$beta_total_assumed_std, rho_corrected = result$rho_total_assumed,
                result$beta_total_assumed, cov_corrected = result$cov_total_assumed,
                result$sib_total$std[sib_map_names], result$sib_total$raw[sib_map_names])
    names(answer) <- row_names
    if (any(!is.finite(answer)))
      stop("This pi/AM assumption is incompatible with the fitted moments")
    answer
  }
  estimate <- evaluate(x)
  J <- matrix(0, length(estimate), n_core + 1L, dimnames = list(names(estimate), joint_labels))
  pi_index <- n_core + 1L
  for (j in seq_along(x)) {
    if (j == pi_index && external$variance == 0 && all(cross == 0)) next
    step <- 1e-6 * max(1, abs(x[j]))
    if (j == pi_index && (pi < step || 1 - pi < step)) {
      direction <- if (pi < step) 1 else -1
      step <- min(step, if (direction == 1) (1 - pi) / 3 else pi / 3)
      near <- far <- x
      near[j] <- near[j] + direction * step
      far[j] <- far[j] + direction * 2 * step
      J[, j] <- direction * (-3 * estimate + 4 * evaluate(near) - evaluate(far)) / (2 * step)
      next
    }
    upper <- lower <- x
    upper[j] <- upper[j] + step
    lower[j] <- lower[j] - step
    J[, j] <- (evaluate(upper) - evaluate(lower)) / (2 * step)
  }
  if (any(!is.finite(J))) stop("Correction derivatives unavailable near an admissibility boundary")
  partial_J <- J
  J[, "pi"] <- J[, "pi"] + drop(J[, labels, drop = FALSE] %*% sensitivity)
  covariance <- J %*% V %*% t(J)
  covariance <- (covariance + t(covariance)) / 2
  if (any(!is.finite(covariance)) || any(diag(covariance) < -1e-12))
    stop("Invalid propagated correction covariance")
  se <- sqrt(pmax(0, diag(covariance)))
  z <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
  statistic <- ifelse(se == 0, ifelse(estimate == 0, 0, Inf), abs(estimate / se))
  out <- data.frame(outcome_variable = packet$outcome$name %||% "Y",
    outcome_scale = packet$outcome$scale %||% if (identical(packet$outcome$type, "binary_liability"))
      "SEM_probit_liability" else "prepared_phenotype_SD",
    residual_model = packet$settings$residual_model %||% "unspecified",
    scale = row_scale,
    parameter = row_parameter,
    estimate = unname(estimate), se = unname(se), ci_low = unname(estimate - z * se),
    ci_high = unname(estimate + z * se), p_value = unname(2 * stats::pnorm(-statistic)),
    inference = if (input$uncertainty == "bootstrap") "bootstrap_covariance_pi_delta_wald" else "lavaan_pi_delta_wald",
    pi = pi, pi_se = external$raw_se, pi_uncertainty_source = external_source,
    pi_ci_scale = if (external_source == "ci") pi_ci_scale else NA_character_,
    am_correction = am_correction,
    am_n_generations = if (am_correction == "n_gen") as.integer(am_n_generations) else NA_integer_,
    n_boot_success = input$n_success, uncertainty_status = "ok", stringsAsFactors = FALSE)
  attr(out, "covariance") <- covariance
  attr(out, "ci_level") <- ci_level
  attr(out, "constraint_spec") <- packet$settings$constraint_spec
  attr(out, "jacobian") <- J
  attr(out, "partial_jacobian") <- partial_J
  attr(out, "joint_input_covariance") <- V
  attr(out, "variance_components") <- list(
    core = J[, labels, drop = FALSE] %*% packet$vcov_core %*% t(J[, labels, drop = FALSE]),
    pi = tcrossprod(J[, "pi"]) * external$variance,
    cross = tcrossprod(drop(J[, labels, drop = FALSE] %*% cross), J[, "pi"]) +
            tcrossprod(J[, "pi"], drop(J[, labels, drop = FALSE] %*% cross)))
  attr(out, "pi_uncertainty") <- list(pi = pi, source = external_source,
    ci = if (external_source == "ci") pi_ci else NULL,
    ci_level = if (external_source == "ci") pi_ci_level else NULL,
    ci_scale = if (external_source == "ci") pi_ci_scale else NULL,
    working_scale = if (external_source == "ci" && pi_ci_scale == "logit") "logit" else "raw probability",
    working_scale_se = external$working_se,
    raw_scale_se = external$raw_se, independent_of_core = independent,
    cross_covariance_scale = "raw pi", method = "first-order joint delta method",
    core_pi_sensitivity = sensitivity,
    core_covariance_conditioning = "core estimator evaluated at fixed pi",
    correction_ci_level = ci_level, core_uncertainty = input$uncertainty,
    fixed_am_history = am_correction != "estimated_state", boundary_truncation = FALSE)
  if (am_correction == "estimated_state") {
    state <- .hidden_estimated_state(packet$theta["rhoMF"], packet$theta["eta_m"],
      packet$theta["eta_f"], pi)
    if (is.null(state)) stop("The inferred assortment state is incompatible with pi")
    state_J6 <- .hidden_estimated_state_jacobian(packet$theta["rhoMF"],
      packet$theta["eta_m"], packet$theta["eta_f"], pi)
    state_J <- matrix(0, nrow(state_J6), n_core + 1L, dimnames = list(rownames(state_J6), joint_labels))
    state_J[, colnames(state_J6)] <- state_J6
    state_covariance <- state_J %*% V %*% t(state_J)
    state_covariance <- (state_covariance + t(state_covariance)) / 2
    if (any(!is.finite(state_covariance)) || any(diag(state_covariance) < -1e-12))
      stop("Invalid propagated assortment-state covariance")
    state_estimate <- c(state$estimated_state$B, state$vG,
      state$estimated_state$vTp, state$rho_T)
    state_se <- sqrt(pmax(0, diag(state_covariance)))
    attr(out, "estimated_state") <- c(state$estimated_state, list(
      summary = data.frame(parameter = rownames(state_J), estimate = unname(state_estimate),
        se = unname(state_se), ci_low = unname(state_estimate - z * state_se),
        ci_high = unname(state_estimate + z * state_se), stringsAsFactors = FALSE),
      covariance = state_covariance, jacobian = state_J))
  }
  out
}

# Validate the optional external-estimate covariance before expensive fitting.
.ravel_pi_validate_cross <- function(core_pi_cov, variance, labels = .report_core_labels) {
  if (is.null(core_pi_cov)) return(invisible(NULL))
  if (!is.numeric(core_pi_cov) || is.complex(core_pi_cov) || length(core_pi_cov) != length(labels) ||
      any(!is.finite(core_pi_cov)) || is.null(names(core_pi_cov)) ||
      anyDuplicated(names(core_pi_cov)) || !setequal(names(core_pi_cov), labels))
    stop(if (length(labels) == 6L) "core_pi_cov must be a finite vector named for all six core parameters, on the raw pi scale" else
      paste0("core_pi_cov must be a finite vector named for all core parameters (",
             paste(labels, collapse = ", "), "), on the raw pi scale"))
  if (variance == 0 && any(core_pi_cov != 0))
    stop("A zero-variance pi cannot have nonzero core_pi_cov")
  invisible(NULL)
}

# Build display blocks without altering the fixed-pi covariance used for fitting,
# bootstrap packets, sensitivity refits, or subsequent external-pi corrections.
.ravel_pi_merge_uncertainty <- function(result, conditional, source) {
  if (!is.null(result$pi_correction_error)) {
    # A failed requested external-pi correction must not fall back to fixed-pi
    # inference. The visible model's uncertainty remains available.
    for (name in c("total_std", "total_raw", "sib_total_std", "sib_total_raw")) {
      if (!is.null(conditional[[name]])) for (field in c("se", "ci", "cov"))
        if (!is.null(conditional[[name]][[field]])) conditional[[name]][[field]][] <- NA_real_
    }
    conditional["pi_correction"] <- list(NULL)
    conditional$pi_correction_error <- result$pi_correction_error
    return(conditional)
  }
  corrected <- result$pi_correction
  if (is.null(corrected)) return(conditional)
  meta <- attr(corrected, "pi_uncertainty")
  if (!identical(meta$core_uncertainty, source)) {
    joint <- attr(corrected, "joint_input_covariance")
    cross <- if (isTRUE(meta$independent_of_core)) NULL else
      joint[setdiff(rownames(joint), "pi"), "pi"]
    corrected <- ravel_correct_pi(result, pi = meta$pi, pi_se = meta$raw_scale_se,
      ci_level = meta$correction_ci_level, uncertainty = source,
      core_pi_cov = cross, core_pi_sensitivity = meta$core_pi_sensitivity,
      covariate_mapping_confirmed = TRUE)
  }
  block <- function(scale, parameters, labels) {
    rows <- which(corrected$scale == scale & corrected$parameter %in% parameters)
    part <- corrected[rows, , drop = FALSE]
    labels <- labels[match(part$parameter, parameters)]
    covariance <- attr(corrected, "covariance")[rows, rows, drop = FALSE]
    dimnames(covariance) <- list(labels, labels)
    list(estimate = setNames(part$estimate, labels), se = setNames(part$se, labels),
      ci = rbind(lower = setNames(part$ci_low, labels), upper = setNames(part$ci_high, labels)),
      cov = covariance, method = unique(part$inference))
  }
  trio_labels <- c("bD", "bm", "bf", "rhoMF")
  conditional$total_std <- block("corrected_std", c("beta_D", "beta_M", "beta_P", "rho_corrected"), trio_labels)
  conditional$total_raw <- block("corrected_raw", c("beta_D", "beta_M", "beta_P", "cov_corrected"), trio_labels)
  if (!is.null(result$parental_sibs)) {
    map <- .parental_sib_coef_map(result$parental_sibs$roles)
    conditional$sib_total_std <- block("corrected_std", unname(map), names(map))
    conditional$sib_total_raw <- block("corrected_raw", unname(map), names(map))
  }
  conditional$pi_correction <- corrected
  conditional
}
