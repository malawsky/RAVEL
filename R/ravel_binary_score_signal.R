# Internal exploratory binary score-signal calculation, retained for research
# and reproducible calibration. Its nominal Wald p-values failed null-size
# calibration. ravel_score_diagnostic() reports them as nominal results with
# an explicit calibration warning; the continuous comparison is unchanged.

.ravel_binary_score_map <- function(theta, pt, scores, latents) {
  L <- matrix(0, length(scores), length(latents), dimnames = list(scores, latents))
  P <- matrix(0, length(latents), length(latents), dimnames = list(latents, latents))
  E <- matrix(0, length(scores), length(scores), dimnames = list(scores, scores))
  value <- pt$est
  free <- which(pt$free > 0L)
  value[free] <- theta[pt$free[free]]
  for (i in seq_len(nrow(pt))) {
    a <- pt$lhs[i]; b <- pt$rhs[i]
    if (pt$op[i] == "=~" && a %in% latents && b %in% scores) L[b, a] <- value[i]
    if (pt$op[i] == "~~" && a %in% latents && b %in% latents) P[a, b] <- P[b, a] <- value[i]
    if (pt$op[i] == "~~" && a %in% scores && b %in% scores) E[a, b] <- E[b, a] <- value[i]
  }
  S <- L %*% P %*% t(L) + E
  A <- solve(S, L %*% P)
  list(A = A, conditional = P - P %*% t(L) %*% A)
}

.ravel_test_binary_score_signal <- function(result, data = NULL) {
  out <- structure(list(status = "not_run", reason = "", statistic = NA_real_,
    df = NA_real_, p_value = NA_real_, raw_statistic = NA_real_,
    test = "Stacked-influence Wald test of unadjusted probit score-signal restrictions",
    N = NA_integer_, roles = character(),
    diagnostics = list(calibrated = FALSE, inference_status = "exploratory_nominal"),
    fits = list()),
    class = "ravel_score_signal_test")
  finish <- function(status, reason) { out$status <- status; out$reason <- reason; out }
  tryCatch({
    if (!is.character(result$model) || length(result$model) != 1L ||
        (!is.null(result$outcome) && !identical(result$outcome, "binary_liability")))
      return(finish("needs_fit", "An unadjusted binary RAVEL model specification is required."))
    original <- if (inherits(result$fit, "lavaan")) result$fit else NULL
    if (length(result$fixed_effects) || length(result$interactions) ||
        isTRUE(result$score_signal_input$dependent_families) ||
        (!is.null(original) && (length(original@Data@cluster) || original@Data@ngroups != 1L)))
      return(finish("unsupported_design", "This auxiliary diagnostic requires independent unadjusted binary observations."))
    if (identical(result$constraint_spec$mode, "none") &&
        !length(result[["fixed_crossperson"]]) && !length(result[["fixed_crossperson_corr"]]))
      return(finish("unsupported_free_measurement", "The unrestricted score-only measurement model is unidentified."))
    .ensure_lavaan_runtime()
    if (is.null(data)) {
      if (is.null(original)) return(finish("needs_data", "Raw data are required when the auxiliary outcome fit is unavailable."))
      data <- as.data.frame(lavaan::lavInspect(original, "data"))
    } else {
      data <- as.data.frame(data)
      if ("FID" %in% names(data) && (anyNA(data$FID) || anyDuplicated(data$FID)))
        return(finish("unsupported_design", "Supplied observations contain repeated or missing extended-family identifiers."))
      mapping <- result$score_signal_input$column_map
      for (name in names(mapping)) if (mapping[[name]] %in% names(data)) data[[name]] <- data[[mapping[[name]]]]
    }
    template <- lavaan::lavaanify(result$model, std.lv = TRUE)
    latents <- unique(template$lhs[template$op == "=~"])
    scores <- unique(template$rhs[template$op == "=~"])
    paths <- template[template$lhs == "Y" & template$op == "~", , drop = FALSE]
    if (!setequal(paths$rhs, latents) || any(paths$free == 0L))
      return(finish("unsupported_fixed_paths", "Every latent outcome path must be free in the shared-factor null; fixed paths are not silently released."))
    if (length(scores) != 2L * length(latents) || !all(c(scores, "Y") %in% names(data)))
      return(finish("unsupported_model", "Exactly two scores per genetic role and a binary outcome are required."))
    data <- data[c(scores, "Y")]
    if (any(!stats::complete.cases(data)))
      return(finish("unsupported_missingness", "Complete scores and outcomes are required."))
    if (any(!vapply(data, is.numeric, logical(1))) || any(!is.finite(as.matrix(data))))
      return(finish("invalid_data", "All score and outcome entries must be finite numeric values."))
    yy <- sort(unique(data$Y))
    if (length(yy) != 2L) return(finish("invalid_data", "Both binary outcome categories are required."))
    data$Y <- as.numeric(data$Y == yy[2L])
    n <- nrow(data); out$N <- n; out$roles <- latents
    # Canonical units prevent arbitrary PGS scales from changing the numerical
    # rank, finite-difference derivatives, or the chosen Wald restriction basis.
    # Reconstruct only a verified generated model, converting fixed covariances
    # and absolute bounds as in the continuous diagnostic.
    if (is.null(result$constraint_spec$mode))
      return(finish("needs_fit_settings", "Genetic constraint metadata are required, including when the auxiliary outcome fit failed."))
    # ravel_fit's failure record may omit its normalized residual settings.
    # Reapply the same validated preparation from the explicit genetic contract.
    # This restores prescribed constrained defaults, and rejects conflicting
    # custom values rather than silently replacing them.
    role_codes <- .parental_sib_roles(result$parental_sibs$mode %||% "none")
    prepared <- .ravel_prepare_constraints(result$constraint_spec$mode,
      result$constraint_spec$h, result$constraint_spec$mode,
      result$constraint_spec$n_generations, TRUE,
      result$crossperson_structure %||% "free",
      .normalize_fixed_crossperson(result$crossperson_structure %||% "free", result[["fixed_crossperson"]], sib_roles = role_codes),
      .normalize_fixed_crossperson(result$crossperson_structure %||% "free", result[["fixed_crossperson_corr"]],
        arg_name = "fixed_crossperson_corr", sib_roles = role_codes), sib_roles = role_codes)
    result$fixed_crossperson <- prepared$fixed_cov[!is.na(prepared$fixed_cov)]
    result$fixed_crossperson_corr <- prepared$fixed_corr[!is.na(prepared$fixed_corr)]
    result$crossperson_structure <- prepared$crossperson
    reconstructed <- .ravel_score_signal_model(result)
    if (!identical(.ravel_score_signal_model_key(result$model), .ravel_score_signal_model_key(reconstructed)))
      return(finish("unsupported_model", "The supplied model does not match its genetic settings; score restrictions cannot be safely rescaled."))
    units <- vapply(data[scores], stats::sd, numeric(1))
    centers <- vapply(data[scores], mean, numeric(1))
    if (any(!is.finite(units) | units <= 0))
      return(finish("invalid_data", "Each score needs positive finite variance."))
    data[scores] <- lapply(scores, function(name) (data[[name]] - centers[[name]]) / units[[name]])
    out$diagnostics$input_standard_deviations <- units
    out$diagnostics$input_score_means <- centers
    out$diagnostics$parameter_units <- "Scores are centred and divided by their sample SDs for this auxiliary diagnostic."
    out$diagnostics$numerical_bounds <- "The source model's absolute numerical bounds are transformed with the score units; a source bound that is active after extreme rescaling is not silently relaxed."
    model <- .ravel_score_signal_model(result, units)
    mode <- result$constraint_spec$mode
    if (mode != "none") model <- .ravel_constraint_starts(model, data, mode,
      result$constraint_spec$h, result$constraint_spec$n_generations,
      .parental_sib_roles(result$parental_sibs$mode %||% "none"))
    score_model <- .ravel_score_signal_syntax(model, scores, latents, "score")
    reference <- lavaan::sem(score_model, data = data[scores],
      std.lv = TRUE, fixed.x = FALSE, conditional.x = FALSE,
      estimator = "MLR", missing = "fiml", meanstructure = TRUE, do.fit = FALSE)
    score <- .ravel_score_signal_fit_checked(score_model, data[scores], reference, "MLR")
    out$diagnostics$score_warnings <- score$warnings
    if (is.null(score$value)) return(finish("fit_error", score$error %||% "Score-only fit failed."))
    sf <- score$value; out$fits$score <- sf
    check <- .ravel_score_signal_validate_fit(sf)
    out$diagnostics$score <- check
    if (!isTRUE(check$valid)) return(finish("invalid_score_model", check$reason))
    if (length(score$warnings)) return(finish("fit_warning", paste(score$warnings, collapse = "; ")))
    X <- cbind(intercept = 1, as.matrix(data[scores]))
    glm_warnings <- character()
    g <- withCallingHandlers(stats::glm.fit(X, data$Y, family = stats::binomial("probit")),
      warning = function(w) { glm_warnings <<- c(glm_warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    out$diagnostics$probit_warnings <- glm_warnings
    if (!isTRUE(g$converged) || g$rank != ncol(X) || any(!is.finite(g$coefficients)) || length(glm_warnings))
      return(finish("invalid_probit_fit", paste(c("The unrestricted probit fit was not regular.", glm_warnings), collapse = " ")))
    eta <- as.vector(X %*% g$coefficients); pr <- stats::pnorm(eta); den <- stats::dnorm(eta)
    if (any(pr <= 1e-10 | pr >= 1 - 1e-10))
      return(finish("invalid_probit_fit", "Extreme fitted probabilities prevent stable probit influence calculations."))
    v <- pr * (1 - pr); residual <- data$Y - pr
    Gscores <- X * as.vector(den * residual / v)
    # Observed rather than Fisher information: differentiate each probit score.
    w <- den^2 / v + eta * den * residual / v + den^2 * residual * (1 - 2 * pr) / v^2
    Ig <- crossprod(X, X * as.vector(w)) / n
    if (min(eigen(Ig, symmetric = TRUE, only.values = TRUE)$values) <= 0)
      return(finish("invalid_probit_fit", "The probit information matrix is not positive definite."))
    theta <- lavaan::coef(sf); p <- length(theta)
    tangent <- sweep(check$rank$tangent, 1, check$rank$parameter_scale, "*")
    information <- lavaan::lavTech(sf, "information")
    B <- tangent %*% solve(crossprod(tangent, information %*% tangent)) %*% t(tangent)
    U <- lavaan::lavScores(sf, ignore.constraints = TRUE, remove.duplicated = FALSE)
    if (!identical(colnames(U), names(theta)) || nrow(U) != n)
      return(finish("invalid_influence", "The measurement influence rows or parameter ordering did not match."))
    If <- U %*% B
    expected <- lavaan::vcov(sf)
    discrepancy <- max(abs(crossprod(If) / n^2 - expected)) / max(abs(expected))
    out$diagnostics$measurement_vcov_relative_error <- discrepancy
    if (!is.finite(discrepancy) || discrepancy > 1e-5)
      return(finish("invalid_influence", "The reconstructed measurement sandwich did not reproduce lavaan uncertainty."))
    influence <- cbind(If, Gscores %*% solve(Ig))
    influence <- sweep(influence, 2, colMeans(influence))
    joint <- crossprod(influence) / n^2
    pt <- lavaan::parameterTable(sf)
    map <- .ravel_binary_score_map(theta, pt, scores, latents)
    k <- length(latents)
    singular <- svd(map$A, nu = 0L, nv = 0L)$d
    out$diagnostics$prediction_singular_values <- singular
    if (length(singular) != k || any(!is.finite(singular)) ||
        min(singular) <= max(singular) * 1e-8)
      return(finish("invalid_score_model", "The genetic prediction matrix is rank deficient or numerically singular."))
    # The two scores for a role can imply almost proportional prediction rows.
    # Default LINPACK QR can retain those near-duplicate first rows and produce
    # an unstable restriction basis. LAPACK uses norm pivoting; its reported
    # rank is not informative, so rank is checked independently above.
    qrA <- qr(t(map$A), LAPACK = TRUE)
    selected <- qrA$pivot[seq_len(k)]; other <- setdiff(seq_along(scores), selected)
    out$diagnostics$restriction_basis_scores <- scores[selected]
    out$diagnostics$restriction_basis_condition <- kappa(map$A[selected, , drop = FALSE])
    q <- g$coefficients[-1L]
    d <- solve(crossprod(map$A), crossprod(map$A, q))
    slack <- 1 - drop(crossprod(d, map$conditional %*% d))
    out$diagnostics$null_conditional_variance_slack <- slack
    if (!is.finite(slack) || slack <= 1e-7)
      return(finish("invalid_null_scale", "Projected score associations imply a nonpositive residual variance under the genetic-factor null."))
    transform <- function(z) {
      a <- .ravel_binary_score_map(z[seq_len(p)], pt, scores, latents)$A
      slope <- z[p + 1L + seq_along(scores)]
      as.vector(slope[other] - a[other, , drop = FALSE] %*% solve(a[selected, , drop = FALSE], slope[selected]))
    }
    z <- c(theta, g$coefficients)
    h <- transform(z); J <- .finite_difference_jacobian(transform, z)
    V <- J %*% joint %*% t(J); V <- (V + t(V)) / 2
    ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(ev)) || min(ev) <= max(ev) * 1e-8)
      return(finish("invalid_influence", "The covariance of the tested restrictions is singular or indefinite."))
    out$statistic <- max(0, drop(crossprod(h, solve(V, h))))
    out$raw_statistic <- out$statistic; out$df <- length(h)
    out$p_value <- stats::pchisq(out$statistic, out$df, lower.tail = FALSE)
    out$diagnostics$restrictions <- h; out$diagnostics$restriction_vcov <- V
    out$diagnostics$joint_parameter_vcov <- joint
    out$diagnostics$probit_coefficients <- g$coefficients
    out$diagnostics$conditional_on_pi <- isTRUE(result$constraint_spec$depends_on_pi)
    out$diagnostics$scope <- "Independent complete binary observations; Gaussian score measurement model; unadjusted conditional-probit association restrictions. This is not the continuous likelihood-ratio test."
    finish("ok", "")
  }, error = function(e) finish("error", conditionMessage(e)))
}
