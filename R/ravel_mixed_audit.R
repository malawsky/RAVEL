# Optional extraction only: no estimation or change to the fitted model.

.ravel_mixed_audit_values <- function(par, spec, n_children) {
  g <- .ravel_mixed_genetic(n_children, par, spec)
  s <- g$state
  labels <- c("mother", "father", paste0("child", seq_len(n_children)))
  sd_total <- sqrt(c(s$parent, s$parent, rep(s$child, n_children)))
  total <- g$T * outer(sd_total, sd_total)
  neutral <- matrix(.5, length(labels), length(labels))
  diag(neutral) <- 1; neutral[1, 2] <- neutral[2, 1] <- 0
  visible <- (1 - spec$pi) * total + spec$pi * neutral
  matrices <- list(total_founder_covariance = total,
    total_correlation = g$T, visible_founder_covariance = visible,
    visible_correlation = stats::cov2cor(visible))
  entries <- lapply(names(matrices), function(kind) {
    m <- matrices[[kind]]
    ix <- which(upper.tri(m, diag = TRUE), arr.ind = TRUE)
    stats::setNames(m[ix], paste(kind, labels[ix[, 1]], labels[ix[, 2]], sep = "|"))
  })
  c(stats::setNames(unlist(s), paste0("state|", names(s), "|")),
    stats::setNames(unlist(entries, use.names = FALSE), unlist(lapply(entries, names))))
}

#' Extract detailed tables from a mixed family fit
#'
#' State and genetic covariance uncertainty is conditional on the supplied pi.
#' This does not alter the primary regression tables or refit the model.
#' @export
ravel_mixed_audit <- function(fit, n_children = 2L, design_group = 1L) {
  if (!inherits(fit, "ravel_mixed_results")) stop("fit must be a ravel_mixed_fit result.")
  if (!is.numeric(n_children) || length(n_children) != 1L || !is.finite(n_children) ||
      n_children < 1 || n_children > 20 || n_children != as.integer(n_children))
    stop("n_children must be an integer from 1 to 20 for the reference full-sibling family.")
  par <- fit$parameters
  transformed <- .ravel_mixed_delta(function(x)
    .ravel_mixed_audit_values(x, fit$settings, as.integer(n_children)), par)
  V <- fit$vcov
  available <- !is.null(V) && all(names(par) %in% rownames(V))
  covariance <- matrix(NA_real_, length(transformed$value), length(transformed$value),
    dimnames = list(names(transformed$value), names(transformed$value)))
  if (available) covariance <- transformed$jacobian %*%
    V[names(par), names(par), drop = FALSE] %*% t(transformed$jacobian)
  components <- strsplit(names(transformed$value), "|", fixed = TRUE)
  tab <- data.frame(block = vapply(components, `[`, character(1), 1L),
    lhs = vapply(components, `[`, character(1), 2L),
    rhs = vapply(components, function(x) if (length(x) > 2L) x[3L] else "", character(1)),
    estimate = unname(transformed$value), se = sqrt(pmax(0, diag(covariance))),
    uncertainty_conditions = "FID sandwich; conditional on supplied pi", row.names = NULL)
  tab$ci_low <- tab$estimate - stats::qnorm(.975) * tab$se
  tab$ci_high <- tab$estimate + stats::qnorm(.975) * tab$se
  all_parameters <- c(par, fit$means[setdiff(names(fit$means), names(par))])
  parameter_se <- rep(NA_real_, length(all_parameters))
  if (!is.null(V)) {
    take <- match(names(all_parameters), rownames(V))
    parameter_se[!is.na(take)] <- sqrt(pmax(0, diag(V)[take[!is.na(take)]]))
  }
  parameter_table <- data.frame(parameter = names(all_parameters),
    estimate = unname(all_parameters), se = parameter_se,
    scale = "internal estimation parameterization",
    uncertainty_conditions = "FID sandwich; conditional on supplied pi", row.names = NULL)
  implied <- .ravel_mixed_audit_implied(fit, design_group)
  list(fit_statistics = data.frame(statistic = names(fit$fit_statistics),
      value = unname(fit$fit_statistics),
      basis = "Internally standardized data; maternal-family working-independence likelihood"),
    counts = data.frame(outcomes = fit$N, individuals = fit$n_subjects,
      maternal_families = fit$n_families, extended_families = fit$n_clusters),
    status = data.frame(converged = fit$converged, optimizer_converged = fit$optimizer_converged,
      inference_status = fit$inference_status),
    state = tab[tab$block == "state", , drop = FALSE],
    genetic_covariances = tab[tab$block != "state", , drop = FALSE],
    state_and_genetic_vcov = covariance, parameters = parameter_table,
    coefficients = fit$table, coefficient_vcov = fit$coefficient_vcov,
    variance_components = fit$variance_components, variance_vcov = fit$variance_vcov,
    parameter_vcov = fit$vcov, model_parameter_vcov = fit$model_vcov,
    one_factor = fit$one_factor, reference_n_children = as.integer(n_children),
    covariate_scaling = data.frame(covariate = names(fit$settings$centers),
      center = unname(fit$settings$centers), scale = unname(fit$settings$scales)),
    score_scaling = fit$score_scaling,
    implied_covariances = implied$table, implied_covariance_vcov = implied$vcov,
    implied_outcome_design = implied$design, implied_status = implied$status,
    notes = c("Genetic tables describe an exchangeable full-sibling reference family, not a particular observed family.",
      "State, genetic covariance and internal parameter SEs condition on pi, even if pi_se was supplied. Coefficient and residual variance tables retain the fit's external-pi adjustment.",
      "Founder covariance uses founder genetic units; correlations use present-day role SDs. Residual random-effect variances are in the separate variance_components table.",
      "Implied covariances describe one retained covariance-design group, in input score and outcome units (binary Y denotes the conditional-probit latent response). They include the genetic and residual random-effect contributions and condition on pi.",
      "Likelihood and AIC use internal numerical units. Compare only fits of the same observations with the same scaling; with dependent maternal families within FID, these are working-independence quantities."))
}

.ravel_mixed_audit_implied <- function(fit, design_group) {
  unavailable <- function(reason) list(table = data.frame(), vcov = NULL,
    design = data.frame(), status = data.frame(status = "not_extracted", reason = reason))
  if (is.null(design_group)) return(unavailable("No covariance-design group requested."))
  if (!is.numeric(design_group) || length(design_group) != 1L || !is.finite(design_group) ||
      design_group < 1L || design_group != as.integer(design_group))
    stop("design_group must be a positive integer or NULL.")
  if (is.null(fit$prepared)) return(unavailable("Prepared design was not retained; fit with keep_data=TRUE."))
  if (design_group > length(fit$prepared$groups)) stop("design_group exceeds the number of retained covariance-design groups.")
  z <- fit$prepared$groups[[design_group]]; spec <- fit$settings
  roles <- c("mother", "father", rep("child", z$nchild))
  persons <- c("mother", "father", paste0("child", seq_len(z$nchild)))
  score_names <- unlist(lapply(1:2, function(k) paste0("PGS", k, "_", persons)))
  columns <- unlist(lapply(1:2, function(k) unname(fit$prepared$pgs_columns[paste0(roles, k)])))
  units <- c(fit$score_scaling$scale[match(columns, fit$score_scaling$column)],
    rep(spec$outcome_scale, length(z$child)))
  outcome_labels <- if (length(z$child)) paste0("Y", seq_along(z$child)) else character()
  labels <- c(score_names, outcome_labels)
  dimension <- length(labels)
  ij <- which(upper.tri(matrix(0, dimension, dimension), diag = TRUE), arr.ind = TRUE)
  fun <- function(par) {
    if (identical(spec$outcome, "binary_liability")) {
      # Extend the linear family-variance contribution across zero for the
      # derivative; do not differentiate log(variance) at a boundary.
      theta <- c(par, residual_sd = 0)
      current <- spec; current$family_random <- FALSE
      S <- .ravel_mixed_covariance(theta, z, current)
      if (spec$family_random && length(z$child)) {
        yi <- 2L * (z$nchild + 2L) + seq_along(z$child)
        S[yi, yi] <- S[yi, yi] + par[["family_variance"]]
      }
    } else S <- .ravel_mixed_covariance(par, z, spec)
    S <- S * outer(units, units)
    stats::setNames(S[ij], paste(labels[ij[, 1]], labels[ij[, 2]], sep = "|"))
  }
  map <- .ravel_mixed_delta(fun, fit$parameters)
  V <- matrix(NA_real_, length(map$value), length(map$value),
    dimnames = list(names(map$value), names(map$value)))
  if (!is.null(fit$vcov)) V <- map$jacobian %*%
    fit$vcov[names(fit$parameters), names(fit$parameters), drop = FALSE] %*% t(map$jacobian)
  tab <- data.frame(design_group = design_group, lhs = labels[ij[, 1]], rhs = labels[ij[, 2]],
    estimate = unname(map$value), se = sqrt(pmax(0, diag(V))),
    units = if (identical(spec$outcome, "binary_liability"))
      "input PGS units; Y is conditional-probit latent response" else "input PGS and outcome units",
    uncertainty_conditions = "FID sandwich; conditional on supplied pi", row.names = NULL)
  design <- data.frame(design_group = rep(design_group, length(z$child)),
    response = outcome_labels, child = if (length(z$child)) paste0("child", z$child) else character())
  for (j in seq_along(spec$interactions)) {
    nm <- spec$interactions[j]
    design[[paste0("moderator:", nm)]] <- z$moderators[, j + 1L] * spec$scales[[nm]] + spec$centers[[nm]]
  }
  if (!is.null(spec$random_slope)) design[[paste0("random_slope:", spec$random_slope)]] <-
    z$slope * spec$scales[[spec$random_slope]] + spec$centers[[spec$random_slope]]
  list(table = tab, vcov = V, design = design,
    status = data.frame(status = "ok", reason = paste("One of", length(fit$prepared$groups),
      "retained covariance-design groups; fixed-effect means do not enter this covariance.")))
}
