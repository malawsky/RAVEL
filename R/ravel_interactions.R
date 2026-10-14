# Genetic-liability x covariate interactions (visible scale).
#
# A moderator x is observed without error, so for each score k and role r
#   PGS_kr * x = lambda_kr (G_r * x) + epsilon_kr * x.
# The products of the centred scores with the centred moderator are therefore
# indicators of the latent interaction G_r * x with the SAME loadings as G_r.
# The latent interaction has a free variance and free covariances with the
# other latents and covariates. Product-indicator residuals have free
# within-score covariances with one another and with the base indicators
# (the latter equal theta * E[x], which is nonzero when x is not centred at
# its mean). Repeated measurements (long format) are handled by
# cluster-robust standard errors.

.ravel_interaction_role_codes <- c(child = "c", mother = "m", father = "f")
.ravel_interaction_role_labels <- c(c = "D", m = "M", f = "F")

#' Validate interaction settings
#' @noRd
.ravel_interaction_spec <- function(interactions, interaction_roles, interaction_center,
                                    dat, outcome) {
  if (is.null(interactions) || !length(interactions)) {
    if (!is.null(interaction_center))
      stop("interaction_center requires interactions.")
    return(NULL)
  }
  if (!is.character(interactions) || anyNA(interactions) || any(!nzchar(interactions)) ||
      anyDuplicated(interactions))
    stop("interactions must be a character vector of distinct covariate column names.")
  if (outcome != "continuous")
    stop("Genetic x covariate interactions currently require a continuous outcome.")
  missing <- setdiff(interactions, names(dat))
  if (length(missing))
    stop("Missing interaction covariate columns: ", paste(missing, collapse = ", "))
  bad_name <- interactions[make.names(interactions) != interactions]
  if (length(bad_name))
    stop("Interaction covariate names must be syntactic R names: ", paste(bad_name, collapse = ", "))
  for (x in interactions) {
    v <- dat[[x]]
    if (!is.numeric(v) || any(!is.na(v) & !is.finite(v)))
      stop("Interaction covariate '", x, "' must be numeric and finite (or NA).")
    if (sum(!is.na(v)) < 50L || stats::var(v, na.rm = TRUE) <= 0)
      stop("Interaction covariate '", x, "' needs at least 50 observed, non-constant values.")
  }
  if (!is.character(interaction_roles) || !length(interaction_roles) ||
      !all(interaction_roles %in% names(.ravel_interaction_role_codes)) ||
      anyDuplicated(interaction_roles))
    stop("interaction_roles must contain distinct values from 'child', 'mother', 'father'.")
  roles <- names(.ravel_interaction_role_codes)[names(.ravel_interaction_role_codes) %in% interaction_roles]

  centers <- vapply(interactions, function(x) mean(dat[[x]], na.rm = TRUE), numeric(1))
  if (!is.null(interaction_center)) {
    if (!is.numeric(interaction_center) || is.null(names(interaction_center)) ||
        any(!is.finite(interaction_center)) || anyDuplicated(names(interaction_center)) ||
        !all(names(interaction_center) %in% interactions))
      stop("interaction_center must be a finite numeric vector named by interaction covariates.")
    centers[names(interaction_center)] <- interaction_center
  }
  # Moderators are also divided by their SD inside the model, purely for
  # numerical conditioning (e.g. squared ages); coefficients are reported per
  # unit of the original covariate.
  scales <- vapply(interactions, function(x) stats::sd(dat[[x]], na.rm = TRUE), numeric(1))
  list(moderators = interactions, roles = roles,
       codes = unname(.ravel_interaction_role_codes[roles]),
       centers = centers, scales = scales)
}

#' Internal name of a product indicator
#' @noRd
.ravel_product_name <- function(k, role, j) sprintf("PGS%d_%s_x%d", k, role, j)

#' Centre the moderators and add product indicators
#'
#' Scores are centred at their sample means before multiplying, so that the
#' product does not carry an intercept times x term.
#' @noRd
.ravel_interaction_data <- function(dat, spec) {
  if (is.null(spec)) return(dat)
  for (j in seq_along(spec$moderators)) {
    x <- spec$moderators[j]
    dat[[x]] <- (dat[[x]] - spec$centers[[x]]) / spec$scales[[x]]
    for (role in spec$roles) for (k in 1:2) {
      pgs <- dat[[sprintf("PGS%d_%s", k, role)]]
      dat[[.ravel_product_name(k, role, j)]] <- (pgs - mean(pgs, na.rm = TRUE)) * dat[[x]]
    }
  }
  dat
}

#' Product-indicator syntax for the latent interactions
#' @noRd
.ravel_interaction_syntax <- function(spec, fixed_effects, base_roles) {
  if (is.null(spec)) return("")
  lines <- "# Genetic x covariate interactions (product indicators share G's loadings)"
  products <- list(`1` = character(0), `2` = character(0))
  for (j in seq_along(spec$moderators)) {
    for (i in seq_along(spec$roles)) {
      role <- spec$roles[i]
      code <- spec$codes[i]
      lv <- sprintf("GX%s_%d", code, j)
      p <- vapply(1:2, function(k) .ravel_product_name(k, role, j), character(1))
      lines <- c(lines,
        sprintf("%s =~ l%s1*%s + l%s2*%s", lv, code, p[1], code, p[2]),
        sprintf("%s ~~ NA*%s", lv, lv),
        sprintf("Y ~ g%s_%d*%s", .ravel_interaction_role_labels[[code]], j, lv))
      if (length(fixed_effects))
        lines <- c(lines, sprintf("%s ~~ %s", lv, paste(fixed_effects, collapse = " + ")))
      products[["1"]] <- c(products[["1"]], p[1])
      products[["2"]] <- c(products[["2"]], p[2])
    }
  }
  for (k in 1:2) {
    pk <- products[[as.character(k)]]
    if (length(pk) > 1L) {
      pairs <- utils::combn(pk, 2L)
      lines <- c(lines, sprintf("%s ~~ %s", pairs[1, ], pairs[2, ]))
    }
    base <- sprintf("PGS%d_%s", k, base_roles)
    lines <- c(lines, sprintf("%s ~~ %s", rep(base, each = length(pk)), rep(pk, times = length(base))))
  }
  paste0("    ", lines, collapse = "\n")
}

#' Visible-scale interaction coefficients with lavaan inference
#' @noRd
.ravel_interaction_table <- function(fit, spec, ci_level) {
  if (is.null(spec)) return(NULL)
  PE <- lavaan::parameterEstimates(fit, level = ci_level)
  rows <- list()
  for (j in seq_along(spec$moderators)) for (i in seq_along(spec$roles)) {
    lab_role <- .ravel_interaction_role_labels[[spec$codes[i]]]
    label <- sprintf("g%s_%d", lab_role, j)
    r <- which(!is.na(PE$label) & PE$label == label & PE$op == "~")[1]
    s <- spec$scales[[spec$moderators[j]]]  # model uses x / SD(x); report per unit of x
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = sprintf("beta_%s:%s", lab_role, spec$moderators[j]),
      label = label, role = spec$roles[i], moderator = spec$moderators[j],
      center = unname(spec$centers[[spec$moderators[j]]]),
      estimate = PE$est[r] / s, se = PE$se[r] / s, ci_low = PE$ci.lower[r] / s,
      ci_high = PE$ci.upper[r] / s, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$p_value <- .uncertainty_p(out$estimate, out$se)
  out$uncertainty_source <- if (identical(lavaan::lavInspect(fit, "options")$se, "robust.cluster"))
    "lavaan cluster-robust" else "lavaan/delta-method"
  rownames(out) <- NULL
  out
}

#' Print interaction coefficients
#' @noRd
.print_interactions <- function(result, digits) {
  tab <- result$interactions$table
  if (is.null(tab)) return(invisible(NULL))
  cat(sprintf("\nGenetic x covariate interactions (visible scale; per SD of G per unit of covariate; %s):\n",
              tab$uncertainty_source[1]))
  for (i in seq_len(nrow(tab))) {
    ci <- if (is.na(tab$se[i])) "" else
      sprintf(" [%.*f, %.*f]", digits, tab$ci_low[i], digits, tab$ci_high[i])
    cat(sprintf("  %-14s %.*f%s\n", tab$parameter[i], digits, tab$estimate[i], ci))
  }
  centers <- result$interactions$centers
  cat("  Main effects are at", paste(sprintf("%s = %s", names(centers),
    format(signif(centers, 4))), collapse = ", "), "\n")
  invisible(NULL)
}
