# Linear contrasts for ordinary PGS regressions, using the fitted joint covariance.
.ravel_pgs_contrast_table <- function(b, V, mother, father, mat_sib = NULL,
                                      pat_sib = NULL, ci_level = .95, df = Inf,
                                      parent_contrasts = TRUE) {
  weights <- list()
  if (isTRUE(parent_contrasts)) {
    weights <- list(beta_mid = stats::setNames(c(1, 1), c(mother, father)),
      delta_M_mid = stats::setNames(c(1, -1), c(mother, father)),
      delta_P_mid = stats::setNames(c(-1, 1), c(mother, father)))
  }
  if (length(mat_sib)) weights$delta_M <- stats::setNames(c(1, -1), c(mother, mat_sib))
  if (length(pat_sib)) weights$delta_F <- stats::setNames(c(1, -1), c(father, pat_sib))
  if (!length(weights)) return(NULL)
  critical <- stats::qt(1 - (1 - ci_level) / 2, df)
  rows <- lapply(names(weights), function(parameter) {
    w <- weights[[parameter]]
    needed <- names(w)
    estimate <- sum(w * b[needed])
    variance <- if (all(needed %in% rownames(V)))
      drop(t(w) %*% V[needed, needed, drop = FALSE] %*% w) else NA_real_
    se <- sqrt(pmax(0, variance))
    data.frame(parameter = parameter, estimate = unname(estimate), se = se,
      ci_low = unname(estimate - critical * se), ci_high = unname(estimate + critical * se),
      p_value = unname(2 * stats::pt(-abs(estimate / se), df)), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Also supports older fit objects that retain the underlying lm/glm fits.
.ravel_vanilla_parent_contrasts <- function(vanilla) {
  if (!is.null(vanilla$parent_contrasts)) return(vanilla$parent_contrasts)
  rows <- lapply(c("PGS1", "PGS2"), function(version) {
    object <- vanilla[[version]]
    fit <- object$fit
    if (is.null(fit)) return(NULL)
    tab <- .ravel_pgs_contrast_table(stats::coef(fit), stats::vcov(fit),
      ".pgs_mother", ".pgs_father", ci_level = vanilla$ci_level %||% .95,
      df = if (inherits(fit, "glm")) Inf else fit$df.residual)
    data.frame(regression = version, tab, family = object$table$family[1],
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.print_vanilla_contrasts <- function(result, parent_contrasts, digits) {
  vanilla <- result$vanilla_regression
  tab <- vanilla$table
  if (is.null(tab)) return(invisible(NULL))
  tab <- tab[tab$parameter %in% c("delta_M", "delta_F"), , drop = FALSE]
  if (isTRUE(parent_contrasts)) tab <- rbind(tab, .ravel_vanilla_parent_contrasts(vanilla))
  if (!nrow(tab)) return(invisible(NULL))
  cat("\nObserved-PGS contrasts (ordinary model SEs):\n")
  for (i in seq_len(nrow(tab))) {
    n <- vanilla[[tab$regression[i]]]$N
    cat(sprintf("  %s %-12s %.*f (SE %.*f) [%.*f, %.*f], p = %s; N = %s\n",
      tab$regression[i], tab$parameter[i], digits, tab$estimate[i], digits, tab$se[i],
      digits, tab$ci_low[i], digits, tab$ci_high[i], format.pval(tab$p_value[i], digits = 3), n))
  }
  invisible(NULL)
}
