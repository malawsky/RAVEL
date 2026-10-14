# Genetic-liability x covariate interactions.
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
#
# Total-liability scale. If x is independent of every family member's visible
# and hidden components, the projection of Y = beta'T + gamma'(T x) + c x + e
# on (G, G x, x) gives b = A^{-1} B beta and g = A^{-1} B gamma for any
# centring, with A = Sigma_GG and B = Cov(G, T). The interaction vector is
# therefore corrected by the main-effect map, provided every role in the model
# has an interaction term (otherwise the omitted terms load on the included
# ones and, off-mean, on the main effects).

.ravel_interaction_role_codes <- c(child = "c", mother = "m", father = "f",
                                   mat_sib = "ms", pat_sib = "fs")
.ravel_interaction_role_labels <- c(c = "D", m = "M", f = "F", ms = "MS", fs = "FS")

#' Validate interaction settings
#' @noRd
.ravel_interaction_spec <- function(interactions, interaction_roles, interaction_center,
                                    dat, outcome, available_roles = c("child", "mother", "father")) {
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
      !all(interaction_roles %in% available_roles) || anyDuplicated(interaction_roles))
    stop("interaction_roles must contain distinct values from ",
         paste(sprintf("'%s'", available_roles), collapse = ", "), ".")
  roles <- available_roles[available_roles %in% interaction_roles]

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

#' lavaan label of an interaction coefficient (model scale)
#' @noRd
.ravel_interaction_label <- function(code, j) sprintf("g%s_%d", .ravel_interaction_role_labels[[code]], j)

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
        sprintf("Y ~ %s*%s", .ravel_interaction_label(code, j), lv))
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

#' Model-scale interaction coefficients (named by lavaan label)
#' @noRd
.ravel_interaction_theta <- function(fit, spec) {
  if (is.null(spec)) return(NULL)
  PE <- lavaan::parameterEstimates(fit)
  labels <- unlist(lapply(seq_along(spec$moderators), function(j)
    vapply(spec$codes, .ravel_interaction_label, character(1), j = j)), use.names = FALSE)
  setNames(vapply(labels, function(l)
    PE$est[!is.na(PE$label) & PE$label == l & PE$op == "~"][1], numeric(1)), labels)
}

#' Variance of the full genetic predictor in an interaction model
#'
#' The genetic predictor is b'G + sum_j g_j'(G x_j), x_j on the internal
#' centred/SD scale. Re-expressed at the moderator means m_j,
#'   b'G + sum_j g_j'GX_j = b_bar'G + sum_j g_j'(GX_j - m_j G),
#'   b_bar = b + sum_j m_j g_j,
#' so each component is invariant to the centre. Uses the model-implied latent
#' covariance matrix (no independence assumption) and model-scale coefficients.
#' @param b Model-scale main coefficients named by role code (c, m, f, ms, fs).
#' @noRd
.ravel_interaction_variance <- function(fit, spec, b) {
  J <- seq_along(spec$moderators)
  theta <- .ravel_interaction_theta(fit, spec)
  g <- lapply(J, function(j) {
    out <- setNames(numeric(length(b)), names(b))
    out[spec$codes] <- theta[vapply(spec$codes, .ravel_interaction_label, character(1), j = j)]
    out
  })
  m <- lavaan::lavInspect(fit, "mean.ov")[spec$moderators]
  shift <- Reduce(`+`, Map(`*`, g, m))
  base_lv <- unname(c(c = "Gc", m = "Gm", f = "Gf", ms = "Gms", fs = "Gfs")[names(b)])
  gx_lv <- unlist(lapply(J, function(j) sprintf("GX%s_%d", spec$codes, j)))
  lv <- c(base_lv, gx_lv)
  Psi <- lavaan::lavInspect(fit, "cov.lv")[lv, lv]
  main <- c(b + shift, numeric(length(gx_lv)))
  inter <- c(-shift, unlist(lapply(g, function(gj) gj[spec$codes])))
  q <- function(u, v) drop(crossprod(u, Psi %*% v))
  list(b_bar = b + shift, g = g, moderator_means = m,
       moderator_cov = lavaan::lavInspect(fit, "cov.ov")[spec$moderators, spec$moderators, drop = FALSE],
       components_vis = c(interaction = q(inter, inter), cross_main_interaction = 2 * q(main, inter)))
}

#' Total-scale interaction variance under moderator-genetic independence:
#' sum_jk Cov(x_j, x_k) gamma_j' Sigma_TT gamma_k, gamma_j = map(g_j)
#' @noRd
.ravel_interaction_total_variance <- function(iv, map, Sigma_TT) {
  gam <- lapply(iv$g, map)
  C <- iv$moderator_cov
  out <- 0
  for (j in seq_along(gam)) for (k in seq_along(gam))
    out <- out + C[j, k] * drop(crossprod(gam[[j]], Sigma_TT %*% gam[[k]]))
  out
}

#' Rows for one moderator: per-role coefficients plus parent-specific contrasts
#' @noRd
.ravel_interaction_rows <- function(spec, j) {
  codes <- spec$codes
  x <- spec$moderators[j]
  rows <- data.frame(parameter = sprintf("beta_%s:%s", .ravel_interaction_role_labels[codes], x),
    label = vapply(codes, .ravel_interaction_label, character(1), j = j),
    contrast_of = NA_character_, stringsAsFactors = FALSE)
  for (pair in list(c(parent = "m", sib = "ms", name = "M"), c(parent = "f", sib = "fs", name = "F")))
    if (all(pair[c("parent", "sib")] %in% codes))
      rows <- rbind(rows, data.frame(parameter = sprintf("delta_%s:%s", pair[["name"]], x),
        label = NA_character_, contrast_of = paste(.ravel_interaction_label(pair[["parent"]], j),
          .ravel_interaction_label(pair[["sib"]], j), sep = "-"), stringsAsFactors = FALSE))
  rows
}

#' Visible-scale interaction coefficients with lavaan or bootstrap inference
#'
#' With \code{draws} (bootstrap replicates of the model-scale interaction
#' coefficients, one column per lavaan label), SEs are bootstrap SDs and
#' intervals are percentile intervals of the same linear combinations.
#' @noRd
.ravel_interaction_table <- function(fit, spec, ci_level, draws = NULL) {
  if (is.null(spec)) return(NULL)
  theta <- .ravel_interaction_theta(fit, spec)
  V <- tryCatch(lavaan::vcov(fit), error = function(e) NULL)
  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  alpha <- 1 - ci_level
  if (!is.null(draws)) {
    draws <- draws[, names(theta), drop = FALSE]
    draws <- draws[apply(is.finite(draws), 1L, all), , drop = FALSE]
  }
  out <- do.call(rbind, lapply(seq_along(spec$moderators), function(j) {
    rows <- .ravel_interaction_rows(spec, j)
    s <- spec$scales[[spec$moderators[j]]]  # model uses x / SD(x); report per unit of x
    weights <- lapply(seq_len(nrow(rows)), function(i) {
      w <- setNames(numeric(length(theta)), names(theta))
      if (is.na(rows$contrast_of[i])) w[rows$label[i]] <- 1 else {
        parts <- strsplit(rows$contrast_of[i], "-", fixed = TRUE)[[1]]
        w[parts[1]] <- 1; w[parts[2]] <- -1
      }
      w / s
    })
    est <- vapply(weights, function(w) sum(w * theta), numeric(1))
    if (!is.null(draws)) {
      d <- draws %*% do.call(cbind, weights)
      se <- apply(d, 2L, stats::sd)
      lo <- apply(d, 2L, stats::quantile, probs = alpha / 2, names = FALSE)
      hi <- apply(d, 2L, stats::quantile, probs = 1 - alpha / 2, names = FALSE)
    } else {
      se <- vapply(weights, function(w) {
        if (is.null(V) || !all(names(theta) %in% rownames(V))) return(NA_real_)
        sqrt(max(0, drop(t(w) %*% V[names(theta), names(theta)] %*% w)))
      }, numeric(1))
      lo <- est - z * se
      hi <- est + z * se
    }
    data.frame(parameter = rows$parameter, role = c(spec$roles, rep(NA_character_, nrow(rows) - length(spec$roles))),
      moderator = spec$moderators[j], center = unname(spec$centers[[spec$moderators[j]]]),
      estimate = est, se = unname(se), ci_low = unname(lo), ci_high = unname(hi), stringsAsFactors = FALSE)
  }))
  out$p_value <- .uncertainty_p(out$estimate, out$se)
  out$uncertainty_source <- if (!is.null(draws)) sprintf("bootstrap (%d reps)", nrow(draws)) else
    if (identical(lavaan::lavInspect(fit, "options")$se, "robust.cluster"))
    "lavaan cluster-robust" else "lavaan/delta-method"
  rownames(out) <- NULL
  out
}

#' Fitted correlation of each visible genetic factor with each moderator
#'
#' The total-liability mapping of interactions assumes the moderator is
#' independent of the visible and hidden genetic components. Independence
#' implies these are zero; zero values do not imply independence.
#' @noRd
.ravel_interaction_independence <- function(fit, spec, latents) {
  if (is.null(spec)) return(NULL)
  PE <- lavaan::parameterEstimates(fit)
  rows <- list()
  for (x in spec$moderators) {
    vx <- PE$est[PE$lhs == x & PE$op == "~~" & PE$rhs == x][1]
    for (lv in latents) {
      i <- which(PE$op == "~~" & ((PE$lhs == lv & PE$rhs == x) | (PE$lhs == x & PE$rhs == lv)))[1]
      if (is.na(i)) next
      rows[[length(rows) + 1L]] <- data.frame(latent = lv, moderator = x,
        correlation = PE$est[i] / sqrt(vx), z = PE$est[i] / PE$se[i],
        p_value = 2 * stats::pnorm(-abs(PE$est[i] / PE$se[i])), stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Total-liability correction of the interaction coefficients
#'
#' Applies the main-effect visible-to-total map to the interaction vector of
#' each moderator, with delta-method uncertainty over the joint SEM covariance
#' of the core parameters and the interaction coefficients, optionally adding
#' external pi uncertainty (and, for pi-dependent fits, the refit sensitivity).
#' External pi must be independent of the whole SEM estimator: core_pi_cov
#' covers only the core estimates, and the interaction/pi covariance is not an
#' input, so correlated-pi interaction inference is reported as unavailable.
#' \code{fixed_labels} names core parameters the model fixes (bm and bf with
#' \code{fix_parent_paths}); they enter as constants with zero variance.
#' @noRd
.ravel_interaction_corrected <- function(fit, estimates, spec, h, am_correction, am_n_generations,
                                         sib_roles, ci_level, pi_se = NULL, core_pi_cov = NULL,
                                         sensitivity = NULL, draws = NULL,
                                         fixed_labels = character(0)) {
  if (is.null(spec)) return(NULL)
  core_labels <- c(.report_core_labels, .sib_core_labels(sib_roles))
  core <- c(bD = unname(estimates$beta_vis["bD"]), bm = unname(estimates$beta_vis["bm"]),
            bf = unname(estimates$beta_vis["bf"]), rhoMF = unname(estimates$latent_cov_vis["rhoMF"]),
            eta_m = unname(estimates$latent_cov_vis["eta_m"]), eta_f = unname(estimates$latent_cov_vis["eta_f"]),
            .sib_theta(estimates$parental_sibs, sib_roles))
  gtheta <- .ravel_interaction_theta(fit, spec)
  V_full <- tryCatch(lavaan::vcov(fit), error = function(e) NULL)
  with_pi <- !is.null(pi_se) && pi_se > 0
  # Zero-filling Cov(interaction, pi) beside a nonzero core/pi block is an
  # unstated assumption and need not give a valid joint covariance.
  correlated_pi <- with_pi && !is.null(core_pi_cov) && any(core_pi_cov[core_labels] != 0)
  if (correlated_pi && is.null(draws))
    warning("Corrected interaction uncertainty is unavailable with a nonzero core_pi_cov: ",
            "the interaction/pi covariance is not supplied. Estimates are retained.", call. = FALSE)
  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  alpha <- 1 - ci_level
  n_draws <- NA_integer_
  sib_code <- c(ms = "ms", fs = "fs")[sib_roles]
  status <- "ok"

  out <- lapply(seq_along(spec$moderators), function(j) {
    g_labels <- vapply(spec$codes, .ravel_interaction_label, character(1), j = j)
    names(g_labels) <- spec$codes
    s <- spec$scales[[spec$moderators[j]]]
    x_labels <- c(core_labels, unname(g_labels))
    x <- c(core[core_labels], gtheta[g_labels], if (with_pi) c(pi = h))
    transform <- function(v) {
      v <- setNames(as.numeric(v), names(x))
      sib <- if (length(sib_roles)) {
        si <- .sib_input_from_theta(v, sib_roles)
        si$b <- setNames(v[g_labels[sib_code]], c(ms = "bms", fs = "bfs")[sib_roles])
        si
      } else NULL
      est <- .compute_scale_estimates(v[g_labels[c("c", "m", "f")]], v["eta_m"], v["eta_f"], v["rhoMF"],
        if (with_pi) v["pi"] else h, am_correction, am_n_generations, sib = sib)
      std <- c(est$beta_total_assumed_std, if (length(sib_roles)) est$sib_total$std)
      raw <- c(est$beta_total_assumed, if (length(sib_roles)) est$sib_total$raw)
      c(std, raw) / s
    }
    estimate <- transform(x)
    rows <- .ravel_interaction_rows(spec, j)
    keys <- c("bD", "bm", "bf", if ("ms" %in% sib_roles) c("bms", "delta_m"), if ("fs" %in% sib_roles) c("bfs", "delta_f"))
    names_std <- c(beta_D = "bD", beta_M = "bm", beta_F = "bf", beta_MS = "bms", delta_M = "delta_m",
                   beta_FS = "bfs", delta_F = "delta_f")
    want <- names_std[sub(":.*$", "", rows$parameter)]
    n_half <- length(estimate) / 2
    idx <- match(want, keys)
    se <- lo <- hi <- rep(NA_real_, length(estimate))
    # Joint bootstrap draws of (core, interaction) on the model scale.
    dj <- NULL
    if (!is.null(draws) && all(x_labels %in% colnames(draws))) {
      dj <- draws[, x_labels, drop = FALSE]
      dj <- dj[apply(is.finite(dj), 1L, all), , drop = FALSE]
    }
    # Never zero-fill an unavailable block: missing SEM covariance, labels or
    # draws, an unsupported cross-covariance, or an invalid joint matrix gives NA.
    Vsrc <- if (correlated_pi) NULL else if (!is.null(draws)) {
      if (!is.null(dj) && nrow(dj) > 1L) stats::cov(dj) else NULL
    } else if (!is.null(V_full) &&
               all(setdiff(x_labels, rownames(V_full)) %in% fixed_labels)) {
      # Core parameters fixed by the model (bm = bf = 0 with fix_parent_paths)
      # are constants with zero sampling variance; only they may be absent.
      free_x <- intersect(x_labels, rownames(V_full))
      Vs <- matrix(0, length(x_labels), length(x_labels), dimnames = list(x_labels, x_labels))
      Vs[free_x, free_x] <- V_full[free_x, free_x, drop = FALSE]
      Vs
    } else NULL
    if (!is.null(dj)) n_draws <<- nrow(dj)
    if (!is.null(Vsrc) && !with_pi && !is.null(dj)) {
      # Percentile intervals of the per-replicate corrected coefficients, as for
      # the main corrected effects.
      td <- t(apply(dj, 1L, function(r) tryCatch(transform(r),
        error = function(e) rep(NA_real_, length(estimate)))))
      td <- td[apply(is.finite(td), 1L, all), , drop = FALSE]
      n_draws <<- nrow(td)
      if (nrow(td) > 1L) {
        se <- apply(td, 2L, stats::sd)
        lo <- apply(td, 2L, stats::quantile, probs = alpha / 2, names = FALSE)
        hi <- apply(td, 2L, stats::quantile, probs = 1 - alpha / 2, names = FALSE)
      } else status <<- "invalid"
    } else if (!is.null(Vsrc)) {
      J <- tryCatch(.finite_difference_jacobian(transform, x), error = function(e) NULL)
      Vx <- matrix(0, length(x), length(x), dimnames = list(names(x), names(x)))
      Vx[x_labels, x_labels] <- Vsrc
      if (with_pi) Vx["pi", "pi"] <- pi_se^2
      Vx <- tryCatch(.ravel_pi_validate_covariance(Vx, names(x), "Joint core/interaction/pi covariance"),
                     error = function(e) NULL)
      if (!is.null(J) && all(is.finite(J)) && !is.null(Vx)) {
        colnames(J) <- names(x)
        if (with_pi && !is.null(sensitivity)) {
          s <- setNames(sensitivity[x_labels], x_labels)
          s[intersect(fixed_labels, x_labels)] <- 0  # constants do not respond to pi
          J[, "pi"] <- J[, "pi"] + drop(J[, x_labels, drop = FALSE] %*% s)
        }
        cov <- J %*% Vx %*% t(J)
        d <- diag((cov + t(cov)) / 2)
        if (all(is.finite(d)) && all(d >= -1e-12)) se <- sqrt(pmax(0, d)) else status <<- "invalid"
      } else status <<- "invalid"
      lo <- estimate - z * se
      hi <- estimate + z * se
    } else if (!correlated_pi) status <<- "invalid"
    do.call(rbind, lapply(c("corrected_std", "corrected_raw"), function(scale) {
      k <- idx + if (scale == "corrected_raw") n_half else 0
      data.frame(scale = scale, parameter = rows$parameter, moderator = spec$moderators[j],
        estimate = unname(estimate[k]), se = unname(se[k]),
        ci_low = unname(lo[k]), ci_high = unname(hi[k]),
        stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, out)
  out$p_value <- .uncertainty_p(out$estimate, out$se)
  cluster <- identical(lavaan::lavInspect(fit, "options")$se, "robust.cluster")
  out$uncertainty_source <- if (correlated_pi)
    "unavailable: correlated external pi (interaction/pi covariance not supplied)" else
    if (status == "invalid") "unavailable: SEM covariance missing or not positive semidefinite" else
    if (!is.null(draws)) {
      if (with_pi) "bootstrap covariance + pi (delta-Wald)" else sprintf("bootstrap (%d reps)", n_draws)
    } else paste0(if (cluster) "lavaan cluster-robust" else "lavaan/delta-method",
      if (with_pi) " + pi (delta-Wald)" else "")
  rownames(out) <- NULL
  out
}

#' Interaction rows for the requested scale and uncertainty source
#'
#' Mirrors \code{.parental_sib_table}: \code{unc = NULL} (uncertainty "none")
#' keeps the point estimates and drops all inference; a bootstrap selection uses
#' the stored bootstrap tables, or falls back to lavaan with an explicit label
#' when no interaction draws exist; corrected rows are unavailable when the
#' requested external-pi correction failed.
#' @noRd
.ravel_interaction_rows_for <- function(result, unc, scale = c("visible", "corrected_std", "corrected_raw")) {
  scale <- match.arg(scale)
  ix <- result$interactions
  if (is.null(ix)) return(NULL)
  boot <- identical(unc$source, "bootstrap")
  pick <- function(lav, bt) if (boot && !is.null(bt)) bt else lav
  tab <- if (scale == "visible") pick(ix$table, ix$bootstrap$table) else {
    ct <- pick(ix$corrected, ix$bootstrap$corrected)
    if (is.null(ct)) NULL else ct[ct$scale == scale, , drop = FALSE]
  }
  if (is.null(tab) || !nrow(tab)) return(NULL)
  rownames(tab) <- NULL
  blank <- function(source, method) {
    tab[c("se", "ci_low", "ci_high", "p_value")] <- NA_real_
    tab$uncertainty_source <- source
    tab$p_value_method <- method
    tab
  }
  if (is.null(unc)) return(blank("none", "none"))
  if (scale != "visible" && (!is.null(result$pi_correction_error) || !is.null(unc$pi_correction_error)))
    return(blank("unavailable (incompatible pi/assortment state)", "unavailable"))
  src <- tab$uncertainty_source[1]
  tab$p_value_method <- if (startsWith(src, "unavailable")) "unavailable" else
    if (grepl("pi (delta-Wald)", src, fixed = TRUE)) "normal_approx_from_joint_pi_delta_se" else
    if (startsWith(src, "bootstrap")) "normal_approx_from_bootstrap_se" else "normal_approx_from_lavaan_se"
  if (boot && !startsWith(src, "bootstrap") && !startsWith(src, "unavailable"))
    tab$uncertainty_source <- paste0(tab$uncertainty_source, " (no bootstrap draws for interactions)")
  tab
}

#' Print interaction coefficients
#' @noRd
.print_interactions <- function(result, unc, digits, scale = c("visible", "corrected_std")) {
  scale <- match.arg(scale)
  tab <- .ravel_interaction_rows_for(result, unc, scale)
  if (is.null(tab) || !nrow(tab)) return(invisible(NULL))
  src <- if (identical(tab$uncertainty_source[1], "none")) "" else paste0("; ", tab$uncertainty_source[1])
  if (scale == "visible") {
    cat(sprintf("\nGenetic x covariate interactions (visible scale; per SD of G per unit of covariate%s):\n", src))
  } else {
    cat(sprintf("  Genetic x covariate interactions (standardized total liability%s):\n", src))
  }
  indent <- if (scale == "visible") "  " else "    "
  for (i in seq_len(nrow(tab))) {
    ci <- if (is.na(tab$se[i])) "" else
      sprintf(" [%.*f, %.*f]", digits, tab$ci_low[i], digits, tab$ci_high[i])
    cat(sprintf("%s%-14s %.*f%s\n", indent, tab$parameter[i], digits, tab$estimate[i], ci))
  }
  if (scale == "visible") {
    centers <- result$interactions$centers
    cat("  Main effects are at", paste(sprintf("%s = %s", names(centers),
      format(signif(centers, 4))), collapse = ", "), "\n")
  }
  invisible(NULL)
}
