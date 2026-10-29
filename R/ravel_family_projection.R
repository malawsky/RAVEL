# Reduced family regressions are projections of the fitted joint moments.
# Auxiliary family members remain in the measurement model; their coefficients
# are nuisance parameters, not zero restrictions on the requested regression.

.ravel_family_role_codes <- c(child = "c", mother = "m", father = "f",
                             mat_sib = "ms", pat_sib = "fs")

# Build a light evaluator from lavaan free-parameter IDs. This avoids modifying
# a fitted lavaan object or relying on lavaan's internal matrix setters.
.ravel_family_lavaan_design <- function(fit, roles, interaction_spec = NULL,
                                        fixed_effects = character(0)) {
  pt <- lavaan::parTable(fit)
  ordered_y <- "Y" %in% lavaan::lavNames(fit, "ov.ord")
  fixed_effects <- fixed_effects %||% character(0)
  lv <- paste0("G", unname(.ravel_family_role_codes[roles]))
  dn <- roles
  dr <- roles
  dm <- rep("", length(roles))
  units <- rep(1, length(roles))
  if (!is.null(interaction_spec)) {
    for (j in seq_along(interaction_spec$moderators)) {
      moderator <- interaction_spec$moderators[j]
      ir <- interaction_spec$roles
      lv <- c(lv, sprintf("GX%s_%d", .ravel_family_role_codes[ir], j))
      dn <- c(dn, paste(ir, moderator, sep = ":"))
      dr <- c(dr, ir)
      dm <- c(dm, rep(moderator, length(ir)))
      units <- c(units, rep(interaction_spec$scales[[moderator]], length(ir)))
    }
  }
  value <- function(theta, row, default = 0) {
    if (!length(row)) return(default)
    row <- row[1L]
    if (pt$free[row] > 0L) theta[pt$free[row]] else pt$est[row]
  }
  row <- function(lhs, op, rhs) which(pt$lhs == lhs & pt$op == op & pt$rhs == rhs)
  covariance <- function(theta, a, b) {
    r <- row(a, "~~", b)
    if (!length(r)) r <- row(b, "~~", a)
    value(theta, r)
  }
  matrix_cov <- function(theta, a, b) {
    out <- matrix(0, length(a), length(b), dimnames = list(a, b))
    for (i in seq_along(a)) for (j in seq_along(b))
      out[i, j] <- covariance(theta, a[i], b[j])
    out
  }
  function(theta) {
    V <- matrix_cov(theta, lv, lv)
    dimnames(V) <- list(dn, dn)
    B <- setNames(vapply(lv, function(x) value(theta, row("Y", "~", x)), numeric(1)), dn)
    CW <- matrix_cov(theta, lv, fixed_effects)
    rownames(CW) <- dn
    WW <- matrix_cov(theta, fixed_effects, fixed_effects)
    bw <- setNames(vapply(fixed_effects, function(x)
      value(theta, row("Y", "~", x)), numeric(1)), fixed_effects)
    residual <- value(theta, row("Y", "~~", "Y"), NA_real_)
    all_cov <- rbind(cbind(V, CW), cbind(t(CW), WW))
    all_b <- c(B, bw)
    vy <- if (ordered_y) 1 else unname(drop(crossprod(all_b, all_cov %*% all_b)) + residual)
    yw <- if (length(fixed_effects)) drop(t(CW) %*% B + WW %*% bw) else numeric(0)
    names(yw) <- fixed_effects
    mods <- interaction_spec$moderators %||% character(0)
    mm <- setNames(vapply(mods, function(x) value(theta, row(x, "~1", "")), numeric(1)), mods)
    mc <- matrix_cov(theta, mods, mods)
    list(beta_visible = B[roles], Sigma_GG = V[roles, roles, drop = FALSE],
      Sigma_GW = CW[roles, , drop = FALSE], Sigma_WW = WW,
      beta_covariate = bw, outcome_variance = vy, outcome_covariate_cov = yw,
      design_cov_visible = V, design_beta_visible = B,
      design_covariate_cov_visible = CW, design_names = dn,
      design_roles = dr, design_moderators = dm, design_units = units,
      moderator_means = mm, moderator_covariance = mc)
  }
}

.ravel_family_legacy_descriptor <- function(object) {
  if (!inherits(object, "ravel_results") || !isTRUE(object$converged) || is.null(object$fit))
    stop("object must be a converged RAVEL fit or a family fit with a family_descriptor.")
  fit <- object$fit
  if (!suppressWarnings(isTRUE(lavaan::lavInspect(fit, "post.check"))))
    stop("The source SEM is inadmissible.")
  if (!.ravel_family_identified(fit))
    stop("The source SEM is not locally identified with regular, finite uncertainty.")
  sib <- object$parental_sibs$roles %||% character(0)
  roles <- c("child", "mother", "father", c(ms = "mat_sib", fs = "pat_sib")[sib])
  roles <- unname(roles)
  theta <- lavaan::coef(fit)
  names(theta) <- make.unique(names(theta))
  V <- lavaan::vcov(fit)
  dimnames(V) <- list(names(theta), names(theta))
  pt <- lavaan::parTable(fit)
  source_row <- function(label) {
    r <- which(pt$label == label & pt$op %in% c("~", "~~"))[1L]
    # Legacy fix_parent_paths writes Y ~ 0*Gm + 0*Gf without the bm/bf
    # labels. Recover the structural row, retaining its fixed value and free
    # ID zero; an absent label is not an absent family member.
    if (is.na(r)) {
      outcome_lv <- c(bD = "Gc", bm = "Gm", bf = "Gf", bms = "Gms", bfs = "Gfs")
      if (label %in% names(outcome_lv))
        r <- which(pt$lhs == "Y" & pt$op == "~" & pt$rhs == outcome_lv[[label]])[1L]
    }
    r
  }
  grab <- function(theta, label) {
    r <- source_row(label)
    if (is.na(r)) stop("Missing source parameter: ", label)
    if (pt$free[r] > 0L) unname(theta[pt$free[r]]) else pt$est[r]
  }
  design <- .ravel_family_lavaan_design(fit, roles, object$interactions, object$fixed_effects)
  moments <- function(theta, pi) {
    out <- design(theta)
    core_names <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f", .sib_core_labels(sib))
    core <- setNames(vapply(core_names, function(x) grab(theta, x), numeric(1)), core_names)
    sc <- suppressWarnings(.compute_scale_estimates(core[c("bD", "bm", "bf")],
      core["eta_m"], core["eta_f"], core["rhoMF"], pi,
      object$am_correction, object$am_n_generations,
      sib = .sib_input_from_theta(core, sib)))
    if (length(sib)) {
      TT <- sc$sib_total$Sigma_TT
      M <- sc$sib_total$M
      br <- c(sc$beta_total_assumed,
        sc$sib_total$raw[c(ms = "bms", fs = "bfs")[sib]])
      reason <- sc$sib_total$invalid_reason
    } else {
      hd <- sc$hidden_correction_details
      a <- sqrt(1 - pi) * sqrt(c(hd$v_Gc, hd$v_G, hd$v_G))
      b <- sqrt(pi) * sqrt(c(hd$v_Hc, hd$v_H, hd$v_H))
      TT <- if (all(is.finite(c(a, b))) && !is.null(hd$Sigma_hh))
        diag(a) %*% hd$Sigma_xx %*% diag(a) +
          diag(a) %*% hd$Sigma_xh %*% diag(b) +
          diag(b) %*% t(hd$Sigma_xh) %*% diag(a) +
          diag(b) %*% hd$Sigma_hh %*% diag(b) else NULL
      M <- hd$M
      br <- sc$beta_total_assumed
      reason <- NULL
    }
    if (!is.null(TT)) dimnames(TT) <- list(roles, roles)
    if (!is.null(M)) dimnames(M) <- list(roles, roles)
    names(br) <- roles
    c(out, list(beta_total_raw = br, Sigma_TT = TT, M = M,
      state = sc$hidden_correction_details, invalid_reason = reason))
  }
  draws <- object$bootstrap$parameter_draws
  if (!is.null(draws) && ncol(draws) == length(theta)) colnames(draws) <- names(theta)
  sensitivity <- NULL
  if (isTRUE(object$constraint_spec$depends_on_pi)) {
    sensitivity <- function(pi) {
      args <- object$constraint_refit
      if (is.null(args)) stop("The source fit does not retain data for pi-sensitivity refitting.")
      evaluate <- function(p) {
        args$h <- p
        z <- suppressWarnings(do.call(ravel_fit, args))
        if (!isTRUE(z$converged)) stop("A pi-sensitivity refit did not converge.")
        b <- lavaan::coef(z$fit)
        names(b) <- make.unique(names(b))
        if (!setequal(names(b), names(theta))) stop("Pi refitting changed the parameter structure.")
        b[names(theta)]
      }
      .ravel_family_pi_sensitivity(evaluate, pi, theta)
    }
  }
  sensitivity <- .ravel_family_cached_sensitivity(sensitivity, names(theta))
  list(theta = theta, vcov = V, roles = roles, pi = object$h,
    pi_se = object$pi_se %||% 0, moments = moments, bootstrap_draws = draws,
    pi_sensitivity = sensitivity, pi_dependent = isTRUE(object$constraint_spec$depends_on_pi),
    has_covariates = length(object$fixed_effects) > 0L,
    has_interactions = !is.null(object$interactions),
    model_source = "legacy_joint_family_SEM", N = object$N,
    outcome_variable = object$outcome_variable, outcome = object$outcome,
    has_legacy_bootstrap = !is.null(object$bootstrap),
    correlated_external_pi = identical(attr(object$pi_correction, "pi_uncertainty")$independent_of_core, FALSE),
    covariance_assumptions = "Existing RAVEL full visible-to-total projection model",
    core_parameter_indices = setNames(vapply(c(.report_core_labels, .sib_core_labels(sib)), function(l) {
      r <- source_row(l)
      if (is.na(r)) 0L else pt$free[r]
    }, integer(1)), c(.report_core_labels, .sib_core_labels(sib))))
}

# Each descriptor owns its cache. Regressions on different subsets share the
# same validated derivative; failed or malformed refits are never cached.
.ravel_family_cached_sensitivity <- function(callback, parameter_names) {
  if (is.null(callback)) return(NULL)
  if (!is.function(callback) || !is.character(parameter_names) ||
      !length(parameter_names) || anyNA(parameter_names) || anyDuplicated(parameter_names))
    stop("Sensitivity caching requires a function and distinct parameter names.")
  cache <- new.env(parent = emptyenv())
  force(callback)
  force(parameter_names)
  function(pi) {
    if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi) || pi < 0 || pi >= 1)
      stop("Sensitivity pi must be a finite scalar in [0, 1).")
    key <- sprintf("%.17g", pi)
    if (exists(key, envir = cache, inherits = FALSE))
      return(get(key, envir = cache, inherits = FALSE))
    value <- callback(pi)
    if (!is.numeric(value) || is.complex(value) || is.null(names(value)) ||
        anyDuplicated(names(value)) || !setequal(names(value), parameter_names) ||
        any(!is.finite(value)))
      stop("The full-parameter pi sensitivity must be finite and uniquely named for every source parameter.")
    value <- value[parameter_names]
    assign(key, value, envir = cache)
    value
  }
}

.ravel_family_pi_sensitivity <- function(evaluate, pi, center) {
  cache <- new.env(parent = emptyenv())
  get <- function(p) {
    key <- format(p, digits = 17)
    if (!exists(key, cache, inherits = FALSE)) assign(key, evaluate(p), cache)
    base::get(key, cache, inherits = FALSE)
  }
  steps <- unique(pmin(c(.001, .002, .004), (1 - pi) / 4))
  for (step in steps) {
    if (step < 100 * .Machine$double.eps) next
    if (pi >= step) {
      coarse <- (get(pi + step) - get(pi - step)) / (2 * step)
      fine <- (get(pi + step / 2) - get(pi - step / 2)) / step
    } else {
      coarse <- (get(pi + step) - center) / step
      fine <- (get(pi + step / 2) - center) / (step / 2)
    }
    if (all(is.finite(fine)) && all(abs(coarse - fine) <= .001 + .02 * abs(fine)))
      return(fine)
  }
  stop("Full-parameter pi-sensitivity refits are numerically unstable.")
}

.ravel_family_conditional_covariance <- function(S, CW, WW) {
  if (!ncol(CW)) return(S)
  .ravel_family_positive_definite(WW, "Covariate covariance")
  S - CW %*% solve(WW, t(CW))
}

.ravel_family_positive_definite <- function(S, what) {
  if (!is.matrix(S) || !nrow(S) || nrow(S) != ncol(S) || any(!is.finite(S)))
    stop(what, " is not a finite square covariance matrix.")
  E <- eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(E) <= 1e-10 * max(1, max(abs(E))))
    stop(what, " is singular or not positive definite.")
  invisible(TRUE)
}

# Evaluate a chosen regression from the full moments. Unselected liabilities
# are marginalized, rather than assigned zero coefficients.
.ravel_family_projection_values <- function(m, roles, selected, pi,
                                             covariate_mapping_confirmed) {
  if (is.null(m$design_cov_visible)) {
    m$design_cov_visible <- m$Sigma_GG
    m$design_beta_visible <- m$beta_visible
    m$design_covariate_cov_visible <- m$Sigma_GW %||% matrix(0, length(roles), 0)
    m$design_names <- m$design_roles <- roles
    m$design_moderators <- rep("", length(roles))
    m$design_units <- rep(1, length(roles))
  }
  dn <- m$design_names
  choose <- which(m$design_roles %in% selected)
  # Return requested roles first in the user's order, with moderator blocks in
  # their existing order. This changes labels, never matrix indexing semantics.
  blocks <- unique(m$design_moderators)
  choose <- unlist(lapply(blocks, function(z) {
    k <- which(m$design_moderators == z & m$design_roles %in% selected)
    k[order(match(m$design_roles[k], selected))]
  }), use.names = FALSE)
  cn <- dn[choose]
  WW <- m$Sigma_WW %||% matrix(0, 0, 0)
  CW <- m$design_covariate_cov_visible
  if (is.null(CW)) CW <- matrix(0, length(dn), 0)
  SV <- m$design_cov_visible
  bv <- m$design_beta_visible
  .ravel_family_positive_definite(SV, "Visible predictor covariance")
  CV <- .ravel_family_conditional_covariance(SV, CW, WW)
  .ravel_family_positive_definite(CV, "Conditional visible predictor covariance")
  visible <- if (length(choose) == length(dn)) bv[choose] else
    solve(CV[choose, choose, drop = FALSE], CV[choose, , drop = FALSE] %*% bv)
  visible <- setNames(as.numeric(visible), cn)
  vy <- unname(m$outcome_variance)
  yw <- m$outcome_covariate_cov %||% numeric(0)
  vy_cond <- vy - if (length(yw)) drop(crossprod(yw, solve(WW, yw))) else 0
  if (!all(is.finite(c(vy, vy_cond))) || min(vy, vy_cond) <= 0)
    stop("The outcome variance, including its covariate-adjusted variance, must be positive.")
  out <- list(visible = visible / m$design_units[choose])
  contrast <- function(beta) {
    ans <- numeric(0)
    for (mod in blocks) for (pr in list(c("mother", "mat_sib", "delta_m"),
                                       c("father", "pat_sib", "delta_f"))) {
      key <- if (!nzchar(mod)) pr[1:2] else paste(pr[1:2], mod, sep = ":")
      if (all(key %in% names(beta))) ans[if (!nzchar(mod)) pr[3] else paste(pr[3], mod, sep = ":")] <-
        beta[key[1]] - beta[key[2]]
    }
    ans
  }
  metrics <- function(beta, S, C) {
    explained <- drop(crossprod(beta, S[choose, choose, drop = FALSE] %*% beta))
    conditional <- drop(crossprod(beta, C[choose, choose, drop = FALSE] %*% beta))
    c(genetic_variance = explained, genetic_R2 = explained / vy,
      conditional_genetic_variance = conditional,
      partial_R2 = conditional / vy_cond)
  }
  result <- c(setNames(c(out$visible, contrast(out$visible)),
                      paste0("visible:", names(c(out$visible, contrast(out$visible))))),
    setNames(metrics(visible, SV, CV), paste0("visible:", names(metrics(visible, SV, CV)))))
  corrected_names <- c(cn, names(contrast(visible)), names(metrics(visible, SV, CV)))
  corrected_empty <- setNames(rep(NA_real_, 2 * length(corrected_names)),
    c(paste0("total_raw:", corrected_names), paste0("total_std:", corrected_names)))
  reason <- m$invalid_reason
  valid <- !is.null(m$Sigma_TT) && !is.null(m$M) &&
    all(is.finite(m$Sigma_TT)) && all(is.finite(m$M)) &&
    all(is.finite(m$beta_total_raw))
  if (!valid) return(list(value = c(result, corrected_empty),
    invalid_reason = reason %||% "incompatible_pi_or_assortment_state", selected_design = cn))
  if (pi > 0 && ncol(CW) && !isTRUE(covariate_mapping_confirmed))
    return(list(value = c(result, corrected_empty),
      invalid_reason = "Corrected covariate-adjusted reduced regressions require covariate_mapping_confirmed = TRUE: Cov(T - M'G, W) = 0 (also for included interaction predictors).",
      selected_design = cn))
  TT <- m$Sigma_TT[roles, roles, drop = FALSE]
  M <- m$M[roles, roles, drop = FALSE]
  fullB <- m$beta_total_raw[roles]
  scale <- sqrt(diag(TT))[match(m$design_roles, roles)]
  if (length(blocks) == 1L) {
    ST <- TT
    CTW <- t(M) %*% CW
    bt <- fullB
  } else if (pi == 0) {
    # At pi=0 the total transformation is diagonal; retain the fitted product
    # covariance, including departures from moderator-genetic independence.
    a <- diag(M)[match(m$design_roles, roles)]
    ST <- SV * outer(a, a)
    CTW <- CW * a
    bt <- bv / a
  } else {
    if (any(vapply(blocks[-1], function(mod)
      !setequal(m$design_roles[m$design_moderators == mod], roles), logical(1))))
      return(list(value = c(result, corrected_empty),
        invalid_reason = "Corrected family interactions require an interaction for every genetic-model role.",
        selected_design = cn))
    mods <- blocks[-1]
    means <- m$moderator_means[mods]
    Cmod <- m$moderator_covariance[mods, mods, drop = FALSE]
    zz <- rbind(c(1, means), cbind(means, Cmod + tcrossprod(means)))
    ST <- kronecker(zz, TT)
    CTW <- matrix(0, length(dn), ncol(CW))
    bt <- numeric(length(dn))
    for (j in seq_along(blocks)) {
      idx <- which(m$design_moderators == blocks[j])
      order_roles <- match(roles, m$design_roles[idx])
      idx <- idx[order_roles]
      CTW[idx, ] <- t(M) %*% CW[idx, , drop = FALSE]
      bt[idx] <- solve(M, bv[idx])
    }
    # kronecker uses canonical role order; design may contain a different
    # interaction-role order, so explicitly align every block.
    kn <- unlist(lapply(blocks, function(z)
      if (!nzchar(z)) roles else paste(roles, z, sep = ":")), use.names = FALSE)
    dimnames(ST) <- list(kn, kn)
    ST <- ST[dn, dn, drop = FALSE]
  }
  dimnames(ST) <- list(dn, dn)
  checked <- tryCatch({
    .ravel_family_positive_definite(ST, "Total-liability predictor covariance")
    CT <- .ravel_family_conditional_covariance(ST, CTW, WW)
    .ravel_family_positive_definite(CT, "Conditional total-liability predictor covariance")
    CT
  }, error = function(e) e)
  if (inherits(checked, "error")) return(list(value = c(result, corrected_empty),
    invalid_reason = conditionMessage(checked), selected_design = cn))
  CT <- checked
  full_explained <- drop(crossprod(bt, CT %*% bt))
  if (!is.finite(full_explained) || full_explained > vy_cond + 1e-8 * max(1, vy_cond))
    return(list(value = c(result, corrected_empty),
      invalid_reason = "The corrected joint model implies a negative outcome residual variance.",
      selected_design = cn))
  total <- setNames(as.numeric(if (length(choose) == length(dn)) bt[choose] else
    solve(CT[choose, choose, drop = FALSE], CT[choose, , drop = FALSE] %*% bt)), cn)
  raw <- total / m$design_units[choose]
  std <- raw * scale[choose]
  met <- metrics(total, ST, CT)
  raw_out <- c(raw, contrast(raw), met)
  std_out <- c(std, contrast(std), met)
  list(value = c(result, setNames(raw_out, paste0("total_raw:", names(raw_out))),
    setNames(std_out, paste0("total_std:", names(std_out)))),
    invalid_reason = NULL, selected_design = cn)
}

#' Regress an outcome on a chosen subset of family liabilities
#'
#' The joint source SEM retains every genetic-model member. This function
#' reports reduced population linear projections by marginalizing unselected
#' liabilities, rather than setting their full-model coefficients to zero.
#' The total-liability result inherits the source fit's hidden-component model.
#'
#' @details
#' If A contains all genetic-model members and R the requested members, the
#' raw total-liability target is
#' \eqn{\Sigma_{TT,RR}^{-1}\Sigma_{TT,RA}\beta_A}.
#' It is not the R subvector of the full-model coefficient vector. The source
#' correction assumes that the full projection residual is uncorrelated with
#' the visible and hidden family genetic components. Thus omitted coefficients
#' need not be zero, but a reduced regression is not a new causal estimand.
#'
#' With covariates, projection uses conditional (Schur-complement) covariance
#' matrices. The visible adjustment uses fitted genetic-covariate moments.
#' Corrected adjustment additionally requires
#' \eqn{\operatorname{Cov}(T-M^\top G,W)=0}, giving
#' \eqn{\Sigma_{TW}=M^\top\Sigma_{GW}}. Set
#' \code{covariate_mapping_confirmed=TRUE} only when this assumption is intended.
#' Without it the visible results remain available. At pi zero this additional
#' assumption is unnecessary. Corrected interactions inherit the source model's
#' independence assumption for moderators and every family genetic component;
#' an interaction must be present for every genetic-model member.
#'
#' Coefficients retain the phenotype units of the source SEM. Total raw units
#' are founder-scaled liability units; total standardized coefficients use the
#' current-generation marginal liability SD. Interaction coefficients are per
#' unit of the original moderator and main effects retain its chosen center.
#' Genetic variance is the variance of the selected genetic linear predictor.
#' Genetic R-squared divides it by marginal outcome variance. Partial R-squared
#' divides covariate-adjusted genetic variance by covariate-adjusted outcome
#' variance. These are joint regression summaries, not per-coefficient causal
#' variance partitions. Binary outcomes retain the source liability scale.
#'
#' Delta uncertainty propagates the complete source parameter covariance.
#' External pi uncertainty includes refitting sensitivity when fitting depends
#' on pi. Visible intervals remain conditional on the supplied pi. Bootstrap
#' inference requires full stored parameter draws; old core-only draws omit
#' uncertainty needed for covariance projections and variance explained.
#'
#'
#' @param object A live RAVEL fit or family fit.
#' @param outcome_roles Distinct role names to retain in the outcome regression.
#' @param pi Hidden founder variance fraction; defaults to the source fit.
#' @param pi_se External standard error of pi; defaults to the source fit or zero.
#' @param ci_level Confidence level.
#' @param uncertainty One of auto, lavaan, bootstrap or none.
#' @param covariate_mapping_confirmed For corrected reduced regressions with
#'   covariates, assert Cov(T - M'G, W)=0, also for included product predictors.
#' @param parameter_pi_cov Optional covariance of every source free SEM
#'   parameter with external pi, named exactly like the source parameter vector.
#' @return A \code{ravel_family_regression} list with a tidy \code{table}, named
#'   \code{coefficients}, \code{contrasts}, \code{variance_explained}, their joint
#'   \code{vcov}, role metadata, correction and inference status, and the
#'   unchanged \code{source_fit}. Table scales are visible, total_raw and
#'   total_std. Parameters include role names, role:moderator interactions,
#'   delta_m/delta_f contrasts where available, and variance summaries.
#' @seealso \code{ravel_family_fit}, \code{ravel_fit}
#' @export
ravel_family_regression <- function(object, outcome_roles, pi = NULL, pi_se = NULL,
                                     ci_level = .95,
                                     uncertainty = c("auto", "lavaan", "bootstrap", "none"),
                                     covariate_mapping_confirmed = FALSE,
                                     parameter_pi_cov = NULL) {
  uncertainty <- match.arg(uncertainty)
  if (inherits(object, "ravel_family_regression")) object <- object$source_fit
  d <- object$family_descriptor
  if (is.null(d)) d <- .ravel_family_legacy_descriptor(object)
  roles <- unname(d$roles)
  if (!is.character(outcome_roles) || !length(outcome_roles) || anyNA(outcome_roles) ||
      anyDuplicated(outcome_roles) || !all(outcome_roles %in% roles))
    stop("outcome_roles must be distinct genetic-model roles: ", paste(roles, collapse = ", "))
  if (!is.numeric(ci_level) || length(ci_level) != 1L || !is.finite(ci_level) ||
      ci_level <= 0 || ci_level >= 1) stop("ci_level must be strictly between zero and one.")
  pi <- pi %||% d$pi %||% 0
  pi_se <- pi_se %||% d$pi_se %||% 0
  .ravel_pi_se_variance(pi, pi_se)
  if (isTRUE(d$pi_dependent) && abs(pi - d$pi) > 1e-12)
    stop("The source fit depends on pi; refit it at the requested pi.")
  theta <- d$theta
  if (!is.numeric(theta) || any(!is.finite(theta)) || is.null(names(theta)) ||
      anyDuplicated(names(theta))) stop("The source descriptor needs uniquely named finite parameters.")
  bootstrap_note <- NULL
  if (uncertainty == "auto") {
    available_draws <- d$bootstrap_draws
    enough <- is.matrix(available_draws) && nrow(available_draws) >= 10L &&
      all(names(theta) %in% colnames(available_draws)) &&
      sum(apply(is.finite(available_draws), 1L, all)) >= max(10L, ceiling(.8 * nrow(available_draws)))
    uncertainty <- if (enough) "bootstrap" else "lavaan"
    if (!enough && (!is.null(available_draws) || isTRUE(d$has_legacy_bootstrap))) {
      bootstrap_note <- "Insufficient stored full bootstrap draws for family projections; auto uses joint lavaan delta inference. Request uncertainty = 'bootstrap' to require bootstrap inference."
      warning(bootstrap_note, call. = FALSE)
    }
  }
  if (pi_se > 0 && isTRUE(d$correlated_external_pi) && is.null(parameter_pi_cov))
    stop("The source used correlated external pi. Supply parameter_pi_cov for every source parameter; the legacy core-only covariance is insufficient for reduced regressions and variance explained.")
  evaluate <- function(th, p = pi) .ravel_family_projection_values(
    d$moments(th, p), roles, outcome_roles, p, covariate_mapping_confirmed)
  point <- evaluate(theta)
  estimate <- point$value
  finite <- is.finite(estimate)
  Vout <- matrix(NA_real_, length(estimate), length(estimate),
    dimnames = list(names(estimate), names(estimate)))
  se <- rep(NA_real_, length(estimate))
  lower <- upper <- estimate * NA_real_
  errors <- character(0)
  ndraw <- 0L
  draws <- NULL
  method <- rep(if (uncertainty == "none") "not_requested" else "unavailable", length(estimate))
  if (uncertainty != "none") {
    V <- d$vcov
    if (uncertainty == "bootstrap") {
      draws <- d$bootstrap_draws
      if (!is.matrix(draws) || !all(names(theta) %in% colnames(draws)))
        stop("Bootstrap projection needs stored full parameter_draws; use uncertainty = 'lavaan' or refit with the updated bootstrap.")
      draws <- draws[, names(theta), drop = FALSE]
      complete <- apply(is.finite(draws), 1L, all)
      ndraw <- sum(complete)
      if (ndraw < max(10L, ceiling(.8 * nrow(draws)))) stop("Insufficient complete joint bootstrap draws.")
      draws <- draws[complete, , drop = FALSE]
      V <- stats::cov(draws)
    }
    V <- .ravel_pi_validate_covariance(V, names(theta), "Source parameter covariance")
    x <- c(theta, .external_pi = pi)
    joint <- matrix(0, length(x), length(x), dimnames = list(names(x), names(x)))
    joint[names(theta), names(theta)] <- V
    joint[".external_pi", ".external_pi"] <- pi_se^2
    if (!is.null(parameter_pi_cov)) {
      if (!is.numeric(parameter_pi_cov) || is.null(names(parameter_pi_cov)) ||
          anyDuplicated(names(parameter_pi_cov)) || !setequal(names(parameter_pi_cov), names(theta)) ||
          any(!is.finite(parameter_pi_cov)))
        stop("parameter_pi_cov must be finite and named for every source parameter.")
      if (pi_se == 0 && any(parameter_pi_cov != 0))
        stop("Nonzero parameter_pi_cov requires positive pi_se.")
      joint[names(theta), ".external_pi"] <- joint[".external_pi", names(theta)] <- parameter_pi_cov[names(theta)]
    }
    joint <- .ravel_pi_validate_covariance(joint, names(x), "Joint parameter/pi covariance")
    transform <- function(x) evaluate(x[names(theta)], unname(x[".external_pi"]))$value[finite]
    # Derivatives at pi=0 use a one-sided difference; negative hidden fractions
    # are not admissible and must never enter the correction.
    J <- tryCatch({
      base <- transform(x)
      jac <- matrix(NA_real_, length(base), length(x), dimnames = list(names(base), names(x)))
      for (j in seq_along(x)) {
        if (j == length(x) && pi_se == 0) {
          jac[, j] <- 0
          next
        }
        step <- 1e-5 * max(1, abs(x[j]))
        plus <- minus <- x
        plus[j] <- x[j] + step
        minus[j] <- x[j] - step
        if (j == length(x)) {
          step <- min(step, (1 - pi) / 4)
          plus[j] <- pi + step
          minus[j] <- pi - step
        }
        jac[, j] <- if (j == length(x) && pi < step)
          (transform(plus) - base) / step else (transform(plus) - transform(minus)) / (2 * step)
      }
      if (pi_se > 0 && isTRUE(d$pi_dependent)) {
        if (!is.function(d$pi_sensitivity)) stop("Pi-dependent source fit lacks a full-parameter sensitivity callback.")
        sensitivity <- d$pi_sensitivity(pi)[names(theta)]
        if (any(!is.finite(sensitivity))) stop("The pi-refit sensitivity is not finite.")
        jac[, ".external_pi"] <- jac[, ".external_pi"] +
          drop(jac[, names(theta), drop = FALSE] %*% sensitivity)
      }
      # Visible estimates retain inference conditional on supplied pi, as in
      # existing RAVEL reporting. Their pi column is exactly zero.
      jac[startsWith(rownames(jac), "visible:"), ".external_pi"] <- 0
      jac
    }, error = function(e) { errors <<- c(errors, conditionMessage(e)); NULL })
    if (!is.null(J)) {
      good <- apply(is.finite(J), 1L, all)
      idx <- which(finite)[good]
      cov <- J[good, , drop = FALSE] %*% joint %*% t(J[good, , drop = FALSE])
      Vout[idx, idx] <- (cov + t(cov)) / 2
      se[idx] <- sqrt(pmax(diag(cov), 0))
      z <- stats::qnorm((1 + ci_level) / 2)
      lower[idx] <- estimate[idx] - z * se[idx]
      upper[idx] <- estimate[idx] + z * se[idx]
      method[idx] <- paste0(if (uncertainty == "bootstrap") "bootstrap_covariance" else "lavaan",
        if (pi_se > 0) "_pi_delta_wald" else "_delta_wald")
      if (pi_se > 0) method[idx[startsWith(names(estimate)[idx], "visible:")]] <-
        paste0(if (uncertainty == "bootstrap") "bootstrap_covariance" else "lavaan", "_conditional_pi_delta_wald")
      if (any(!good)) errors <- c(errors, "Some projection derivatives were non-finite; their uncertainty is unavailable.")
    }
    if (uncertainty == "bootstrap" && pi_se == 0) {
      transformed <- t(vapply(seq_len(nrow(draws)), function(i) {
        tryCatch(evaluate(draws[i, ])$value, error = function(e) estimate * NA_real_)
      }, estimate))
      success <- apply(is.finite(transformed[, finite, drop = FALSE]), 1L, all)
      if (sum(success) >= max(10L, ceiling(.8 * nrow(draws)))) {
        vv <- stats::cov(transformed[success, finite, drop = FALSE])
        Vout[finite, finite] <- vv
        se[finite] <- sqrt(pmax(diag(vv), 0))
        ci <- apply(transformed[success, finite, drop = FALSE], 2, stats::quantile,
          probs = c((1 - ci_level) / 2, (1 + ci_level) / 2), names = FALSE)
        lower[finite] <- ci[1, ]
        upper[finite] <- ci[2, ]
        method[finite] <- "bootstrap_percentile"
      } else errors <- c(errors, "Too few compatible transformed bootstrap draws for percentile intervals; delta intervals retained.")
    }
  }
  scale <- sub(":.*$", "", names(estimate))
  parameter <- sub("^[^:]+:", "", names(estimate))
  metrics <- c("genetic_variance", "genetic_R2", "conditional_genetic_variance", "partial_R2")
  status <- ifelse(is.finite(estimate), "ok", point$invalid_reason %||% "unavailable")
  table <- data.frame(scale = scale, parameter = parameter, estimate = unname(estimate),
    se = se, ci_low = unname(lower), ci_high = unname(upper),
    p_value = ifelse(is.finite(se) & se > 0, 2 * stats::pnorm(-abs(estimate / se)), NA_real_),
    kind = ifelse(parameter %in% metrics, "variance", ifelse(grepl("^delta_[mf](:|$)", parameter), "contrast", "coefficient")),
    status = status, corrected_status = ifelse(scale == "visible", "not_applicable", status),
    inference_method = method,
    inference_status = ifelse(uncertainty == "none", "not_requested", ifelse(is.finite(se), "ok", "unavailable")),
    stringsAsFactors = FALSE, row.names = NULL)
  table$p_value[table$kind == "variance"] <- NA_real_
  extract <- function(kind) setNames(lapply(c("visible", "total_raw", "total_std"), function(sc) {
    k <- table$scale == sc & table$kind == kind
    setNames(table$estimate[k], table$parameter[k])
  }), c("visible", "total_raw", "total_std"))
  structure(list(table = table, coefficients = extract("coefficient"),
    contrasts = extract("contrast"), variance_explained = extract("variance"),
    vcov = Vout, genetic_roles = roles, outcome_roles = outcome_roles,
    N = d$N %||% object$N, outcome_variable = d$outcome_variable %||% object$outcome_variable,
    pi = pi, pi_se = pi_se, uncertainty = uncertainty, ci_level = ci_level,
    n_bootstrap = ndraw, correction_status = point$invalid_reason %||% "ok",
    uncertainty_errors = unique(errors), uncertainty_note = bootstrap_note,
    estimand = "Reduced linear projection on selected liabilities, retaining supplied covariates and interactions",
    covariance_assumptions = c(d$covariance_assumptions,
      if (pi > 0 && isTRUE(d$has_covariates) && isTRUE(covariate_mapping_confirmed))
        "Cov(T - M'G, W) = 0, including product predictors",
      if (pi > 0 && isTRUE(d$has_interactions))
        "Moderators independent of every family visible and hidden genetic component"),
    source_fit = object), class = c("ravel_family_regression", "list"))
}

#' Print a reduced family regression
#' @param x A \code{ravel_family_regression} object.
#' @param ... Additional arguments, currently unused.
#' @return The input object, invisibly.
#' @export
print.ravel_family_regression <- function(x, ...) {
  cat("RAVEL family regression: ", x$outcome_variable, "\n", sep = "")
  cat("Genetic model: ", paste(x$genetic_roles, collapse = ", "), "\n", sep = "")
  cat("Outcome predictors: ", paste(x$outcome_roles, collapse = ", "), "\n", sep = "")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Convert a reduced family regression to a data frame
#' @param x A \code{ravel_family_regression} object.
#' @param row.names Optional row names.
#' @param optional Passed to \code{as.data.frame}.
#' @param ... Additional arguments passed to \code{as.data.frame}.
#' @return The regression's tidy coefficient, contrast and variance table.
#' @export
as.data.frame.ravel_family_regression <- function(x, row.names = NULL,
                                                 optional = FALSE, ...) {
  as.data.frame(x$table, row.names = row.names, optional = optional, ...)
}
