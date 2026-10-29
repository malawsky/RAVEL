# Family-member selection is separate from outcome-predictor selection.  The
# legacy trio fitter remains the implementation for every complete-trio model.

.ravel_family_roles <- function(x) {
  allowed <- c("child", "mother", "father", "mat_sib", "pat_sib")
  if (!is.character(x) || !length(x) || anyNA(x) || anyDuplicated(x) ||
      any(!x %in% allowed))
    stop("genetic_roles must contain distinct roles from: ", paste(allowed, collapse = ", "))
  allowed[allowed %in% x]
}

# Validate resampling before dispatching to either fitter. PLINK files commonly
# contain a constant FID; it must never produce zero-width bootstrap intervals.
.ravel_family_resampling <- function(dat, args) {
  cluster <- args[["cluster"]]
  validate_unit <- function(name) {
    if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name) ||
        !name %in% names(dat)) stop("The declared resampling/cluster column is missing or invalid.")
    ids <- dat[[name]]
    if (anyNA(ids) || any(as.character(ids) == "") ||
        (is.numeric(ids) && any(!is.finite(ids))))
      stop("Resampling/cluster IDs must be observed and finite: ", name)
    if (length(unique(ids)) < 2L)
      stop("Resampling requires at least two independent units in ", name,
        ". For independent rows with a dummy FID, explicitly use cluster_var = NULL and cluster = NULL.")
  }
  if (!is.null(cluster)) validate_unit(cluster)
  n <- args[["n_boot"]] %||% 0L
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != as.integer(n))
    stop("n_boot must be a nonnegative integer.")
  if (n > 0L) {
    explicit <- "cluster_var" %in% names(args)
    unit <- if (explicit) args[["cluster_var"]] else
      cluster %||% if ("FID" %in% names(dat)) "FID" else NULL
    if (explicit && is.null(unit)) {
      if (!is.null(cluster)) stop("Row bootstrap cannot be combined with clustered inference; supply the cluster resampling unit.")
      dat <- as.data.frame(dat)
      unit <- ".ravel_bootstrap_row"
      while (unit %in% names(dat)) unit <- paste0(unit, "_")
      dat[[unit]] <- seq_len(nrow(dat))
    }
    if (!is.null(unit)) {
      validate_unit(unit)
      args[["cluster_var"]] <- unit
    }
  }
  list(dat = dat, args = args)
}

.ravel_adult_pair <- function(a, b) {
  pair <- paste(sort(c(a, b)), collapse = ":")
  unname(c("f:m" = "rhoMF", "m:ms" = "kap_m", "f:fs" = "kap_f",
    "f:ms" = "chi_ms", "fs:m" = "chi_fs", "fs:ms" = "zeta")[pair])
}

.ravel_adult_history_variance <- function(r, n) {
  if (n == 1L) return(1)
  # Recursion avoids cancellation around r = 0 or r = 1.
  v <- 1
  for (j in seq_len(n - 1L)) v <- .5 * (1 + r) * v + .5
  v
}

.ravel_adult_state <- function(Sigma, roles, pi, mode, n_generations,
                               parent_visible_variance = NULL) {
  pair <- all(c("mother", "father") %in% roles)
  sib_pairs <- list(c("mother", "mat_sib"), c("father", "pat_sib"))
  sib_pairs <- Filter(function(x) all(x %in% roles), sib_pairs)
  rho <- if (pair) Sigma["mother", "father"] else NA_real_
  kappa <- if (length(sib_pairs)) mean(vapply(sib_pairs,
    function(x) Sigma[x[1], x[2]], numeric(1))) else NA_real_
  fail <- function(reason) list(valid = FALSE, invalid_reason = reason)
  if (!is.finite(pi) || pi < 0 || pi >= 1) return(fail("invalid_pi"))

  if (!is.null(parent_visible_variance)) {
    vG <- parent_visible_variance
    source <- "supplied_parent_visible_variance"
  } else if (mode == "estimated_state" || !pair) {
    if (!is.finite(kappa)) return(fail("variance_state_not_identified_without_child_or_sibling"))
    if (mode == "n_gen" && n_generations == 1L) {
      vG <- 1
      source <- "founder_adults"
    } else {
      vG <- 1 / (2 * (1 - kappa))
      source <- "sibling_correlation"
    }
  } else {
    # The adult generation is n - 1, matching the existing trio convention.
    visible_rho <- function(r) {
      V <- if (mode == "equilibrium") 1 / (1 - r) else
        .ravel_adult_history_variance(r, n_generations)
      r * (1 - pi) * V / (pi + (1 - pi) * V)
    }
    interval <- c(-1 + 1e-8, 1 - 1e-8)
    values <- vapply(interval, visible_rho, numeric(1)) - rho
    if (any(!is.finite(values)) || prod(sign(values)) > 0)
      return(fail("spouse_correlation_incompatible_with_pi_and_history"))
    r <- stats::uniroot(function(r) visible_rho(r) - rho, interval,
      tol = 1e-11)$root
    V <- if (mode == "equilibrium") 1 / (1 - r) else
      .ravel_adult_history_variance(r, n_generations)
    vG <- pi + (1 - pi) * V
    source <- paste0(mode, "_spouse_correlation")
  }
  V <- (vG - pi) / (1 - pi)
  if (!is.finite(vG) || !is.finite(V) || vG <= .5 || V <= .5)
    return(fail("sibling_variance_state_incompatible_with_pi"))
  r <- if (pair) rho * vG / (vG - pi) else NA_real_
  if (pair && (!is.finite(r) || abs(r) >= 1))
    return(fail("spouse_correlation_incompatible_with_variance_state"))
  if (!pair && mode == "equilibrium") r <- 1 - 1 / V
  if (!pair && mode == "n_gen" && n_generations > 1L) {
    interval <- c(-1 + 1e-8, 1 - 1e-8)
    values <- vapply(interval, .ravel_adult_history_variance, numeric(1),
      n = n_generations) - V
    if (prod(sign(values)) > 0)
      return(fail("sibling_variance_state_incompatible_with_generation_history"))
    r <- stats::uniroot(function(r)
      .ravel_adult_history_variance(r, n_generations) - V,
      interval, tol = 1e-11)$root
  }
  list(valid = TRUE, visible_parent_variance = vG, total_parent_variance = V,
    rho_total = r, rho_visible = rho, kappa_visible = kappa,
    state_source = source, current_assortment_observed = pair,
    history_imposed = mode != "estimated_state",
    generation = if (mode == "n_gen") n_generations - 1L else NA_integer_,
    state_fixed = !is.null(parent_visible_variance))
}

.ravel_adult_correction <- function(beta, Sigma, roles, pi, mode, n_generations,
                                    parent_visible_variance = NULL) {
  n <- length(roles)
  state <- .ravel_adult_state(Sigma, roles, pi, mode, n_generations,
    parent_visible_variance)
  bad <- list(beta_total_raw = setNames(rep(NA_real_, n), roles),
    Sigma_TT = matrix(NA_real_, n, n, dimnames = list(roles, roles)),
    M = matrix(NA_real_, n, n, dimnames = list(roles, roles)), state = state,
    invalid_reason = state$invalid_reason)
  if (!isTRUE(state$valid)) return(bad)
  w <- c(sqrt(1 - pi), sqrt(pi))
  V <- state$total_parent_variance
  S <- diag(2L) + (V - 1) * tcrossprod(w)
  A <- S - diag(2L) / 2
  r <- state$rho_total
  # A same-lineage sibling pair does not require a current spouse parameter.
  K <- if (is.finite(r)) r * V * tcrossprod(w) else matrix(0, 2L, 2L)
  k <- 1 - 1 / (2 * V)
  line <- c(mother = "m", mat_sib = "m", father = "f", pat_sib = "f")
  issib <- roles %in% c("mat_sib", "pat_sib")
  GH <- HH <- matrix(0, n, n, dimnames = list(roles, roles))
  sc <- sqrt(diag(S))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    B <- if (i == j) S else if (line[roles[i]] == line[roles[j]]) A else
      K * k^(as.integer(issib[i]) + as.integer(issib[j]))
    B <- B / outer(sc, sc)
    GH[i, j] <- B[1, 2]
    HH[i, j] <- B[2, 2]
  }
  joint <- rbind(cbind(Sigma, GH), cbind(t(GH), HH))
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  ej <- eigen(joint, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-10 * max(1, max(ev)) ||
      (pi > 0 && min(ej) < -1e-8 * max(1, max(ej)))) {
    bad$invalid_reason <- "incompatible_family_hidden_covariance"
    return(bad)
  }
  a <- sqrt(1 - pi) * sc[1]
  b <- sqrt(pi) * sc[2]
  M <- diag(a, n) + solve(Sigma, GH) * b
  dimnames(M) <- list(roles, roles)
  TT <- a^2 * Sigma + a * b * (GH + t(GH)) + b^2 * HH
  dimnames(TT) <- list(roles, roles)
  list(beta_total_raw = setNames(drop(solve(M, beta)), roles),
    Sigma_TT = TT, M = M, state = state, invalid_reason = NULL)
}

.ravel_family_identified <- function(fit, tolerance = 1e-7) {
  tryCatch({
    if (!isTRUE(lavaan::lavInspect(fit, "converged")) ||
        !suppressWarnings(isTRUE(lavaan::lavInspect(fit, "post.check")))) return(FALSE)
    V <- as.matrix(lavaan::vcov(fit))
    if (any(!is.finite(V)) || any(diag(V) < -1e-8)) return(FALSE)
    scales <- sqrt(pmax(abs(diag(V)), max(abs(diag(V))) * .Machine$double.eps))
    if (any(!is.finite(scales)) || any(scales <= 0)) return(FALSE)
    ev <- eigen((V + t(V)) / (2 * outer(scales, scales)),
      symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-7 * max(1, max(abs(ev)))) return(FALSE)
    D <- lavaan::lavInspect(fit, "delta")
    if (is.list(D)) D <- do.call(rbind, D)
    D <- sweep(as.matrix(D), 2, scales, "*")
    C <- fit@Model@con.jac
    eq <- attr(C, "ceq.idx")
    inequality <- attr(C, "cin.idx")
    inactive <- attr(C, "inactive.idx")
    if (length(setdiff(inequality, inactive))) return(FALSE)
    normalize <- function(x) {
      s <- sqrt(rowSums(x^2)); x[s > 0, , drop = FALSE] / s[s > 0]
    }
    tangent <- diag(ncol(D))
    if (is.matrix(C) && length(eq)) {
      C <- normalize(sweep(C[eq, , drop = FALSE], 2, scales, "*"))
      sv <- svd(C, nu = 0, nv = ncol(C))
      rank <- sum(sv$d > tolerance * max(sv$d))
      if (rank == ncol(C)) return(FALSE)
      tangent <- sv$v[, seq.int(rank + 1L, ncol(C)), drop = FALSE]
    }
    singular <- svd(normalize(D) %*% tangent, nu = 0, nv = 0)$d
    length(singular) == ncol(tangent) && max(singular) > 0 &&
      sum(singular > tolerance * max(singular)) == ncol(tangent)
  }, error = function(e) FALSE)
}

.ravel_adult_model <- function(dat, roles, fixed_effects, mode, pi, n_generations,
                                fixed_crossperson, fixed_crossperson_corr,
                                parent_visible_variance) {
  codes <- unname(.ravel_interaction_role_codes[roles])
  lv <- paste0("G", codes)
  lines <- character(0)
  # Outcome covariances identify the loading-ratio direction in free models.
  reliability <- matrix(.5, length(roles), 2L)
  for (i in seq_along(roles)) {
    p <- paste0("PGS", 1:2, "_", roles[i])
    cv <- stats::cov(dat[c(p, "Y")], use = "pairwise.complete.obs")
    cand <- c(cv[1, 2] * cv[1, 3] / cv[2, 3] / cv[1, 1],
      cv[1, 2] * cv[2, 3] / cv[1, 3] / cv[2, 2])
    if (all(is.finite(cand))) reliability[i, ] <- pmin(.95, pmax(.05, cand))
  }
  if (mode != "none") reliability[,] <- rep(apply(reliability, 2, stats::median),
    each = length(roles))
  for (i in seq_along(roles)) {
    r <- roles[i]; c <- codes[i]
    for (j in 1:2) {
      col <- paste0("PGS", j, "_", r)
      v <- stats::var(dat[[col]], na.rm = TRUE)
      loading <- sqrt(v * reliability[i, j])
      error <- v * (1 - reliability[i, j])
      lines <- c(lines,
        sprintf("G%s =~ start(%.12g)*l%s%d*%s", c, loading, c, j, col),
        sprintf("%s ~~ start(%.12g)*th%d_%s*%s", col, error, j, c, col))
      if (mode != "none") lines <- c(lines,
        sprintf("l%s%d > 0.00000001", c, j), sprintf("th%d_%s > 0.00000001", j, c))
    }
    lines <- c(lines, sprintf("Y ~ b%s*G%s", c, c))
  }
  pairs <- utils::combn(seq_along(roles), 2L)
  valid_residual <- character(0)
  for (p in seq_len(ncol(pairs))) {
    i <- pairs[1, p]; j <- pairs[2, p]
    label <- .ravel_adult_pair(codes[i], codes[j])
    start <- if (grepl("^kap", label)) .55 else .1
    lines <- c(lines, sprintf("%s ~~ start(%g)*%s*%s", lv[i], start, label, lv[j]))
    for (v in 1:2) {
      residual <- if (label == "rhoMF") paste0("e", v, "_mf") else
        paste0("e", v, "_", codes[i], "_", codes[j])
      valid_residual <- c(valid_residual, residual)
      covariance <- fixed_crossperson[residual]
      correlation <- fixed_crossperson_corr[residual]
      if (mode != "none") {
        benchmark <- if (grepl("^kap", label)) .5 else 0
        if (length(correlation) && !is.na(correlation) && abs(correlation - benchmark) > 1e-10)
          stop("Constrained residual correlation conflicts with the pedigree benchmark: ", residual)
        correlation <- benchmark
      }
      has_cov <- length(covariance) && !is.na(covariance)
      has_cor <- length(correlation) && !is.na(correlation)
      if (has_cov && has_cor) stop("Residual specified as both covariance and correlation: ", residual)
      val <- if ((has_cov && covariance == 0) || (has_cor && correlation == 0)) "0" else residual
      lhs <- paste0("PGS", v, "_", roles[i]); rhs <- paste0("PGS", v, "_", roles[j])
      lines <- c(lines, sprintf("%s ~~ %s*%s", lhs, val, rhs))
      if (val != "0" && has_cov) lines <- c(lines,
        sprintf("%s == %.15g", residual, covariance))
      if (val != "0" && has_cor) lines <- c(lines,
        sprintf("%s == %.15g*sqrt(th%d_%s*th%d_%s)", residual, correlation,
          v, codes[i], v, codes[j]))
    }
  }
  for (values in list(fixed_crossperson, fixed_crossperson_corr)) {
    if (length(values) && (is.null(names(values)) || anyDuplicated(names(values)) ||
        any(!names(values) %in% valid_residual)))
      stop("Fixed residual names must identify a pair present in genetic_roles: ",
        paste(valid_residual, collapse = ", "))
  }
  if (mode != "none") {
    for (i in seq.int(2L, length(roles))) for (j in 1:2)
      lines <- c(lines, sprintf("l%s%d^2*th%d_%s == l%s%d^2*th%d_%s",
        codes[1], j, j, codes[i], codes[i], j, j, codes[1]))
    pair <- all(c("mother", "father") %in% roles)
    own <- c(if (all(c("mother", "mat_sib") %in% roles)) "kap_m",
      if (all(c("father", "pat_sib") %in% roles)) "kap_f")
    if (pair) lines <- c(lines, "rhoMF > -0.999999", "rhoMF < 0.999999")
    if (length(own)) {
      for (k in own) lines <- c(lines, paste(k, "> 0.00000001"), paste(k, "< 0.99999999"))
      if (length(own) == 2L) lines <- c(lines, "kap_m == kap_f", "chi_ms == chi_fs")
      vG <- paste0("(1/(2*(1-", own[1], ")))" )
      if (!is.null(parent_visible_variance)) {
        lines <- c(lines, sprintf("%s == %.15g", vG, parent_visible_variance))
      } else if (mode == "equilibrium" && pair) {
        lines <- c(lines, paste0(own[1], " == (1+rhoMF)/2"))
      } else if (mode == "n_gen" && n_generations == 1L) {
        lines <- c(lines, paste0(own[1], " == 0.5"))
      } else if (mode == "n_gen" && pair) {
        lines <- c(lines,
          sprintf("family_vG := %s", vG),
          sprintf("family_rT := rhoMF*%s/(%s-%.15g)", vG, vG, pi),
          sprintf("family_V := (%s-%.15g)/(1-%.15g)", vG, pi, pi),
          sprintf("family_V*(1-family_rT) == 1-family_rT*((1+family_rT)/2)^%d", n_generations - 1L),
          "family_rT > -0.999999", "family_rT < 0.999999")
      }
    }
  }
  if (length(fixed_effects)) {
    for (i in seq_along(fixed_effects)) lines <- c(lines,
      sprintf("Y ~ c_%d*%s", i, fixed_effects[i]))
    for (g in lv) lines <- c(lines, paste(g, "~~", paste(fixed_effects, collapse = " + ")))
    if (length(fixed_effects) > 1L) {
      pp <- utils::combn(fixed_effects, 2L)
      lines <- c(lines, paste(pp[1, ], "~~", pp[2, ]))
    }
  }
  paste(lines, collapse = "\n")
}

#' Fit RAVEL using selected genetic family members
#'
#' The full outcome covariance model is retained over genetic_roles. Requested
#' outcome_roles define a reduced total-liability projection, not zero paths in
#' the visible-score model. Complete trios use ravel_fit() unchanged.
#'
#' @param dat Data frame with two independent PGS versions for each requested
#'   genetic member, the outcome, and any supplied covariates.
#' @param genetic_roles Character vector drawn from \code{child}, \code{mother},
#'   \code{father}, \code{mat_sib}, and \code{pat_sib}. A child requires both
#'   parents. A parental sibling requires that parent's scores. Child-free fits
#'   require a full-sibling pair or both partners.
#' @param outcome_roles Nonempty subset of \code{genetic_roles} defining the
#'   reported regression. Other members remain in the joint genetic and outcome
#'   covariance model as auxiliary variables; they are marginalized from the
#'   requested regression.
#' @param parent_visible_variance Optional fixed, founder-scaled visible genetic
#'   variance of the adult generation, greater than 0.5. Only available for
#'   child-free \code{estimated_state} fits. It is required for partner-only
#'   estimated-state analyses and is treated as known, without a supplied-state
#'   standard error.
#' @param ... Named arguments accepted by \code{ravel_fit}, including PGS column
#'   mappings, \code{col_outcome}, \code{fit_constraints}, \code{am_correction},
#'   \code{pi}, \code{pi_se}, covariates, interactions, clustering and bootstrap
#'   settings. Child-free models do not support \code{fix_parent_paths},
#'   \code{crossperson = "equal_pc"}, or trio-only \code{core_pi_cov}.
#' @details
#' The two PGS versions must measure the same visible liability and have
#' independent cross-version errors, as in RAVEL. In child-free estimated-state
#' fits, a full-sibling visible correlation \eqn{\kappa} identifies adult visible
#' variance \eqn{V_G=1/[2(1-\kappa)]}; the supplied founder hidden fraction gives
#' \eqn{V_T=(V_G-\pi)/(1-\pi)}. Current spouse assortment is estimated only when
#' both partners are observed. Sibling-only correction does not need current
#' spouse assortment. Constant-history or equilibrium assortment inferred without
#' partners is an assumption-derived quantity, not an observed current mating
#' correlation. Generation counts retain the existing convention: adults are
#' generation \code{am_n_generations - 1} relative to random-mating founders.
#'
#' Free measurement models use the outcome jointly to identify score reliabilities;
#' a score-only first-stage fit is not substituted. Local identification and the
#' positive semidefiniteness of corrected covariance blocks are checked. Continuous
#' models use FIML for missing observations; binary models use listwise deletion.
#' The existing restrictions on binary clustering and interactions are retained.
#' @return Complete-trio calls return the existing \code{ravel_results} object with
#'   added \code{genetic_roles}, \code{outcome_roles}, and \code{regression} fields.
#'   Child-free calls return a \code{ravel_family_fit} object containing
#'   \code{converged}, the lavaan \code{fit}, a reusable \code{family_descriptor},
#'   and the requested \code{regression}. A failed fit has \code{converged = FALSE}
#'   and an explanatory \code{error}. Use \code{ravel_family_regression} to obtain
#'   further outcome-predictor subsets without refitting the same joint model.
#' @seealso \code{\link{ravel_fit}}, \code{\link{ravel_family_regression}}
#' @export
ravel_family_fit <- function(dat, genetic_roles = c("child", "mother", "father"), ...,
                             outcome_roles = genetic_roles,
                             parent_visible_variance = NULL) {
  roles <- .ravel_family_roles(genetic_roles)
  if (!is.character(outcome_roles) || !length(outcome_roles) || anyNA(outcome_roles) ||
      anyDuplicated(outcome_roles) || !all(outcome_roles %in% roles))
    stop("outcome_roles must be a nonempty, distinct subset of genetic_roles.")
  args <- list(...)
  if (length(args) && (is.null(names(args)) || anyNA(names(args)) ||
      any(!nzchar(names(args))) || anyDuplicated(names(args))))
    stop("All fitting arguments in ... must have distinct, nonempty names.")
  unknown <- setdiff(names(args), names(formals(ravel_fit)))
  if (length(unknown)) stop("Unknown fitting arguments: ", paste(unknown, collapse = ", "))
  if ("dat" %in% names(args)) stop("Supply dat only once.")
  sampling <- .ravel_family_resampling(dat, args)
  dat <- sampling$dat
  args <- sampling$args
  if (!is.null(parent_visible_variance) && (!is.numeric(parent_visible_variance) ||
      length(parent_visible_variance) != 1L || !is.finite(parent_visible_variance) ||
      parent_visible_variance <= .5))
    stop("parent_visible_variance must be a finite founder-scaled adult visible variance greater than 0.5.")
  if ("child" %in% roles) {
    if (!all(c("mother", "father") %in% roles))
      stop("Models containing a child currently require both parents in genetic_roles.")
    if (!is.null(parent_visible_variance))
      stop("parent_visible_variance is for child-free models; a complete trio supplies its state through the existing correction.")
    siblings <- roles[roles %in% c("mat_sib", "pat_sib")]
    mode <- if (length(siblings) == 2L) "both" else if (identical(siblings, "mat_sib"))
      "maternal" else if (identical(siblings, "pat_sib")) "paternal" else "none"
    if (!is.null(args[["parental_sibs"]]) && args[["parental_sibs"]] != mode)
      stop("parental_sibs contradicts genetic_roles.")
    args[["parental_sibs"]] <- mode
    result <- do.call(ravel_fit, c(list(dat = dat), args))
    result$genetic_roles <- roles
    result$outcome_roles <- outcome_roles
    if (isTRUE(result$converged)) {
      result$family_descriptor <- .ravel_family_legacy_descriptor(result)
      result$regression <- ravel_family_regression(result,
        outcome_roles, ci_level = args[["ci_level"]] %||% .95,
        covariate_mapping_confirmed = args[["covariate_mapping_confirmed"]] %||% FALSE)
    }
    class(result) <- c("ravel_family_fit", class(result))
    return(result)
  }
  .ravel_adult_fit(dat, roles, outcome_roles, parent_visible_variance, args)
}

.ravel_adult_fit <- function(dat, roles, outcome_roles, parent_visible_variance, args) {
  pair <- all(c("mother", "father") %in% roles)
  own_pair <- all(c("mother", "mat_sib") %in% roles) ||
    all(c("father", "pat_sib") %in% roles)
  if (length(roles) < 2L || (!pair && !own_pair))
    stop("Child-free genetic_roles require a full-sibling pair or both partners.")
  if ("mat_sib" %in% roles && !"mother" %in% roles ||
      "pat_sib" %in% roles && !"father" %in% roles)
    stop("Include each parental sibling's own parent in genetic_roles; outcome_roles may omit that parent.")
  pi <- args[["pi"]] %||% args[["h"]] %||% 0
  if (!is.null(args[["pi_se"]]) && is.null(args[["pi"]]) && is.null(args[["h"]]))
    stop("Supply pi (or h) explicitly when using pi_se.")
  .ravel_pi_se_variance(pi, args[["pi_se"]] %||% 0)
  if (!is.null(args[["pi"]]) && !is.null(args[["h"]]) && args[["pi"]] != args[["h"]])
    stop("h and pi specify the same hidden fraction; supply matching values.")
  if (!is.null(args[["core_pi_cov"]]))
    stop("For child-free fits supply full parameter_pi_cov to ravel_family_regression(), rather than trio core_pi_cov.")
  fit_mode <- match.arg(args[["fit_constraints"]] %||% "none",
    c("none", "equilibrium", "n_gen", "estimated_state"))
  correction <- match.arg(args[["am_correction"]] %||% if (fit_mode == "none") "n_gen" else fit_mode,
    c("n_gen", "equilibrium", "estimated_state"))
  if (fit_mode != "none" && correction != fit_mode)
    stop("fit_constraints and am_correction must match when fitting constrained child-free models.")
  if (correction == "estimated_state" && fit_mode != "estimated_state")
    stop("am_correction = 'estimated_state' requires fit_constraints = 'estimated_state'.")
  if (correction == "estimated_state" && !own_pair && is.null(parent_visible_variance))
    stop("Partner-only estimated_state needs parent_visible_variance: without child or sibling scores the accumulated state is not identified.")
  if (!is.null(parent_visible_variance) && correction != "estimated_state")
    stop("parent_visible_variance is supported only with estimated_state; it cannot override a generation or equilibrium history.")
  n_generations <- args[["am_n_generations"]] %||% 1L
  if (!is.numeric(n_generations) || length(n_generations) != 1L ||
      !is.finite(n_generations) || n_generations < 1 || n_generations != as.integer(n_generations))
    stop("am_n_generations must be a positive integer.")
  outcome <- match.arg(args[["outcome"]] %||% "continuous", c("continuous", "binary_liability"))
  if (isTRUE(args[["fix_parent_paths"]]))
    stop("fix_parent_paths applies to legacy child-outcome models; choose outcome_roles for adult regressions.")
  if (!is.null(args[["measurement"]]) && args[["measurement"]] != "free_loadings")
    stop("Only measurement = 'free_loadings' is supported.")
  if (!is.null(args[["crossperson"]]) && !args[["crossperson"]] %in% c("free", "equal_pc"))
    stop("crossperson must be 'free' or 'equal_pc'.")
  if (identical(args[["crossperson"]], "equal_pc"))
    stop("crossperson = 'equal_pc' requires a child and is not applicable to child-free models.")
  if (!is.null(args[["parental_sibs"]])) {
    expected <- if (all(c("mat_sib", "pat_sib") %in% roles)) "both" else if
      ("mat_sib" %in% roles) "maternal" else if ("pat_sib" %in% roles) "paternal" else "none"
    if (args[["parental_sibs"]] != expected) stop("parental_sibs contradicts genetic_roles.")
  }
  dat <- as.data.frame(dat)
  colmap <- setNames(vapply(roles, function(r) args[[paste0("col_pgs1_", r)]] %||%
    paste0("PGS1_", r), character(1)), paste0("PGS1_", roles))
  colmap <- c(colmap, setNames(vapply(roles, function(r) args[[paste0("col_pgs2_", r)]] %||%
    paste0("PGS2_", r), character(1)), paste0("PGS2_", roles)),
    Y = args[["col_outcome"]] %||% "Y")
  missing_cols <- setdiff(unname(colmap), names(dat))
  if (length(missing_cols)) stop("Missing model columns: ", paste(missing_cols, collapse = ", "))
  original <- dat
  for (nm in names(colmap)) dat[[nm]] <- original[[colmap[[nm]]]]
  fixed_effects <- args[["fixed_effects"]] %||% character(0)
  if (any(!fixed_effects %in% names(dat))) stop("Missing fixed-effect columns.")
  if (any(make.names(fixed_effects) != fixed_effects)) stop("Covariate names must be syntactic R names.")
  for (nm in unique(c(names(colmap), fixed_effects))) {
    if (!is.numeric(dat[[nm]]) || any(!is.na(dat[[nm]]) & !is.finite(dat[[nm]])))
      stop("Model columns must be numeric and finite or NA: ", nm)
    if (sum(!is.na(dat[[nm]])) < 2L || stats::var(dat[[nm]], na.rm = TRUE) <= 0)
      stop("Model columns must have observed variation: ", nm)
  }
  if (nrow(dat) < 50L) stop("At least 50 families are required.")
  if (outcome == "binary_liability") {
    if (!setequal(unique(stats::na.omit(dat$Y)), c(0, 1))) stop("Binary Y must contain 0 and 1.")
    keep <- stats::complete.cases(dat[unique(c(names(colmap), fixed_effects))])
    if (any(!keep)) warning("Binary family models use listwise deletion; dropping ", sum(!keep), " rows.")
    dat <- dat[keep, , drop = FALSE]
    if (nrow(dat) < 50L) stop("At least 50 complete families are required for binary models.")
  }
  interaction_roles <- args[["interaction_roles"]] %||% roles
  if (!is.null(args[["interactions"]]) && pi > 0 && !all(roles %in% interaction_roles))
    stop("With pi > 0, interaction_roles must include every genetic role.")
  interaction_spec <- .ravel_interaction_spec(args[["interactions"]], interaction_roles,
    args[["interaction_center"]], dat, outcome, available_roles = roles)
  if (!is.null(interaction_spec)) {
    fixed_effects <- unique(c(fixed_effects, interaction_spec$moderators))
    dat <- .ravel_interaction_data(dat, interaction_spec)
  }
  cluster <- args[["cluster"]]
  cluster_var <- args[["cluster_var"]] %||% cluster %||% "FID"
  if (!is.null(cluster)) {
    if (outcome != "continuous") stop("Cluster-robust fitting currently requires a continuous outcome.")
    if (!is.character(cluster) || length(cluster) != 1L || !cluster %in% names(dat) || anyNA(dat[[cluster]]))
      stop("cluster must name a complete identifier column.")
    if (!cluster_var %in% names(dat) || anyNA(dat[[cluster_var]])) stop("Invalid cluster_var.")
    units <- unique(data.frame(cluster = dat[[cluster]], unit = dat[[cluster_var]]))
    if (anyDuplicated(units$cluster)) stop("The bootstrap unit must equal or nest the cluster.")
  }
  n_boot <- args[["n_boot"]] %||% 0L
  if (length(n_boot) != 1L || !is.finite(n_boot) || n_boot < 0 || n_boot != as.integer(n_boot))
    stop("n_boot must be a nonnegative integer.")
  fixed_cov <- args[["fixed_crossperson"]] %||% numeric(0)
  fixed_cor <- args[["fixed_crossperson_corr"]] %||% numeric(0)
  if (any(!is.finite(fixed_cov) & !is.na(fixed_cov)) || any(!is.finite(fixed_cor) & !is.na(fixed_cor)) ||
      any(abs(fixed_cor) > 1, na.rm = TRUE)) stop("Invalid fixed residual covariance/correlation.")
  if (fit_mode != "none" && any(!is.na(fixed_cov)))
    stop("Constrained family models fix residual correlations; do not also fix residual covariances.")
  model_at <- function(h) {
    syntax <- .ravel_adult_model(dat, roles, fixed_effects, fit_mode, h, n_generations,
      fixed_cov, fixed_cor, parent_visible_variance)
    if (!is.null(interaction_spec)) syntax <- paste(syntax,
      .ravel_interaction_syntax(interaction_spec, fixed_effects, roles), sep = "\n")
    syntax
  }
  .ensure_lavaan_runtime()
  fit_once <- function(data, h = pi, start = NULL) {
    model <- model_at(h)
    input <- if (outcome == "binary_liability" && length(fixed_effects))
      .ravel_binary_moments(data, model) else list(data = data)
    call <- c(list(model = model, std.lv = TRUE, fixed.x = FALSE, conditional.x = FALSE,
      estimator = if (outcome == "continuous") "MLR" else "WLSMV",
      missing = if (outcome == "continuous") "fiml" else "listwise",
      ordered = if (outcome == "continuous") NULL else "Y", check.gradient = FALSE,
      warn = FALSE), input)
    if (!is.null(cluster)) call$cluster <- cluster
    if (!is.null(start)) call$start <- start
    suppressWarnings(tryCatch(do.call(lavaan::sem, call), error = function(e) NULL))
  }
  if (isTRUE(args[["verbose"]])) cat(model_at(pi), "\n")
  if (!is.null(args[["seed"]])) set.seed(args[["seed"]])
  fit <- fit_once(dat)
  valid <- function(f) !is.null(f) && .ravel_family_identified(f) &&
    (fit_mode == "none" || .ravel_constraint_fit_valid(f))
  if (!valid(fit)) fit <- fit_once(dat, start = "simple")
  if (!valid(fit)) return(structure(list(converged = FALSE, fit = fit,
    genetic_roles = roles, outcome_roles = outcome_roles,
    error = "Family SEM did not produce an admissible, locally identified fit with finite uncertainty."),
    class = "ravel_family_fit"))
  theta <- lavaan::coef(fit)
  names(theta) <- make.unique(names(theta))
  Vtheta <- as.matrix(lavaan::vcov(fit)); dimnames(Vtheta) <- list(names(theta), names(theta))
  pt <- lavaan::parameterTable(fit)
  codes <- unname(.ravel_interaction_role_codes[roles])
  moment_base <- .ravel_family_lavaan_design(fit, roles, interaction_spec, fixed_effects)
  moments <- function(theta, pi) {
    base <- moment_base(theta)
    correction <- .ravel_adult_correction(base$beta_visible, base$Sigma_GG, roles,
      pi, correction, n_generations, parent_visible_variance)
    # Distinct names are essential: correction's M is a coefficient map.
    for (nm in names(correction)) base[[nm]] <- correction[[nm]]
    if (!is.null(base$interaction_beta_visible)) {
      base$interaction_beta_total_raw <- if (is.null(correction$invalid_reason))
        solve(correction$M, base$interaction_beta_visible) else
        base$interaction_beta_visible * NA_real_
    }
    base
  }
  pi_dependent <- fit_mode == "n_gen" && n_generations > 1L && pair && own_pair &&
    is.null(parent_visible_variance)
  sensitivity <- if (pi_dependent) function(h = pi) {
    evaluate <- function(p) {
      fp <- fit_once(dat, p, start = unname(theta))
      if (!valid(fp)) stop("Pi sensitivity refits were not admissible.")
      setNames(lavaan::coef(fp), names(theta))
    }
    .ravel_family_pi_sensitivity(evaluate, h, theta)
  } else NULL
  sensitivity <- .ravel_family_cached_sensitivity(sensitivity, names(theta))
  draws <- NULL
  if (n_boot > 0L) {
    draws <- matrix(NA_real_, n_boot, length(theta), dimnames = list(NULL, names(theta)))
    units <- if (cluster_var %in% names(dat)) unique(dat[[cluster_var]]) else seq_len(nrow(dat))
    for (b in seq_len(n_boot)) {
      selected <- sample(seq_along(units), length(units), replace = TRUE)
      if (cluster_var %in% names(dat)) {
        pieces <- lapply(seq_along(selected), function(j) {
          d <- dat[dat[[cluster_var]] == units[selected[j]], , drop = FALSE]
          if (!is.null(cluster)) d[[cluster]] <- paste0(j, ":", d[[cluster]])
          d
        })
        db <- do.call(rbind, pieces)
      } else db <- dat[selected, , drop = FALSE]
      fb <- fit_once(db)
      if (valid(fb)) draws[b, ] <- lavaan::coef(fb)
    }
  }
  descriptor <- list(theta = theta, vcov = Vtheta, roles = roles, pi = pi,
    pi_se = args[["pi_se"]], ci_level = args[["ci_level"]] %||% .95,
    model_source = "child_free_family_sem", has_covariates = length(fixed_effects) > 0L,
    has_interactions = !is.null(interaction_spec), moments = moments,
    bootstrap_draws = draws, pi_dependent = pi_dependent, pi_sensitivity = sensitivity,
    parameter_pi_cov = NULL, N = as.integer(lavaan::lavInspect(fit, "nobs")),
    outcome_variable = args[["col_outcome"]] %||% "Y", outcome = outcome,
    covariance_assumptions = "RAVEL two-component family covariance and full visible-to-total projection model")
  result <- structure(list(converged = TRUE, fit = fit, family_descriptor = descriptor,
    genetic_roles = roles, outcome_roles = outcome_roles, N = as.integer(lavaan::lavInspect(fit, "nobs")),
    outcome = outcome, outcome_variable = args[["col_outcome"]] %||% "Y",
    fixed_effects = fixed_effects, interactions = interaction_spec,
    cluster = cluster, n_clusters = if (is.null(cluster)) NULL else length(unique(dat[[cluster]])),
    h = pi, pi = pi, am_correction = correction, am_n_generations = n_generations,
    constraint_spec = list(mode = fit_mode, depends_on_pi = pi_dependent),
    parent_visible_variance = parent_visible_variance, model = model_at(pi)),
    class = "ravel_family_fit")
  result$regression <- ravel_family_regression(result, outcome_roles,
    ci_level = descriptor$ci_level,
    covariate_mapping_confirmed = args[["covariate_mapping_confirmed"]] %||% FALSE)
  result
}

#' @export
print.ravel_family_fit <- function(x, ...) {
  if (!isTRUE(x$converged)) {
    cat("RAVEL family fit did not converge: ", x$error, "\n", sep = "")
  } else print(x$regression, ...)
  invisible(x)
}

#' @export
as.data.frame.ravel_family_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (!isTRUE(x$converged)) stop("A failed family fit has no regression table.")
  as.data.frame(x$regression$table, row.names = row.names, optional = optional, ...)
}
