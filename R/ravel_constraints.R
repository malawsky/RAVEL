# Optional, strongly constrained measurement and inheritance models.
# The unconstrained ravel_fit() path does not call these helpers.

#' Build constraints on reliability and Mendelian variance dynamics
#' @noRd
.ravel_constraint_syntax <- function(mode, h, n_generations) {
  if (identical(mode, "none")) return("")
  out <- c(
    "# Distinct PGS reliabilities; unrestricted observed marginal scales",
    "eta_m == eta_f",
    "rhoMF > -0.999999", "rhoMF < 0.999999",
    "eta_m > 0.000001", "eta_m < 0.999999"
  )
  for (j in 1:2) {
    out <- c(out,
      sprintf("lf%d^2*th%d_m == lm%d^2*th%d_f", j,j,j,j),
      sprintf("lc%d^2*th%d_m == am_growth*lm%d^2*th%d_c", j,j,j,j))
    for (role in c("m", "f", "c")) {
      # Orient both scores toward the same genetic liability. Squared-odds
      # constraints alone would also admit opposite signs across family roles.
      out <- c(out, sprintf("l%s%d > 0.00000001", role, j),
        sprintf("th%d_%s > 0.00000001", j, role))
    }
  }
  if (identical(mode, "equilibrium")) {
    out <- c(out, "am_growth := 1 + 0*rhoMF",
      "eta_m == (1+rhoMF)/2")
  } else if (identical(mode, "estimated_state")) {
    # Estimate the accumulated state and current assortment separately. This
    # preserves Mendelian segregation but does not impose a constant-strength
    # integer-generation history or use the external hidden fraction in fitting.
    out <- c(out,
      "am_growth := ((1+rhoMF)/(2*eta_m))^2",
      "am_inverse_vG_parent := 2*am_growth-1-rhoMF",
      "am_vG_parent := 1/am_inverse_vG_parent",
      "am_inverse_vG_parent > 0.00000001",
      "am_inverse_vG_parent < 1.99999999")
  } else if (n_generations == 1L) {
    # Parents are founders; the first child generation has vGc/vG=1+rhoG/2.
    # The hidden fraction cancels, unlike in later generations.
    out <- c(out, "am_growth := 1+rhoMF/2",
      "4*eta_m^2*am_growth == (1+rhoMF)^2")
  } else {
    h_text <- .format_lavaan_numeric(h)
    # Constant total-liability spouse correlation and random-mating founders:
    # V_T,parent=(1-rT*((1+rT)/2)^(n-1))/(1-rT).
    # Eliminate rT and V_G,parent using the two fitted latent correlations.
    # This avoids an unidentified auxiliary/phantom factor in lavaan.
    out <- c(out,
      "am_growth := ((1+rhoMF)/(2*eta_m))^2",
      "am_vG_parent := 1/(2*am_growth-1-rhoMF)",
      sprintf("am_rho_total := rhoMF*am_vG_parent/(am_vG_parent-%s)", h_text),
      sprintf("am_vT_parent := (am_vG_parent-%s)/(1-%s)", h_text, h_text),
      sprintf(paste0("am_vT_parent*(1-am_rho_total) == ",
        "1-am_rho_total*((1+am_rho_total)/2)^%d"), n_generations-1L),
      "am_rho_total > -0.999999", "am_rho_total < 0.999999",
      "am_vG_parent > 0.00000001", "am_vT_parent > 0.00000001")
  }
  paste(out, collapse = "\n")
}

#' Starting values inside the finite-generation constraint manifold
#' @noRd
.ravel_constraint_starts <- function(model, dat, mode, h, n_generations) {
  # Starts use observed score moments only; no true reliability is supplied.
  score_names <- unlist(lapply(1:2, function(j)
    paste0("PGS", j, "_", c("mother", "father", "child"))))
  S <- stats::cov(dat[, score_names, drop = FALSE], use = "pairwise.complete.obs")
  sd_score <- sqrt(diag(S))
  R <- stats::cov2cor(S)
  product <- mean(c(R[1, 4], R[2, 5]), na.rm=TRUE)
  if (!is.finite(product)) product <- .5
  product <- max(.03, min(.95, product))
  rho <- mean(c(R[1, 5], R[2, 4]), na.rm=TRUE)/product
  if (!is.finite(rho)) rho <- .1
  rho <- max(-.8, min(.8, rho))
  # Match the observed spouse correlation by a one-dimensional starting-value
  # inversion only. Neither rho nor either reliability is fixed during fitting.
  if (mode == "equilibrium") {
    eta <- (1+rho)/2
    growth <- 1
  } else if (mode == "estimated_state") {
    # Cross-score covariances identify eta without knowing either reliability.
    # This initialization is also independent of supplied h and generations.
    child_product <- R[3,6]
    cross_child_parent <- mean(c(R[3,4],R[3,5]),na.rm=TRUE)
    cross_parent_child <- mean(c(R[1,6],R[2,6]),na.rm=TRUE)
    eta_squared <- cross_child_parent*cross_parent_child/(product*child_product)
    inverse_variance <- if (is.finite(eta_squared) && eta_squared > 0)
      (1+rho)^2/(2*eta_squared)-1-rho else NA_real_
    # A midpoint between the founder and equilibrium growth factors is a
    # stable interior fallback when sampling noise spoils the moment inversion.
    if (!is.finite(inverse_variance)) inverse_variance <- 1-rho/2
    inverse_variance <- max(.05,min(1.95,inverse_variance))
    growth <- (inverse_variance+1+rho)/2
    eta <- (1+rho)/(2*sqrt(growth))
  } else if (n_generations == 1L) {
    growth <- 1+rho/2
    eta <- (1+rho)/(2*sqrt(growth))
  } else {
    lower <- .am_state_components(-.95,h,mode,n_generations)$rho_G
    upper <- .am_state_components(.95,h,mode,n_generations)$rho_G
    rho <- max(lower+1e-8,min(upper-1e-8,rho))
    f <- function(r) .am_state_components(r,h,mode,n_generations)$rho_G-rho
    total_rho <- tryCatch(stats::uniroot(f,c(-.95,.95),tol=1e-10)$root,
      error=function(e) max(-.7,min(.7,rho/(1-h))))
    state <- .am_state_components(total_rho,h,mode,n_generations)
    rho <- state$rho_G
    eta <- .5*(1+rho)*sqrt(state$vG/state$vGc)
    growth <- state$vGc/state$vG
  }
  starting_y <- dat$Y
  if (is.factor(starting_y) || is.logical(starting_y)) starting_y <- as.numeric(starting_y)
  outcome_cov <- stats::cor(dat[,score_names,drop=FALSE],starting_y,use="pairwise.complete.obs")
  ratio <- sqrt(sum(outcome_cov[1:2]^2)/max(sum(outcome_cov[4:5]^2),1e-12))
  if (!is.finite(ratio)) ratio <- 1
  ratio <- max(product/.95,min(.95/product,ratio))
  rel_parent <- c(product*ratio,product/ratio)
  starts <- c(rhoMF=rho, eta_m=eta, eta_f=eta)
  for (j in 1:2) for (role in c("m", "f", "c")) {
    long_role <- c(m="mother", f="father", c="child")[[role]]
    rel <- if (role == "c") growth*rel_parent[j]/(1-rel_parent[j]+growth*rel_parent[j]) else rel_parent[j]
    sd <- sd_score[[paste0("PGS", j, "_", long_role)]]
    starts[[paste0("l", role, j)]] <- sqrt(rel)*sd
    starts[[paste0("th", j, "_", role)]] <- (1-rel)*sd^2
  }
  for (j in 1:2) for (parent in c("m","f")) {
    starts[[paste0("e",j,"_",parent,"c")]] <- .5*sqrt(
      starts[[paste0("th",j,"_",parent)]]*starts[[paste0("th",j,"_c")]])
  }
  for (label in names(starts)) {
    if (is.finite(starts[[label]])) {
      model <- gsub(paste0("\\b", label, "\\*(?=[A-Z])"),
        paste0("start(", .format_lavaan_numeric(starts[[label]]), ")*", label, "*"), model, perl = TRUE)
    }
  }
  model
}

#' Validate the opt-in constraint mode without changing legacy defaults
#' @noRd
.ravel_prepare_constraints <- function(mode, h, am_correction, n_generations,
                                      correction_missing, crossperson,
                                      fixed_cov, fixed_corr) {
  if (mode == "none") return(list(mode=mode, am_correction=am_correction,
    crossperson=crossperson, fixed_cov=fixed_cov, fixed_corr=fixed_corr))
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h < 0 || h >= 1)
    stop("h must be a single finite number in [0, 1).")
  if (!correction_missing && am_correction != mode)
    stop("fit_constraints must agree with am_correction; use '", mode, "' for both.")
  if (mode == "n_gen" && (!is.numeric(n_generations) || length(n_generations) != 1L ||
      !is.finite(n_generations) || n_generations < 1 ||
      n_generations > .Machine$integer.max ||
      n_generations != as.integer(n_generations)))
    stop("am_n_generations must be a single integer >= 1 for constrained n_gen fitting.")
  if (crossperson != "free")
    stop("Constrained fitting requires crossperson = 'free'; its residual correlations are fixed internally.")
  intended_cov <- c(e1_mf=0, e2_mf=0)
  intended_corr <- c(e1_mc=.5, e1_fc=.5, e2_mc=.5, e2_fc=.5)
  cov_used <- fixed_cov[!is.na(fixed_cov)]
  corr_used <- fixed_corr[!is.na(fixed_corr)]
  # Allow matching explicit settings, but never silently override a conflicting
  # user constraint or equate residual covariance with residual correlation.
  if (any(!names(cov_used) %in% names(intended_cov)) ||
      any(cov_used != intended_cov[names(cov_used)]) ||
      any(!names(corr_used) %in% c(names(intended_corr), names(intended_cov))) ||
      any(corr_used != c(intended_corr,intended_cov)[names(corr_used)]))
    stop("fixed_crossperson/fixed_crossperson_corr conflict with constrained fitting (spouse 0; parent-child correlation 0.5).")
  list(mode=mode, am_correction=mode, crossperson="free",
    fixed_cov=.normalize_fixed_crossperson("free", intended_cov),
    fixed_corr=.normalize_fixed_crossperson("free", intended_corr))
}

#' Total derivative of the constrained visible fit with respect to external pi
#' @noRd
.ravel_constraint_core_h_sensitivity <- function(object) {
  labels <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f")
  spec <- object$constraint_spec
  if (is.null(spec) || !isTRUE(spec$depends_on_pi))
    return(setNames(rep(0, 6), labels))
  args <- object$constraint_refit
  if (is.null(args)) stop("This finite-generation fit needs retained refit data to propagate pi uncertainty.")
  cache <- new.env(parent=emptyenv())
  get_core <- function(h) {
    key <- format(h,digits=17,scientific=TRUE)
    if (exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache))
    args$h <- h
    fit <- suppressWarnings(do.call(ravel_fit, args))
    if (!isTRUE(fit$converged)) stop("A finite-generation pi-sensitivity refit did not converge.")
    out <- c(fit$beta_vis, fit$latent_cov_vis)
    out <- out[labels]
    if (any(!is.finite(out))) stop("A finite-generation pi-sensitivity refit was not finite.")
    assign(key,out,envir=cache)
    out
  }
  h <- spec$h
  steps <- unique(pmin(c(.001,.002,.004),(1-h)/4))
  if (steps[1] < 100*.Machine$double.eps)
    stop("h is too close to 1 for a stable pi-sensitivity refit.")
  # Small finite differences can amplify an otherwise negligible optimizer
  # error. Retry at larger local steps under the SAME accuracy criterion;
  # cached fits mean each larger step needs only two additional refits.
  for (step in steps) {
    if (h >= step) {
      coarse <- (get_core(h+step)-get_core(h-step))/(2*step)
      fine <- (get_core(h+step/2)-get_core(h-step/2))/step
    } else {
      center <- c(object$beta_vis,object$latent_cov_vis)[labels]
      coarse <- (get_core(h+step)-center)/step
      fine <- (get_core(h+step/2)-center)/(step/2)
    }
    if (all(abs(coarse-fine) <= .001+.02*abs(fine)))
      return(setNames(as.numeric(fine),labels))
  }
  stop("Finite-generation pi-sensitivity refits are numerically unstable; supply an independently validated core_pi_sensitivity.")
}

#' Check that correction assumptions match those used during constrained fitting
#' @noRd
.ravel_check_constraint_correction <- function(spec, h, mode, n_generations) {
  if (is.null(spec) || is.null(spec$mode) || spec$mode == "none") return(invisible(TRUE))
  if (mode != spec$mode || (mode == "n_gen" && n_generations != spec$n_generations))
    stop("Correction regime/generation must match the constrained fitting assumptions.")
  if (isTRUE(spec$depends_on_pi) && !isTRUE(all.equal(as.numeric(h), as.numeric(spec$h), tolerance=1e-12)))
    stop("Finite-generation constrained fitting depends on h/pi; refit the model at the requested value before correction.")
  invisible(TRUE)
}

#' Check local identification and regular uncertainty in the estimated-state model
#' @noRd
.ravel_estimated_state_fit_valid <- function(fit, rank_tolerance = 1e-7) {
  # Convergence alone can conceal the unidentified reliability direction at
  # random mating with no outcome/covariate signal. Check all free parameters,
  # not only the six parameters retained in a correction packet.
  tryCatch({
    if (is.null(fit) || !isTRUE(lavaan::lavInspect(fit, "converged")) ||
        !suppressWarnings(isTRUE(lavaan::lavInspect(fit, "post.check")))) return(FALSE)
    V <- as.matrix(lavaan::vcov(fit))
    if (!nrow(V) || nrow(V) != ncol(V) || any(!is.finite(V))) return(FALSE)
    V <- (V + t(V)) / 2
    variance_scale <- max(abs(diag(V)))
    if (variance_scale <= 0 || any(diag(V) < -1e-8 * variance_scale)) return(FALSE)
    # Work in uncertainty-scaled coordinates to avoid treating observed score
    # units as evidence of nonidentification.
    scales <- sqrt(pmax(abs(diag(V)), variance_scale * .Machine$double.eps))
    covariance <- V / outer(scales, scales)
    ev <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-7 * max(1, max(abs(ev)))) return(FALSE)
    delta <- lavaan::lavInspect(fit, "delta")
    if (is.list(delta)) delta <- do.call(rbind, delta)
    delta <- as.matrix(delta)
    if (ncol(delta) != ncol(V) || any(!is.finite(delta))) return(FALSE)
    # lavaan stores its FINAL constraint Jacobian here, including nonlinear
    # constraints. Its ceq.JAC slot alone need not contain the final derivative.
    # Fail closed if a future lavaan version changes this representation.
    constraint <- fit@Model@con.jac
    equality <- attr(constraint, "ceq.idx")
    inequality <- attr(constraint, "cin.idx")
    inactive <- attr(constraint, "inactive.idx")
    if (!is.matrix(constraint) || ncol(constraint) != ncol(V) ||
        is.null(equality) || !length(equality) || any(!is.finite(constraint))) return(FALSE)
    # Active inequality bounds are nonregular; do not certify ordinary Wald SEs.
    if (length(setdiff(inequality, inactive))) return(FALSE)
    normalize_rows <- function(A) {
      norms <- sqrt(rowSums(A^2))
      keep <- is.finite(norms) & norms > 0
      A[keep, , drop = FALSE] / norms[keep]
    }
    C <- normalize_rows(sweep(constraint[equality, , drop = FALSE], 2, scales, "*"))
    D <- normalize_rows(sweep(delta, 2, scales, "*"))
    decomposition <- svd(C, nu = 0, nv = ncol(C))
    constraint_rank <- sum(decomposition$d > rank_tolerance * max(decomposition$d))
    dimension <- ncol(C) - constraint_rank
    if (dimension <= 0) return(FALSE)
    tangent <- decomposition$v[, seq.int(constraint_rank + 1L, ncol(C)), drop = FALSE]
    singular <- svd(D %*% tangent, nu = 0, nv = 0)$d
    if (length(singular) < dimension || !length(singular) || max(singular) <= 0) return(FALSE)
    sum(singular > rank_tolerance * max(singular)) == dimension
  }, error = function(e) FALSE)
}

#' Verify actual constraint residuals, not only the optimizer convergence flag
#' @noRd
.ravel_constraint_fit_valid <- function(fit, tolerance=1e-5) {
  if (is.null(fit) || !isTRUE(lavaan::lavInspect(fit,"converged"))) return(FALSE)
  if (!suppressWarnings(isTRUE(lavaan::lavInspect(fit,"post.check")))) return(FALSE)
  pe <- lavaan::parameterEstimates(fit)
  if ("am_inverse_vG_parent" %in% pe$label && !.ravel_estimated_state_fit_valid(fit))
    return(FALSE)
  values <- new.env(parent=baseenv())
  for (i in which(nzchar(pe$label))) assign(pe$label[i],pe$est[i],envir=values)
  pt <- lavaan::parameterTable(fit)
  rows <- which(pt$op %in% c("==",">","<"))
  checks <- vapply(rows,function(i) {
    lhs <- tryCatch(eval(parse(text=pt$lhs[i]),envir=values),error=function(e) NA_real_)
    rhs <- tryCatch(eval(parse(text=pt$rhs[i]),envir=values),error=function(e) NA_real_)
    if (length(lhs)!=1L || length(rhs)!=1L || !all(is.finite(c(lhs,rhs)))) return(FALSE)
    tol <- tolerance*(1+max(abs(c(lhs,rhs))))
    bound_tol <- 1e-10*(1+max(abs(c(lhs,rhs))))
    switch(pt$op[i],"=="=abs(lhs-rhs)<=tol,">"=lhs>=rhs-bound_tol,"<"=lhs<=rhs+bound_tol)
  },logical(1))
  all(checks)
}
