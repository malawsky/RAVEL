# =============================================================================
# ravel_sem.R
# =============================================================================
# RAVEL Structural Equation Model for Trio Genetic Analysis
#
# Author: Daniel Malawsky (d.malawsky@gmail.com)
# =============================================================================

# =============================================================================
# MODEL BUILDING FUNCTIONS (internal)
# =============================================================================

#' Build measurement model syntax (free loadings only)
#' @noRd
.build_measurement <- function() {
  "
  # Measurement model: PGS as indicators of a latent (visible) genetic liability
  # Free loadings allow role-specific reliability
  Gm =~ lm1*PGS1_mother + lm2*PGS2_mother
  Gf =~ lf1*PGS1_father + lf2*PGS2_father
  Gc =~ lc1*PGS1_child  + lc2*PGS2_child

  # Label indicator residual variances so cross-person residual correlations
  # can be constrained on the correlation scale when requested.
  PGS1_mother ~~ th1_m*PGS1_mother
  PGS1_father ~~ th1_f*PGS1_father
  PGS1_child  ~~ th1_c*PGS1_child
  PGS2_mother ~~ th2_m*PGS2_mother
  PGS2_father ~~ th2_f*PGS2_father
  PGS2_child  ~~ th2_c*PGS2_child
  "
}

#' Build latent genetic covariance structure
#' @noRd
.build_genetic_structure <- function() {
  "
  # Latent genetic covariances (visible scale)
  # rhoMF: assortative mating covariance between parental latent genetics
  # eta_m, eta_f: parent-child genetic covariances (Mendelian + assortment)
  Gm ~~ rhoMF*Gf
  Gc ~~ eta_m*Gm
  Gc ~~ eta_f*Gf
  "
}

#' Format a numeric value for direct use in lavaan syntax
#' @noRd
.format_lavaan_numeric <- function(x) {
  format(signif(x, 15), trim = TRUE, scientific = FALSE)
}

#' Normalize optional fixed cross-person residual specifications
#' @noRd
.normalize_fixed_crossperson <- function(crossperson,
                                         values = NULL,
                                         arg_name = "fixed_crossperson",
                                         sib_roles = character(0)) {
  allowed <- if (crossperson == "free") {
    c("e1_mf", "e1_mc", "e1_fc", "e2_mf", "e2_mc", "e2_fc")
  } else {
    c("e1_mf", "e1_pc", "e2_mf", "e2_pc")
  }

  allowed <- c(allowed, names(.parental_sib_residual_benchmarks(sib_roles)))
  out <- setNames(rep(NA_real_, length(allowed)), allowed)

  if (is.null(values)) {
    return(out)
  }

  if (!is.numeric(values) || is.complex(values) || is.null(names(values)) ||
      anyNA(names(values)) || any(!nzchar(names(values))) || anyDuplicated(names(values))) {
    stop(arg_name, " must be a numeric vector with unique parameter names.")
  }

  bad_names <- setdiff(names(values), allowed)
  if (length(bad_names) > 0) {
    stop(
      "Invalid ", arg_name, " names for crossperson = '", crossperson, "': ",
      paste(bad_names, collapse = ", ")
    )
  }

  if (any(!is.na(values) & !is.finite(values))) {
    stop(arg_name, " values must be finite numerics.")
  }

  if (arg_name == "fixed_crossperson_corr" && any(abs(values) > 1, na.rm = TRUE))
    stop("fixed_crossperson_corr values must lie in [-1, 1].")
  out[names(values)] <- unname(values)
  out
}

#' Build cross-person residual covariances
#'
#' @param crossperson Character. Cross-person residual covariance constraints:
#'   - 'free': separate e_mc and e_fc parameters (default; most flexible)
#'   - 'equal_pc': constrain mother-child and father-child residual covariances
#'                to be equal within each PGS version (e_mc == e_fc)
#' @noRd
.build_crossperson <- function(crossperson = c("free", "equal_pc"),
                               fixed_crossperson = NULL,
                               fixed_crossperson_corr = NULL) {
  crossperson <- match.arg(crossperson)
  fixed_crossperson <- .normalize_fixed_crossperson(
    crossperson,
    fixed_crossperson,
    arg_name = "fixed_crossperson"
  )
  fixed_crossperson_corr <- .normalize_fixed_crossperson(
    crossperson,
    fixed_crossperson_corr,
    arg_name = "fixed_crossperson_corr"
  )

  if (crossperson != "free" && any(!is.na(fixed_crossperson_corr))) {
    stop("fixed_crossperson_corr is only supported when crossperson = 'free'.")
  }

  overlap <- names(fixed_crossperson)[
    !is.na(fixed_crossperson) & !is.na(fixed_crossperson_corr)
  ]
  if (length(overlap) > 0) {
    stop(
      "The same cross-person residual cannot be fixed on both the covariance and correlation scales: ",
      paste(overlap, collapse = ", ")
    )
  }

  cross_line <- function(lhs, rhs, label) {
    multiplier <- if (!is.na(fixed_crossperson[[label]])) {
      .format_lavaan_numeric(fixed_crossperson[[label]])
    } else {
      label
    }
    sprintf("    %s ~~ %s*%s", lhs, multiplier, rhs)
  }

  corr_constraints <- character(0)
  if (crossperson == "free") {
    theta_map <- list(
      e1_mf = c("th1_m", "th1_f"),
      e1_mc = c("th1_m", "th1_c"),
      e1_fc = c("th1_f", "th1_c"),
      e2_mf = c("th2_m", "th2_f"),
      e2_mc = c("th2_m", "th2_c"),
      e2_fc = c("th2_f", "th2_c")
    )

    corr_labels <- names(fixed_crossperson_corr)[!is.na(fixed_crossperson_corr)]
    corr_constraints <- vapply(corr_labels, function(label) {
      theta_labels <- theta_map[[label]]
      sprintf(
        "    %s == %s*sqrt(%s*%s)",
        label,
        .format_lavaan_numeric(fixed_crossperson_corr[[label]]),
        theta_labels[1],
        theta_labels[2]
      )
    }, character(1))
  }

  if (crossperson == "free") {
    out <- paste(
      "    # Cross-person residual covariances",
      "    # Same GWAS weights create correlated measurement errors across relatives",
      cross_line("PGS1_mother", "PGS1_father", "e1_mf"),
      cross_line("PGS1_mother", "PGS1_child",  "e1_mc"),
      cross_line("PGS1_father", "PGS1_child",  "e1_fc"),
      cross_line("PGS2_mother", "PGS2_father", "e2_mf"),
      cross_line("PGS2_mother", "PGS2_child",  "e2_mc"),
      cross_line("PGS2_father", "PGS2_child",  "e2_fc"),
      sep = "\n"
    )
  } else {
    out <- paste(
      "    # Cross-person residual covariances (with equality constraints)",
      "    # Constrain parent-child residual covariances to be equal:",
      "    #   Cov(e_m, e_c) = Cov(e_f, e_c) within each PGS version",
      cross_line("PGS1_mother", "PGS1_father", "e1_mf"),
      cross_line("PGS1_mother", "PGS1_child",  "e1_pc"),
      cross_line("PGS1_father", "PGS1_child",  "e1_pc"),
      cross_line("PGS2_mother", "PGS2_father", "e2_mf"),
      cross_line("PGS2_mother", "PGS2_child",  "e2_pc"),
      cross_line("PGS2_father", "PGS2_child",  "e2_pc"),
      sep = "\n"
    )
  }

  if (length(corr_constraints) > 0) {
    out <- paste(
      out,
      "",
      "    # Fixed cross-person residual correlations",
      paste(corr_constraints, collapse = "\n"),
      sep = "\n"
    )
  }

  out
}

#' Parental-sibling roles included in the model
#'
#' Role codes: "ms" is a full sibling of the mother, "fs" a full sibling of the
#' father. Relative to the child these are avuncular relatives.
#' @noRd
.parental_sib_roles <- function(parental_sibs = c("none", "maternal", "paternal", "both")) {
  parental_sibs <- match.arg(parental_sibs)
  switch(parental_sibs, none = character(0), maternal = "ms", paternal = "fs",
    both = c("ms", "fs"))
}

#' Internal names of the parental-sibling score columns
#' @noRd
.parental_sib_columns <- function(roles) {
  cols <- list(ms = c("PGS1_mat_sib", "PGS2_mat_sib"), fs = c("PGS1_pat_sib", "PGS2_pat_sib"))
  unlist(cols[roles], use.names = FALSE)
}

#' Pedigree residual correlations for included parental siblings
#' @noRd
.parental_sib_residual_benchmarks <- function(roles) {
  out <- numeric(0)
  own <- c(ms = "m", fs = "f")
  for (r in roles) for (j in 1:2) for (p in c(own[[r]], setdiff(c("m", "f"), own[[r]]), "c"))
    out[[sprintf("e%d_%s_%s", j, p, r)]] <- if (p == own[[r]]) .5 else if (p == "c") .25 else 0
  if (length(roles) == 2L) for (j in 1:2) out[[paste0("e", j, "_ms_fs")]] <- 0
  out
}

#' Build the parental-sibling extension of the SEM
#'
#' @details
#' Each parental sibling has the same one-factor measurement model as the trio
#' members. Its visible liability correlates freely with every other family
#' member: \code{kap_*} with the sibling's own parent, \code{chi_*} with that
#' parent's spouse, \code{eta_*} with the child, and \code{zeta} between the two
#' parental siblings. Within each score version, residual covariances with all
#' other family members are free, as for the trio. The outcome regression gains
#' one coefficient per sibling, and the parent-specific contrasts
#' \code{delta_m := bm - bms} and \code{delta_f := bf - bfs} are defined.
#' @noRd
.build_parental_sibs <- function(parental_sibs = "none", constrained = FALSE,
                                fixed_crossperson = NULL, fixed_crossperson_corr = NULL) {
  roles <- .parental_sib_roles(parental_sibs)
  if (!length(roles)) return("")
  # Constrained modes fix within-score residual correlations to the relatedness
  # benchmarks: 1/2 with the sibling's own parent, 1/4 with the child, and 0
  # with the parent's spouse and between the two parental siblings.
  benchmarks <- .parental_sib_residual_benchmarks(roles)
  sibling_residual <- function(lhs, rhs, label, variance_labels) {
    covariance <- fixed_crossperson[label]
    correlation <- if (constrained) benchmarks[[label]] else fixed_crossperson_corr[label]
    has_cov <- length(covariance) && !is.na(covariance)
    has_cor <- length(correlation) && !is.na(correlation)
    if (has_cov && has_cor && !constrained)
      stop("The same cross-person residual cannot be fixed on both the covariance and correlation scales: ", label)
    # Label every nonzero fixed covariance as well, so extraction and all
    # downstream covariance calculations retain the specified value.
    zero <- (has_cov && covariance == 0) || (has_cor && correlation == 0)
    out <- sprintf("%s ~~ %s*%s", lhs, if (zero) "0" else label, rhs)
    if (!zero && has_cov)
      out <- c(out, sprintf("%s == %s", label, .format_lavaan_numeric(covariance)))
    if (!zero && has_cor)
      out <- c(out, sprintf("%s == %s*sqrt(%s*%s)", label,
        .format_lavaan_numeric(correlation), variance_labels[1], variance_labels[2]))
    out
  }

  spec <- list(
    ms = list(lv = "Gms", cols = c("PGS1_mat_sib", "PGS2_mat_sib"), own = "m",
      spouse = "f", label = "maternal sibling (full sibling of the mother)"),
    fs = list(lv = "Gfs", cols = c("PGS1_pat_sib", "PGS2_pat_sib"), own = "f",
      spouse = "m", label = "paternal sibling (full sibling of the father)")
  )
  trio_lv <- c(m = "Gm", f = "Gf", c = "Gc")
  trio_col <- c(m = "mother", f = "father", c = "child")

  lines <- character(0)
  for (r in roles) {
    s <- spec[[r]]
    lines <- c(lines,
      "",
      sprintf("# Parental sibling: %s", s$label),
      sprintf("%s =~ l%s1*%s + l%s2*%s", s$lv, r, s$cols[1], r, s$cols[2]),
      sprintf("%s ~~ th1_%s*%s", s$cols[1], r, s$cols[1]),
      sprintf("%s ~~ th2_%s*%s", s$cols[2], r, s$cols[2]),
      sprintf("%s ~~ kap_%s*%s", s$lv, s$own, trio_lv[[s$own]]),
      sprintf("%s ~~ chi_%s*%s", s$lv, r, trio_lv[[s$spouse]]),
      sprintf("%s ~~ eta_%s*Gc", s$lv, r))
    for (j in 1:2) for (p in c(s$own, s$spouse, "c")) {
      label <- sprintf("e%d_%s_%s", j, p, r)
      lines <- c(lines, sibling_residual(s$cols[j], sprintf("PGS%d_%s", j, trio_col[[p]]),
        label, c(sprintf("th%d_%s", j, p), sprintf("th%d_%s", j, r))))
    }
  }
  if (length(roles) == 2L) {
    lines <- c(lines, "", "# Covariances between the two parental siblings", "Gms ~~ zeta*Gfs")
    for (j in 1:2) lines <- c(lines, sibling_residual(paste0("PGS", j, "_mat_sib"),
      paste0("PGS", j, "_pat_sib"), paste0("e", j, "_ms_fs"),
      c(paste0("th", j, "_ms"), paste0("th", j, "_fs"))))
  }

  lines <- c(lines, "", "# Parental-sibling (lineage) coefficients and parent-specific contrasts")
  if ("ms" %in% roles) lines <- c(lines, "Y ~ bms*Gms", "delta_m := bm - bms")
  if ("fs" %in% roles) lines <- c(lines, "Y ~ bfs*Gfs", "delta_f := bf - bfs")
  paste0("    ", lines, collapse = "\n")
}

#' Build outcome regression
#' @noRd
.build_outcome <- function(fix_parent_paths = FALSE, fixed_effects = NULL,
                           latents = c("Gc", "Gm", "Gf")) {

  if (fix_parent_paths) {
    y_syntax <- "
    # Outcome regression (parental paths fixed to 0)
    Y ~ bD*Gc + 0*Gm + 0*Gf
    "
  } else {
    y_syntax <- "
    # Outcome regression
    # bD: direct genetic effect (child genetics -> own phenotype)
    # bm: maternal visible-liability coefficient
    # bf: paternal visible-liability coefficient (legacy father-role name)
    Y ~ bD*Gc + bm*Gm + bf*Gf
    "
  }

  # Add fixed effects if specified
  if (!is.null(fixed_effects) && length(fixed_effects) > 0) {
    fe_lines <- paste(sprintf("Y ~ c_%d*%s", seq_along(fixed_effects), fixed_effects),
                      collapse = "\n    ")
    y_syntax <- paste0(y_syntax, "\n    # Fixed effects (covariates)\n    ", fe_lines)
    for (g in latents) {
      y_syntax <- paste(y_syntax, paste(g, "~~", paste(fixed_effects, collapse = " + ")),
                        sep = "\n    ")
    }
    if (length(fixed_effects) > 1L) {
      pairs <- utils::combn(fixed_effects, 2L)
      y_syntax <- paste(y_syntax,
                        paste(pairs[1, ], "~~", pairs[2, ], collapse = "\n    "),
                        sep = "\n    ")
    }
  }

  y_syntax
}

#' Build complete RAVEL SEM model
#' @noRd
.build_model <- function(fix_parent_paths = FALSE,
                         fixed_effects = NULL,
                         crossperson = "free",
                         fixed_crossperson = NULL,
                         fixed_crossperson_corr = NULL,
                         parental_sibs = "none",
                         constrained = FALSE) {
  roles <- .parental_sib_roles(parental_sibs)
  sib_latents <- c(ms = "Gms", fs = "Gfs")[roles]
  fixed_crossperson <- .normalize_fixed_crossperson(crossperson, fixed_crossperson,
    sib_roles = roles)
  fixed_crossperson_corr <- .normalize_fixed_crossperson(crossperson, fixed_crossperson_corr,
    arg_name = "fixed_crossperson_corr", sib_roles = roles)
  if (crossperson != "free" && any(!is.na(fixed_crossperson_corr)))
    stop("fixed_crossperson_corr is only supported when crossperson = 'free'.")
  overlap <- names(fixed_crossperson)[!is.na(fixed_crossperson) & !is.na(fixed_crossperson_corr)]
  if (length(overlap)) stop("The same cross-person residual cannot be fixed on both the covariance and correlation scales: ",
    paste(overlap, collapse = ", "))
  trio_labels <- names(.normalize_fixed_crossperson(crossperson))
  trio_only <- function(x) if (is.null(x)) NULL else x[intersect(names(x), trio_labels)]
  paste0(
    "# =============================================================================",
    "\n# RAVEL Structural Equation Model",
    "\n# =============================================================================",
    "\n", .build_measurement(),
    "\n", .build_genetic_structure(),
    "\n", .build_crossperson(
      crossperson,
      fixed_crossperson = trio_only(fixed_crossperson),
      fixed_crossperson_corr = trio_only(fixed_crossperson_corr)
    ),
    "\n", .build_outcome(fix_parent_paths, fixed_effects,
      latents = c("Gc", "Gm", "Gf", unname(sib_latents))),
    if (length(sib_latents)) paste0("\n", .build_parental_sibs(parental_sibs, constrained,
      fixed_crossperson, fixed_crossperson_corr)) else ""
  )
}

#' Ensure lavaan can initialize on systems with invalid detectCores()
#' @noRd
.ensure_lavaan_runtime <- function(verbose = FALSE) {
  detected <- suppressWarnings(parallel::detectCores())

  if (is.numeric(detected) && length(detected) == 1L &&
      !is.na(detected) && detected >= 1L) {
    return(invisible(TRUE))
  }

  ns_parallel <- asNamespace("parallel")
  was_locked <- bindingIsLocked("detectCores", ns_parallel)

  if (was_locked) {
    unlockBinding("detectCores", ns_parallel)
  }

  assign(
    "detectCores",
    function(all.tests = FALSE, logical = TRUE) 1L,
    envir = ns_parallel
  )

  if (was_locked) {
    lockBinding("detectCores", ns_parallel)
  }

  if (isTRUE(verbose)) {
    message("parallel::detectCores() returned an invalid value; using 1 core for lavaan.")
  }

  invisible(TRUE)
}

# =============================================================================
# POST-ESTIMATION HELPERS (internal)
# =============================================================================

#' Safe linear solve (avoids hard failure on singular matrices)
#' @noRd
.safe_solve <- function(A, B = NULL) {
  out <- tryCatch({
    if (is.null(B)) solve(A) else solve(A, B)
  }, error = function(e) {
    tryCatch({
      # Fall back to Moore-Penrose generalized inverse if available
      MASS::ginv(A) %*% B
    }, error = function(e2) {
      matrix(NA_real_, nrow = nrow(A), ncol = if (is.null(B)) ncol(A) else ncol(as.matrix(B)))
    })
  })
  out
}

#' @noRd
.diag3 <- function(v) {
  M <- matrix(0, nrow = 3L, ncol = 3L)
  diag(M) <- as.numeric(v)
  M
}

#' Component-variance state under n generations (or equilibrium) of assortment
#'
#' @description
#' Tracks the raw variances of the visible and hidden components relative to
#' their random-mating values, alongside the standardized within-person
#' visible-hidden covariance d. Generalises the fixed-unit component variance
#' recursion, which is the special case v_G = v_H = 1.
#'
#' @details
#' The variance recursion uses raw Mendelian transmission (one half) and a
#' segregation variance equal to half the random-mating component variance.
#' The fitted transmission coefficients are used for covariance blocks, scaled
#' to the hidden side by `t_ratio`.
#'
#' @noRd
.am_state_components <- function(rho_T, h, am_correction, n_generations,
                                 tol = 1e-13, max_iter = 10000L) {

  a0 <- sqrt(1 - h)
  b0 <- sqrt(h)

  evaluate <- function(vG, vH, d) {
    a <- a0 * sqrt(vG)
    b <- b0 * sqrt(vH)
    varT <- a * a + b * b + 2 * a * b * d
    if (!is.finite(varT) || varT <= 0) return(NULL)
    cGT <- a + b * d
    cHT <- b + a * d
    list(
      a = a, b = b, varT = varT,
      rho_G = rho_T * cGT * cGT / varT,
      q     = rho_T * cGT * cHT / varT,
      rho_H = rho_T * cHT * cHT / varT
    )
  }

  advance <- function(vG, vH, d, s) {
    vGn <- 0.5 * vG * (1 + s$rho_G) + 0.5
    vHn <- 0.5 * vH * (1 + s$rho_H) + 0.5
    if (!is.finite(vGn) || !is.finite(vHn) || vGn <= 0 || vHn <= 0) return(NULL)
    dn <- 0.5 * sqrt(vG * vH) * (d + s$q) / sqrt(vGn * vHn)
    list(vG = vGn, vH = vHn, d = dn)
  }

  vG <- 1
  vH <- 1
  d  <- 0

  if (identical(am_correction, "equilibrium")) {
    converged <- FALSE
    for (i in seq_len(max_iter)) {
      s <- evaluate(vG, vH, d)
      if (is.null(s)) return(NULL)
      nx <- advance(vG, vH, d, s)
      if (is.null(nx)) return(NULL)
      delta <- max(abs(nx$vG - vG), abs(nx$vH - vH), abs(nx$d - d))
      vG <- nx$vG
      vH <- nx$vH
      d <- nx$d
      if (delta < tol) {
        converged <- TRUE
        break
      }
    }
    if (!converged) return(NULL)
  } else {
    if (n_generations > 1L) {
      for (g in seq_len(n_generations - 1L)) {
        s <- evaluate(vG, vH, d)
        if (is.null(s)) return(NULL)
        nx <- advance(vG, vH, d, s)
        if (is.null(nx)) return(NULL)
        vG <- nx$vG
        vH <- nx$vH
        d <- nx$d
      }
    }
  }

  s <- evaluate(vG, vH, d)
  if (is.null(s)) return(NULL)
  ch <- advance(vG, vH, d, s)
  if (is.null(ch)) return(NULL)

  a_par <- a0 * sqrt(vG)
  b_par <- b0 * sqrt(vH)

  list(
    vG = vG, vH = vH, d = d,
    vGc = ch$vG, vHc = ch$vH,
    q = s$q, rho_G = s$rho_G, rho_H = s$rho_H,
    t_ratio = sqrt((vH / ch$vH) / (vG / ch$vG)),
    a_par = a_par, b_par = b_par,
    a_ch  = a0 * sqrt(ch$vG),
    b_ch  = b0 * sqrt(ch$vH),
    pi_present = b_par^2 / (a_par^2 + b_par^2)
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Total-liability variance recursion independent of the visible/hidden split
#' @noRd
.total_liability_state <- function(rho_T_corr, tm, tf, am_correction, am_n_generations) {
  S <- tm^2 + tf^2
  u <- tm * tf
  persistence <- S + 2 * u * rho_T_corr
  segregation_var <- 1 - S

  if (am_correction == "n_gen") {
    var_parent <- 1
    if (am_n_generations > 1L) {
      for (g in seq_len(am_n_generations - 1L)) {
        var_parent <- persistence * var_parent + segregation_var
      }
    }
    var_child <- persistence * var_parent + segregation_var
  } else {
    denom <- 1 - persistence
    if (denom <= 0) {
      return(list(
        var_total = c(child = NA_real_, mother = NA_real_, father = NA_real_),
        cov = c(rhoMF = NA_real_, eta_m = NA_real_, eta_f = NA_real_),
        corr = c(rhoMF = NA_real_, eta_m = NA_real_, eta_f = NA_real_),
        persistence = persistence,
        segregation_var = segregation_var
      ))
    }
    var_parent <- segregation_var / denom
    var_child <- var_parent
  }

  var_total <- c(child = var_child, mother = var_parent, father = var_parent)
  cov_mf <- rho_T_corr * var_parent
  cov_cm <- (tm + tf * rho_T_corr) * var_parent
  cov_cf <- (tf + tm * rho_T_corr) * var_parent

  list(
    var_total = var_total,
    cov = c(rhoMF = cov_mf, eta_m = cov_cm, eta_f = cov_cf),
    corr = .latent_corr_from_cov(
      var_total,
      cov_cm = cov_cm,
      cov_cf = cov_cf,
      cov_mf = cov_mf
    ),
    persistence = persistence,
    segregation_var = segregation_var
  )
}

#' @noRd
.hidden_correction_na_result <- function(beta_vis, Sigma_xx, t_m, t_f, h,
                                         n_generations = NULL) {
  na_beta <- c(bD = NA_real_, bm = NA_real_, bf = NA_real_)
  out <- list(
    beta_h_only_assumed = beta_vis / sqrt(1 - h),
    beta_total_assumed = na_beta,
    rho_total_assumed = NA_real_,
    rho_hidden_assumed = NA_real_,
    eta_total_assumed = c(eta_m = NA_real_, eta_f = NA_real_),
    eta_hidden_assumed = c(eta_m = NA_real_, eta_f = NA_real_),
    var_total_assumed = c(child = NA_real_, mother = NA_real_, father = NA_real_),
    rho_T_corr_assumed = NA_real_,
    kappa_T_assumed = NA_real_,
    pi_present = NA_real_,
    t_ratio = NA_real_,
    tm = t_m, tf = t_f,
    d = NA_real_, d_child = NA_real_, q = NA_real_,
    Sigma_xx = Sigma_xx,
    Sigma_xh = matrix(NA_real_, 3, 3,
                      dimnames = list(c("Gc", "Gm", "Gf"), c("Hc", "Hm", "Hf"))),
    Sigma_hh = matrix(NA_real_, 3, 3,
                      dimnames = list(c("Hc", "Hm", "Hf"), c("Hc", "Hm", "Hf"))),
    M = matrix(NA_real_, 3, 3)
  )
  if (!is.null(n_generations)) out$n_generations <- n_generations
  out
}

#' Trio hidden-variant correction with component-specific AM inflation
#'
#' @details
#' Implements the visible-to-total map
#'   b_vis = [ diag(a) + Sigma_xx^{-1} Sigma_xh diag(b) ] beta
#' with role-specific loadings
#'   a_r = sqrt(1-h) sqrt(v_G,r),   b_r = sqrt(h) sqrt(v_H,r),
#' and role-specific total-liability variances
#'   Var(T_r^*) = a_r^2 + b_r^2 + 2 a_r b_r d_r.
#'
#' `h` is the hidden fraction of additive genetic variance on the founder or
#' random-mating reference scale.
#'
#' @noRd
.hidden_correction_trio_core <- function(beta_vis, cov_cm, cov_cf, cov_mf, h,
                                         am_correction = c("n_gen", "equilibrium", "estimated_state"),
                                         n_generations = 1L) {

  am_correction <- match.arg(am_correction)

  beta_vis <- as.numeric(beta_vis)
  names(beta_vis) <- c("bD", "bm", "bf")
  cov_cm <- unname(cov_cm)
  cov_cf <- unname(cov_cf)
  cov_mf <- unname(cov_mf)

  Sigma_xx <- matrix(
    c(1, cov_cm, cov_cf,
      cov_cm, 1, cov_mf,
      cov_cf, cov_mf, 1),
    nrow = 3, byrow = TRUE
  )
  colnames(Sigma_xx) <- rownames(Sigma_xx) <- c("Gc", "Gm", "Gf")

  den <- 1 - cov_mf^2
  if (abs(den) < 1e-8) {
    t_m <- t_f <- NA_real_
  } else {
    t_m <- (cov_cm - cov_mf * cov_cf) / den
    t_f <- (cov_cf - cov_mf * cov_cm) / den
  }

  n_gen_field <- if (identical(am_correction, "n_gen")) as.integer(n_generations) else NULL

  if (any(!is.finite(c(t_m, t_f)))) {
    warning("Hidden correction: non-finite transmission coefficients; returning NAs.")
    return(.hidden_correction_na_result(beta_vis, Sigma_xx, t_m, t_f, h, n_gen_field))
  }

  if (identical(am_correction, "estimated_state")) {
    st <- .hidden_estimated_state(cov_mf, cov_cm, cov_cf, h)
    if (is.null(st)) {
      warning("Hidden correction: estimated variance state is incompatible with the supplied founder hidden fraction; returning NAs.")
      return(.hidden_correction_na_result(beta_vis, Sigma_xx, t_m, t_f, h, n_gen_field))
    }
    rho_T <- st$rho_T
  } else if (h == 0) {
    rho_T <- cov_mf
  } else if (abs(cov_mf) < 1e-12) {
    rho_T <- 0
  } else if (identical(am_correction, "n_gen") && as.integer(n_generations) == 1L) {
    rho_T <- cov_mf / (1 - h)
    if (!is.finite(rho_T) || abs(rho_T) >= 1) {
      warning("Hidden correction: implied total-liability spouse correlation is ",
              "not in (-1, 1); the visible spouse correlation is too large for ",
              "the assumed hidden fraction. Returning NAs.")
      return(.hidden_correction_na_result(beta_vis, Sigma_xx, t_m, t_f, h, n_gen_field))
    }
  } else {
    f_root <- function(r) {
      st <- .am_state_components(r, h, am_correction, as.integer(n_generations))
      if (is.null(st)) return(NA_real_)
      st$rho_G - cov_mf
    }
    f0 <- f_root(0)
    rho_T <- NA_real_
    if (is.finite(f0) && abs(f0) < 1e-12) {
      rho_T <- 0
    } else if (is.finite(f0)) {
      search <- c(0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.98, 0.99, 0.995, 0.999)
      search <- if (f0 < 0) search else -search
      prev_x <- 0
      prev_y <- f0
      for (x in search) {
        y <- f_root(x)
        if (!is.finite(y)) next
        if (abs(y) < 1e-12) {
          rho_T <- x
          break
        }
        if (prev_y * y < 0) {
          rho_T <- stats::uniroot(f_root, interval = sort(c(prev_x, x)), tol = 1e-12)$root
          break
        }
        prev_x <- x
        prev_y <- y
      }
    }
    if (!is.finite(rho_T)) {
      grid <- seq(-0.999, 0.999, length.out = 201L)
      vals <- vapply(grid, f_root, numeric(1))
      finite <- is.finite(vals)
      exact <- which(finite & abs(vals) < 1e-12)
      sign_change <- which(
        finite[-length(finite)] & finite[-1L] &
          vals[-length(vals)] * vals[-1L] < 0
      )
      if (length(exact) > 0L) {
        rho_T <- grid[exact[1L]]
      } else if (length(sign_change) > 0L) {
        i <- sign_change[1L]
        rho_T <- stats::uniroot(f_root, interval = c(grid[i], grid[i + 1L]), tol = 1e-12)$root
      }
    }
    if (!is.finite(rho_T)) {
      warning("Hidden correction: could not bracket a root for the implied ",
              "total-liability spouse correlation; returning NAs.")
      return(.hidden_correction_na_result(beta_vis, Sigma_xx, t_m, t_f, h, n_gen_field))
    }
  }

  if (!identical(am_correction, "estimated_state"))
    st <- .am_state_components(rho_T, h, am_correction, as.integer(n_generations))
  if (is.null(st)) {
    warning("Hidden correction: component-variance state is not well defined; returning NAs.")
    return(.hidden_correction_na_result(beta_vis, Sigma_xx, t_m, t_f, h, n_gen_field))
  }

  d <- st$d
  q <- st$q

  tH_m <- st$t_ratio * t_m
  tH_f <- st$t_ratio * t_f

  Cov_Gc_Hm <- t_m * d + t_f * q
  Cov_Gc_Hf <- t_f * d + t_m * q
  Cov_Gm_Hc <- tH_m * d + tH_f * q
  Cov_Gf_Hc <- tH_f * d + tH_m * q
  Cov_Gc_Hc <- t_m * tH_m * d + t_m * tH_f * q +
               t_f * tH_m * q + t_f * tH_f * d

  Sigma_xh <- matrix(
    c(Cov_Gc_Hc, Cov_Gc_Hm, Cov_Gc_Hf,
      Cov_Gm_Hc, d,         q,
      Cov_Gf_Hc, q,         d),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("Gc", "Gm", "Gf"), c("Hc", "Hm", "Hf"))
  )

  a_vec <- c(st$a_ch, st$a_par, st$a_par)
  b_vec <- c(st$b_ch, st$b_par, st$b_par)

  K <- .safe_solve(Sigma_xx, Sigma_xh)
  M <- .diag3(a_vec) + K %*% .diag3(b_vec)

  beta_total <- as.numeric(.safe_solve(M, beta_vis))
  names(beta_total) <- names(beta_vis)

  d_vec <- c(Cov_Gc_Hc, d, d)
  var_total <- a_vec^2 + b_vec^2 + 2 * a_vec * b_vec * d_vec
  names(var_total) <- c("child", "mother", "father")

  rho_H   <- st$rho_H
  eta_H_m <- tH_m + tH_f * rho_H
  eta_H_f <- tH_f + tH_m * rho_H

  cov_Tmf <- a_vec[2] * a_vec[3] * cov_mf +
             a_vec[2] * b_vec[3] * q +
             b_vec[2] * a_vec[3] * q +
             b_vec[2] * b_vec[3] * rho_H

  cov_Tcm <- a_vec[1] * a_vec[2] * cov_cm +
             a_vec[1] * b_vec[2] * Cov_Gc_Hm +
             b_vec[1] * a_vec[2] * Cov_Gm_Hc +
             b_vec[1] * b_vec[2] * eta_H_m

  cov_Tcf <- a_vec[1] * a_vec[3] * cov_cf +
             a_vec[1] * b_vec[3] * Cov_Gc_Hf +
             b_vec[1] * a_vec[3] * Cov_Gf_Hc +
             b_vec[1] * b_vec[3] * eta_H_f

  Sigma_hh <- matrix(
    c(1,       eta_H_m, eta_H_f,
      eta_H_m, 1,       rho_H,
      eta_H_f, rho_H,   1),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("Hc", "Hm", "Hf"), c("Hc", "Hm", "Hf"))
  )

  out <- list(
    beta_h_only_assumed = beta_vis / sqrt(1 - h),
    beta_total_assumed = beta_total,
    rho_total_assumed = unname(cov_Tmf),
    rho_hidden_assumed = rho_H,
    eta_total_assumed = c(eta_m = unname(cov_Tcm), eta_f = unname(cov_Tcf)),
    eta_hidden_assumed = c(eta_m = eta_H_m, eta_f = eta_H_f),
    var_total_assumed = var_total,
    rho_T_corr_assumed = rho_T,
    kappa_T_assumed = unname(cov_Tmf),
    pi_present = st$pi_present,
    t_ratio = st$t_ratio,
    tm = t_m, tf = t_f,
    tH_m = tH_m, tH_f = tH_f,
    d = d, d_child = Cov_Gc_Hc, q = q,
    v_G = st$vG, v_H = st$vH, v_Gc = st$vGc, v_Hc = st$vHc,
    Sigma_xx = Sigma_xx,
    Sigma_xh = Sigma_xh,
    Sigma_hh = Sigma_hh,
    M_component = M,
    M = M
  )
  if (!is.null(st$estimated_state)) out$estimated_state <- st$estimated_state
  if (!is.null(n_gen_field)) out$n_generations <- n_gen_field
  out
}

#' @noRd
.hidden_correction_trio_n_gen <- function(beta_vis, cov_cm, cov_cf, cov_mf, h,
                                          n_generations = 1L) {
  if (!is.numeric(n_generations) || length(n_generations) != 1 ||
      is.na(n_generations) || n_generations < 1) {
    stop("n_generations must be a single integer >= 1.")
  }
  .hidden_correction_trio_core(
    beta_vis, cov_cm, cov_cf, cov_mf, h,
    am_correction = "n_gen",
    n_generations = as.integer(n_generations)
  )
}

#' @noRd
.hidden_correction_trio_equilibrium <- function(beta_vis, cov_cm, cov_cf, cov_mf, h) {
  .hidden_correction_trio_core(
    beta_vis, cov_cm, cov_cf, cov_mf, h,
    am_correction = "equilibrium",
    n_generations = NA_integer_
  )
}

#' Extend a trio hidden correction to parental siblings
#'
#' @details
#' Uses the assortative-mating state of the trio correction \code{hc} (founder
#' scale) to build the visible-hidden cross-block for the child, both parents
#' and the requested parental siblings, then applies the same projection map
#'   b_vis = [ diag(a) + Sigma_xx^{-1} Sigma_xh diag(b) ] beta_raw.
#' With parental raw component covariance S, spouse block K and total weights
#' w (assortment only through T*, A3), the raw person-to-person blocks are
#'   parent and own sibling:          S - I/2   (Var of the grandparental midparent)
#'   parent and spouse's sibling:     X = rho_T S w w'(S - I/2) / V_T
#'   maternal and paternal siblings:  (S - I/2)' S^{-1} X
#' so the same-lineage visible-hidden entry equals the within-person d, and the
#' cross-lineage entries are k q and k^2 q with k = 1 - 1/(2 V_T). Child rows
#' use the fitted visible and hidden transmission coefficients, as in the trio.
#' Siblings share the parents' generation, variances and loadings a, b.
#' @noRd
.hidden_correction_sib_extend <- function(hc, beta_vis, Sigma_xx, h, roles) {
  all_roles <- c("c", "m", "f", roles)
  n <- length(all_roles)
  sib_names <- c(ms = "mat_sib", fs = "pat_sib")[roles]
  var_names <- c("child", "mother", "father", unname(sib_names))
  dimnames(Sigma_xx) <- list(paste0("G", all_roles), paste0("G", all_roles))
  na_result <- list(
    beta_total = setNames(rep(NA_real_, n), names(beta_vis)),
    beta_total_std = setNames(rep(NA_real_, n), names(beta_vis)),
    var_total = setNames(rep(NA_real_, n), var_names),
    Sigma_xx = Sigma_xx, Sigma_xh = NULL, Sigma_hh = NULL, Sigma_TT = NULL,
    lineage_factor = NA_real_, M = NULL, invalid_reason = NULL)
  state <- c(hc$d, hc$q, hc$v_G, hc$v_H, hc$v_Gc, hc$v_Hc, hc$tm, hc$tf,
             hc$tH_m, hc$tH_f, hc$rho_T_corr_assumed, hc$var_total_assumed)
  if (length(state) != 14L || any(!is.finite(state)) || is.null(hc$Sigma_xh) ||
      any(!is.finite(hc$Sigma_xh)) || any(!is.finite(hc$Sigma_hh)) ||
      any(!is.finite(Sigma_xx)))
    return(na_result)

  vG <- hc$v_G
  vH <- hc$v_H
  w <- c(sqrt(1 - h), sqrt(h))
  S <- matrix(c(vG, hc$d * sqrt(vG * vH), hc$d * sqrt(vG * vH), vH), 2L, 2L)
  V_T <- drop(t(w) %*% S %*% w)
  rho_T <- hc$rho_T_corr_assumed
  Sw <- drop(S %*% w)
  K <- rho_T * tcrossprod(Sw) / V_T
  A <- S - diag(2L) / 2
  X <- rho_T * tcrossprod(Sw, drop(t(w) %*% A)) / V_T
  Z <- t(A) %*% .safe_solve(S, X)
  sc <- sqrt(c(vG, vH))
  lineage <- c(m = "m", f = "f", ms = "m", fs = "f")
  is_sib <- c(m = FALSE, f = FALSE, ms = TRUE, fs = TRUE)
  block <- function(r, s) {
    B <- if (r == s) S else if (lineage[[r]] == lineage[[s]]) A else
      if (!is_sib[[r]] && !is_sib[[s]]) K else
      if (is_sib[[r]] && is_sib[[s]]) (if (r == "ms") Z else t(Z)) else
      if (!is_sib[[r]]) X else t(X)
    B / tcrossprod(sc)
  }

  GH <- HH <- matrix(0, n, n, dimnames = list(all_roles, all_roles))
  adults <- c("m", "f", roles)
  for (r in adults) for (s in adults) {
    B <- block(r, s)
    GH[r, s] <- B[1L, 2L]
    HH[r, s] <- B[2L, 2L]
  }
  GH[1:3, 1:3] <- hc$Sigma_xh
  HH[1:3, 1:3] <- hc$Sigma_hh
  for (s in roles) {
    GH["c", s] <- hc$tm * GH["m", s] + hc$tf * GH["f", s]
    GH[s, "c"] <- hc$tH_m * GH[s, "m"] + hc$tH_f * GH[s, "f"]
    HH["c", s] <- HH[s, "c"] <- hc$tH_m * HH["m", s] + hc$tH_f * HH["f", s]
  }
  dimnames(GH) <- list(paste0("G", all_roles), paste0("H", all_roles))
  dimnames(HH) <- list(paste0("H", all_roles), paste0("H", all_roles))

  # A valid visible SEM alone does not guarantee that its fitted sibling
  # covariances can coexist with the hidden blocks implied by pi and history.
  # Do not turn an impossible joint covariance into finite corrected effects.
  joint <- rbind(cbind(Sigma_xx, GH), cbind(t(GH), HH))
  visible_eigen <- eigen(Sigma_xx, symmetric = TRUE, only.values = TRUE)$values
  joint_eigen <- eigen(joint, symmetric = TRUE, only.values = TRUE)$values
  if (min(visible_eigen) <= 1e-10 * max(1, max(visible_eigen)) ||
      (h > 0 && min(joint_eigen) < -1e-8 * max(1, max(joint_eigen)))) {
    warning("Hidden correction: parental-sibling moments are incompatible with the supplied pi/assortment state; returning NAs.")
    na_result$invalid_reason <- "incompatible_parental_sibling_covariance"
    return(na_result)
  }

  a_par <- sqrt(1 - h) * sqrt(vG)
  b_par <- sqrt(h) * sqrt(vH)
  a <- c(sqrt(1 - h) * sqrt(hc$v_Gc), rep(a_par, n - 1L))
  b <- c(sqrt(h) * sqrt(hc$v_Hc), rep(b_par, n - 1L))
  M <- diag(a, n) + .safe_solve(Sigma_xx, GH) %*% diag(b, n)
  beta_total <- setNames(as.numeric(.safe_solve(M, beta_vis)), names(beta_vis))
  var_total <- setNames(c(unname(hc$var_total_assumed[c("child", "mother", "father")]),
    rep(unname(hc$var_total_assumed["mother"]), length(roles))), var_names)
  TT <- diag(a, n) %*% Sigma_xx %*% diag(a, n) + diag(a, n) %*% GH %*% diag(b, n) +
    diag(b, n) %*% t(GH) %*% diag(a, n) + diag(b, n) %*% HH %*% diag(b, n)
  dimnames(TT) <- list(var_names, var_names)
  list(beta_total = beta_total,
    beta_total_std = setNames(sqrt(var_total) * beta_total, names(beta_vis)),
    var_total = var_total, Sigma_xx = Sigma_xx, Sigma_xh = GH, Sigma_hh = HH,
    Sigma_TT = TT, lineage_factor = 1 - 1 / (2 * V_T), M = M)
}

#' Visible latent covariance matrix for the trio and any parental siblings
#' @noRd
.sib_latent_matrix <- function(cov_cm, cov_cf, cov_mf, sib) {
  roles <- c("c", "m", "f", sib$roles)
  S <- diag(length(roles))
  dimnames(S) <- list(roles, roles)
  S["c", "m"] <- S["m", "c"] <- cov_cm
  S["c", "f"] <- S["f", "c"] <- cov_cf
  S["m", "f"] <- S["f", "m"] <- cov_mf
  cv <- sib$cov
  if ("ms" %in% sib$roles) {
    S["m", "ms"] <- S["ms", "m"] <- cv[["kap_m"]]
    S["f", "ms"] <- S["ms", "f"] <- cv[["chi_ms"]]
    S["c", "ms"] <- S["ms", "c"] <- cv[["eta_ms"]]
  }
  if ("fs" %in% sib$roles) {
    S["f", "fs"] <- S["fs", "f"] <- cv[["kap_f"]]
    S["m", "fs"] <- S["fs", "m"] <- cv[["chi_fs"]]
    S["c", "fs"] <- S["fs", "c"] <- cv[["eta_fs"]]
  }
  if (all(c("ms", "fs") %in% sib$roles)) S["ms", "fs"] <- S["fs", "ms"] <- cv[["zeta"]]
  S
}

#' Quadratic-form variance decomposition for any set of roles
#' @noRd
.variance_components_general <- function(beta, Sigma, names) {
  beta <- unname(as.numeric(beta))
  n <- length(beta)
  out <- c()
  for (i in seq_len(n)) out[[names[i]]] <- beta[i]^2 * Sigma[i, i]
  for (i in seq_len(n - 1L)) for (j in seq.int(i + 1L, n))
    out[[paste0("cross_", names[i], "_", names[j])]] <- 2 * beta[i] * beta[j] * Sigma[i, j]
  unlist(out)
}

#' Compute genetic variance components given coefficients and latent covariances
#' @noRd
.variance_components_trio <- function(beta, cov_cm, cov_cf, cov_mf) {
  bD <- beta["bD"]; bm <- beta["bm"]; bf <- beta["bf"]
  comps <- c(
    child_direct = bD^2,
    maternal = bm^2,
    paternal = bf^2,
    cross_child_mother = 2 * bD * bm * cov_cm,
    cross_child_father = 2 * bD * bf * cov_cf,
    cross_mother_father = 2 * bm * bf * cov_mf
  )
  comps
}

#' Compute genetic variance components for possibly non-unit-variance trio factors
#' @noRd
.variance_components_trio_general <- function(beta,
                                              var_c, var_m, var_f,
                                              cov_cm, cov_cf, cov_mf) {
  bD <- beta["bD"]; bm <- beta["bm"]; bf <- beta["bf"]
  comps <- c(
    child_direct        = bD^2 * var_c,
    maternal            = bm^2 * var_m,
    paternal            = bf^2 * var_f,
    cross_child_mother  = 2 * bD * bm * cov_cm,
    cross_child_father  = 2 * bD * bf * cov_cf,
    cross_mother_father = 2 * bm * bf * cov_mf
  )
  comps
}

#' Convert raw total-scale trio covariances to standardized correlations
#' given role-specific variances.
#' @noRd
.latent_corr_from_cov <- function(var_total, cov_cm, cov_cf, cov_mf) {
  sd_c <- sqrt(unname(var_total["child"]))
  sd_m <- sqrt(unname(var_total["mother"]))
  sd_f <- sqrt(unname(var_total["father"]))
  c(
    rhoMF = unname(cov_mf) / (sd_m * sd_f),
    eta_m = unname(cov_cm) / (sd_c * sd_m),
    eta_f = unname(cov_cf) / (sd_c * sd_f)
  )
}

#' Standardize raw total-liability coefficients by role-specific total variances
#' @noRd
.standardize_total_betas <- function(beta_total, var_total) {
  c(
    bD = unname(beta_total["bD"]) * sqrt(unname(var_total["child"])),
    bm = unname(beta_total["bm"]) * sqrt(unname(var_total["mother"])),
    bf = unname(beta_total["bf"]) * sqrt(unname(var_total["father"]))
  )
}

#' Build symmetric normal-approximation confidence intervals from estimates and SEs
#' @noRd
.normal_ci <- function(est, se, ci_level) {
  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  ci <- rbind(
    lower = est - z * se,
    upper = est + z * se
  )
  colnames(ci) <- names(est)
  ci
}

#' Summarize bootstrap/posterior-like samples into SEs and empirical CIs
#' @noRd
.summarize_sample_matrix <- function(samples, ci_level) {
  samples <- as.matrix(samples)
  alpha <- 1 - ci_level
  se <- apply(samples, 2, stats::sd, na.rm = TRUE)
  ci <- apply(
    samples,
    2,
    stats::quantile,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE,
    names = FALSE
  )
  if (is.null(dim(ci))) {
    ci <- matrix(ci, ncol = ncol(samples))
  }
  rownames(ci) <- c("lower", "upper")
  colnames(ci) <- colnames(samples)
  list(
    se = setNames(as.numeric(se), colnames(samples)),
    ci = ci
  )
}

#' Numerical Jacobian via central finite differences
#' @noRd
.finite_difference_jacobian <- function(fun, theta, rel_step = 1e-6) {
  theta <- as.numeric(theta)
  names(theta) <- names(theta) %||% seq_along(theta)
  f0 <- fun(theta)
  J <- matrix(
    NA_real_,
    nrow = length(f0),
    ncol = length(theta),
    dimnames = list(names(f0), names(theta))
  )

  for (j in seq_along(theta)) {
    step <- rel_step * max(1, abs(theta[j]))
    theta_hi <- theta
    theta_lo <- theta
    theta_hi[j] <- theta_hi[j] + step
    theta_lo[j] <- theta_lo[j] - step

    f_hi <- fun(theta_hi)
    f_lo <- fun(theta_lo)
    J[, j] <- (f_hi - f_lo) / (2 * step)
  }

  J
}

#' Compute visible and hidden-corrected genetic-scale estimates from core SEM parameters
#' @noRd
.compute_scale_estimates <- function(beta_vis, cov_cm, cov_cf, cov_mf,
                                     h, am_correction, am_n_generations,
                                     sib = NULL) {
  beta_vis <- as.numeric(beta_vis)
  names(beta_vis) <- c("bD", "bm", "bf")
  cov_cm <- unname(cov_cm)
  cov_cf <- unname(cov_cf)
  cov_mf <- unname(cov_mf)

  latent_cov_vis <- c(rhoMF = cov_mf, eta_m = cov_cm, eta_f = cov_cf)

  hc <- switch(
    am_correction,
    "n_gen" = .hidden_correction_trio_n_gen(
      beta_vis, cov_cm, cov_cf, cov_mf, h,
      n_generations = am_n_generations
    ),
    "equilibrium" = .hidden_correction_trio_equilibrium(beta_vis, cov_cm, cov_cf, cov_mf, h),
    "estimated_state" = .hidden_correction_trio_core(beta_vis, cov_cm, cov_cf, cov_mf, h,
      am_correction = "estimated_state")
  )

  beta_h_only_assumed <- hc$beta_h_only_assumed
  beta_total_assumed <- hc$beta_total_assumed

  latent_cov_total_assumed <- c(
    rhoMF = unname(hc$rho_total_assumed),
    eta_m = unname(hc$eta_total_assumed["eta_m"]),
    eta_f = unname(hc$eta_total_assumed["eta_f"])
  )

  var_total_assumed <- if (is.null(hc$var_total_assumed)) {
    c(child = 1, mother = 1, father = 1)
  } else if (length(hc$var_total_assumed) == 1) {
    c(
      child = hc$var_total_assumed,
      mother = hc$var_total_assumed,
      father = hc$var_total_assumed
    )
  } else {
    hc$var_total_assumed
  }

  latent_corr_total_assumed <- .latent_corr_from_cov(
    var_total_assumed,
    cov_cm = latent_cov_total_assumed["eta_m"],
    cov_cf = latent_cov_total_assumed["eta_f"],
    cov_mf = latent_cov_total_assumed["rhoMF"]
  )

  beta_total_assumed_std <- .standardize_total_betas(
    beta_total_assumed,
    var_total_assumed
  )

  # With parental siblings in the model, bm and bf are conditional on the
  # siblings' liabilities, so all coefficients are mapped jointly.
  sib_total <- NULL
  if (!is.null(sib) && length(sib$roles)) {
    sib_labels <- c(ms = "bms", fs = "bfs")[sib$roles]
    ext <- .hidden_correction_sib_extend(hc,
      c(beta_vis, setNames(unname(sib$b[sib_labels]), sib_labels)),
      .sib_latent_matrix(cov_cm, cov_cf, cov_mf, sib), h, sib$roles)
    beta_total_assumed <- ext$beta_total[c("bD", "bm", "bf")]
    beta_total_assumed_std <- ext$beta_total_std[c("bD", "bm", "bf")]
    contrast <- function(beta) c(
      if ("ms" %in% sib$roles) c(bms = unname(beta["bms"]), delta_m = unname(beta["bm"] - beta["bms"])),
      if ("fs" %in% sib$roles) c(bfs = unname(beta["bfs"]), delta_f = unname(beta["bf"] - beta["bfs"])))
    sib_total <- list(raw = contrast(ext$beta_total), std = contrast(ext$beta_total_std),
      var_total = ext$var_total, Sigma_TT = ext$Sigma_TT, Sigma_xx = ext$Sigma_xx,
      Sigma_xh = ext$Sigma_xh, Sigma_hh = ext$Sigma_hh,
      lineage_factor = ext$lineage_factor, M = ext$M, invalid_reason = ext$invalid_reason)
  }

  list(
    beta_vis = beta_vis,
    beta_h_only_assumed = beta_h_only_assumed,
    beta_total_assumed = beta_total_assumed,
    beta_total_assumed_std = beta_total_assumed_std,
    rho_vis = unname(latent_cov_vis["rhoMF"]),
    rho_total_assumed = unname(latent_corr_total_assumed["rhoMF"]),
    cov_total_assumed = unname(latent_cov_total_assumed["rhoMF"]),
    latent_cov_vis = latent_cov_vis,
    latent_cov_total_assumed = latent_cov_total_assumed,
    latent_corr_total_assumed = latent_corr_total_assumed,
    var_total_assumed = var_total_assumed,
    sib_total = sib_total,
    hidden_correction_details = c(
      list(
        am_correction = am_correction,
        am_n_generations = if (am_correction == "n_gen") hc$n_generations else NULL,
        tm = hc$tm,
        tf = hc$tf,
        beta_total_assumed_std = beta_total_assumed_std,
        rho_total_assumed = unname(latent_corr_total_assumed["rhoMF"]),
        cov_total_assumed = unname(latent_cov_total_assumed["rhoMF"]),
        rho_hidden_assumed = hc$rho_hidden_assumed,
        eta_total_assumed = hc$eta_total_assumed,
        eta_hidden_assumed = hc$eta_hidden_assumed,
        var_total_assumed = var_total_assumed,
        latent_corr_total_assumed = latent_corr_total_assumed,
        rho_T_corr_assumed = if (is.null(hc$rho_T_corr_assumed)) NULL else hc$rho_T_corr_assumed,
        kappa_T_assumed = if (is.null(hc$kappa_T_assumed)) NULL else hc$kappa_T_assumed,
        pi_present = if (is.null(hc$pi_present)) NULL else hc$pi_present,
        t_ratio = if (is.null(hc$t_ratio)) NULL else hc$t_ratio,
        tH_m = if (is.null(hc$tH_m)) NULL else hc$tH_m,
        tH_f = if (is.null(hc$tH_f)) NULL else hc$tH_f,
        v_G = if (is.null(hc$v_G)) NULL else hc$v_G,
        v_H = if (is.null(hc$v_H)) NULL else hc$v_H,
        v_Gc = if (is.null(hc$v_Gc)) NULL else hc$v_Gc,
        v_Hc = if (is.null(hc$v_Hc)) NULL else hc$v_Hc,
        Sigma_xx = hc$Sigma_xx,
        Sigma_xh = hc$Sigma_xh,
        Sigma_hh = if (is.null(hc$Sigma_hh)) NULL else hc$Sigma_hh,
        M = hc$M
      ),
      if (!is.null(sib_total)) list(parental_sibs = sib_total),
      if (!is.null(hc$estimated_state)) list(estimated_state = hc$estimated_state),
      list(d = hc$d, d_child = if (is.null(hc$d_child)) hc$d else hc$d_child, q = hc$q)
    )
  )
}

#' Build lavaan-based uncertainty for visible and corrected estimates
#' @noRd
.build_lavaan_uncertainty <- function(fit, estimates, h, am_correction, am_n_generations,
                                      ci_level, sib_roles = character(0)) {
  V_full <- tryCatch(lavaan::vcov(fit), error = function(e) NULL)
  if (is.null(V_full)) {
    return(NULL)
  }

  sib_labels <- .sib_core_labels(sib_roles)
  core_labels <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f", sib_labels)
  theta_hat <- c(
    bD = unname(estimates$beta_vis["bD"]),
    bm = unname(estimates$beta_vis["bm"]),
    bf = unname(estimates$beta_vis["bf"]),
    rhoMF = unname(estimates$latent_cov_vis["rhoMF"]),
    eta_m = unname(estimates$latent_cov_vis["eta_m"]),
    eta_f = unname(estimates$latent_cov_vis["eta_f"]),
    .sib_theta(estimates$parental_sibs, sib_roles)
  )

  V <- matrix(0, nrow = length(core_labels), ncol = length(core_labels),
              dimnames = list(core_labels, core_labels))
  shared <- intersect(core_labels, rownames(V_full))
  if (length(shared) > 0) {
    V[shared, shared] <- V_full[shared, shared, drop = FALSE]
  }

  transform_blocks <- function(theta) {
    theta <- setNames(as.numeric(theta), core_labels)
    scale_est <- .compute_scale_estimates(
      beta_vis = theta[c("bD", "bm", "bf")],
      cov_cm = theta["eta_m"],
      cov_cf = theta["eta_f"],
      cov_mf = theta["rhoMF"],
      h = h,
      am_correction = am_correction,
      am_n_generations = am_n_generations,
      sib = .sib_input_from_theta(theta, sib_roles)
    )

    c(list(
      visible = c(
        bD = unname(scale_est$beta_vis["bD"]),
        bm = unname(scale_est$beta_vis["bm"]),
        bf = unname(scale_est$beta_vis["bf"]),
        rhoMF = unname(scale_est$rho_vis)
      ),
      total_raw = c(
        bD = unname(scale_est$beta_total_assumed["bD"]),
        bm = unname(scale_est$beta_total_assumed["bm"]),
        bf = unname(scale_est$beta_total_assumed["bf"]),
        rhoMF = unname(scale_est$cov_total_assumed)
      ),
      total_std = c(
        bD = unname(scale_est$beta_total_assumed_std["bD"]),
        bm = unname(scale_est$beta_total_assumed_std["bm"]),
        bf = unname(scale_est$beta_total_assumed_std["bf"]),
        rhoMF = unname(scale_est$rho_total_assumed)
      )
    ), if (length(sib_roles)) list(
      sib_total_raw = scale_est$sib_total$raw,
      sib_total_std = scale_est$sib_total$std
    ))
  }

  block_estimates <- transform_blocks(theta_hat)

  summarize_block <- function(block_name) {
    est <- block_estimates[[block_name]]
    if (is.null(est)) {
      return(NULL)
    }

    f_block <- function(theta) {
      transform_blocks(theta)[[block_name]]
    }

    J <- tryCatch(.finite_difference_jacobian(f_block, theta_hat), error = function(e) NULL)
    if (is.null(J) || any(!is.finite(J))) {
      se <- setNames(rep(NA_real_, length(est)), names(est))
      ci <- .normal_ci(est, se, ci_level)
      return(list(estimate = est, se = se, ci = ci, cov = matrix(NA_real_, length(est), length(est))))
    }

    cov_block <- J %*% V %*% t(J)
    cov_block <- (cov_block + t(cov_block)) / 2
    se <- sqrt(pmax(diag(cov_block), 0))
    names(se) <- names(est)

    list(
      estimate = est,
      se = se,
      ci = .normal_ci(est, se, ci_level),
      cov = cov_block
    )
  }

  out <- list(
    source = "lavaan",
    method = "lavaan_delta",
    ci_level = ci_level,
    visible = summarize_block("visible"),
    total_raw = summarize_block("total_raw"),
    total_std = summarize_block("total_std")
  )
  if (length(sib_roles)) {
    out$sib_total_raw <- summarize_block("sib_total_raw")
    out$sib_total_std <- summarize_block("sib_total_std")
  }
  out
}

#' Named sibling core parameters from the extracted sibling estimates
#' @noRd
.sib_theta <- function(sib_estimates, roles) {
  if (is.null(sib_estimates) || !length(roles)) return(NULL)
  vals <- c(sib_estimates$coefficients, sib_estimates$latent_cov)
  vals[.sib_core_labels(roles)]
}

#' Choose which uncertainty summary to display
#' @noRd
.resolve_uncertainty <- function(result, uncertainty = c("auto", "bootstrap", "lavaan", "none")) {
  uncertainty <- match.arg(uncertainty)

  if (uncertainty == "auto") {
    if (!is.null(result$bootstrap)) {
      uncertainty <- "bootstrap"
    } else if (!is.null(result$lavaan_inference)) {
      uncertainty <- "lavaan"
    } else {
      uncertainty <- "none"
    }
  }

  if (uncertainty == "bootstrap") {
    if (is.null(result$bootstrap)) {
      warning("Bootstrap uncertainty requested, but no bootstrap results are available; reporting point estimates only.")
      return(NULL)
    }
    bootstrap <- result$bootstrap
    bootstrap$visible$estimate <- c(result$beta_vis, rhoMF = unname(result$rho_vis))
    bootstrap$total_raw$estimate <- c(result$beta_total_assumed, rhoMF = unname(result$cov_total_assumed))
    bootstrap$total_std$estimate <- c(result$beta_total_assumed_std, rhoMF = unname(result$rho_total_assumed))
    for (sc in c("std", "raw")) {
      blk <- bootstrap[[paste0("parental_sibs_total_", sc)]]
      if (!is.null(blk)) {
        blk$estimate <- result$parental_sibs$corrected[[sc]][colnames(blk$ci)]
        bootstrap[[paste0("sib_total_", sc)]] <- blk
      }
    }
    return(.ravel_pi_merge_uncertainty(result, bootstrap, uncertainty))
  }

  if (uncertainty == "lavaan") {
    if (is.null(result$lavaan_inference)) {
      warning("lavaan uncertainty requested, but no lavaan-based uncertainty results are available; reporting point estimates only.")
      return(NULL)
    }
    return(.ravel_pi_merge_uncertainty(result, result$lavaan_inference, uncertainty))
  }

  NULL
}

#' Factor-score determinacy for each latent genetic factor (role-specific)
#'
#' @details
#' For each role r, with two indicators y = (PGS1,PGS2) and loadings λ_r,
#' residual variances θ_r, and Var(G_r)=1, determinacy is:
#'
#'   det_r^2 = λ_r' Σ_r^{-1} λ_r,   where Σ_r = λ_r λ_r' + Θ_r.
#'
#' det_r is the correlation between the latent factor and its regression
#' factor score (under normal-theory assumptions).
#'
#' @noRd
.compute_determinacy <- function(PE) {

  grab_load <- function(lv, ind) {
    x <- PE[PE$lhs == lv & PE$op == "=~" & PE$rhs == ind, "est"]
    if (length(x) == 0) NA_real_ else x[1]
  }
  grab_theta <- function(v) {
    x <- PE[PE$lhs == v & PE$op == "~~" & PE$rhs == v, "est"]
    if (length(x) == 0) NA_real_ else x[1]
  }

  roles <- list(
    mother = list(lv = "Gm", p1 = "PGS1_mother", p2 = "PGS2_mother"),
    father = list(lv = "Gf", p1 = "PGS1_father", p2 = "PGS2_father"),
    child  = list(lv = "Gc", p1 = "PGS1_child",  p2 = "PGS2_child")
  )

  out <- sapply(names(roles), function(r) {
    lv <- roles[[r]]$lv
    p1 <- roles[[r]]$p1
    p2 <- roles[[r]]$p2

    lam <- c(grab_load(lv, p1), grab_load(lv, p2))
    th  <- c(grab_theta(p1), grab_theta(p2))

    if (any(is.na(lam)) || any(is.na(th))) return(NA_real_)

    Sigma <- matrix(c(lam[1]^2 + th[1], lam[1]*lam[2],
                      lam[1]*lam[2], lam[2]^2 + th[2]), 2, 2, byrow = TRUE)

    det2 <- as.numeric(t(lam) %*% .safe_solve(Sigma, lam))
    det <- if (is.finite(det2) && det2 >= 0) sqrt(det2) else NA_real_
    det
  })

  out
}

#' Sibling inputs to the hidden correction (coefficients and latent covariances)
#' @noRd
.sib_correction_input <- function(sib_estimates, roles) {
  if (is.null(sib_estimates) || !length(roles)) return(NULL)
  b_labels <- c(ms = "bms", fs = "bfs")[roles]
  list(roles = roles, b = sib_estimates$coefficients[b_labels], cov = sib_estimates$latent_cov)
}

#' Labels of the sibling parameters entering the hidden correction
#' @noRd
.sib_core_labels <- function(roles) {
  c(if ("ms" %in% roles) c("bms", "kap_m", "chi_ms", "eta_ms"),
    if ("fs" %in% roles) c("bfs", "kap_f", "chi_fs", "eta_fs"),
    if (length(roles) == 2L) "zeta")
}

#' Rebuild the correction input from a named core vector
#' @noRd
.sib_input_from_theta <- function(theta, roles) {
  if (!length(roles)) return(NULL)
  b_labels <- c(ms = "bms", fs = "bfs")[roles]
  cov_labels <- setdiff(.sib_core_labels(roles), b_labels)
  list(roles = roles, b = theta[b_labels], cov = theta[cov_labels])
}

#' Labels and reporting names of the parental-sibling coefficients
#' @noRd
.parental_sib_coef_map <- function(roles) {
  map <- c(bms = "beta_MS", delta_m = "delta_M", bfs = "beta_FS", delta_f = "delta_F")
  keep <- c(if ("ms" %in% roles) c("bms", "delta_m"), if ("fs" %in% roles) c("bfs", "delta_f"))
  map[keep]
}

#' Point estimates for the parental-sibling extension
#' @noRd
.parental_sib_estimates <- function(PE, roles) {
  if (!length(roles)) return(NULL)
  grab <- function(lab) {
    x <- PE[!is.na(PE$label) & PE$label == lab, "est"]
    if (length(x) == 0) NA_real_ else x[1]
  }
  coef_labels <- names(.parental_sib_coef_map(roles))
  coefficients <- vapply(coef_labels, grab, numeric(1))

  own <- c(ms = "m", fs = "f")
  cov_labels <- c(
    if ("ms" %in% roles) c("kap_m", "chi_ms", "eta_ms"),
    if ("fs" %in% roles) c("kap_f", "chi_fs", "eta_fs"),
    if (length(roles) == 2L) "zeta")
  latent_cov <- vapply(cov_labels, grab, numeric(1))

  theta_label <- function(j, p) paste0("th", j, "_", p)
  pairs <- list()
  for (r in roles) for (p in c(own[[r]], setdiff(c("m", "f"), own[[r]]), "c"))
    pairs[[length(pairs) + 1L]] <- c(p, r)
  if (length(roles) == 2L) pairs[[length(pairs) + 1L]] <- c("ms", "fs")
  resid_cov <- resid_corr <- numeric(0)
  for (j in 1:2) for (pr in pairs) {
    lab <- sprintf("e%d_%s_%s", j, pr[1], pr[2])
    est <- grab(lab)
    if (is.na(est)) est <- 0  # fixed at zero in constrained modes
    resid_cov[[lab]] <- est
    resid_corr[[lab]] <- est / sqrt(grab(theta_label(j, pr[1])) * grab(theta_label(j, pr[2])))
  }

  role_names <- c(ms = "mat_sib", fs = "pat_sib")
  reliability <- do.call(rbind, lapply(roles, function(r) data.frame(
    role = role_names[[r]], PGS = c("PGS1", "PGS2"),
    reliability = vapply(1:2, function(j) {
      l <- grab(sprintf("l%s%d", r, j))
      l^2 / (l^2 + grab(theta_label(j, r)))
    }, numeric(1)),
    stringsAsFactors = FALSE)))

  determinacy <- vapply(roles, function(r) {
    lam <- c(grab(sprintf("l%s1", r)), grab(sprintf("l%s2", r)))
    th <- c(grab(theta_label(1, r)), grab(theta_label(2, r)))
    if (anyNA(c(lam, th))) return(NA_real_)
    Sigma <- tcrossprod(lam) + diag(th)
    det2 <- as.numeric(t(lam) %*% .safe_solve(Sigma, lam))
    if (is.finite(det2) && det2 >= 0) sqrt(det2) else NA_real_
  }, numeric(1))
  names(determinacy) <- role_names[roles]

  list(coefficients = coefficients, latent_cov = latent_cov,
    resid_cov = resid_cov, resid_corr = resid_corr,
    reliability = reliability, determinacy = determinacy)
}

#' lavaan (delta-method) inference for the parental-sibling coefficients
#' @noRd
.parental_sib_lavaan_table <- function(fit, roles, ci_level) {
  map <- .parental_sib_coef_map(roles)
  PE <- tryCatch(lavaan::parameterEstimates(fit, level = ci_level),
    error = function(e) NULL)
  if (is.null(PE)) return(NULL)
  rows <- lapply(names(map), function(lab) {
    i <- which(!is.na(PE$label) & PE$label == lab & PE$op %in% c("~", ":="))[1]
    data.frame(parameter = map[[lab]], label = lab,
      estimate = PE$est[i], se = PE$se[i], ci_low = PE$ci.lower[i], ci_high = PE$ci.upper[i],
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$p_value <- .uncertainty_p(out$estimate, out$se)
  out$uncertainty_source <- "lavaan/delta-method"
  rownames(out) <- NULL
  out
}

#' Parental-sibling inference table for the requested scale and uncertainty source
#'
#' Visible rows use lavaan's delta method for the defined contrasts (or the
#' bootstrap); corrected rows use the matching \code{sib_total_std} /
#' \code{sib_total_raw} uncertainty blocks.
#' @noRd
.parental_sib_table <- function(result, unc, scale = c("visible", "corrected_std", "corrected_raw")) {
  scale <- match.arg(scale)
  sib <- result$parental_sibs
  if (is.null(sib)) return(NULL)
  map <- .parental_sib_coef_map(sib$roles)
  make <- function(est, se, lo, hi, source) {
    out <- data.frame(parameter = unname(map), label = names(map), estimate = unname(est),
      se = unname(se), ci_low = unname(lo), ci_high = unname(hi), stringsAsFactors = FALSE)
    out$p_value <- .uncertainty_p(out$estimate, out$se)
    out$uncertainty_source <- source
    rownames(out) <- NULL
    out
  }
  none <- rep(NA_real_, length(map))
  if (scale != "visible") {
    est <- sib$corrected[[if (scale == "corrected_std") "std" else "raw"]][names(map)]
    block <- if (is.null(unc)) NULL else unc[[if (scale == "corrected_std") "sib_total_std" else "sib_total_raw"]]
    if (is.null(block)) return(make(est, none, none, none, "none"))
    return(make(est, block$se[names(map)], block$ci["lower", names(map)],
      block$ci["upper", names(map)], .uncertainty_label(unc, "total_std")))
  }
  if (!is.null(unc) && identical(unc$source, "bootstrap") && !is.null(result$bootstrap$parental_sibs)) {
    block <- result$bootstrap$parental_sibs
    return(make(sib$coefficients[names(map)], block$se[names(map)], block$ci["lower", names(map)],
      block$ci["upper", names(map)], sprintf("bootstrap (%d reps)", result$bootstrap$n_successful)))
  }
  if (is.null(unc)) {
    out <- sib$lavaan_table
    out[c("se", "ci_low", "ci_high", "p_value")] <- NA_real_
    out$uncertainty_source <- "none"
    return(out)
  }
  sib$lavaan_table
}

# =============================================================================
# MAIN FITTING FUNCTION
# =============================================================================

#' Fit RAVEL SEM to Trio Data
#'
#' @description
#' Fits a structural equation model to family trio data with two polygenic scores
#' per individual, estimating direct-offspring and parental genetic-liability-
#' correlated coefficients while accounting for measurement error, assortative
#' mating, and correlated cross-person PGS errors due to shared GWAS weights.
#'
#' @details
#' The SEM has three components:
#'
#' **(1) Measurement model (ORIV-as-SEM):**
#' For each mother, father, or child role and each of the two PGS versions,
#' \deqn{PGS_{k,r} = \lambda_{k,r} G_r + \varepsilon_{k,r},}
#' with \eqn{\operatorname{Var}(G_r)=1} (\code{std.lv = TRUE}) and free
#' residual variances.
#'
#' **(2) Trio genetic covariance structure:**
#'   Cov(G_m, G_f) = ρ_MF (assortative mating)
#'   Cov(G_c, G_m) = η_m
#'   Cov(G_c, G_f) = η_f
#'
#' **(3) Outcome model:**
#'   Y = β_D G_c + β_M G_m + β_F G_f + Γ'W + ε_Y
#'
#' **Total-liability correction (assumption-based):**
#' RAVEL returns post-estimation coefficients on raw founder-scaled total
#' liabilities \eqn{T^*} and on role-standardized total liabilities \eqn{T}.
#' It recovers standardized child-on-parent projection coefficients from the
#' visible trio covariance block and uses them to build the cross-component block
#' \eqn{\Sigma_{GH}} induced by assortative mating. See the
#' "beta_total_assumed" and "beta_total_assumed_std" outputs and the Methods
#' write-up for details.
#'
#' @param dat Data frame with PGS columns for mother/father/child (two PGS each),
#'   an outcome column, and optionally fixed-effect columns. Column names can be
#'   specified via the mapping parameters below.
#' @param h Numeric [0,1). Hidden additive-variance fraction \eqn{\pi} on the
#'   founder or random-mating reference scale. At \code{h = 0} there is no
#'   hidden component.
#' @param outcome Character. "continuous" (uses MLR) or "binary_liability"
#'   (uses WLSMV with probit link for a liability-threshold model).
#' @param measurement Character. Only "free_loadings" is supported. The previous
#'   "equal_loadings" option has been removed (it was redundant and can hinder
#'   identification with cross-person residuals).
#' @param crossperson Character. Cross-person residual covariance constraints:
#'   "free" (default) or "equal_pc".
#' @param fixed_crossperson Optional named numeric vector fixing selected
#'   observed cross-person residual covariances to user-specified values.
#'   Names must match the active residual labels:
#'   \code{e1_mf}, \code{e1_mc}, \code{e1_fc}, \code{e2_mf}, \code{e2_mc},
#'   \code{e2_fc} when \code{crossperson = "free"}, or
#'   \code{e1_mf}, \code{e1_pc}, \code{e2_mf}, \code{e2_pc} when
#'   \code{crossperson = "equal_pc"}.
#' @param fixed_crossperson_corr Optional named numeric vector fixing selected
#'   observed cross-person residual correlations to user-specified values.
#'   Currently supported only when \code{crossperson = "free"} because the
#'   equal-parent-child structure is defined on the covariance scale. Names
#'   must be drawn from \code{e1_mf}, \code{e1_mc}, \code{e1_fc},
#'   \code{e2_mf}, \code{e2_mc}, \code{e2_fc}. The underlying residual
#'   covariances are then constrained as
#'   \code{cov = corr * sqrt(var_1 * var_2)} using the estimated indicator
#'   residual variances.
#'   With parental siblings, both residual arguments also accept labels
#'   \code{e1_m_ms}, \code{e1_f_ms}, \code{e1_c_ms} for the maternal
#'   sibling, \code{e1_f_fs}, \code{e1_m_fs}, \code{e1_c_fs} for the
#'   paternal sibling, and \code{e1_ms_fs} when both are included, plus the
#'   corresponding \code{e2_*} labels. Unspecified sibling residuals are free
#'   unless a constrained fitting mode supplies the pedigree restrictions.
#' @param am_correction Character. Raw assortative-mating state used for the
#'   total-liability mapping:
#'   \code{"n_gen"} (default) assumes a finite user-specified number of generations of
#'   assortative mating and recursively propagates the covariance state forward; setting
#'   \code{am_n_generations = 1} reproduces the one-generation correction.
#'   \code{"equilibrium"} assumes assortative mating has reached equilibrium over
#'   many generations. This choice also determines raw-scale assortative-mating
#'   inflation at \code{h = 0}, although standardized total coefficients then
#'   equal visible coefficients.
#'   \code{"estimated_state"} recovers the founder-scaled covariance state
#'   from fitted spouse and parent-child correlations, without a generation
#'   count. Requires matching \code{fit_constraints = "estimated_state"}.
#' @param am_n_generations Integer. Number of assortative-mating generations to
#'   use when \code{am_correction = "n_gen"}. Must be at least 1. The value
#'   \code{1} means the mating parents are on the founder scale and the child is
#'   produced by the first mating round. Ignored otherwise.
#' @param fix_parent_paths Logical. If TRUE, constrain β_M = β_F = 0.
#' @param fixed_effects Character vector of numeric covariate column names to
#'   include jointly. Encode categorical covariates as treatment dummies, or use
#'   \code{ravel_dataframe_report}. With a binary outcome, conditional-probit
#'   moments from all included scores and these covariates are mapped to a marginal
#'   variance-one liability before SEM fitting. Positive working weights and a
#'   full joint sandwich covariance are used for robust SEs. This assumes a
#'   linear, homoscedastic conditional-probit liability model.
#' @param fit_constraints Character. Optional strong measurement/inheritance
#'   restrictions: \code{"none"} (default), \code{"equilibrium"},
#'   \code{"n_gen"}, or \code{"estimated_state"}. Scores retain separate unknown reliabilities and free
#'   marginal scales, independent cross-score residuals, zero spouse residual
#'   correlation and parent-child residual correlation 0.5. At equilibrium,
#'   each score's reliability is shared across roles and latent parent-child
#'   correlation is \eqn{(1+\rho)/2}. Finite-generation constraints follow
#'   random-mating founders and constant total-liability spouse correlation;
#'   child reliability evolves with genetic variance. For two or more
#'   generations, these fitting restrictions are conditional on supplied
#'   \code{h}. The matching \code{am_correction} is selected automatically
#'   if omitted; an explicit conflicting value is rejected. Loadings must have
#'   the same positive orientation. See the package constrained-fitting notes.
#'   Estimated-state fitting estimates variance growth independently of
#'   \code{h}, retaining Mendelian segregation and a shared genetic target.
#'   It relaxes exact constant-strength integer-generation history; a
#'   compatibility diagnostic is returned, not an estimated integer count.
#' @param n_boot Integer. Number of bootstrap replicates (0 = no bootstrap).
#' @param cluster_var Character. Column name for cluster bootstrap sampling.
#' @param ci_level Numeric. Confidence level for bootstrap CIs.
#' @param seed Integer. Random seed for reproducibility.
#' @param verbose Logical. Print model syntax?
#' @param col_pgs1_mother Character. Column name for mother's PGS1.
#' @param col_pgs2_mother Character. Column name for mother's PGS2.
#' @param col_pgs1_father Character. Column name for father's PGS1.
#' @param col_pgs2_father Character. Column name for father's PGS2.
#' @param col_pgs1_child Character. Column name for child's PGS1.
#' @param col_pgs2_child Character. Column name for child's PGS2.
#' @param col_outcome Character. Column name for the outcome variable.
#' @param pi Optional preferred name for the external hidden genetic variance
#'   fraction, in [0, 1). It is the same quantity as \code{h}. If both are
#'   supplied explicitly, they must agree.
#' @param pi_se Optional nonnegative standard error of \code{pi} on its raw
#'   probability scale. When supplied, corrected estimates automatically include
#'   this uncertainty by the joint delta method; visible estimates retain their
#'   uncertainty conditional on pi. Omit it to treat pi as fixed. Supply
#'   \code{pi} or \code{h} explicitly when using this argument.
#' @param core_pi_cov Optional named vector of covariances between all
#'   conditional core estimates and external pi, on its raw probability scale.
#'   Requires \code{pi_se}; omitted means independent external estimation.
#' @param covariate_mapping_confirmed Logical. Set TRUE with \code{pi_se} and
#'   jointly fitted covariates only after confirming the visible-to-total
#'   correction remains valid conditional on those covariates.

#' @param parental_sibs Character. Add full siblings of the parents:
#'   \code{"none"} (default), \code{"maternal"} (a sibling of the mother),
#'   \code{"paternal"} (a sibling of the father), or \code{"both"}. Each
#'   sibling contributes the same two scores, built with the same weights as
#'   the trio. The outcome is then also regressed on each sibling's visible
#'   liability. The sibling coefficient can absorb shared lineage background;
#'   \code{delta_M = bm - bMS} (or \code{delta_F = bf - bFS}) is the
#'   parent-specific contrast conditional on the child's and other parent's
#'   genetic liabilities. Its causal interpretation requires equal parent/sibling
#'   liability distributions, equal conditional lineage contributions, no
#'   effect of the sibling on the child, and no remaining parent-specific
#'   confounding or selection bias. Continuous outcomes retain families
#'   with missing sibling scores through FIML; binary outcomes use complete
#'   cases, and the reported sample size counts retained families. Works with
#'   every \code{fit_constraints} mode; with \code{pi > 0}, corrected sibling
#'   coefficients and contrasts are reported. Incompatible joint visible/hidden
#'   covariances produce unavailable corrected estimates with a warning. If
#'   \code{pi_se} was requested, the visible fit is retained, corrected
#'   uncertainty is unavailable, and \code{pi_correction_error} records why.
#' @param col_pgs1_mat_sib,col_pgs2_mat_sib Column names for the maternal
#'   sibling's PGS1 and PGS2.
#' @param col_pgs1_pat_sib,col_pgs2_pat_sib Column names for the paternal
#'   sibling's PGS1 and PGS2.
#'
#' @return A list with components (selected):
#'   \describe{
#'     \item{fit}{lavaan fit object}
#'     \item{model}{Character string with model syntax}
#'     \item{N}{Sample size used for fitting}
#'     \item{beta_vis}{Named vector \code{bD}, \code{bm}, \code{bf} on the role-standardized visible-liability scale}
#'     \item{beta_h_only_assumed}{Assumption-based h-only scaling correction}
#'     \item{beta_total_assumed}{Coefficients \eqn{\beta^{raw}} on raw founder-scaled total liabilities \eqn{T^*}}
#'     \item{beta_total_assumed_std}{Coefficients \eqn{\beta=D_T\beta^{raw}} on role-standardized total liabilities \eqn{T}}
#'     \item{rho_vis}{Visible-scale spouse correlation}
#'     \item{rho_total_assumed}{Standardized spouse correlation on the assumed total-liability scale}
#'     \item{cov_total_assumed}{Raw spouse covariance on the assumed total-liability scale}
#'     \item{latent_cov_vis}{Entries \code{rhoMF}, \code{eta_m}, and \code{eta_f}; the latter two are offspring-parent visible covariance parameters}
#'     \item{latent_cov_total_assumed}{Raw covariances among \eqn{T_c^*,T_m^*,T_p^*}; the spouse entry is \eqn{\kappa_T} despite its legacy \code{rhoMF} name}
#'     \item{latent_corr_total_assumed}{Correlations among role-standardized total liabilities; the spouse entry is \eqn{\rho_T}}
#'     \item{var_total_assumed}{Role-specific variances \eqn{\operatorname{Var}(T_r^*)}}
#'     \item{crossperson_resid}{Estimated observed cross-person residual covariances}
#'     \item{crossperson_resid_corr}{Estimated observed cross-person residual correlations}
#'     \item{components_vis}{Quadratic variance decomposition over all included family roles on the visible scale}
#'     \item{components_total_assumed}{Variance decomposition on assumed total scale}
#'     \item{R2_model}{Variance explained by the visible latent linear predictor}
#'     \item{R2_total_assumed}{Variance explained by the assumed total-liability linear predictor, conditional on the supplied hidden-variance and assortment model}
#'     \item{reliability}{Indicator-level reliability estimates}
#'     \item{determinacy}{Factor-score determinacy per role}
#'     \item{within_person_corr}{Within-person Cor(PGS1,PGS2) per role}
#'     \item{vanilla_regression}{Plain standardized trio regressions on observed PGS1 and PGS2 using \code{lm()} for continuous outcomes or \code{glm(link = "probit")} for binary outcomes, with the same fixed effects if supplied}
#'     \item{lavaan_inference}{lavaan-based SEs and normal-approximation CIs for visible and corrected estimates}
#'     \item{bootstrap}{Bootstrap results (if n_boot > 0), conditional on supplied pi}
#'     \item{pi_correction}{When pi_se is supplied, the external-pi corrected table
#'       with full propagation metadata. Print and data-frame methods use its
#'       corrected uncertainty automatically. Core covariance and the stored
#'       lavaan/bootstrap summaries remain conditional on pi.}
#'     \item{parental_sibs}{When \code{parental_sibs} is not \code{"none"}: the
#'       sibling coefficients (\code{bms}, \code{bfs}) and parent-specific
#'       contrasts (\code{delta_m}, \code{delta_f}) with a lavaan inference table,
#'       the sibling latent covariances, within-score residual covariances and
#'       correlations, reliabilities, determinacy, and the number of families
#'       with each sibling observed.}
#'   }
#'
#' @export
#' @importFrom lavaan sem inspect parameterEstimates
#' @importFrom stats var sd quantile cor
ravel_fit <- function(dat,
                      h = 0,
                      outcome = c("continuous", "binary_liability"),
                      measurement = c("free_loadings"),
                      crossperson = c("free", "equal_pc"),
                      fixed_crossperson = NULL,
                      fixed_crossperson_corr = NULL,
                      am_correction = c("n_gen", "equilibrium", "estimated_state"),
                      am_n_generations = 1L,
                      fix_parent_paths = FALSE,
                      fixed_effects = NULL,
                      n_boot = 0,
                      cluster_var = "FID",
                      ci_level = 0.95,
                      seed = NULL,
                      verbose = FALSE,
                      col_pgs1_mother = "PGS1_mother",
                      col_pgs2_mother = "PGS2_mother",
                      col_pgs1_father = "PGS1_father",
                      col_pgs2_father = "PGS2_father",
                      col_pgs1_child  = "PGS1_child",
                      col_pgs2_child  = "PGS2_child",
                      col_outcome     = "Y",
                      fit_constraints = c("none", "equilibrium", "n_gen", "estimated_state"),
                      pi = NULL, pi_se = NULL, core_pi_cov = NULL,
                      covariate_mapping_confirmed = FALSE,
                      parental_sibs = c("none", "maternal", "paternal", "both"),
                      col_pgs1_mat_sib = "PGS1_mat_sib",
                      col_pgs2_mat_sib = "PGS2_mat_sib",
                      col_pgs1_pat_sib = "PGS1_pat_sib",
                      col_pgs2_pat_sib = "PGS2_pat_sib") {

  # Resolve the preferred pi name before constructing any fitting constraints.
  h_supplied <- !missing(h)
  if (!is.null(pi)) {
    .ravel_pi_se_variance(pi, 0)
    if (h_supplied && !isTRUE(all.equal(as.numeric(h), as.numeric(pi), tolerance = 1e-12)))
      stop("h and pi specify the same hidden fraction; supply only pi or matching values")
    h <- pi
  }
  if (!is.null(pi_se)) {
    if (is.null(pi) && !h_supplied) stop("Supply pi (or h) explicitly when using pi_se")
    .ravel_pi_se_variance(h, pi_se)
    if (!is.logical(covariate_mapping_confirmed) || length(covariate_mapping_confirmed) != 1L ||
        is.na(covariate_mapping_confirmed))
      stop("covariate_mapping_confirmed must be TRUE or FALSE")
    if (length(fixed_effects) > 0L && !covariate_mapping_confirmed)
      stop("Confirm that the standard visible-to-total mapping applies conditional on the jointly fitted covariates")
    .ravel_pi_validate_cross(core_pi_cov, pi_se^2,
      c(.report_core_labels, .sib_core_labels(.parental_sib_roles(match.arg(parental_sibs)))))
  } else if (!is.null(core_pi_cov)) {
    stop("core_pi_cov requires pi_se")
  }

  # Input validation
  correction_missing <- missing(am_correction)
  fit_constraints <- match.arg(fit_constraints)
  outcome <- match.arg(outcome)
  measurement <- match.arg(measurement)
  parental_sibs <- match.arg(parental_sibs)
  sib_roles <- .parental_sib_roles(parental_sibs)
  crossperson <- match.arg(crossperson)
  fixed_crossperson <- .normalize_fixed_crossperson(
    crossperson,
    fixed_crossperson,
    arg_name = "fixed_crossperson", sib_roles = sib_roles
  )
  fixed_crossperson_corr <- .normalize_fixed_crossperson(
    crossperson,
    fixed_crossperson_corr,
    arg_name = "fixed_crossperson_corr", sib_roles = sib_roles
  )
  residual_overlap <- names(fixed_crossperson)[!is.na(fixed_crossperson) & !is.na(fixed_crossperson_corr)]
  if (length(residual_overlap))
    stop("The same cross-person residual cannot be fixed on both the covariance and correlation scales: ",
      paste(residual_overlap, collapse = ", "))
  am_correction <- match.arg(am_correction)
  constraint_options <- .ravel_prepare_constraints(
    fit_constraints, h, am_correction, am_n_generations, correction_missing,
    crossperson, fixed_crossperson, fixed_crossperson_corr, sib_roles = sib_roles)
  am_correction <- constraint_options$am_correction
  if (am_correction == "estimated_state" && fit_constraints != "estimated_state")
    stop("am_correction = 'estimated_state' requires fit_constraints = 'estimated_state'.")
  crossperson <- constraint_options$crossperson
  fixed_crossperson <- constraint_options$fixed_cov
  fixed_crossperson_corr <- constraint_options$fixed_corr

  if (am_correction == "n_gen") {
    if (!is.numeric(am_n_generations) || length(am_n_generations) != 1 ||
        is.na(am_n_generations) || am_n_generations < 1) {
      stop("am_n_generations must be a single integer >= 1 when am_correction = 'n_gen'.")
    }
    am_n_generations <- as.integer(am_n_generations)
  }

  stopifnot(h >= 0 && h < 1)

  parental_sibs <- match.arg(parental_sibs)
  sib_roles <- .parental_sib_roles(parental_sibs)
  if (length(sib_roles) && isTRUE(fix_parent_paths))
    stop("fix_parent_paths cannot be combined with parental_sibs.")

  if (!is.null(seed)) set.seed(seed)

  # Some R installations report NA from parallel::detectCores(), which breaks
  # lavaan option initialization before fitting starts. Patch to a safe single-
  # core fallback for the current R session when needed.
  .ensure_lavaan_runtime(verbose = verbose)

  # --- Column name mapping ---
  # Rename user-specified columns to internal names so the rest of the code
  # (model syntax, extraction, etc.) can use a fixed column vocabulary.
  col_map <- c(
    PGS1_mother = col_pgs1_mother,
    PGS2_mother = col_pgs2_mother,
    PGS1_father = col_pgs1_father,
    PGS2_father = col_pgs2_father,
    PGS1_child  = col_pgs1_child,
    PGS2_child  = col_pgs2_child,
    Y           = col_outcome
  )
  sib_col_map <- c(
    PGS1_mat_sib = col_pgs1_mat_sib,
    PGS2_mat_sib = col_pgs2_mat_sib,
    PGS1_pat_sib = col_pgs1_pat_sib,
    PGS2_pat_sib = col_pgs2_pat_sib
  )
  col_map <- c(col_map, sib_col_map[.parental_sib_columns(sib_roles)])

  # Check that user-specified columns exist
  missing_user <- col_map[!col_map %in% names(dat)]
  if (length(missing_user) > 0) {
    stop("Missing columns in data: ",
         paste(sprintf("%s (expected '%s')", names(missing_user), missing_user),
               collapse = ", "))
  }

  # Rename to internal names (only if they differ from defaults)
  if (any(names(col_map) != col_map)) {
    for (internal_name in names(col_map)) {
      user_name <- col_map[[internal_name]]
      if (user_name != internal_name) {
        dat[[internal_name]] <- dat[[user_name]]
      }
    }
  }

  # Parental-sibling scores may be missing for families without that sibling.
  n_with_sib <- NULL
  if (length(sib_roles)) {
    sib_cols <- .parental_sib_columns(sib_roles)
    if (!all(vapply(dat[sib_cols], is.numeric, logical(1))))
      stop("Parental-sibling score columns must be numeric.")
    if (any(vapply(dat[sib_cols], function(x) any(!is.na(x) & !is.finite(x)), logical(1))))
      stop("Parental-sibling score columns must be finite or missing (NA).")
    role_cols <- list(mat_sib = c("PGS1_mat_sib", "PGS2_mat_sib"),
                      pat_sib = c("PGS1_pat_sib", "PGS2_pat_sib"))[c(ms = "mat_sib", fs = "pat_sib")[sib_roles]]
    n_with_sib <- vapply(role_cols, function(cols)
      as.numeric(sum(stats::complete.cases(dat[cols]))), numeric(1))
    if (any(n_with_sib < 50))
      stop("Each requested parental sibling needs both scores observed in at least 50 families (observed: ",
           paste(sprintf("%s = %d", names(n_with_sib), n_with_sib), collapse = ", "), ").")
  }

  # Validate fixed effects
  if (!is.null(fixed_effects)) {
    missing_fe <- setdiff(fixed_effects, names(dat))
    if (length(missing_fe) > 0) {
      stop("Missing fixed effect columns: ", paste(missing_fe, collapse = ", "))
    }
  } else {
    fixed_effects <- character(0)
  }

  # Use the same retained sample for binary sibling fitting, diagnostics and
  # comparators. Continuous fits keep partially observed families through FIML.
  if (length(sib_roles) && outcome == "binary_liability") {
    analysis_columns <- unique(c(names(col_map), fixed_effects))
    complete <- stats::complete.cases(dat[analysis_columns])
    if (any(!complete)) {
      warning("Binary outcomes are fitted with listwise deletion; dropping ",
        sum(!complete), " families with missing model variables, including parental-sibling scores.")
      dat <- dat[complete, , drop = FALSE]
      n_with_sib[] <- nrow(dat)
      if (nrow(dat) < 50L)
        stop("Binary parental-sibling fitting requires at least 50 complete families after listwise deletion.")
    }
  }

  # Estimation settings
  estimator <- if (outcome == "binary_liability") "WLSMV" else "MLR"
  ordered_Y <- if (outcome == "binary_liability") "Y" else NULL

  # Build model syntax (std.lv=TRUE is required for free-loadings identification)
  model <- .build_model(
    fix_parent_paths,
    fixed_effects,
    crossperson,
    fixed_crossperson = fixed_crossperson,
    fixed_crossperson_corr = fixed_crossperson_corr,
    parental_sibs = parental_sibs,
    constrained = fit_constraints != "none"
  )

  if (fit_constraints != "none") {
    model <- paste(model, .ravel_constraint_syntax(fit_constraints, h, am_n_generations, sib_roles), sep = "\n")
    model <- .ravel_constraint_starts(model, dat, fit_constraints, h, am_n_generations, sib_roles)
  } else if (length(sib_roles) && any(!is.na(
      fixed_crossperson_corr[names(.parental_sib_residual_benchmarks(sib_roles))]))) {
    # Starting values only: explicit residual correlations do not impose a
    # variance history or reliability constraints on an unconstrained fit.
    model <- .ravel_constraint_starts(model, dat, "n_gen", h, 1L, sib_roles)
  }

  if (verbose) {
    cat("=== MODEL SYNTAX ===\n")
    cat(model)
    cat("\n====================\n\n")
  }

  # ---------------------------------------------------------------------------
  # Fitting function
  # ---------------------------------------------------------------------------
  .fit_once <- function(data) {
    input <- if (outcome == "binary_liability" && length(fixed_effects) > 0L) {
      .ravel_binary_moments(data, model)
    } else {
      list(data = data)
    }
    fit_model <- function(...) do.call(lavaan::sem, c(list(
      model = model,
      std.lv = TRUE,
      fixed.x = FALSE,
      conditional.x = FALSE,
      estimator = estimator,
      missing = if (outcome == "continuous") "fiml" else "listwise",
      ordered = ordered_Y,
      check.gradient = FALSE,
      warn = FALSE,
      ...
    ), input))

    if (length(fixed_effects) == 0L) return(fit_model())

    valid_fit <- function(fit) {
      if (is.null(fit) || !isTRUE(lavaan::lavInspect(fit, "converged"))) return(FALSE)
      if (!suppressWarnings(isTRUE(lavaan::lavInspect(fit, "post.check")))) return(FALSE)
      V <- tryCatch(lavaan::vcov(fit), error = function(e) NULL)
      core <- c("bD", if (!fix_parent_paths) c("bm", "bf"), "rhoMF", "eta_m", "eta_f")
      is.matrix(V) && all(core %in% rownames(V)) &&
        all(is.finite(V)) && all(diag(V)[core] > 0)
    }

    fit <- suppressWarnings(tryCatch(fit_model(), error = function(e) NULL))
    if (!valid_fit(fit)) {
      fit <- suppressWarnings(tryCatch(fit_model(start = "simple"), error = function(e) NULL))
    }
    if (!valid_fit(fit)) {
      stop("Covariate SEM did not produce an admissible converged fit with finite SEs.")
    }
    fit
  }

  # ---------------------------------------------------------------------------
  # Extraction function
  # ---------------------------------------------------------------------------
  .extract <- function(fit, data_for_r2 = NULL) {

    PE <- lavaan::parameterEstimates(fit)

    # Helper to get parameter estimate by label
    grab <- function(lab) {
      x <- PE[PE$label == lab, "est"]
      if (length(x) == 0) NA_real_ else x[1]
    }

    grab_or_fixed <- function(lab) {
      if (lab %in% names(fixed_crossperson) && !is.na(fixed_crossperson[[lab]])) {
        return(unname(fixed_crossperson[[lab]]))
      }
      grab(lab)
    }

    # Effect estimates (visible latent factors)
    bD <- grab("bD")
    bm <- if (fix_parent_paths) 0 else grab("bm")
    bf <- if (fix_parent_paths) 0 else grab("bf")

    beta_vis <- c(bD = bD, bm = bm, bf = bf)

    # Latent covariances (visible scale; Var(G)=1 by construction)
    cov_mf <- grab("rhoMF")
    cov_cm <- grab("eta_m")
    cov_cf <- grab("eta_f")

    latent_cov_vis <- c(rhoMF = cov_mf, eta_m = cov_cm, eta_f = cov_cf)

    # Cross-person residual covariances
    crossperson_est <- if (crossperson == "free") {
      c(
        e1_mf = grab_or_fixed("e1_mf"),
        e1_mc = grab_or_fixed("e1_mc"),
        e1_fc = grab_or_fixed("e1_fc"),
        e2_mf = grab_or_fixed("e2_mf"),
        e2_mc = grab_or_fixed("e2_mc"),
        e2_fc = grab_or_fixed("e2_fc")
      )
    } else {
      c(
        e1_mf = grab_or_fixed("e1_mf"),
        e1_pc = grab_or_fixed("e1_pc"),
        e2_mf = grab_or_fixed("e2_mf"),
        e2_pc = grab_or_fixed("e2_pc")
      )
    }

    theta_obs <- c(
      PGS1_mother = grab("th1_m"),
      PGS1_father = grab("th1_f"),
      PGS1_child  = grab("th1_c"),
      PGS2_mother = grab("th2_m"),
      PGS2_father = grab("th2_f"),
      PGS2_child  = grab("th2_c")
    )

    crossperson_corr <- if (crossperson == "free") {
      setNames(c(
        unname(crossperson_est["e1_mf"]) /
          sqrt(unname(theta_obs["PGS1_mother"]) * unname(theta_obs["PGS1_father"])),
        unname(crossperson_est["e1_mc"]) /
          sqrt(unname(theta_obs["PGS1_mother"]) * unname(theta_obs["PGS1_child"])),
        unname(crossperson_est["e1_fc"]) /
          sqrt(unname(theta_obs["PGS1_father"]) * unname(theta_obs["PGS1_child"])),
        unname(crossperson_est["e2_mf"]) /
          sqrt(unname(theta_obs["PGS2_mother"]) * unname(theta_obs["PGS2_father"])),
        unname(crossperson_est["e2_mc"]) /
          sqrt(unname(theta_obs["PGS2_mother"]) * unname(theta_obs["PGS2_child"])),
        unname(crossperson_est["e2_fc"]) /
          sqrt(unname(theta_obs["PGS2_father"]) * unname(theta_obs["PGS2_child"]))
      ), c("e1_mf", "e1_mc", "e1_fc", "e2_mf", "e2_mc", "e2_fc"))
    } else {
      setNames(c(
        unname(crossperson_est["e1_mf"]) /
          sqrt(unname(theta_obs["PGS1_mother"]) * unname(theta_obs["PGS1_father"])),
        unname(crossperson_est["e1_pc"]) /
          sqrt(unname(theta_obs["PGS1_mother"]) * unname(theta_obs["PGS1_child"])),
        unname(crossperson_est["e1_pc"]) /
          sqrt(unname(theta_obs["PGS1_father"]) * unname(theta_obs["PGS1_child"])),
        unname(crossperson_est["e2_mf"]) /
          sqrt(unname(theta_obs["PGS2_mother"]) * unname(theta_obs["PGS2_father"])),
        unname(crossperson_est["e2_pc"]) /
          sqrt(unname(theta_obs["PGS2_mother"]) * unname(theta_obs["PGS2_child"])),
        unname(crossperson_est["e2_pc"]) /
          sqrt(unname(theta_obs["PGS2_father"]) * unname(theta_obs["PGS2_child"]))
      ), c("e1_mf", "e1_pc_mother", "e1_pc_father",
           "e2_mf", "e2_pc_mother", "e2_pc_father"))
    }

    sib_estimates <- .parental_sib_estimates(PE, sib_roles)
    sib_input <- .sib_correction_input(sib_estimates, sib_roles)

    scale_estimates <- .compute_scale_estimates(
      beta_vis = beta_vis,
      cov_cm = cov_cm,
      cov_cf = cov_cf,
      cov_mf = cov_mf,
      h = h,
      am_correction = am_correction,
      am_n_generations = am_n_generations,
      sib = sib_input
    )

    beta_h_only_assumed <- scale_estimates$beta_h_only_assumed
    beta_total_assumed <- scale_estimates$beta_total_assumed
    beta_total_assumed_std <- scale_estimates$beta_total_assumed_std
    rho_vis <- scale_estimates$rho_vis
    rho_total_assumed <- scale_estimates$rho_total_assumed
    cov_total_assumed <- scale_estimates$cov_total_assumed
    latent_cov_total_assumed <- scale_estimates$latent_cov_total_assumed
    latent_corr_total_assumed <- scale_estimates$latent_corr_total_assumed
    var_total_assumed <- scale_estimates$var_total_assumed

    # SEM-implied R2 for Y explained by visible latent factors
    # Note: when h>0, the SEM residual includes omitted (hidden) genetic variance.
    theta_Y <- PE[PE$lhs == "Y" & PE$op == "~~" & PE$rhs == "Y", "est"]
    theta_Y <- if (length(theta_Y) == 0) NA_real_ else theta_Y[1]

    components_vis <- if (is.null(sib_input)) {
      .variance_components_trio(beta_vis, cov_cm, cov_cf, cov_mf)
    } else {
      .variance_components_general(c(beta_vis, sib_input$b),
        .sib_latent_matrix(cov_cm, cov_cf, cov_mf, sib_input),
        c("child", "mother", "father", c(ms = "mat_sib", fs = "pat_sib")[sib_roles]))
    }
    s2_gen_vis <- sum(components_vis)
    varY_model <- if (length(fixed_effects) > 0L) {
      lavaan::lavInspect(fit, "cov.ov")["Y", "Y"]
    } else {
      s2_gen_vis + theta_Y
    }
    R2_model <- s2_gen_vis / varY_model

    components_total_assumed <- .variance_components_trio_general(
      beta_total_assumed,
      var_c = var_total_assumed["child"],
      var_m = var_total_assumed["mother"],
      var_f = var_total_assumed["father"],
      cov_cm = latent_cov_total_assumed["eta_m"],
      cov_cf = latent_cov_total_assumed["eta_f"],
      cov_mf = latent_cov_total_assumed["rhoMF"]
    )
    if (!is.null(sib_input)) {
      st <- scale_estimates$sib_total
      components_total_assumed <- .variance_components_general(
        c(beta_total_assumed, st$raw[names(sib_input$b)]), st$Sigma_TT, rownames(st$Sigma_TT))
    }
    s2_gen_total_assumed <- sum(components_total_assumed)

    # "Total genetic R2" is only directly interpretable if Y is standardized.
    # For convenience, compute:
    #   - continuous: divide by empirical Var(Y)
    #   - binary liability: divide by the SEM-implied liability variance
    if (outcome == "continuous") {
      varY <- if (is.null(data_for_r2)) NA_real_ else stats::var(data_for_r2$Y, na.rm = TRUE)
      R2_total_assumed <- s2_gen_total_assumed / varY
    } else {
      R2_total_assumed <- s2_gen_total_assumed / varY_model
    }

    list(
      beta_vis = beta_vis,
      beta_h_only_assumed = beta_h_only_assumed,
      beta_total_assumed = beta_total_assumed,
      beta_total_assumed_std = beta_total_assumed_std,

      rho_vis = rho_vis,
      rho_total_assumed = rho_total_assumed,
      cov_total_assumed = cov_total_assumed,

      latent_cov_vis = latent_cov_vis,
      latent_cov_total_assumed = latent_cov_total_assumed,
      latent_corr_total_assumed = latent_corr_total_assumed,
      var_total_assumed = var_total_assumed,

      crossperson = crossperson_est,
      crossperson_corr = crossperson_corr,

      components_vis = components_vis,
      s2_genetic_vis = s2_gen_vis,
      R2_model = R2_model,

      components_total_assumed = components_total_assumed,
      s2_genetic_total_assumed = s2_gen_total_assumed,
      R2_total_assumed = R2_total_assumed,

      hidden_correction_details = scale_estimates$hidden_correction_details,
      parental_sibs = sib_estimates,
      parental_sibs_total = scale_estimates$sib_total
    )
  }

  # ---------------------------------------------------------------------------
  # Fit main model
  # ---------------------------------------------------------------------------
  fit <- tryCatch(.fit_once(dat), error = function(e) {
    warning("SEM fitting failed: ", e$message)
    NULL
  })

  if (is.null(fit) || !lavaan::inspect(fit, "converged") ||
      (fit_constraints != "none" && !.ravel_constraint_fit_valid(fit))) {
    return(structure(list(
      converged = FALSE,
      N = nrow(dat),
      model = model,
      fit = fit,
      error = if (fit_constraints == "none") "Model did not converge" else if (fit_constraints == "estimated_state")
        "Estimated-state fit failed convergence, constraint, local-identification, or full-covariance checks; inference is unavailable" else
        "Constrained model did not produce an admissible converged fit satisfying its restrictions",
      outcome_variable = col_outcome
    ), class = c("ravel_results", "list")))
  }

  # Extract estimates
  estimates <- .extract(fit, data_for_r2 = dat)

  lavaan_inference <- .build_lavaan_uncertainty(
    fit = fit,
    estimates = estimates,
    h = h,
    am_correction = am_correction,
    am_n_generations = am_n_generations,
    ci_level = ci_level,
    sib_roles = sib_roles
  )

  # ---------------------------------------------------------------------------
  # Bootstrap (if requested)
  # ---------------------------------------------------------------------------
  boot_results <- NULL

  if (n_boot > 0) {
    if (cluster_var %in% names(dat)) {
      uids <- unique(dat[[cluster_var]])
    } else {
      uids <- seq_len(nrow(dat))
      cluster_var <- NULL
    }

    # Store visible, raw-total, and standardized-total uncertainty for all fits.
    boot_vis_bD <- boot_vis_bM <- boot_vis_bF <- boot_vis_rho <- rep(NA_real_, n_boot)
    boot_raw_bD <- boot_raw_bM <- boot_raw_bF <- rep(NA_real_, n_boot)
    boot_cov_Tmf <- boot_corr_Tmf <- rep(NA_real_, n_boot)
    boot_std_bD <- boot_std_bM <- boot_std_bF <- rep(NA_real_, n_boot)
    core_cols <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f", .sib_core_labels(sib_roles))
    boot_core <- matrix(NA_real_, n_boot, length(core_cols), dimnames = list(NULL, core_cols))
    sib_labels <- names(.parental_sib_coef_map(sib_roles))
    boot_sib <- boot_sib_std <- boot_sib_raw <-
      matrix(NA_real_, n_boot, length(sib_labels), dimnames = list(NULL, sib_labels))
    n_ok <- 0

    for (b in seq_len(n_boot)) {
      # Sample with replacement
      if (!is.null(cluster_var)) {
        samp_ids <- sample(uids, length(uids), replace = TRUE)
        idx <- unlist(lapply(samp_ids, function(id) which(dat[[cluster_var]] == id)))
        dat_b <- dat[idx, , drop = FALSE]
      } else {
        dat_b <- dat[sample.int(nrow(dat), nrow(dat), replace = TRUE), , drop = FALSE]
      }

      fit_b <- tryCatch(.fit_once(dat_b), error = function(e) NULL)
      if (is.null(fit_b) || !lavaan::inspect(fit_b, "converged")) next
      if (fit_constraints != "none" && !.ravel_constraint_fit_valid(fit_b)) next

      ext_b <- tryCatch(.extract(fit_b, data_for_r2 = dat_b), error = function(e) NULL)
      if (is.null(ext_b)) next

      boot_vis_bD[b] <- ext_b$beta_vis["bD"]
      boot_vis_bM[b] <- ext_b$beta_vis["bm"]
      boot_vis_bF[b] <- ext_b$beta_vis["bf"]
      boot_vis_rho[b] <- ext_b$rho_vis
      boot_core[b, ] <- c(ext_b$beta_vis, ext_b$latent_cov_vis,
        .sib_theta(ext_b$parental_sibs, sib_roles))[colnames(boot_core)]
      if (length(sib_labels)) {
        boot_sib[b, ] <- ext_b$parental_sibs$coefficients[sib_labels]
        boot_sib_std[b, ] <- ext_b$parental_sibs_total$std[sib_labels]
        boot_sib_raw[b, ] <- ext_b$parental_sibs_total$raw[sib_labels]
      }

      boot_raw_bD[b] <- ext_b$beta_total_assumed["bD"]
      boot_raw_bM[b] <- ext_b$beta_total_assumed["bm"]
      boot_raw_bF[b] <- ext_b$beta_total_assumed["bf"]
      boot_cov_Tmf[b] <- ext_b$cov_total_assumed
      boot_corr_Tmf[b] <- ext_b$rho_total_assumed
      boot_std_bD[b] <- ext_b$beta_total_assumed_std["bD"]
      boot_std_bM[b] <- ext_b$beta_total_assumed_std["bm"]
      boot_std_bF[b] <- ext_b$beta_total_assumed_std["bf"]
      n_ok <- n_ok + 1
    }

    if (n_ok >= 10) {
      visible_block <- .summarize_sample_matrix(
        cbind(
          bD = boot_vis_bD,
          bm = boot_vis_bM,
          bf = boot_vis_bF,
          rhoMF = boot_vis_rho
        ),
        ci_level = ci_level
      )

      total_raw_block <- .summarize_sample_matrix(
        cbind(
          bD = boot_raw_bD,
          bm = boot_raw_bM,
          bf = boot_raw_bF,
          rhoMF = boot_cov_Tmf
        ),
        ci_level = ci_level
      )

      total_std_block <- .summarize_sample_matrix(
        cbind(
          bD = boot_std_bD,
          bm = boot_std_bM,
          bf = boot_std_bF,
          rhoMF = boot_corr_Tmf
        ),
        ci_level = ci_level
      )

      boot_results <- list(
        source = "bootstrap",
        n_successful = n_ok,
        core_draws = boot_core,
        minimum_success = max(10L, ceiling(.8 * n_boot)),
        ci_level = ci_level,
        visible = visible_block,
        total_raw = total_raw_block,
        total_std = total_std_block,
        parental_sibs = if (length(sib_labels))
          .summarize_sample_matrix(boot_sib, ci_level = ci_level) else NULL,
        parental_sibs_total_std = if (length(sib_labels))
          .summarize_sample_matrix(boot_sib_std, ci_level = ci_level) else NULL,
        parental_sibs_total_raw = if (length(sib_labels))
          .summarize_sample_matrix(boot_sib_raw, ci_level = ci_level) else NULL,
        visible_bD_se = visible_block$se["bD"],
        visible_bD_ci = visible_block$ci[, "bD"],
        visible_bM_se = visible_block$se["bm"],
        visible_bM_ci = visible_block$ci[, "bm"],
        visible_bF_se = visible_block$se["bf"],
        visible_bF_ci = visible_block$ci[, "bf"],
        visible_rho_se = visible_block$se["rhoMF"],
        visible_rho_ci = visible_block$ci[, "rhoMF"],
        bD_se = total_raw_block$se["bD"],
        bD_ci = total_raw_block$ci[, "bD"],
        bM_se = total_raw_block$se["bm"],
        bM_ci = total_raw_block$ci[, "bm"],
        bF_se = total_raw_block$se["bf"],
        bF_ci = total_raw_block$ci[, "bf"],
        total_std_bD_se = total_std_block$se["bD"],
        total_std_bD_ci = total_std_block$ci[, "bD"],
        total_std_bM_se = total_std_block$se["bm"],
        total_std_bM_ci = total_std_block$ci[, "bm"],
        total_std_bF_se = total_std_block$se["bf"],
        total_std_bF_ci = total_std_block$ci[, "bf"],
        # Raw total-scale spouse covariance (kappa_T). At h = 0 its
        # standardized correlation, rather than its raw covariance, coincides
        # with the visible-scale spouse correlation.
        cov_Tmf_se  = total_raw_block$se["rhoMF"],
        cov_Tmf_ci  = total_raw_block$ci[, "rhoMF"],
        # Standardized total-scale spouse correlation (rho_T).
        corr_Tmf_se = total_std_block$se["rhoMF"],
        corr_Tmf_ci = total_std_block$ci[, "rhoMF"]
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Compute indicator reliability and determinacy diagnostics
  # ---------------------------------------------------------------------------
  PE <- lavaan::parameterEstimates(fit)

  get_load <- function(lv, ind) {
    x <- PE[PE$lhs == lv & PE$op == "=~" & PE$rhs == ind, "est"]
    if (length(x) == 0) NA_real_ else x[1]
  }
  get_theta <- function(v) {
    x <- PE[PE$lhs == v & PE$op == "~~" & PE$rhs == v, "est"]
    if (length(x) == 0) NA_real_ else x[1]
  }

  lm1 <- get_load("Gm", "PGS1_mother"); lm2 <- get_load("Gm", "PGS2_mother")
  lf1 <- get_load("Gf", "PGS1_father"); lf2 <- get_load("Gf", "PGS2_father")
  lc1 <- get_load("Gc", "PGS1_child");  lc2 <- get_load("Gc", "PGS2_child")

  reliability <- data.frame(
    role = rep(c("mother", "father", "child"), each = 2),
    PGS = rep(c("PGS1", "PGS2"), 3),
    reliability = c(
      lm1^2 / (lm1^2 + get_theta("PGS1_mother")),
      lm2^2 / (lm2^2 + get_theta("PGS2_mother")),
      lf1^2 / (lf1^2 + get_theta("PGS1_father")),
      lf2^2 / (lf2^2 + get_theta("PGS2_father")),
      lc1^2 / (lc1^2 + get_theta("PGS1_child")),
      lc2^2 / (lc2^2 + get_theta("PGS2_child"))
    ),
    stringsAsFactors = FALSE
  )

  determinacy <- .compute_determinacy(PE)

  within_person_corr <- c(
    mother = cor(dat$PGS1_mother, dat$PGS2_mother, use = "pairwise.complete.obs"),
    father = cor(dat$PGS1_father, dat$PGS2_father, use = "pairwise.complete.obs"),
    child  = cor(dat$PGS1_child,  dat$PGS2_child,  use = "pairwise.complete.obs")
  )

  sib_estimates <- estimates$parental_sibs
  if (!is.null(sib_estimates)) {
    reliability <- rbind(reliability, sib_estimates$reliability)
    determinacy <- c(determinacy, sib_estimates$determinacy)
    sib_name <- c(ms = "mat_sib", fs = "pat_sib")[sib_roles]
    within_person_corr <- c(within_person_corr, setNames(vapply(sib_name, function(nm)
      cor(dat[[paste0("PGS1_", nm)]], dat[[paste0("PGS2_", nm)]], use = "pairwise.complete.obs"),
      numeric(1)), sib_name))
  }

  vanilla_regression <- .fit_vanilla_pgs_regression(
    dat = dat,
    outcome = outcome,
    fixed_effects = fixed_effects,
    ci_level = ci_level
  )

  # Ratios of spousal residual covariance to latent rho (descriptive; no constraint)
  rho_vis <- estimates$latent_cov_vis["rhoMF"]
  crossperson_ratio <- NULL
  if (!is.na(rho_vis) && abs(rho_vis) > 1e-12) {
    crossperson_ratio <- c(
      e1_mf_over_rho = unname(estimates$crossperson["e1_mf"]) / unname(rho_vis),
      e2_mf_over_rho = unname(estimates$crossperson["e2_mf"]) / unname(rho_vis)
    )
  }

  # ---------------------------------------------------------------------------
  # Return results
  # ---------------------------------------------------------------------------
  result <- structure(list(
    converged = TRUE,
    N = if (length(sib_roles)) as.integer(lavaan::lavInspect(fit, "nobs")) else nrow(dat),
    fit = fit,
    model = model,

    # Effect estimates
    beta_vis = estimates$beta_vis,
    beta_h_only_assumed = estimates$beta_h_only_assumed,
    beta_total_assumed = estimates$beta_total_assumed,
    beta_total_assumed_std = estimates$beta_total_assumed_std,

    # Latent structure
    rho_vis = estimates$rho_vis,
    rho_total_assumed = estimates$rho_total_assumed,
    cov_total_assumed = estimates$cov_total_assumed,
    latent_cov_vis = estimates$latent_cov_vis,
    latent_cov_total_assumed = estimates$latent_cov_total_assumed,
    latent_corr_total_assumed = estimates$latent_corr_total_assumed,
    var_total_assumed = estimates$var_total_assumed,
    crossperson_resid = estimates$crossperson,
    crossperson_resid_corr = estimates$crossperson_corr,
    crossperson_ratio = crossperson_ratio,

    # Variance decomposition
    components_vis = estimates$components_vis,
    s2_genetic_vis = estimates$s2_genetic_vis,
    R2_model = estimates$R2_model,

    components_total_assumed = estimates$components_total_assumed,
    s2_genetic_total_assumed = estimates$s2_genetic_total_assumed,
    R2_total_assumed = estimates$R2_total_assumed,

    # Diagnostics
    reliability = reliability,
    determinacy = determinacy,
    within_person_corr = within_person_corr,
    vanilla_regression = vanilla_regression,

    # Fixed effects
    fixed_effects = fixed_effects,
    fixed_crossperson = {
      fixed_used <- fixed_crossperson[!is.na(fixed_crossperson)]
      if (length(fixed_used) == 0) NULL else fixed_used
    },
    fixed_crossperson_corr = {
      fixed_used <- fixed_crossperson_corr[!is.na(fixed_crossperson_corr)]
      if (length(fixed_used) == 0) NULL else fixed_used
    },

    # Bootstrap
    lavaan_inference = lavaan_inference,
    bootstrap = boot_results,

    # Opt-in constraint metadata and finite-generation sensitivity refit inputs.
    # The latter is retained only when fitting itself depends on external h.
    constraint_spec = list(mode = fit_constraints, h = h,
      n_generations = if (fit_constraints == "n_gen") am_n_generations else NA_integer_,
      depends_on_pi = fit_constraints == "n_gen" && am_n_generations > 1L),
    constraint_refit = if (fit_constraints == "n_gen" && am_n_generations > 1L) list(
      dat = dat, h = h, outcome = outcome, measurement = measurement,
      crossperson = crossperson, am_correction = am_correction,
      am_n_generations = am_n_generations, fix_parent_paths = fix_parent_paths,
      fixed_effects = fixed_effects, n_boot = 0, ci_level = ci_level,
      fit_constraints = fit_constraints, parental_sibs = parental_sibs) else NULL,

    # Settings
    h = h,
    am_n_generations = if (am_correction == "n_gen") am_n_generations else NA_integer_,
    measurement = measurement,
    crossperson_structure = crossperson,
    outcome = outcome,
    outcome_variable = col_outcome,
    am_correction = am_correction,

    # Internals for transparency
    hidden_correction_details = estimates$hidden_correction_details,

    # Parental-sibling extension (NULL when parental_sibs = "none")
    parental_sibs = if (is.null(sib_estimates)) NULL else c(
      list(mode = parental_sibs, roles = sib_roles, n_with_sib = n_with_sib),
      sib_estimates,
      list(lavaan_table = .parental_sib_lavaan_table(fit, sib_roles, ci_level),
           corrected = estimates$parental_sibs_total))
  ), class = c("ravel_results", "list"))
  if (!is.null(pi)) result$pi <- h
  if (!is.null(pi_se)) {
    result$pi_se <- pi_se
    if (identical(result$parental_sibs$corrected$invalid_reason,
        "incompatible_parental_sibling_covariance")) {
      # Preserve the fitted visible model when only the requested hidden state
      # is impossible. Invalid user covariance inputs still fail explicitly.
      input <- .ravel_pi_extract(result, "auto", covariate_mapping_confirmed)
      joint <- matrix(0, length(input$labels) + 1L, length(input$labels) + 1L,
        dimnames = list(c(input$labels, "pi"), c(input$labels, "pi")))
      joint[input$labels, input$labels] <- input$packet$vcov_core
      joint["pi", "pi"] <- pi_se^2
      if (!is.null(core_pi_cov))
        joint[input$labels, "pi"] <- joint["pi", input$labels] <- core_pi_cov[input$labels]
      .ravel_pi_validate_covariance(joint, c(input$labels, "pi"), "Joint core/pi covariance")
      result["pi_correction"] <- list(NULL)
      result$pi_correction_error <- paste("External-pi corrected inference is unavailable because",
        "the parental-sibling moments are incompatible with the supplied pi/assortment state.")
      warning(result$pi_correction_error, " The visible fit is retained.")
    } else {
      result$pi_correction <- ravel_correct_pi(result, pi = h, pi_se = pi_se,
        ci_level = ci_level, core_pi_cov = core_pi_cov,
        covariate_mapping_confirmed = covariate_mapping_confirmed)
    }
  }
  result
}

# =============================================================================
# OUTPUT FUNCTIONS
# =============================================================================

.ravel_print_assumptions <- function(result, digits) {
  cat("Fit constraints:", result$constraint_spec$mode %||% "none", "\n")
  cat("AM correction:", result$am_correction %||% "unspecified")
  if (identical(result$am_correction, "n_gen"))
    cat("; generations =", result$am_n_generations)
  if (identical(result$am_correction, "estimated_state")) cat("; generation count not specified")
  cat("\n")
  if (!is.null(result$pi_correction)) {
    meta <- attr(result$pi_correction, "pi_uncertainty")
    cat(sprintf("External pi = %.*f; SE = %.*f; stored corrected inference includes pi uncertainty.\n",
      digits, meta$pi, digits, meta$raw_scale_se))
    cat("Stored visible inference is conditional on pi.\n")
  }
  if (!is.null(result$pi_correction_error))
    cat("Corrected inference unavailable:", result$pi_correction_error, "\n")
  state <- result$hidden_correction_details$estimated_state
  if (!is.null(state)) {
    cat("Estimated variance state (founder reference):\n")
    for (nm in intersect(c("B", "vGp", "vTp", "rT"), names(state)))
      cat(sprintf("  %s = %.*f\n", if (nm == "B") "growth (B)" else nm, digits, state[[nm]]))
    if (!is.null(state$history_status)) cat("  History diagnostic:", state$history_status, "\n")
  }
  invisible(NULL)
}

.ravel_print_correction_table <- function(result, digits) {
  if (!all(c("scale", "parameter", "estimate", "se", "ci_low", "ci_high", "am_correction") %in% names(result)))
    stop("Supply a fitted RAVEL result or a table from ravel_correct_packet, ravel_correct_pi, or ravel_correct_grid.")
  cat("RAVEL corrected estimates\n")
  cat("AM correction:", paste(unique(result$am_correction), collapse = ", "), "\n")
  columns <- intersect(c("outcome_variable", "residual_model", "scale", "parameter", "estimate", "se",
    "ci_low", "ci_high", "scenario_id", "pi", "pi_se", "h", "hidden_fraction", "am_correction",
    "am_n_generations", "n_generations", "inference", "status", "correction_status", "error"), names(result))
  print.data.frame(result[columns], digits = digits, row.names = FALSE)
  state <- attr(result, "estimated_state")
  if (!is.null(state$summary)) {
    cat("\nEstimated variance state (same uncertainty propagation):\n")
    print.data.frame(state$summary, digits = digits, row.names = FALSE)
  } else if (!is.null(state)) {
    cat("\nEstimated variance state (point estimates):\n")
    for (nm in intersect(c("B", "vGp", "vTp", "rT"), names(state)))
      cat(sprintf("  %s = %.*f\n", if (nm == "B") "growth (B)" else nm, digits, state[[nm]]))
    if (!is.null(state$history_status)) cat("  History diagnostic:", state$history_status, "\n")
  }
  invisible(result)
}

#' Fetch a named CI vector from an uncertainty block
#' @noRd
.uncertainty_ci <- function(uncertainty, block_name, param_name) {
  if (is.null(uncertainty) || is.null(uncertainty[[block_name]])) {
    return(NULL)
  }
  block <- uncertainty[[block_name]]
  if (is.null(block$ci) || !(param_name %in% colnames(block$ci))) {
    return(NULL)
  }
  unname(block$ci[, param_name])
}

#' Human-readable uncertainty source label
#' @noRd
.uncertainty_label <- function(uncertainty, block = NULL) {
  if (is.null(uncertainty)) return(NULL)
  visible <- if (identical(uncertainty$source, "bootstrap")) {
    sprintf("bootstrap (%d reps)", uncertainty$n_successful)
  } else if (identical(uncertainty$source, "lavaan")) "lavaan/delta-method" else "uncertainty"
  if (!is.null(uncertainty$pi_correction_error)) {
    if (identical(block, "visible")) return(visible)
    if (!is.null(block)) return("unavailable (incompatible pi/assortment state)")
    return(paste0(visible, "; corrected inference unavailable"))
  }
  if (is.null(uncertainty$pi_correction)) return(visible)
  corrected <- if (identical(uncertainty$source, "bootstrap"))
    "bootstrap covariance + pi (delta-Wald)" else "lavaan + pi (delta-Wald)"
  if (identical(block, "visible")) return(paste0(visible, "; conditional on pi"))
  if (!is.null(block)) return(corrected)
  paste0("visible: ", visible, "; corrected: ", corrected)
}

#' Wald-style two-sided p-values from estimates and SEs
#' @noRd
.uncertainty_p <- function(est, se) {
  ifelse(is.na(est) | is.na(se) | se <= 0,
         NA_real_,
         2 * stats::pnorm(abs(est / se), lower.tail = FALSE))
}

#' Convert one uncertainty block to a tidy data.frame
#' @noRd
.results_block_to_df <- function(block, scale_label, rho_label, source_label, p_value_method,
                                 regression = "RAVEL", family = NA_character_) {
  if (is.null(block)) {
    return(NULL)
  }

  estimate <- c(
    bD = unname(block$estimate["bD"]),
    bm = unname(block$estimate["bm"]),
    bf = unname(block$estimate["bf"]),
    rhoMF = unname(block$estimate["rhoMF"])
  )
  se <- c(
    bD = unname(block$se["bD"]),
    bm = unname(block$se["bm"]),
    bf = unname(block$se["bf"]),
    rhoMF = unname(block$se["rhoMF"])
  )
  ci_low <- c(
    bD = unname(block$ci["lower", "bD"]),
    bm = unname(block$ci["lower", "bm"]),
    bf = unname(block$ci["lower", "bf"]),
    rhoMF = unname(block$ci["lower", "rhoMF"])
  )
  ci_high <- c(
    bD = unname(block$ci["upper", "bD"]),
    bm = unname(block$ci["upper", "bm"]),
    bf = unname(block$ci["upper", "bf"]),
    rhoMF = unname(block$ci["upper", "rhoMF"])
  )

  data.frame(
    regression = regression,
    family = family,
    scale = scale_label,
    parameter = c("beta_D", "beta_M", "beta_F", rho_label),
    estimate = unname(estimate),
    se = unname(se),
    ci_low = unname(ci_low),
    ci_high = unname(ci_high),
    p_value = unname(.uncertainty_p(estimate, se)),
    uncertainty_source = source_label,
    p_value_method = p_value_method,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Fit vanilla trio regressions on observed standardized PGS
#' @noRd
.fit_vanilla_pgs_regression <- function(dat, outcome, fixed_effects, ci_level) {
  z <- stats::qnorm(1 - (1 - ci_level) / 2)

  fit_one <- function(version) {
    cols <- c(
      child = paste0(version, "_child"),
      mother = paste0(version, "_mother"),
      father = paste0(version, "_father")
    )

    reg_dat <- data.frame(
      .Y = if (outcome == "continuous") as.numeric(scale(dat$Y)) else dat$Y,
      .pgs_child = as.numeric(scale(dat[[cols["child"]]])),
      .pgs_mother = as.numeric(scale(dat[[cols["mother"]]])),
      .pgs_father = as.numeric(scale(dat[[cols["father"]]])),
      stringsAsFactors = FALSE
    )
    if (length(fixed_effects) > 0) {
      reg_dat <- cbind(reg_dat, dat[fixed_effects], stringsAsFactors = FALSE)
    }

    rhs <- c(".pgs_child", ".pgs_mother", ".pgs_father", fixed_effects)
    form <- stats::reformulate(rhs, response = ".Y")

    fit <- tryCatch(
      if (outcome == "continuous") {
        stats::lm(form, data = reg_dat)
      } else {
        stats::glm(form, data = reg_dat, family = stats::binomial(link = "probit"))
      },
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(list(
        converged = FALSE,
        model = NULL,
        fit = NULL,
        N = NA_integer_,
        coefficients = c(child = NA_real_, mother = NA_real_, father = NA_real_),
        se = c(child = NA_real_, mother = NA_real_, father = NA_real_),
        ci = matrix(NA_real_, nrow = 2, ncol = 3,
                    dimnames = list(c("lower", "upper"), c("child", "mother", "father"))),
        p_value = c(child = NA_real_, mother = NA_real_, father = NA_real_),
        table = data.frame(
          regression = version,
          parameter = c("beta_child", "beta_mother", "beta_father"),
          estimate = NA_real_,
          se = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_,
          p_value = NA_real_,
          family = if (outcome == "continuous") "lm_scaled_outcome" else "glm_probit_binary",
          stringsAsFactors = FALSE
        ),
        error = fit$message
      ))
    }

    coef_tab <- summary(fit)$coefficients
    target_rows <- c(".pgs_child", ".pgs_mother", ".pgs_father")
    param_names <- c("child", "mother", "father")

    extract_col <- function(j) {
      out <- sapply(target_rows, function(r) {
        if (r %in% rownames(coef_tab)) coef_tab[r, j] else NA_real_
      })
      names(out) <- param_names
      out
    }

    est <- extract_col(1)
    se <- extract_col(2)
    p_val <- extract_col(4)
    ci <- rbind(
      lower = est - z * se,
      upper = est + z * se
    )

    table <- data.frame(
      regression = version,
      parameter = c("beta_child", "beta_mother", "beta_father"),
      estimate = unname(est),
      se = unname(se),
      ci_low = unname(ci["lower", ]),
      ci_high = unname(ci["upper", ]),
      p_value = unname(p_val),
      family = if (outcome == "continuous") "lm_scaled_outcome" else "glm_probit_binary",
      stringsAsFactors = FALSE
    )

    list(
      converged = TRUE,
      model = form,
      fit = fit,
      N = stats::nobs(fit),
      coefficients = est,
      se = se,
      ci = ci,
      p_value = p_val,
      table = table,
      error = NULL
    )
  }

  pgs1 <- fit_one("PGS1")
  pgs2 <- fit_one("PGS2")

  list(
    outcome = outcome,
    ci_level = ci_level,
    fixed_effects = fixed_effects,
    PGS1 = pgs1,
    PGS2 = pgs2,
    table = rbind(pgs1$table, pgs2$table)
  )
}

#' Convert a fitted RAVEL result to a tidy data frame
#'
#' @description
#' Returns one row per reported parameter for the visible scale and, when
#' \code{h > 0}, the corrected standardized total-liability scale. Each row
#' contains the point estimate, standard error, confidence interval, and
#' two-sided p-value.
#'
#' @details
#' The returned table always includes the visible-scale rows
#' \code{beta_D}, \code{beta_M}, \code{beta_F}, and \code{rho_visible}. If
#' \code{h > 0}, it additionally includes the corrected standardized total-scale
#' rows \code{beta_D}, \code{beta_M}, \code{beta_F}, and
#' \code{rho_corrected}. Setting \code{include_raw = TRUE} also appends the raw
#' total-liability rows with \code{cov_corrected} in place of the spouse
#' correlation.
#'
#' Parental-sibling fits also include \code{beta_MS}/\code{delta_M} for a
#' maternal sibling and \code{beta_FS}/\code{delta_F} for a paternal sibling
#' on each requested scale. Each contrast subtracts the sibling coefficient
#' from the corresponding parent coefficient on the same liability scale;
#' it is not an effect per standard deviation of the parent-sibling difference.
#'
#' When \code{uncertainty = "lavaan"}, p-values are standard Wald
#' normal-approximation p-values based on the lavaan/delta-method SEs. When
#' \code{uncertainty = "bootstrap"}, the confidence intervals are empirical
#' bootstrap intervals, while p-values are normal-approximation p-values based
#' on the bootstrap SEs. When \code{pi_se} was supplied to \code{ravel_fit},
#' corrected rows instead use joint delta-method Wald intervals, including pi
#' uncertainty, even when the conditional core covariance comes from bootstrap
#' draws. Visible rows keep their conditional intervals.
#'
#' @param result A fitted RAVEL result object from \code{\link{ravel_fit}}.
#' @param uncertainty Character. One of \code{"auto"}, \code{"bootstrap"},
#'   \code{"lavaan"}, or \code{"none"}.
#' @param include_raw Logical. If \code{TRUE}, also return the
#'   corrected raw total-liability rows, including at \code{h = 0}.
#' @param include_vanilla Logical. If \code{TRUE}, append the plain observed-PGS
#'   trio regression rows returned by \code{ravel_fit()}.
#'
#' @return A tidy \code{data.frame} with columns \code{regression},
#'   \code{family}, \code{outcome_variable}, \code{scale}, \code{parameter},
#'   \code{estimate}, \code{se}, \code{ci_low}, \code{ci_high},
#'   \code{p_value}, \code{uncertainty_source}, and \code{p_value_method}.
#'
#' @export
ravel_to_dataframe <- function(result,
                               uncertainty = c("auto", "bootstrap", "lavaan", "none"),
                               include_raw = FALSE,
                               include_vanilla = TRUE) {
  if (!isTRUE(result$converged)) {
    stop("Cannot convert a non-converged RAVEL fit to a parameter table.")
  }

  uncertainty <- match.arg(uncertainty)
  unc <- .resolve_uncertainty(result, uncertainty)

  if (is.null(unc)) {
    stop("Requested uncertainty results are not available in this fit object.")
  }

  source_label <- .uncertainty_label(unc, "visible")
  corrected_source_label <- .uncertainty_label(unc, "total_std")
  p_value_method <- if (identical(unc$source, "bootstrap")) {
    "normal_approx_from_bootstrap_se"
  } else {
    "normal_approx_from_lavaan_se"
  }

  corrected_p_value_method <- if (is.null(unc$pi_correction)) p_value_method else
    "normal_approx_from_joint_pi_delta_se"

  out <- .results_block_to_df(
    unc$visible, "visible", "rho_visible", source_label, p_value_method
  )

  sib_tbl <- .parental_sib_table(result, unc)
  if (!is.null(sib_tbl)) {
    out <- rbind(out, data.frame(
      regression = "RAVEL", family = NA_character_, scale = "visible",
      parameter = sib_tbl$parameter, estimate = sib_tbl$estimate, se = sib_tbl$se,
      ci_low = sib_tbl$ci_low, ci_high = sib_tbl$ci_high, p_value = sib_tbl$p_value,
      uncertainty_source = source_label, p_value_method = p_value_method,
      stringsAsFactors = FALSE))
  }

  sib_rows <- function(scale) {
    tbl <- .parental_sib_table(result, unc, scale)
    if (is.null(tbl)) return(NULL)
    data.frame(regression = "RAVEL", family = NA_character_, scale = scale,
      parameter = tbl$parameter, estimate = tbl$estimate, se = tbl$se,
      ci_low = tbl$ci_low, ci_high = tbl$ci_high, p_value = tbl$p_value,
      uncertainty_source = corrected_source_label, p_value_method = corrected_p_value_method,
      stringsAsFactors = FALSE)
  }

  if (result$h > 0) {
    out <- rbind(
      out,
      .results_block_to_df(
        unc$total_std, "corrected_std", "rho_corrected",
        corrected_source_label, corrected_p_value_method
      ),
      sib_rows("corrected_std")
    )
  }

  if (isTRUE(include_raw)) {
    out <- rbind(
      out,
      .results_block_to_df(unc$total_raw, "corrected_raw", "cov_corrected", corrected_source_label, corrected_p_value_method),
      sib_rows("corrected_raw")
    )
  }

  if (isTRUE(include_vanilla) && !is.null(result$vanilla_regression) && !is.null(result$vanilla_regression$table)) {
    vanilla_tbl <- result$vanilla_regression$table
    vanilla_tbl$scale <- "vanilla"
    vanilla_tbl$uncertainty_source <- "model_summary"
    vanilla_tbl$p_value_method <- ifelse(
      vanilla_tbl$family == "glm_probit_binary",
      "glm_summary_z_test",
      "lm_summary_t_test"
    )
    vanilla_tbl <- vanilla_tbl[, c(
      "regression", "family", "scale", "parameter", "estimate", "se",
      "ci_low", "ci_high", "p_value", "uncertainty_source", "p_value_method"
    )]
    out <- rbind(out, vanilla_tbl)
  }

  out$outcome_variable <- if (!is.null(result$outcome_variable)) result$outcome_variable else "Y"
  out <- out[, c(
    "regression", "family", "outcome_variable", "scale", "parameter",
    "estimate", "se", "ci_low", "ci_high", "p_value",
    "uncertainty_source", "p_value_method"
  )]

  attr(out, "settings") <- list(fit_constraints = result$constraint_spec$mode %||% "none",
    am_correction = result$am_correction, am_n_generations = result$am_n_generations,
    h = result$h, ci_level = unc$ci_level, outcome_type = result$outcome,
    parental_sibs = result$parental_sibs$mode %||% "none")
  attr(out, "constraint_spec") <- result$constraint_spec
  attr(out, "estimated_state") <- result$hidden_correction_details$estimated_state
  attr(out, "ci_level") <- unc$ci_level
  if (!is.null(result$pi_correction_error)) {
    attr(out, "pi_correction_error") <- result$pi_correction_error
    attr(out, "settings")$pi <- result$h
    attr(out, "settings")$pi_se <- result$pi_se
    unavailable <- out$scale %in% c("corrected_std", "corrected_raw")
    out$uncertainty_source[unavailable] <- "unavailable (incompatible pi/assortment state)"
    out$p_value_method[unavailable] <- "unavailable"
  }
  if (!is.null(unc$pi_correction)) {
    meta <- attr(unc$pi_correction, "pi_uncertainty")
    out$inference <- ifelse(out$scale %in% c("corrected_std", "corrected_raw"),
      unique(unc$pi_correction$inference),
      ifelse(out$scale == "vanilla", "model_summary", if (identical(unc$source, "bootstrap"))
        "bootstrap_percentile_CI_normal_p" else "lavaan_delta_wald"))
    out$pi <- result$h
    out$pi_se <- ifelse(out$scale %in% c("corrected_std", "corrected_raw"), meta$raw_scale_se, NA_real_)
    attr(out, "settings")$pi <- result$h
    attr(out, "settings")$pi_se <- meta$raw_scale_se
    attr(out, "pi_uncertainty") <- meta
    attr(out, "pi_correction") <- unc$pi_correction
    attr(out, "estimated_state") <- attr(unc$pi_correction, "estimated_state") %||%
      result$hidden_correction_details$estimated_state
  }
  out
}

#' Coerce fitted RAVEL results to a data frame
#'
#' @param x A fitted RAVEL result object from \code{\link{ravel_fit}}.
#' @param row.names Optional row names.
#' @param optional Ignored. Included for compatibility with the generic.
#' @param uncertainty Character. One of \code{"auto"}, \code{"bootstrap"},
#'   \code{"lavaan"}, or \code{"none"}.
#' @param include_raw Logical. If \code{TRUE}, also return the
#'   corrected raw total-liability rows, including at \code{h = 0}.
#' @param include_vanilla Logical. If \code{TRUE}, append the plain observed-PGS
#'   trio regression rows returned by \code{ravel_fit()}.
#' @param ... Unused.
#'
#' @return A tidy \code{data.frame} of fitted coefficients and uncertainty.
#'
#' @method as.data.frame ravel_results
#' @export
as.data.frame.ravel_results <- function(x, row.names = NULL, optional = FALSE,
                                        uncertainty = c("auto", "bootstrap", "lavaan", "none"),
                                        include_raw = FALSE, include_vanilla = TRUE, ...) {
  out <- ravel_to_dataframe(
    result = x,
    uncertainty = uncertainty,
    include_raw = include_raw,
    include_vanilla = include_vanilla
  )

  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }

  out
}

#' Print the parental-sibling coefficients (and, in full mode, their structure)
#' @noRd
.print_parental_sibs <- function(result, unc, digits, compact = TRUE,
                                 scale = c("visible", "corrected_std")) {
  scale <- match.arg(scale)
  sib <- result$parental_sibs
  if (is.null(sib)) return(invisible(NULL))
  tbl <- .parental_sib_table(result, unc, scale)
  desc <- c(beta_MS = "mother's sibling (lineage)", delta_M = "mother-specific (bm - bMS)",
            beta_FS = "father's sibling (lineage)", delta_F = "father-specific (bf - bFS)")
  if (scale == "visible") {
    cat(if (compact) "\n" else "", "Parental siblings (visible scale; families with sibling: ",
        paste(sprintf("%s = %d", names(sib$n_with_sib), as.integer(sib$n_with_sib)), collapse = ", "),
        "):\n", sep = "")
  } else {
    cat("  Parental siblings (standardized total liability):\n")
  }
  indent <- if (scale == "visible") "  " else "    "
  for (i in seq_len(nrow(tbl))) {
    ci <- if (is.na(tbl$se[i])) "" else
      sprintf(" [%.*f, %.*f]", digits, tbl$ci_low[i], digits, tbl$ci_high[i])
    cat(sprintf("%s%-8s %.*f%s  %s\n", indent, tbl$parameter[i], digits, tbl$estimate[i], ci,
                desc[[tbl$parameter[i]]]))
  }
  if (compact || scale != "visible") return(invisible(NULL))
  labels <- c(kap_m = "kappa_m (Gm, G_mat_sib)", chi_ms = "chi (Gf, G_mat_sib)",
              eta_ms = "eta (Gc, G_mat_sib)", kap_f = "kappa_f (Gf, G_pat_sib)",
              chi_fs = "chi (Gm, G_pat_sib)", eta_fs = "eta (Gc, G_pat_sib)",
              zeta = "zeta (G_mat_sib, G_pat_sib)")
  cat("  Latent covariances (visible):\n")
  for (nm in names(sib$latent_cov))
    cat(sprintf("    %-28s %.*f\n", labels[[nm]], digits, sib$latent_cov[[nm]]))
  cat("  Within-score residual correlations (benchmarks: own parent 1/2, child 1/4, others 0):\n")
  for (nm in names(sib$resid_corr))
    cat(sprintf("    %-10s %.*f\n", nm, digits, sib$resid_corr[[nm]]))
  cat("\n")
  invisible(NULL)
}

#' Print RAVEL Results
#'
#' @param result Output from \code{\link{ravel_fit}}, or a corrected coefficient table.
#' @param gold Optional named vector with gold standard values for comparison
#' @param digits Number of decimal places
#' @param uncertainty Character. One of \code{"auto"}, \code{"bootstrap"},
#'   \code{"lavaan"}, or \code{"none"}.
#'
#' @export
#' @importFrom stats reshape
print_ravel_results <- function(result, gold = NULL, digits = 4,
                                uncertainty = c("auto", "bootstrap", "lavaan", "none")) {

  if (is.data.frame(result)) {
    if (match.arg(uncertainty) != "auto")
      stop("A correction table already contains its chosen uncertainty; recompute it to change the source.")
    return(.ravel_print_correction_table(result, digits))
  }
  if (!result$converged) {
    cat(result$error %||% "Model did not converge.", "\n")
    return(invisible(result))
  }

  uncertainty <- match.arg(uncertainty)
  unc <- .resolve_uncertainty(result, uncertainty)

  cat("================================================================================\n")
  cat("RAVEL RESULTS\n")
  cat("================================================================================\n\n")

  cat("Model specification:\n")
  cat("  Measurement:", result$measurement, "\n")
  cat("  Cross-person residuals:", result$crossperson_structure, "\n")
  if (!is.null(result$fixed_crossperson)) {
    fixed_text <- paste(
      sprintf("%s=%s", names(result$fixed_crossperson), result$fixed_crossperson),
      collapse = ", "
    )
    cat("  Fixed cross-person residual covariances:", fixed_text, "\n")
  }
  if (!is.null(result$fixed_crossperson_corr)) {
    fixed_text <- paste(
      sprintf("%s=%s", names(result$fixed_crossperson_corr), result$fixed_crossperson_corr),
      collapse = ", "
    )
    cat("  Fixed cross-person residual correlations:", fixed_text, "\n")
  }
  cat("  Outcome:", result$outcome, "\n")
  cat("  Hidden variance fraction (h, assumed):", result$h, "\n")
  .ravel_print_assumptions(result, digits)
  cat("\n")

  # Effect estimates
  cat("Effect estimates (visible latent genetics):\n")
  cat(sprintf("  bD (direct):   %.*f\n", digits, result$beta_vis["bD"]))
  cat(sprintf("  bm (maternal): %.*f\n", digits, result$beta_vis["bm"]))
  cat(sprintf("  bf (paternal): %.*f\n", digits, result$beta_vis["bf"]))
  if (!is.null(result$parental_sibs))
    cat("  (bm and bf are conditional on the parental siblings' liabilities)\n")
  cat("\n")
  .print_parental_sibs(result, unc, digits, compact = FALSE)

  if (result$h > 0 || (result$constraint_spec$mode %||% "none") != "none") {
    cat("Assumption-based post-estimation corrections:\n")
    cat("  (1) h-only scaling (assumes no cross-component tagging):\n")
    cat(sprintf("      beta_D: %.*f\n", digits, result$beta_h_only_assumed["bD"]))
    cat(sprintf("      beta_M: %.*f\n", digits, result$beta_h_only_assumed["bm"]))
    cat(sprintf("      beta_F: %.*f\n", digits, result$beta_h_only_assumed["bf"]))

    cat("  (2) transmission-calibrated matrix correction (h + assortative tagging):\n")
    cat(sprintf("      beta_D: %.*f\n", digits, result$beta_total_assumed["bD"]))
    cat(sprintf("      beta_M: %.*f\n", digits, result$beta_total_assumed["bm"]))
    cat(sprintf("      beta_F: %.*f\n", digits, result$beta_total_assumed["bf"]))
    cat("      standardized total-liability coefficients:\n")
    cat(sprintf("        beta_D_std: %.*f\n", digits, result$beta_total_assumed_std["bD"]))
    cat(sprintf("        beta_M_std: %.*f\n", digits, result$beta_total_assumed_std["bm"]))
    cat(sprintf("        beta_F_std: %.*f\n\n", digits, result$beta_total_assumed_std["bf"]))
  }

  # Uncertainty summaries
  if (!is.null(unc)) {
    cat("Reported uncertainty: ", .uncertainty_label(unc), "\n", sep = "")
    cat("Visible-scale intervals:\n")
    cat(sprintf("  beta_D: [%.*f, %.*f]\n", digits, .uncertainty_ci(unc, "visible", "bD")[1],
                digits, .uncertainty_ci(unc, "visible", "bD")[2]))
    cat(sprintf("  beta_M: [%.*f, %.*f]\n", digits, .uncertainty_ci(unc, "visible", "bm")[1],
                digits, .uncertainty_ci(unc, "visible", "bm")[2]))
    cat(sprintf("  beta_F: [%.*f, %.*f]\n", digits, .uncertainty_ci(unc, "visible", "bf")[1],
                digits, .uncertainty_ci(unc, "visible", "bf")[2]))
    cat(sprintf("  rho_MF: [%.*f, %.*f]\n", digits, .uncertainty_ci(unc, "visible", "rhoMF")[1],
                digits, .uncertainty_ci(unc, "visible", "rhoMF")[2]))
    if (result$h > 0 || (result$constraint_spec$mode %||% "none") != "none") {
      cat("Corrected raw total-liability intervals:\n")
      cat(sprintf("  beta_D: [%.*f, %.*f]\n", digits, .uncertainty_ci(unc, "total_raw", "bD")[1],
                  digits, .uncertainty_ci(unc, "total_raw", "bD")[2]))
      cat(sprintf("  beta_M: [%.*f, %.*f]\n", digits, .uncertainty_ci(unc, "total_raw", "bm")[1],
                  digits, .uncertainty_ci(unc, "total_raw", "bm")[2]))
      cat(sprintf("  beta_F: [%.*f, %.*f]\n", digits, .uncertainty_ci(unc, "total_raw", "bf")[1],
                  digits, .uncertainty_ci(unc, "total_raw", "bf")[2]))
      cat(sprintf("  Cov(Tm*,Tf*): [%.*f, %.*f]\n", digits, .uncertainty_ci(unc, "total_raw", "rhoMF")[1],
                  digits, .uncertainty_ci(unc, "total_raw", "rhoMF")[2]))
      cat("Corrected standardized total-liability intervals:\n")
      cat(sprintf("  beta_D_std: [%.*f, %.*f]\n", digits, .uncertainty_ci(unc, "total_std", "bD")[1],
                  digits, .uncertainty_ci(unc, "total_std", "bD")[2]))
      cat(sprintf("  beta_M_std: [%.*f, %.*f]\n", digits, .uncertainty_ci(unc, "total_std", "bm")[1],
                  digits, .uncertainty_ci(unc, "total_std", "bm")[2]))
      cat(sprintf("  beta_F_std: [%.*f, %.*f]\n", digits, .uncertainty_ci(unc, "total_std", "bf")[1],
                  digits, .uncertainty_ci(unc, "total_std", "bf")[2]))
      cat(sprintf("  Corr(Tm,Tf): [%.*f, %.*f]\n\n", digits, .uncertainty_ci(unc, "total_std", "rhoMF")[1],
                  digits, .uncertainty_ci(unc, "total_std", "rhoMF")[2]))
      if (!is.null(result$parental_sibs)) {
        .print_parental_sibs(result, unc, digits, scale = "corrected_std")
        cat("\n")
      }
    } else {
      cat("\n")
    }
  }

  # Variance decomposition (visible)
  cat("Variance decomposition (visible genetic latent factors):\n")
  for (nm in names(result$components_vis)) {
    cat(sprintf("  %-20s: %.*f\n", nm, digits, result$components_vis[nm]))
  }
  cat(sprintf("  %-20s: %.*f\n", "Total genetic (s2)", digits, result$s2_genetic_vis))
  cat(sprintf("\nR2 (SEM / visible): %.*f\n\n", digits, result$R2_model))

  # Variance decomposition (assumed total)
  if (result$h > 0 || (result$constraint_spec$mode %||% "none") != "none") {
    cat("Variance decomposition (assumed total genetic liability):\n")
    for (nm in names(result$components_total_assumed)) {
      cat(sprintf("  %-20s: %.*f\n", nm, digits, result$components_total_assumed[nm]))
    }
    cat(sprintf("  %-20s: %.*f\n", "Total genetic (s2)", digits, result$s2_genetic_total_assumed))
    cat(sprintf("\nR2 (assumed total): %.*f\n\n", digits, result$R2_total_assumed))
  }

  # Latent covariances
  cat("Latent genetic covariances (visible):\n")
  cat(sprintf("  rho (Gm, Gf):   %.*f\n", digits, result$latent_cov_vis["rhoMF"]))
  cat(sprintf("  eta_m (Gc, Gm): %.*f\n", digits, result$latent_cov_vis["eta_m"]))
  cat(sprintf("  eta_f (Gc, Gf): %.*f\n\n", digits, result$latent_cov_vis["eta_f"]))

  if (result$h > 0 || (result$constraint_spec$mode %||% "none") != "none") {
    cat("Latent covariances (assumed total; raw covariances):\n")
    cat(sprintf("  Cov(Tm*, Tf*): %.*f\n", digits, result$latent_cov_total_assumed["rhoMF"]))
    cat(sprintf("  Cov(Tc*, Tm*): %.*f\n", digits, result$latent_cov_total_assumed["eta_m"]))
    cat(sprintf("  Cov(Tc*, Tf*): %.*f\n", digits, result$latent_cov_total_assumed["eta_f"]))
    cat("  Role-specific Var(T*):\n")
    cat(sprintf("    child:  %.*f\n", digits, result$var_total_assumed["child"]))
    cat(sprintf("    mother: %.*f\n", digits, result$var_total_assumed["mother"]))
    cat(sprintf("    father: %.*f\n", digits, result$var_total_assumed["father"]))
    cat("  Standardized total-scale correlations:\n")
    cat(sprintf("    Corr(Tm, Tf): %.*f\n", digits, result$latent_corr_total_assumed["rhoMF"]))
    cat(sprintf("    Corr(Tc, Tm): %.*f\n", digits, result$latent_corr_total_assumed["eta_m"]))
    cat(sprintf("    Corr(Tc, Tf): %.*f\n\n", digits, result$latent_corr_total_assumed["eta_f"]))
  }

  # Cross-person residuals
  cat("Cross-person residual covariances:\n")
  if (result$crossperson_structure == "free") {
    cat(sprintf("  e1_mf: %.*f, e1_mc: %.*f, e1_fc: %.*f\n",
                digits, result$crossperson_resid["e1_mf"],
                digits, result$crossperson_resid["e1_mc"],
                digits, result$crossperson_resid["e1_fc"]))
    cat(sprintf("  e2_mf: %.*f, e2_mc: %.*f, e2_fc: %.*f\n\n",
                digits, result$crossperson_resid["e2_mf"],
                digits, result$crossperson_resid["e2_mc"],
                digits, result$crossperson_resid["e2_fc"]))
  } else {
    cat(sprintf("  e1_mf: %.*f, e1_pc: %.*f\n",
                digits, result$crossperson_resid["e1_mf"],
                digits, result$crossperson_resid["e1_pc"]))
    cat(sprintf("  e2_mf: %.*f, e2_pc: %.*f\n\n",
                digits, result$crossperson_resid["e2_mf"],
                digits, result$crossperson_resid["e2_pc"]))
  }

  cat("Cross-person residual correlations:\n")
  if (result$crossperson_structure == "free") {
    cat(sprintf("  e1_mf: %.*f, e1_mc: %.*f, e1_fc: %.*f\n",
                digits, result$crossperson_resid_corr["e1_mf"],
                digits, result$crossperson_resid_corr["e1_mc"],
                digits, result$crossperson_resid_corr["e1_fc"]))
    cat(sprintf("  e2_mf: %.*f, e2_mc: %.*f, e2_fc: %.*f\n\n",
                digits, result$crossperson_resid_corr["e2_mf"],
                digits, result$crossperson_resid_corr["e2_mc"],
                digits, result$crossperson_resid_corr["e2_fc"]))
  } else {
    cat(sprintf("  e1_mf: %.*f, e1_pc(mother): %.*f, e1_pc(father): %.*f\n",
                digits, result$crossperson_resid_corr["e1_mf"],
                digits, result$crossperson_resid_corr["e1_pc_mother"],
                digits, result$crossperson_resid_corr["e1_pc_father"]))
    cat(sprintf("  e2_mf: %.*f, e2_pc(mother): %.*f, e2_pc(father): %.*f\n\n",
                digits, result$crossperson_resid_corr["e2_mf"],
                digits, result$crossperson_resid_corr["e2_pc_mother"],
                digits, result$crossperson_resid_corr["e2_pc_father"]))
  }

  # Within-person correlations

  # Descriptive ratios (spousal residual covariance / latent rho)
  if (!is.null(result$crossperson_ratio)) {
    cat("Descriptive ratios (spousal residual cov / rho):\n")
    cat(sprintf("  e1_mf / rho: %.*f\n", digits, result$crossperson_ratio["e1_mf_over_rho"]))
    cat(sprintf("  e2_mf / rho: %.*f\n\n", digits, result$crossperson_ratio["e2_mf_over_rho"]))
  }

cat("Within-person Cor(PGS1,PGS2):\n")
  cat(sprintf("  mother: %.*f\n", digits, result$within_person_corr["mother"]))
  cat(sprintf("  father: %.*f\n", digits, result$within_person_corr["father"]))
  cat(sprintf("  child:  %.*f\n", digits, result$within_person_corr["child"]))
  for (nm in setdiff(names(result$within_person_corr), c("mother", "father", "child")))
    cat(sprintf("  %s: %.*f\n", nm, digits, result$within_person_corr[[nm]]))
  cat("\n")

  # Reliability and determinacy
  cat("PGS indicator reliability (Var explained by latent factor):\n")
  rel_wide <- reshape(result$reliability, idvar = "role", timevar = "PGS",
                      direction = "wide")
  print(rel_wide, row.names = FALSE)

  cat("\nFactor-score determinacy (Corr(latent, factor score)):\n")
  cat(sprintf("  mother: %.*f\n", digits, result$determinacy["mother"]))
  cat(sprintf("  father: %.*f\n", digits, result$determinacy["father"]))
  cat(sprintf("  child:  %.*f\n", digits, result$determinacy["child"]))
  for (nm in setdiff(names(result$determinacy), c("mother", "father", "child")))
    cat(sprintf("  %s: %.*f\n", nm, digits, result$determinacy[[nm]]))

  # Gold standard comparison (optional)
  if (!is.null(gold)) {
    cat("\n\nComparison to gold standard:\n")

    est <- if (result$h > 0) result$beta_total_assumed_std else result$beta_vis
    rho_est <- if (result$h > 0) result$latent_corr_total_assumed["rhoMF"] else result$latent_cov_vis["rhoMF"]

    gold_bD <- if ("bD" %in% names(gold)) gold["bD"] else NA
    gold_bM <- if ("bm" %in% names(gold)) gold["bm"] else NA
    gold_bF <- if ("bf" %in% names(gold)) gold["bf"] else NA
    gold_rho <- if ("rho" %in% names(gold)) gold["rho"] else NA

    cat(sprintf("  %-8s %12s %12s %10s\n", "", "Estimate", "Gold", "Bias"))

    if (!is.na(gold_bD) && abs(gold_bD) > 1e-8) {
      bias <- (est["bD"] - gold_bD) / gold_bD * 100
      cat(sprintf("  %-8s %12.*f %12.*f %9.1f%%\n", "beta_D", digits, est["bD"],
                  digits, gold_bD, bias))
    }
    if (!is.na(gold_bM) && abs(gold_bM) > 1e-8) {
      bias <- (est["bm"] - gold_bM) / gold_bM * 100
      cat(sprintf("  %-8s %12.*f %12.*f %9.1f%%\n", "beta_M", digits, est["bm"],
                  digits, gold_bM, bias))
    }
    if (!is.na(gold_bF) && abs(gold_bF) > 1e-8) {
      bias <- (est["bf"] - gold_bF) / gold_bF * 100
      cat(sprintf("  %-8s %12.*f %12.*f %9.1f%%\n", "beta_F", digits, est["bf"],
                  digits, gold_bF, bias))
    }
    if (!is.na(gold_rho) && abs(gold_rho) > 1e-8) {
      bias <- (rho_est - gold_rho) / gold_rho * 100
      cat(sprintf("  %-8s %12.*f %12.*f %9.1f%%\n", "rho", digits, rho_est,
                  digits, gold_rho, bias))
    }
  }

  cat("\n================================================================================\n")
  invisible(result)
}


# =============================================================================
# Compact summary print
# =============================================================================

#' Print compact RAVEL results table
#'
#' Prints a concise summary with point estimates plus either bootstrap
#' or lavaan-based uncertainty for both visible and hidden-corrected scales.
#'
#' @param result A fitted RAVEL result object from \code{ravel_fit}, or a corrected coefficient table.
#' @param digits Integer. Number of decimal places.
#' @param uncertainty Character. One of \code{"auto"}, \code{"bootstrap"},
#'   \code{"lavaan"}, or \code{"none"}.
#'
#' @export
print_simple <- function(result, digits = 4,
                         uncertainty = c("auto", "bootstrap", "lavaan", "none")) {

  if (is.data.frame(result)) {
    if (match.arg(uncertainty) != "auto")
      stop("A correction table already contains its chosen uncertainty; recompute it to change the source.")
    return(.ravel_print_correction_table(result, digits))
  }
  if (!result$converged) {
    cat(result$error %||% "Model did not converge.", "\n")
    return(invisible(result))
  }

  uncertainty <- match.arg(uncertainty)
  unc <- .resolve_uncertainty(result, uncertainty)

  fmt_est <- function(x) sprintf("%.*f", digits, x)

  fmt_ci <- function(est, ci) {
    if (is.null(ci)) return(fmt_est(est))
    sprintf("%.*f [%.*f, %.*f]", digits, est, digits, ci[1], digits, ci[2])
  }

  has_boot <- !is.null(result$bootstrap)
  has_h    <- result$h > 0 || (result$constraint_spec$mode %||% "none") != "none"

  cat(sprintf("RAVEL (h = %s, N = %s", fmt_est(result$h),
              if (!is.null(result$N)) result$N else "?"))
  if (!is.null(unc)) cat(sprintf(", %s", .uncertainty_label(unc)))
  cat(")\n")
  .ravel_print_assumptions(result, digits)
  if (!is.null(result$fixed_crossperson)) {
    cat("Fixed cross-person residual covariances:",
        paste(sprintf("%s=%s", names(result$fixed_crossperson), result$fixed_crossperson),
              collapse = ", "),
        "\n")
  }
  if (!is.null(result$fixed_crossperson_corr)) {
    cat("Fixed cross-person residual correlations:",
        paste(sprintf("%s=%s", names(result$fixed_crossperson_corr), result$fixed_crossperson_corr),
              collapse = ", "),
        "\n")
  }
  cat(strrep("-", 70), "\n\n")

  # Visible scale
  cat("Visible scale:\n")
  cat(sprintf("  beta_D:  %s\n", fmt_ci(result$beta_vis["bD"], .uncertainty_ci(unc, "visible", "bD"))))
  cat(sprintf("  beta_M:  %s\n", fmt_ci(result$beta_vis["bm"], .uncertainty_ci(unc, "visible", "bm"))))
  cat(sprintf("  beta_F:  %s\n", fmt_ci(result$beta_vis["bf"], .uncertainty_ci(unc, "visible", "bf"))))
  cat(sprintf("  rho_MF:  %s\n", fmt_ci(result$latent_cov_vis["rhoMF"], .uncertainty_ci(unc, "visible", "rhoMF"))))
  cat(sprintf("  R2:      %s\n", fmt_est(result$R2_model)))
  .print_parental_sibs(result, unc, digits, compact = TRUE)

  # Hidden-corrected
  if (has_h) {
    cat(sprintf("\nHidden-corrected (h = %s):\n", fmt_est(result$h)))
    cat(sprintf("  beta_D:  %s\n", fmt_ci(result$beta_total_assumed["bD"], .uncertainty_ci(unc, "total_raw", "bD"))))
    cat(sprintf("  beta_M:  %s\n", fmt_ci(result$beta_total_assumed["bm"], .uncertainty_ci(unc, "total_raw", "bm"))))
    cat(sprintf("  beta_F:  %s\n", fmt_ci(result$beta_total_assumed["bf"], .uncertainty_ci(unc, "total_raw", "bf"))))
    cat(sprintf("  Cov(Tm*,Tf*):  %s\n",
                fmt_ci(result$latent_cov_total_assumed["rhoMF"],  .uncertainty_ci(unc, "total_raw", "rhoMF"))))
    cat(sprintf("  Corr(Tm,Tf): %s\n",
                fmt_ci(result$latent_corr_total_assumed["rhoMF"], .uncertainty_ci(unc, "total_std", "rhoMF"))))
    cat("  Standardized total-liability coefficients:\n")
    cat(sprintf("    beta_D_std:  %s\n", fmt_ci(result$beta_total_assumed_std["bD"], .uncertainty_ci(unc, "total_std", "bD"))))
    cat(sprintf("    beta_M_std:  %s\n", fmt_ci(result$beta_total_assumed_std["bm"], .uncertainty_ci(unc, "total_std", "bm"))))
    cat(sprintf("    beta_F_std:  %s\n", fmt_ci(result$beta_total_assumed_std["bf"], .uncertainty_ci(unc, "total_std", "bf"))))
    if (!is.null(result$parental_sibs)) .print_parental_sibs(result, unc, digits, scale = "corrected_std")
    if (!is.null(result$R2_total_assumed) && !is.na(result$R2_total_assumed)) {
      cat(sprintf("  R2:      %s\n", fmt_est(result$R2_total_assumed)))
    }
  }

  # Reliabilities (compact)
  rel <- result$reliability
  cat(sprintf("\nReliabilities:  PGS1 = %.2f / %.2f / %.2f   PGS2 = %.2f / %.2f / %.2f  (m/f/c)\n",
              rel$reliability[rel$role == "mother" & rel$PGS == "PGS1"],
              rel$reliability[rel$role == "father" & rel$PGS == "PGS1"],
              rel$reliability[rel$role == "child"  & rel$PGS == "PGS1"],
              rel$reliability[rel$role == "mother" & rel$PGS == "PGS2"],
              rel$reliability[rel$role == "father" & rel$PGS == "PGS2"],
              rel$reliability[rel$role == "child"  & rel$PGS == "PGS2"]))

  cat(sprintf("Determinacy:    %.3f / %.3f / %.3f  (m/f/c)\n",
              result$determinacy["mother"],
              result$determinacy["father"],
              result$determinacy["child"]))

  cat(strrep("-", 70), "\n")
  invisible(result)
}
