# =============================================================================
# ravel_expected_values.R
# =============================================================================
# Analytic expected-value solver for the RAVEL hidden-variant trio model.
#
# This file is standalone: it defines ravel_expected_values() together with the
# two small internal helpers it depends on (.safe_solve, .latent_corr_from_cov),
# so it can be sourced independently of the main ravel_sem.R file. When sourced
# alongside ravel_sem.R, the internal helpers defined here are identical to the
# ones used by the fitting code, so there is no conflict.
#
# Author: Daniel Malawsky (dm22@sanger.ac.uk)
# =============================================================================


# -----------------------------------------------------------------------------
# Internal helpers (copied from ravel_sem.R so this file is standalone)
# -----------------------------------------------------------------------------

#' Safe linear solve (avoids hard failure on singular matrices)
#' @noRd
.safe_solve <- function(A, B = NULL) {
  out <- tryCatch({
    if (is.null(B)) solve(A) else solve(A, B)
  }, error = function(e) {
    tryCatch({
      MASS::ginv(A) %*% B
    }, error = function(e2) {
      matrix(NA_real_, nrow = nrow(A), ncol = if (is.null(B)) ncol(A) else ncol(as.matrix(B)))
    })
  })
  out
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
      stop("The supplied rho_T_corr and transmission coefficients do not yield a finite total-liability equilibrium variance.")
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

#' Component-variance state for analytical expectations
#'
#' Tracks AM-induced inflation of the visible and hidden component variances,
#' together with their standardized within-person correlation. This is the
#' forward-model analogue of `.am_state_components()` in `ravel_sem.R`.
#' Unlike the fitted correction, the transmission coefficients are known here,
#' so the recursion permits general `tm` and `tf`; `tm = tf = 0.5` reproduces
#' additive Mendelian transmission used by the simulations.
#'
#' @noRd
.expected_am_state_components <- function(rho_T_corr,
                                          h,
                                          tm,
                                          tf,
                                          am_correction,
                                          am_n_generations,
                                          tol = 1e-13,
                                          max_iter = 10000L) {
  S <- tm^2 + tf^2
  u <- tm * tf
  segregation_var <- 1 - S
  a0 <- sqrt(1 - h)
  b0 <- sqrt(h)

  evaluate <- function(vG, vH, d) {
    a <- a0 * sqrt(vG)
    b <- b0 * sqrt(vH)
    varT <- a^2 + b^2 + 2 * a * b * d
    if (!is.finite(varT) || varT <= 0) return(NULL)

    covGT <- a + b * d
    covHT <- b + a * d
    list(
      varT = varT,
      rho_G = rho_T_corr * covGT^2 / varT,
      q = rho_T_corr * covGT * covHT / varT,
      rho_H = rho_T_corr * covHT^2 / varT
    )
  }

  advance <- function(vG, vH, d, state) {
    vGc <- (S + 2 * u * state$rho_G) * vG + segregation_var
    vHc <- (S + 2 * u * state$rho_H) * vH + segregation_var
    if (any(!is.finite(c(vGc, vHc))) || vGc <= 0 || vHc <= 0) return(NULL)

    d_child <- sqrt(vG * vH) * (S * d + 2 * u * state$q) /
      sqrt(vGc * vHc)
    list(vG = vGc, vH = vHc, d = d_child)
  }

  vG <- 1
  vH <- 1
  d <- 0

  if (am_correction == "equilibrium") {
    converged <- FALSE
    for (iteration in seq_len(max_iter)) {
      state <- evaluate(vG, vH, d)
      if (is.null(state)) return(NULL)
      next_state <- advance(vG, vH, d, state)
      if (is.null(next_state)) return(NULL)

      delta <- max(abs(next_state$vG - vG),
                   abs(next_state$vH - vH),
                   abs(next_state$d - d))
      vG <- next_state$vG
      vH <- next_state$vH
      d <- next_state$d
      if (delta < tol) {
        converged <- TRUE
        break
      }
    }
    if (!converged) return(NULL)
  } else if (am_n_generations > 1L) {
    for (generation in seq_len(am_n_generations - 1L)) {
      state <- evaluate(vG, vH, d)
      if (is.null(state)) return(NULL)
      next_state <- advance(vG, vH, d, state)
      if (is.null(next_state)) return(NULL)
      vG <- next_state$vG
      vH <- next_state$vH
      d <- next_state$d
    }
  }

  spouse_state <- evaluate(vG, vH, d)
  if (is.null(spouse_state)) return(NULL)
  child_state <- advance(vG, vH, d, spouse_state)
  if (is.null(child_state)) return(NULL)

  a_parent <- a0 * sqrt(vG)
  b_parent <- b0 * sqrt(vH)

  list(
    vG = vG,
    vH = vH,
    d = d,
    vGc = child_state$vG,
    vHc = child_state$vH,
    d_child = child_state$d,
    rho_G = spouse_state$rho_G,
    q = spouse_state$q,
    rho_H = spouse_state$rho_H,
    a_parent = a_parent,
    b_parent = b_parent,
    a_child = a0 * sqrt(child_state$vG),
    b_child = b0 * sqrt(child_state$vH),
    pi_present = b_parent^2 / (a_parent^2 + b_parent^2)
  )
}

#' Convert standardized total-liability coefficients to raw total-liability coefficients
#' @noRd
.unstandardize_total_betas <- function(beta_total_std, var_total) {
  c(
    bD = unname(beta_total_std["bD"]) / sqrt(unname(var_total["child"])),
    bm = unname(beta_total_std["bm"]) / sqrt(unname(var_total["mother"])),
    bf = unname(beta_total_std["bf"]) / sqrt(unname(var_total["father"]))
  )
}

#' Normalize a beta_D/beta_M/beta_F vector
#' @noRd
.normalize_beta_trio <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 3L || anyNA(x)) {
    stop(arg_name, " must be a numeric vector of length 3.")
  }

  roles <- c("bD", "bm", "bf")
  if (is.null(names(x))) {
    out <- as.numeric(x)
    names(out) <- roles
    return(out)
  }

  alias_map <- c(
    bD = "bD", beta_D = "bD", child = "bD", c = "bD", Gc = "bD",
    bm = "bm", bM = "bm", beta_M = "bm", mother = "bm", m = "bm", Gm = "bm",
    bf = "bf", bF = "bf", beta_F = "bf", father = "bf", f = "bf", Gf = "bf"
  )
  mapped <- unname(alias_map[names(x)])
  if (any(is.na(mapped))) {
    stop(arg_name, " names must identify bD/bm/bf (or child/mother/father).")
  }
  if (length(unique(mapped)) != length(roles)) {
    stop(arg_name, " must supply exactly one value for each of bD, bm, and bf.")
  }

  out <- setNames(rep(NA_real_, length(roles)), roles)
  out[mapped] <- as.numeric(x)
  out
}

#' Normalize a role-indexed numeric vector to child/mother/father
#' @noRd
.normalize_role_vector <- function(x, arg_name) {
  roles <- c("child", "mother", "father")

  if (!is.numeric(x) || anyNA(x)) {
    stop(arg_name, " must be numeric and non-missing.")
  }

  if (length(x) == 1L) {
    out <- rep(as.numeric(x), length(roles))
    names(out) <- roles
    return(out)
  }

  if (length(x) != length(roles)) {
    stop(arg_name, " must have length 1 or 3.")
  }

  if (is.null(names(x))) {
    out <- as.numeric(x)
    names(out) <- roles
    return(out)
  }

  alias_map <- c(
    child = "child", c = "child", Gc = "child",
    mother = "mother", m = "mother", Gm = "mother",
    father = "father", f = "father", Gf = "father"
  )
  mapped <- unname(alias_map[names(x)])
  if (any(is.na(mapped))) {
    stop(arg_name, " names must be drawn from child/mother/father (or c/m/f, Gc/Gm/Gf).")
  }
  if (length(unique(mapped)) != length(roles)) {
    stop(arg_name, " must supply exactly one value for each of child, mother, and father.")
  }

  out <- setNames(rep(NA_real_, length(roles)), roles)
  out[mapped] <- as.numeric(x)
  out
}

#' Normalize cross-person residual correlations to mf/mc/fc
#' @noRd
.normalize_pgs_resid_corr <- function(x) {
  labs <- c("mf", "mc", "fc")

  if (!is.numeric(x) || anyNA(x)) {
    stop("residual_corr must be numeric and non-missing.")
  }

  if (length(x) == 1L) {
    stop("residual_corr must have length 3 or be a named vector with mf/mc/fc entries.")
  }

  if (length(x) != length(labs)) {
    stop("residual_corr must have length 3.")
  }

  if (is.null(names(x))) {
    out <- as.numeric(x)
    names(out) <- labs
  } else {
    alias_map <- c(
      mf = "mf", fm = "mf",
      mc = "mc", cm = "mc",
      fc = "fc", cf = "fc"
    )
    mapped <- unname(alias_map[names(x)])
    if (any(is.na(mapped))) {
      stop("residual_corr names must be mf, mc, fc (or fm/cm/cf).")
    }
    if (length(unique(mapped)) != length(labs)) {
      stop("residual_corr must supply exactly one value for each of mf, mc, and fc.")
    }
    out <- setNames(rep(NA_real_, length(labs)), labs)
    out[mapped] <- as.numeric(x)
  }

  if (any(abs(out) >= 1)) {
    stop("residual_corr entries must satisfy abs(value) < 1.")
  }

  out
}


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

#' Analytical expectations under the RAVEL hidden-variant model
#'
#' @description
#' Computes the population quantities analytically implied by a specified set of
#' total-liability parameters under finite-\eqn{n} or equilibrium assortative
#' mating. This is useful for checking simulations, understanding scale
#' transformations, and verifying what RAVEL should recover in expectation under
#' the maintained assumptions.
#'
#' @details
#' The function assumes
#' \deqn{
#'   T^* = \sqrt{1-h}\,G^* + \sqrt{h}\,H^*,
#' }
#' where \eqn{G^*} is the raw visible component, \eqn{H^*} is the raw hidden
#' component, and \eqn{h} corresponds to the theory parameter \eqn{\pi} on the
#' founder or random-mating reference scale. Role-standardized total liabilities
#' are denoted \eqn{T}.
#' The total-liability variance used to standardize total-scale coefficients is
#' computed from the assortative-mating history and transmission coefficients,
#' independent of \eqn{h}; the hidden fraction changes the projection onto
#' visible and observed-PGS scales rather than the total-liability target.
#'
#' The user supplies:
#' \itemize{
#'   \item the total-liability effect vector \eqn{(\beta_D^{(T)}, \beta_M^{(T)}, \beta_F^{(T)})},
#'   \item the hidden-variance fraction \eqn{h},
#'   \item the standardized spouse correlation on the total-liability scale,
#'         \eqn{\rho_T = \operatorname{Corr}(T_m, T_f)},
#'   \item the raw maternal and paternal transmission coefficients \eqn{t_m, t_p},
#'   \item and either a finite number of generations of assortative mating or the
#'         equilibrium limit.
#' }
#'
#' It returns the analytically expected:
#' \itemize{
#'   \item visible-scale trio covariance structure,
#'   \item total-scale raw covariances and standardized correlations,
#'   \item cross-block \eqn{\Sigma_{GH}} (called \code{Sigma_xh} in code),
#'   \item visible-scale regression coefficients implied by the supplied total
#'         coefficients,
#'   \item and the corrected coefficients recovered by applying the RAVEL matrix
#'         correction back to those implied visible coefficients.
#' }
#'
#' @param beta_total Named numeric vector of length 3 giving the total-liability
#'   coefficients. Names should be \code{c("bD","bm","bf")}. By default these
#'   are interpreted as coefficients \eqn{\beta^{raw}} on \eqn{T^*}, matching the native
#'   simulator. Set \code{beta_total_scale = "standardized"} to interpret them
#'   as coefficients \eqn{\beta} on role-standardized \eqn{T} instead.
#' @param beta_total_scale Character. \code{"raw"} means \code{beta_total}
#'   multiplies the raw total liability
#'   \eqn{T^*=\sqrt{1-h}G^*+\sqrt{h}H^*}; \code{"standardized"} means it
#'   multiplies role-standardized \eqn{T}. Returned objects include both scales.
#' @param h Numeric in \code{[0,1)}. Hidden additive-variance fraction \eqn{\pi}
#'   on the founder or random-mating reference scale.
#' @param rho_T_corr Numeric. Standardized spouse correlation on the total-liability
#'   scale, \eqn{\rho_T = \operatorname{Corr}(T_m,T_f)}.
#' @param tm,tf Numeric. Raw maternal and paternal transmission coefficients.
#'   The legacy argument \code{tf} corresponds to paternal \eqn{t_p} notation.
#'   Defaults to \code{0.5, 0.5}, the maintained additive trio model.
#' @param am_correction Character. Either \code{"n_gen"} for finite-\eqn{n}
#'   assortative mating or \code{"equilibrium"} for the fixed-point limit.
#' @param am_n_generations Integer. Number of assortative-mating generations when
#'   \code{am_correction = "n_gen"}. The value \code{1} reproduces the
#'   one-generation correction.
#'
#' @return A list containing the analytically implied covariance blocks, visible-
#'   and total-scale parameters, and expected visible and corrected coefficients.
#'
#' @export
ravel_expected_values <- function(beta_total,
                                  h = 0,
                                  rho_T_corr = 0,
                                  tm = 0.5,
                                  tf = 0.5,
                                  am_correction = c("n_gen", "equilibrium"),
                                  am_n_generations = 1L,
                                  beta_total_scale = c("raw", "standardized")) {

  am_correction <- match.arg(am_correction)
  beta_total_scale <- match.arg(beta_total_scale)

  beta_total_supplied <- .normalize_beta_trio(beta_total, "beta_total")

  if (!is.numeric(h) || length(h) != 1L || is.na(h) || h < 0 || h >= 1) {
    stop("h must be a single numeric value in [0,1).")
  }
  if (!is.numeric(rho_T_corr) || length(rho_T_corr) != 1L || is.na(rho_T_corr) ||
      abs(rho_T_corr) >= 1) {
    stop("rho_T_corr must be a single numeric value with abs(rho_T_corr) < 1.")
  }
  if (!is.numeric(tm) || !is.numeric(tf) || length(tm) != 1L || length(tf) != 1L ||
      anyNA(c(tm, tf))) {
    stop("tm and tf must be single finite numeric values.")
  }
  if (am_correction == "n_gen") {
    if (!is.numeric(am_n_generations) || length(am_n_generations) != 1L ||
        is.na(am_n_generations) || am_n_generations < 1) {
      stop("am_n_generations must be a single integer >= 1 when am_correction = 'n_gen'.")
    }
    am_n_generations <- as.integer(am_n_generations)
  } else {
    am_n_generations <- NA_integer_
  }

  total_state <- .total_liability_state(
    rho_T_corr = rho_T_corr,
    tm = tm,
    tf = tf,
    am_correction = am_correction,
    am_n_generations = am_n_generations
  )

  state <- .expected_am_state_components(
    rho_T_corr = rho_T_corr,
    h = h,
    tm = tm,
    tf = tf,
    am_correction = am_correction,
    am_n_generations = am_n_generations
  )
  if (is.null(state)) {
    stop("The supplied parameters do not yield a finite component-variance state.")
  }

  d_parent <- state$d
  d_child <- state$d_child
  q <- state$q
  rho_vis <- state$rho_G
  rho_hidden <- state$rho_H

  # Transmission coefficients between generation-standardized components.
  tG_m <- tm * sqrt(state$vG / state$vGc)
  tG_f <- tf * sqrt(state$vG / state$vGc)
  tH_m <- tm * sqrt(state$vH / state$vHc)
  tH_f <- tf * sqrt(state$vH / state$vHc)

  eta_m_vis <- tG_m + rho_vis * tG_f
  eta_f_vis <- tG_f + rho_vis * tG_m
  eta_m_hid <- tH_m + rho_hidden * tH_f
  eta_f_hid <- tH_f + rho_hidden * tH_m

  Cov_Gc_Hm <- tG_m * d_parent + tG_f * q
  Cov_Gc_Hf <- tG_f * d_parent + tG_m * q
  Cov_Gm_Hc <- tH_m * d_parent + tH_f * q
  Cov_Gf_Hc <- tH_f * d_parent + tH_m * q
  Cov_Gc_Hc <- tG_m * tH_m * d_parent +
    tG_m * tH_f * q +
    tG_f * tH_m * q +
    tG_f * tH_f * d_parent

  Sigma_xx <- matrix(
    c(1,         eta_m_vis, eta_f_vis,
      eta_m_vis, 1,         rho_vis,
      eta_f_vis, rho_vis,   1),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("Gc", "Gm", "Gf"), c("Gc", "Gm", "Gf"))
  )

  Sigma_xh <- matrix(
    c(Cov_Gc_Hc, Cov_Gc_Hm, Cov_Gc_Hf,
      Cov_Gm_Hc, d_parent,  q,
      Cov_Gf_Hc, q,         d_parent),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("Gc", "Gm", "Gf"), c("Hc", "Hm", "Hf"))
  )

  Sigma_hh <- matrix(
    c(1,         eta_m_hid, eta_f_hid,
      eta_m_hid, 1,         rho_hidden,
      eta_f_hid, rho_hidden, 1),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("Hc", "Hm", "Hf"), c("Hc", "Hm", "Hf"))
  )

  a_vec <- c(child = state$a_child,
             mother = state$a_parent,
             father = state$a_parent)
  b_vec <- c(child = state$b_child,
             mother = state$b_parent,
             father = state$b_parent)
  d_vec <- c(child = Cov_Gc_Hc, mother = d_parent, father = d_parent)
  var_total_component <- a_vec^2 + b_vec^2 + 2 * a_vec * b_vec * d_vec
  var_total <- total_state$var_total

  if (beta_total_scale == "raw") {
    beta_total_raw <- beta_total_supplied
    beta_total_std <- .standardize_total_betas(beta_total_raw, var_total)
  } else {
    beta_total_std <- beta_total_supplied
    beta_total_raw <- .unstandardize_total_betas(beta_total_std, var_total)
  }

  # Reconstruct raw total covariances from the component state. These should
  # equal the pi-invariant values returned by `.total_liability_state()`.
  kappa_T_mf <- a_vec["mother"] * a_vec["father"] * rho_vis +
    a_vec["mother"] * b_vec["father"] * q +
    b_vec["mother"] * a_vec["father"] * q +
    b_vec["mother"] * b_vec["father"] * rho_hidden
  kappa_T_cm <- a_vec["child"] * a_vec["mother"] * eta_m_vis +
    a_vec["child"] * b_vec["mother"] * Cov_Gc_Hm +
    b_vec["child"] * a_vec["mother"] * Cov_Gm_Hc +
    b_vec["child"] * b_vec["mother"] * eta_m_hid
  kappa_T_cf <- a_vec["child"] * a_vec["father"] * eta_f_vis +
    a_vec["child"] * b_vec["father"] * Cov_Gc_Hf +
    b_vec["child"] * a_vec["father"] * Cov_Gf_Hc +
    b_vec["child"] * b_vec["father"] * eta_f_hid

  latent_cov_vis <- c(rhoMF = rho_vis, eta_m = eta_m_vis, eta_f = eta_f_vis)
  latent_cov_total_component <- c(
    rhoMF = unname(kappa_T_mf),
    eta_m = unname(kappa_T_cm),
    eta_f = unname(kappa_T_cf)
  )
  latent_cov_total <- total_state$cov

  latent_corr_total <- total_state$corr

  M <- diag(as.numeric(a_vec), nrow = 3L, ncol = 3L) +
    .safe_solve(Sigma_xx, Sigma_xh) %*%
    diag(as.numeric(b_vec), nrow = 3L, ncol = 3L)
  dimnames(M) <- list(c("bD", "bm", "bf"), c("bD", "bm", "bf"))

  beta_vis_expected <- as.numeric(M %*% beta_total_raw)
  names(beta_vis_expected) <- names(beta_total_raw)

  beta_total_recovered <- as.numeric(.safe_solve(M, beta_vis_expected))
  names(beta_total_recovered) <- names(beta_total_raw)
  beta_total_recovered_std <- .standardize_total_betas(beta_total_recovered, var_total)
  visible_to_total_matrix <- .safe_solve(M)
  dimnames(visible_to_total_matrix) <- list(c("bD", "bm", "bf"), c("bD", "bm", "bf"))

  list(
    beta_total_supplied = beta_total_supplied,
    beta_total_scale = beta_total_scale,
    beta_total_input = beta_total_raw,
    beta_total_input_raw = beta_total_raw,
    beta_total_input_std = beta_total_std,
    beta_vis_expected = beta_vis_expected,
    beta_total_recovered = beta_total_recovered,
    beta_total_recovered_std = beta_total_recovered_std,
    h = h,
    rho_T_corr = rho_T_corr,
    tm = tm,
    tf = tf,
    d_parent = d_parent,
    d_child = d_child,
    q = q,
    Sigma_xx = Sigma_xx,
    Sigma_xh = Sigma_xh,
    Sigma_hh = Sigma_hh,
    total_to_visible_matrix_component = M,
    total_to_visible_matrix = M,
    visible_to_total_matrix = visible_to_total_matrix,
    latent_cov_vis_expected    = latent_cov_vis,
    latent_cov_total_expected  = latent_cov_total,
    latent_cov_total_component_expected = latent_cov_total_component,
    latent_corr_total_expected = latent_corr_total,
    var_total_expected = var_total,
    var_total_component_expected = var_total_component,
    var_visible_raw_expected = c(child = state$vGc,
                                 mother = state$vG,
                                 father = state$vG),
    var_hidden_raw_expected = c(child = state$vHc,
                                mother = state$vH,
                                father = state$vH),
    pi_present_expected = state$pi_present,
    tG = c(mother = tG_m, father = tG_f),
    tH = c(mother = tH_m, father = tH_f),
    kappa_T_expected = unname(total_state$cov["rhoMF"]),
    am_correction = am_correction,
    am_n_generations = am_n_generations
  )
}

#' Analytical expected trio OLS on observed PGS
#'
#' @description
#' Computes the analytically expected trio OLS coefficients from a regression of
#' a phenotype on one observed polygenic score per trio member. This extends
#' \code{\link{ravel_expected_values}} from the latent visible-genetic scale to
#' the observed standardized-PGS scale.
#'
#' @details
#' The function first uses \code{\link{ravel_expected_values}} to derive the
#' expected visible-scale trio covariance block and the implied visible-scale
#' regression coefficients under the hidden-variant model. It then maps those
#' latent expectations to one observed PGS per role using
#' \deqn{
#'   P_r = \lambda_r G_r + \varepsilon_r,
#' }
#' with \eqn{\operatorname{Var}(P_r)=1}. If
#' \code{reliability_scale = "observed"}, the loading and reliability coincide
#' via
#' \deqn{
#'   \operatorname{Rel}(P_r \mid G_r) = \lambda_r^2,
#' }
#' so if \code{loadings} is omitted the function sets
#' \eqn{\lambda_r = \sqrt{\mathrm{reliability}_r}} and
#' \eqn{\theta_r = 1 - \lambda_r^2}. If
#' \code{reliability_scale = "random_mating"}, the supplied reliability is
#' first converted to the final role-specific reliability implied by the AM
#' history.
#'
#' Cross-person score-residual dependence is specified on the residual
#' correlation scale via \code{residual_corr}, with defaults
#' \code{mf = 0}, \code{mc = 0.5}, and \code{fc = 0.5}. These are converted to
#' observed-score residual covariances using the role-specific residual
#' variances implied by the standardized-score assumption.
#'
#' The expected trio OLS coefficients are then
#' \deqn{
#'   \beta_{\mathrm{PGS}} = \Sigma_{PP}^{-1} \operatorname{Cov}(P, Y),
#' }
#' where \eqn{P = (P_c, P_m, P_f)^\top}. The reported coefficients are on the
#' phenotype scale implied by \code{beta_total}; if \code{outcome_sd = 1}
#' (default), these coincide with the coefficients from a regression on a
#' standardized outcome.
#'
#' @inheritParams ravel_expected_values
#' @param reliability Numeric vector of length 1 or 3 giving the PGS
#'   reliabilities for child, mother, and father. By default these are
#'   interpreted on the random-mating/base visible-variance scale, matching the
#'   simulation convention in which PGS noise is calibrated
#'   before AM-induced visible-variance inflation. Set
#'   \code{reliability_scale = "observed"} to interpret them as final
#'   standardized-score reliabilities, in which case they equal
#'   \eqn{\lambda_r^2}. Ignored when \code{loadings} is supplied.
#' @param reliability_scale Character. \code{"random_mating"} means
#'   \code{reliability} is calibrated against the random-mating/base visible
#'   variance and is converted internally to role-specific final observed-score
#'   reliability under the supplied AM history. \code{"observed"} means
#'   \code{reliability} is already the final reliability of the standardized PGS.
#' @param loadings Optional numeric vector of length 1 or 3 giving the observed
#'   score loadings on the visible latent factors for child, mother, and father.
#'   Under standardized observed scores, reliability is \code{loadings^2}.
#' @param residual_corr Numeric vector giving the cross-person residual
#'   correlations for mother-father (\code{mf}), mother-child (\code{mc}), and
#'   father-child (\code{fc}). Defaults to \code{c(mf = 0, mc = 0.5, fc = 0.5)}.
#' @param outcome_sd Numeric. Standard deviation of the phenotype scale on which
#'   \code{beta_total} is defined. Defaults to \code{1}, corresponding to a
#'   standardized outcome.
#'
#' @return A list containing the expected observed-score covariance block,
#'   the implied \code{Cov(P, Y)} vector, and the expected trio OLS
#'   coefficients on the observed-PGS scale.
#'
#' @export
ravel_expected_pgs_ols <- function(beta_total,
                                   h = 0,
                                   rho_T_corr = 0,
                                   tm = 0.5,
                                   tf = 0.5,
                                   am_correction = c("n_gen", "equilibrium"),
                                   am_n_generations = 1L,
                                   reliability = c(child = 0.5, mother = 0.5, father = 0.5),
                                   reliability_scale = c("random_mating", "observed"),
                                   loadings = NULL,
                                   residual_corr = c(mf = 0, mc = 0.5, fc = 0.5),
                                   outcome_sd = 1,
                                   beta_total_scale = c("raw", "standardized")) {

  beta_total_scale <- match.arg(beta_total_scale)
  reliability_scale <- match.arg(reliability_scale)

  latent <- ravel_expected_values(
    beta_total = beta_total,
    beta_total_scale = beta_total_scale,
    h = h,
    rho_T_corr = rho_T_corr,
    tm = tm,
    tf = tf,
    am_correction = am_correction,
    am_n_generations = am_n_generations
  )

  if (!is.numeric(outcome_sd) || length(outcome_sd) != 1L || is.na(outcome_sd) || outcome_sd <= 0) {
    stop("outcome_sd must be a single positive numeric value.")
  }

  roles <- c("child", "mother", "father")
  var_visible_raw <- latent$var_visible_raw_expected
  reliability_input <- .normalize_role_vector(reliability, "reliability")
  if (any(reliability_input < 0 | reliability_input > 1)) {
    stop("reliability entries must lie in [0, 1].")
  }

  if (is.null(loadings)) {
    if (reliability_scale == "random_mating") {
      noise_var <- ifelse(
        reliability_input == 0,
        Inf,
        (1 - reliability_input) / reliability_input
      )
      reliability <- var_visible_raw / (var_visible_raw + noise_var)
      reliability[is.nan(reliability)] <- 0
      names(reliability) <- roles
    } else {
      reliability <- reliability_input
    }
    loadings <- sqrt(reliability)
  } else {
    loadings <- .normalize_role_vector(loadings, "loadings")
    if (any(abs(loadings) > 1)) {
      stop("Under standardized observed scores, abs(loadings) must be <= 1.")
    }
    reliability <- loadings^2
    names(reliability) <- roles
    reliability_scale <- "observed"
  }

  theta <- 1 - loadings^2
  names(theta) <- roles

  resid_corr <- .normalize_pgs_resid_corr(residual_corr)
  psi_cov <- setNames(
    c(
      unname(resid_corr["mf"] * sqrt(theta["mother"] * theta["father"])),
      unname(resid_corr["mc"] * sqrt(theta["mother"] * theta["child"])),
      unname(resid_corr["fc"] * sqrt(theta["father"] * theta["child"]))
    ),
    c("mf", "mc", "fc")
  )

  Lambda <- matrix(
    0,
    nrow = 3,
    ncol = 3,
    dimnames = list(c("Pc", "Pm", "Pf"), c("Gc", "Gm", "Gf"))
  )
  diag(Lambda) <- c(loadings["child"], loadings["mother"], loadings["father"])

  Psi <- matrix(
    c(
      theta["child"],  psi_cov["mc"], psi_cov["fc"],
      psi_cov["mc"],   theta["mother"], psi_cov["mf"],
      psi_cov["fc"],   psi_cov["mf"], theta["father"]
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("Pc", "Pm", "Pf"), c("Pc", "Pm", "Pf"))
  )

  Sigma_xx_for_pgs <- latent$Sigma_xx
  Cov_gstd_traw <- latent$Sigma_xx %*% latent$total_to_visible_matrix
  dimnames(Cov_gstd_traw) <- list(c("Gc", "Gm", "Gf"), c("bD", "bm", "bf"))
  cov_GY <- as.numeric(Cov_gstd_traw %*% latent$beta_total_input_raw)
  names(cov_GY) <- c("Gc", "Gm", "Gf")
  beta_vis_for_pgs <- latent$beta_vis_expected

  Sigma_pp <- Lambda %*% Sigma_xx_for_pgs %*% t(Lambda) + Psi
  cov_PY <- as.numeric(Lambda %*% cov_GY)
  names(cov_PY) <- c("child", "mother", "father")

  visible_to_pgs_matrix <- .safe_solve(Sigma_pp, Lambda %*% Sigma_xx_for_pgs)
  dimnames(visible_to_pgs_matrix) <- list(
    c("child", "mother", "father"),
    c("bD", "bm", "bf")
  )
  total_to_pgs_matrix <- .safe_solve(Sigma_pp, Lambda %*% Cov_gstd_traw)
  dimnames(total_to_pgs_matrix) <- list(
    c("child", "mother", "father"),
    c("bD", "bm", "bf")
  )

  beta_pgs_raw <- as.numeric(total_to_pgs_matrix %*% latent$beta_total_input_raw)
  names(beta_pgs_raw) <- c("child", "mother", "father")
  beta_pgs_expected <- beta_pgs_raw / outcome_sd

  list(
    beta_total_input = latent$beta_total_input,
    beta_total_input_raw = latent$beta_total_input_raw,
    beta_total_input_std = latent$beta_total_input_std,
    beta_total_scale = latent$beta_total_scale,
    beta_vis_expected = latent$beta_vis_expected,
    beta_vis_expected_for_pgs = beta_vis_for_pgs,
    beta_pgs_expected = beta_pgs_expected,
    beta_pgs_expected_stdY = beta_pgs_expected,
    beta_pgs_expected_rawY = beta_pgs_raw,
    outcome_sd = outcome_sd,
    loadings = loadings,
    reliability = reliability,
    reliability_input = reliability_input,
    reliability_scale = reliability_scale,
    var_visible_raw = var_visible_raw,
    theta = theta,
    residual_corr = resid_corr,
    psi_cov = psi_cov,
    Sigma_pp = Sigma_pp,
    Sigma_xx_for_pgs = Sigma_xx_for_pgs,
    Cov_gstd_traw_for_pgs = Cov_gstd_traw,
    cov_PY = cov_PY,
    visible_to_pgs_matrix = visible_to_pgs_matrix,
    total_to_pgs_matrix = total_to_pgs_matrix,
    latent_expectations = latent
  )
}

#' Infer total-liability coefficients from visible-scale coefficients
#'
#' @description
#' Inverts the RAVEL expected-value projection. Given coefficients from a
#' regression on the visible latent liabilities \eqn{G_c,G_m,G_f}, this function
#' returns the implied coefficients on raw total liabilities and standardized
#' total liabilities under the supplied hidden-variant and assortative-mating
#' assumptions.
#'
#' @inheritParams ravel_expected_values
#' @param beta_visible Numeric vector of length 3 giving visible-scale
#'   coefficients for child, mother, and father. Names may be \code{bD/bm/bf} or
#'   \code{child/mother/father}.
#' @param beta_visible_scale Character. \code{"raw_outcome"} means the supplied
#'   coefficients are from a regression on the original phenotype scale.
#'   \code{"standardized_outcome"} means they are from a regression on a
#'   phenotype divided by \code{outcome_sd}.
#' @param outcome_sd Numeric. Standard deviation of the original phenotype. This
#'   is only used to convert between raw-outcome and standardized-outcome
#'   coefficient scales.
#'
#' @return A list with visible coefficients and implied total-liability
#'   coefficients on raw-outcome and standardized-outcome scales.
#'
#' @export
ravel_expected_from_visible <- function(beta_visible,
                                        h = 0,
                                        rho_T_corr = 0,
                                        tm = 0.5,
                                        tf = 0.5,
                                        am_correction = c("n_gen", "equilibrium"),
                                        am_n_generations = 1L,
                                        beta_visible_scale = c("raw_outcome", "standardized_outcome"),
                                        outcome_sd = 1) {

  am_correction <- match.arg(am_correction)
  beta_visible_scale <- match.arg(beta_visible_scale)

  if (!is.numeric(outcome_sd) || length(outcome_sd) != 1L || is.na(outcome_sd) || outcome_sd <= 0) {
    stop("outcome_sd must be a single positive numeric value.")
  }

  beta_visible_input <- .normalize_beta_trio(beta_visible, "beta_visible")
  beta_visible_rawY <- if (beta_visible_scale == "standardized_outcome") {
    beta_visible_input * outcome_sd
  } else {
    beta_visible_input
  }
  beta_visible_stdY <- beta_visible_rawY / outcome_sd

  latent <- ravel_expected_values(
    beta_total = c(bD = 0, bm = 0, bf = 0),
    beta_total_scale = "raw",
    h = h,
    rho_T_corr = rho_T_corr,
    tm = tm,
    tf = tf,
    am_correction = am_correction,
    am_n_generations = am_n_generations
  )

  beta_total_rawT_rawY <- as.numeric(.safe_solve(latent$total_to_visible_matrix, beta_visible_rawY))
  names(beta_total_rawT_rawY) <- c("bD", "bm", "bf")
  beta_total_stdT_rawY <- .standardize_total_betas(beta_total_rawT_rawY, latent$var_total_expected)
  beta_total_rawT_stdY <- beta_total_rawT_rawY / outcome_sd
  beta_total_stdT_stdY <- beta_total_stdT_rawY / outcome_sd

  list(
    beta_visible_input = beta_visible_input,
    beta_visible_scale = beta_visible_scale,
    beta_visible_rawY = beta_visible_rawY,
    beta_visible_stdY = beta_visible_stdY,
    beta_total_rawT_rawY = beta_total_rawT_rawY,
    beta_total_stdT_rawY = beta_total_stdT_rawY,
    beta_total_rawT_stdY = beta_total_rawT_stdY,
    beta_total_stdT_stdY = beta_total_stdT_stdY,
    outcome_sd = outcome_sd,
    total_to_visible_matrix = latent$total_to_visible_matrix,
    visible_to_total_matrix = latent$visible_to_total_matrix,
    latent_expectations = latent
  )
}

#' Infer visible and total-liability coefficients from trio PGS OLS coefficients
#'
#' @description
#' Inverts the expected observed-PGS trio regression. Given coefficients from
#' \code{lm(Y ~ PGS_child + PGS_mother + PGS_father)} with standardized PGSs,
#' this function returns the implied visible-scale and total-liability
#' coefficients under the supplied reliability, residual-correlation, hidden-
#' variant, and assortative-mating assumptions.
#'
#' @inheritParams ravel_expected_pgs_ols
#' @param beta_pgs Numeric vector of length 3 giving observed-PGS trio
#'   coefficients for child, mother, and father.
#' @param beta_pgs_scale Character. \code{"raw_outcome"} means the supplied PGS
#'   coefficients are from a regression on the original phenotype scale.
#'   \code{"standardized_outcome"} means they are from a regression on a
#'   phenotype divided by \code{outcome_sd}.
#'
#' @return A list with PGS, visible, and total-liability coefficients on
#'   raw-outcome and standardized-outcome scales.
#'
#' @export
ravel_expected_from_pgs_ols <- function(beta_pgs,
                                        h = 0,
                                        rho_T_corr = 0,
                                        tm = 0.5,
                                        tf = 0.5,
                                        am_correction = c("n_gen", "equilibrium"),
                                        am_n_generations = 1L,
                                        reliability = c(child = 0.5, mother = 0.5, father = 0.5),
                                        reliability_scale = c("random_mating", "observed"),
                                        loadings = NULL,
                                        residual_corr = c(mf = 0, mc = 0.5, fc = 0.5),
                                        beta_pgs_scale = c("raw_outcome", "standardized_outcome"),
                                        outcome_sd = 1) {

  am_correction <- match.arg(am_correction)
  reliability_scale <- match.arg(reliability_scale)
  beta_pgs_scale <- match.arg(beta_pgs_scale)

  if (!is.numeric(outcome_sd) || length(outcome_sd) != 1L || is.na(outcome_sd) || outcome_sd <= 0) {
    stop("outcome_sd must be a single positive numeric value.")
  }

  beta_pgs_input <- .normalize_role_vector(beta_pgs, "beta_pgs")
  beta_pgs_rawY <- if (beta_pgs_scale == "standardized_outcome") {
    beta_pgs_input * outcome_sd
  } else {
    beta_pgs_input
  }
  beta_pgs_stdY <- beta_pgs_rawY / outcome_sd

  pgs_zero <- ravel_expected_pgs_ols(
    beta_total = c(bD = 0, bm = 0, bf = 0),
    beta_total_scale = "raw",
    h = h,
    rho_T_corr = rho_T_corr,
    tm = tm,
    tf = tf,
    am_correction = am_correction,
    am_n_generations = am_n_generations,
    reliability = reliability,
    reliability_scale = reliability_scale,
    loadings = loadings,
    residual_corr = residual_corr,
    outcome_sd = outcome_sd
  )

  beta_visible_rawY <- as.numeric(.safe_solve(pgs_zero$visible_to_pgs_matrix, beta_pgs_rawY))
  names(beta_visible_rawY) <- c("bD", "bm", "bf")
  beta_visible_stdY <- beta_visible_rawY / outcome_sd

  visible_inverse <- ravel_expected_from_visible(
    beta_visible = beta_visible_rawY,
    h = h,
    rho_T_corr = rho_T_corr,
    tm = tm,
    tf = tf,
    am_correction = am_correction,
    am_n_generations = am_n_generations,
    beta_visible_scale = "raw_outcome",
    outcome_sd = outcome_sd
  )

  beta_total_rawT_rawY <- as.numeric(.safe_solve(pgs_zero$total_to_pgs_matrix, beta_pgs_rawY))
  names(beta_total_rawT_rawY) <- c("bD", "bm", "bf")
  total_expectation <- ravel_expected_values(
    beta_total = beta_total_rawT_rawY,
    beta_total_scale = "raw",
    h = h,
    rho_T_corr = rho_T_corr,
    tm = tm,
    tf = tf,
    am_correction = am_correction,
    am_n_generations = am_n_generations
  )
  beta_total_stdT_rawY <- total_expectation$beta_total_input_std
  beta_total_rawT_stdY <- beta_total_rawT_rawY / outcome_sd
  beta_total_stdT_stdY <- beta_total_stdT_rawY / outcome_sd

  list(
    beta_pgs_input = beta_pgs_input,
    beta_pgs_scale = beta_pgs_scale,
    beta_pgs_rawY = beta_pgs_rawY,
    beta_pgs_stdY = beta_pgs_stdY,
    beta_visible_rawY = beta_visible_rawY,
    beta_visible_stdY = beta_visible_stdY,
    beta_total_rawT_rawY = beta_total_rawT_rawY,
    beta_total_stdT_rawY = beta_total_stdT_rawY,
    beta_total_rawT_stdY = beta_total_rawT_stdY,
    beta_total_stdT_stdY = beta_total_stdT_stdY,
    outcome_sd = outcome_sd,
    loadings = pgs_zero$loadings,
    reliability = pgs_zero$reliability,
    reliability_input = pgs_zero$reliability_input,
    reliability_scale = pgs_zero$reliability_scale,
    residual_corr = pgs_zero$residual_corr,
    Sigma_pp = pgs_zero$Sigma_pp,
    visible_to_pgs_matrix = pgs_zero$visible_to_pgs_matrix,
    total_to_pgs_matrix = pgs_zero$total_to_pgs_matrix,
    total_expectation = total_expectation,
    visible_inverse = visible_inverse,
    pgs_expectations_at_zero = pgs_zero
  )
}
