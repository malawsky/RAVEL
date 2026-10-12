# Reconstruct an accumulated assortative-mating state without assigning an
# integer generation count. The component covariance is I + u w w', where
# w = (sqrt(1-h), sqrt(h)); current assortment need not equal past assortment.
.hidden_estimated_state <- function(rho, eta_m, eta_f, h) {
  values <- c(rho, eta_m, eta_f, h)
  if (!is.numeric(values) || is.complex(values) || length(values) != 4L ||
      any(!is.finite(values)) || h < 0 || h >= 1 || abs(rho) >= 1) return(NULL)
  rho <- unname(rho)
  eta_m <- unname(eta_m)
  eta_f <- unname(eta_f)
  h <- unname(h)
  # Symmetry is checked at public entry points. A smooth average here extends
  # the map off the equality manifold for the joint six-parameter Jacobian.
  eta <- mean(c(eta_m, eta_f))
  if (eta <= 0) return(NULL)
  B <- ((1 + rho) / (2 * eta))^2
  D <- 2 * B - 1 - rho
  if (!is.finite(B) || !is.finite(D) || B <= 0 || D <= 0) return(NULL)
  vGp <- 1 / D
  vTp <- (vGp - h) / (1 - h)
  # Random-mating founders have variance 1, and each subsequent generation
  # adds segregation variance 1/2 to a strictly positive inherited variance.
  if (!is.finite(vGp) || !is.finite(vTp) || vGp <= 0 || vTp <= .5) return(NULL)
  rho_T <- rho * vGp / (vGp - h)
  if (!is.finite(rho_T) || abs(rho_T) >= 1) return(NULL)
  w <- c(sqrt(1 - h), sqrt(h))
  parent <- diag(2L) + (vTp - 1) * tcrossprod(w)
  spouse <- rho_T * vTp * tcrossprod(w)
  parent_child <- .5 * (parent + spouse)
  child <- parent_child + .5 * diag(2L)
  full <- rbind(cbind(child, parent_child, parent_child),
                cbind(parent_child, parent, spouse),
                cbind(parent_child, spouse, parent))
  if (any(!is.finite(full)) ||
      min(eigen(full, symmetric = TRUE, only.values = TRUE)$values) <= 0) return(NULL)
  vG <- parent[1L, 1L]
  vH <- parent[2L, 2L]
  vGc <- child[1L, 1L]
  vHc <- child[2L, 2L]
  vTc <- .5 * (1 + rho_T) * vTp + .5
  # A diagnostic only: an arbitrary admissible accumulated state can lie
  # outside the trajectory from random-mating founders under today's r_T.
  # Such a state is not rejected or coerced onto that trajectory.
  effective <- NA_real_
  history_status <- "incompatible_with_constant_current_assortment"
  history_compatible <- FALSE
  numerical <- 64 * .Machine$double.eps
  if (rho_T == 0) {
    if (abs(vTp - 1) <= numerical * max(1, abs(vTp))) {
      history_status <- "duration_unidentified_at_random_mating"
      history_compatible <- TRUE
    }
  } else if (abs(rho_T) < sqrt(.Machine$double.eps)) {
    # Dividing a near-cancelled variance residual by a tiny correlation would
    # create an enormous numerical tolerance and could falsely report n=Inf.
    # Check the undivided founder-to-equilibrium variance interval instead;
    # leave duration unresolved when current assortment is this small.
    equilibrium_change <- rho_T / (1 - rho_T)
    observed_change <- vTp - 1
    variance_roundoff <- numerical * max(1, abs(vTp), abs(1 + equilibrium_change))
    if (observed_change >= min(0, equilibrium_change) - variance_roundoff &&
        observed_change <= max(0, equilibrium_change) + variance_roundoff) {
      history_status <- "duration_numerically_unresolved_near_random_mating"
      history_compatible <- NA
    }
  } else {
    persistence <- (1 + rho_T) / 2
    remaining <- (1 - (1 - rho_T) * vTp) / rho_T
    roundoff <- numerical * max(1, abs(vTp), 1 / abs(rho_T))
    if (abs(remaining) <= roundoff) {
      effective <- Inf
      history_status <- "equilibrium_limit"
      history_compatible <- TRUE
    } else if (remaining > 0 && remaining <= 1 + roundoff) {
      effective <- if (abs(remaining - 1) <= roundoff) 1 else
        1 + log(remaining) / log(persistence)
      history_status <- "continuous_duration_compatible"
      history_compatible <- TRUE
    }
  }
  a_par <- sqrt(1 - h) * sqrt(vG)
  b_par <- sqrt(h) * sqrt(vH)
  list(vG = vG, vH = vH, d = parent[1L, 2L] / sqrt(vG * vH),
    vGc = vGc, vHc = vHc,
    q = spouse[1L, 2L] / sqrt(vG * vH),
    rho_G = spouse[1L, 1L] / vG, rho_H = spouse[2L, 2L] / vH,
    t_ratio = sqrt((vH / vHc) / (vG / vGc)),
    a_par = a_par, b_par = b_par,
    a_ch = sqrt(1 - h) * sqrt(vGc), b_ch = sqrt(h) * sqrt(vHc),
    pi_present = b_par^2 / (a_par^2 + b_par^2), rho_T = rho_T,
    estimated_state = list(B = B, growth = B, D = D,
      vGp = vG, vGc = vGc, vTp = vTp, vTc = vTc, rT = rho_T,
      effective_generations = effective,
      constant_current_assortment_compatible = history_compatible,
      history_status = history_status, history_imposed = FALSE))
}

# Analytic derivatives on the same core-plus-pi coordinates used by the effect
# correction. This also retains uncertainty in the inferred accumulated state.
.hidden_estimated_state_jacobian <- function(rho, eta_m, eta_f, h) {
  state <- .hidden_estimated_state(rho, eta_m, eta_f, h)
  if (is.null(state)) stop("The inferred assortment state is incompatible with pi")
  labels <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f", "pi")
  J <- matrix(0, 4L, 7L, dimnames = list(c("growth", "visible_parent_variance",
    "total_parent_variance", "rho_total"), labels))
  B <- state$estimated_state$B
  vG <- state$vG
  eta <- mean(c(eta_m, eta_f))
  J[1L, "rhoMF"] <- 2 * B / (1 + rho)
  J[1L, c("eta_m", "eta_f")] <- -B / eta
  grad_D <- 2 * J[1L, ]
  grad_D["rhoMF"] <- grad_D["rhoMF"] - 1
  J[2L, ] <- -vG^2 * grad_D
  J[3L, ] <- J[2L, ] / (1 - h)
  J[3L, "pi"] <- J[3L, "pi"] + (vG - 1) / (1 - h)^2
  denominator <- vG - h
  J[4L, ] <- -rho * h / denominator^2 * J[2L, ]
  J[4L, "rhoMF"] <- J[4L, "rhoMF"] + vG / denominator
  J[4L, "pi"] <- J[4L, "pi"] + rho * vG / denominator^2
  J
}
