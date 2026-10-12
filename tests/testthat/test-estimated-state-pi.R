estimated_state_reference <- function(history, h, equilibrium = FALSE) {
  w <- c(sqrt(1 - h), sqrt(h))
  P <- diag(2L)
  if (equilibrium) P <- P + history[1] / (1 - history[1]) * tcrossprod(w)
  advance <- function(P, r) {
    cov_total <- drop(P %*% w)
    spouse <- r * tcrossprod(cov_total) / sum(w * cov_total)
    list(child = .5 * (P + spouse) + .5 * diag(2L), spouse = spouse)
  }
  if (length(history) > 1L) for (r in history[-length(history)]) P <- advance(P, r)$child
  next_state <- advance(P, tail(history, 1L))
  C <- next_state$child
  K <- next_state$spouse
  PC <- .5 * (P + K)
  rho <- K[1, 1] / P[1, 1]
  eta <- PC[1, 1] / sqrt(P[1, 1] * C[1, 1])
  list(parent = P, child = C, spouse = K, parent_child = PC, rho = rho, eta = eta,
    vTp = drop(crossprod(w, P %*% w)), vTc = drop(crossprod(w, C %*% w)),
    rT = tail(history, 1L))
}

estimated_state_packet <- function(h = .4, history = rep(.35, 3L)) {
  reference <- estimated_state_reference(history, h)
  labels <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f")
  theta <- setNames(c(.3, .08, .04, reference$rho, reference$eta, reference$eta), labels)
  # eta_m and eta_f are the same fitted parameter; retain that singularity.
  A <- rbind(diag(5L), c(0, 0, 0, 0, 1))
  V <- A %*% diag(1e-5, 5L) %*% t(A)
  dimnames(V) <- list(labels, labels)
  list(schema = "ravel_visible_transfer_v1", diagnostics = list(admissible = TRUE,
    estimated_state_identified = TRUE),
    settings = list(uncertainty = "lavaan", residual_model = "fixed",
      constraint_spec = list(mode = "estimated_state", h = h, n_generations = NA_integer_,
        depends_on_pi = FALSE)), outcome = list(name = "trait"), theta = theta,
    vcov_core = V, software = list(transform_fingerprint = RAVEL:::.report_transform_fingerprint()))
}

test_that("closed estimated state recovers finite and equilibrium component states", {
  for (h in c(0, .4, .8)) for (r in c(-.5, .35)) for (n in c(1L, 2L, 10L, Inf)) {
    ref <- estimated_state_reference(if (is.infinite(n)) r else rep(r, n), h,
                                     equilibrium = is.infinite(n))
    got <- RAVEL:::.hidden_estimated_state(ref$rho, ref$eta, ref$eta, h)
    expect_false(is.null(got))
    expect_equal(c(got$vG, got$vH), diag(ref$parent), tolerance = 2e-10)
    expect_equal(c(got$vGc, got$vHc), diag(ref$child), tolerance = 2e-10)
    expect_equal(got$d, ref$parent[1, 2] / sqrt(prod(diag(ref$parent))), tolerance = 2e-10)
    expect_equal(got$q, ref$spouse[1, 2] / sqrt(prod(diag(ref$parent))), tolerance = 2e-10)
    expect_equal(got$rho_T, ref$rT, tolerance = 2e-10)
    expect_equal(got$estimated_state$vTp, ref$vTp, tolerance = 2e-10)
    expect_equal(got$estimated_state$vTc, ref$vTc, tolerance = 2e-10)
    expect_equal(got$estimated_state$effective_generations, n, tolerance = 2e-7)
  }
})

test_that("estimated state admits changing history without inventing a duration", {
  ref <- estimated_state_reference(c(.7, .7, 0), .4)
  got <- RAVEL:::.hidden_estimated_state(ref$rho, ref$eta, ref$eta, .4)
  expect_false(is.null(got))
  expect_equal(got$estimated_state$vTp, ref$vTp, tolerance = 1e-10)
  expect_equal(got$rho_T, 0)
  expect_false(got$estimated_state$constant_current_assortment_compatible)
  expect_true(is.na(got$estimated_state$effective_generations))
  expect_identical(got$estimated_state$history_status, "incompatible_with_constant_current_assortment")
  random <- RAVEL:::.hidden_estimated_state(0, .5, .5, .4)
  expect_identical(random$estimated_state$history_status, "duration_unidentified_at_random_mating")
  expect_true(is.na(random$estimated_state$effective_generations))
})

test_that("tiny current assortment cannot fabricate an equilibrium duration", {
  for (r in c(-1e-12, 1e-12, -1e-9, 1e-9)) {
    ref <- estimated_state_reference(rep(r, 10L), .4)
    got <- RAVEL:::.hidden_estimated_state(ref$rho, ref$eta, ref$eta, .4)
    expect_true(is.na(got$estimated_state$effective_generations))
    expect_identical(got$estimated_state$history_status,
      "duration_numerically_unresolved_near_random_mating")
    expect_true(is.na(got$estimated_state$constant_current_assortment_compatible))
  }
  previous <- estimated_state_reference(c(.7, .7, 1e-12), .4)
  changed <- RAVEL:::.hidden_estimated_state(previous$rho, previous$eta, previous$eta, .4)
  expect_false(changed$estimated_state$constant_current_assortment_compatible)
  expect_identical(changed$estimated_state$history_status,
    "incompatible_with_constant_current_assortment")
  expect_true(is.na(changed$estimated_state$effective_generations))
})

test_that("estimated state rejects physically impossible founder states", {
  expect_null(RAVEL:::.hidden_estimated_state(.5, .95, .95, .4))
  eta_first <- (1 + .5) / (2 * sqrt(1 + .5 / 2))
  expect_null(RAVEL:::.hidden_estimated_state(.5, eta_first, eta_first, .8))
  # Positive-definite current covariances are insufficient: founders plus
  # segregation variance 1/2 cannot give parental total variance <= 1/2.
  h <- .4
  impossible_vTp <- .3
  vGp <- h + (1 - h) * impossible_vTp
  rho <- 0
  B <- .5 * (1 + rho + 1 / vGp)
  eta <- (1 + rho) / (2 * sqrt(B))
  expect_null(RAVEL:::.hidden_estimated_state(rho, eta, eta, h))
  expect_null(RAVEL:::.hidden_estimated_state(0, .5, .5, 1))
})

test_that("estimated-state analytic Jacobian agrees with independent numerical derivatives", {
  packet <- estimated_state_packet()
  x <- c(packet$theta, pi = .4)
  evaluate <- function(z) {
    rho <- z["rhoMF"]
    eta <- mean(z[c("eta_m", "eta_f")])
    B <- ((1 + rho) / (2 * eta))^2
    vG <- 1 / (2 * B - 1 - rho)
    c(B, vG, (vG - z["pi"]) / (1 - z["pi"]), rho * vG / (vG - z["pi"]))
  }
  J <- matrix(0, 4L, 7L)
  for (j in seq_along(x)) {
    plus <- minus <- x
    plus[j] <- plus[j] + 1e-6
    minus[j] <- minus[j] - 1e-6
    J[, j] <- (evaluate(plus) - evaluate(minus)) / 2e-6
  }
  got <- RAVEL:::.hidden_estimated_state_jacobian(x["rhoMF"], x["eta_m"], x["eta_f"], x["pi"])
  expect_equal(unname(got), J, tolerance = 2e-8)
})

test_that("estimated-state pi correction inherits mode and retains joint state uncertainty", {
  packet <- estimated_state_packet()
  result <- ravel_correct_pi(packet, .4, pi_se = .02)
  expect_true(all(result$am_correction == "estimated_state"))
  expect_true(all(is.na(result$am_n_generations)))
  expect_false(attr(result, "pi_uncertainty")$fixed_am_history)
  expect_equal(unname(attr(result, "pi_uncertainty")$core_pi_sensitivity), rep(0, 6))
  state <- attr(result, "estimated_state")
  expect_equal(state$rT, .35, tolerance = 1e-8)
  expect_equal(nrow(state$summary), 4L)
  expect_equal(state$covariance, state$jacobian %*% attr(result, "joint_input_covariance") %*%
                 t(state$jacobian), tolerance = 1e-12)
  expect_true(all(is.finite(state$summary$se)))
  changed <- ravel_correct_pi(packet, .3, pi_se = .02)
  expect_true(all(is.finite(changed$se)))
  expect_false(isTRUE(all.equal(result$estimate, changed$estimate)))
  fixed <- ravel_correct_pi(packet, .4, pi_se = 0)
  expect_true(all(state$summary$se >= attr(fixed, "estimated_state")$summary$se - 1e-12))
  expect_gt(state$summary$se[4], attr(fixed, "estimated_state")$summary$se[4])
})

test_that("estimated-state public correction enforces fit and symmetry prerequisites", {
  packet <- estimated_state_packet()
  missing <- packet
  missing$settings$constraint_spec <- NULL
  expect_error(ravel_correct_pi(missing, .4, pi_se = .02, am_correction = "estimated_state"), "matching fit_constraints")
  expect_error(ravel_correct_pi(packet, .4, pi_se = .02, am_correction = "equilibrium"), "mating history")
  packet$theta["eta_f"] <- packet$theta["eta_f"] + .01
  expect_error(ravel_correct_pi(packet, .4, pi_se = .02), "equal maternal")
})

test_that("estimated-state packets must attest a successful full-fit identification check", {
  packet <- estimated_state_packet()
  for (flag in list(NULL, FALSE, NA, 1)) {
    packet$diagnostics$estimated_state_identified <- flag
    expect_error(ravel_correct_pi(packet, .4, pi_se = .02), "full-fit identification check")
  }
  packet$diagnostics$estimated_state_identified <- TRUE
  expect_true(all(is.finite(ravel_correct_pi(packet, .4, pi_se = .02)$se)))
})
