# Independent 14-coordinate covariance map. It does not use the lavaan syntax,
# fitted loadings, covariance derivatives, or package reconstruction helper.
esa_covariance <- function(p) {
  rho <- p[1]
  eta <- p[2]
  parent_reliability <- p[3:4]^2
  growth <- ((1 + rho) / (2 * eta))^2
  child_reliability <- growth * parent_reliability /
    (1 - parent_reliability + growth * parent_reliability)
  R <- matrix(c(1, eta, eta, eta, 1, rho, eta, rho, 1), 3L)
  sd_score <- exp(p[5:10])
  L <- rbind(diag(sqrt(c(child_reliability[1], rep(parent_reliability[1], 2L)))),
             diag(sqrt(c(child_reliability[2], rep(parent_reliability[2], 2L))))) * sd_score
  S <- L %*% R %*% t(L)
  E <- matrix(c(1, .5, .5, .5, 1, 0, .5, 0, 1), 3L)
  for (j in 1:2) {
    ix <- (3L * j - 2L):(3L * j)
    residual_sd <- sqrt(1 - c(child_reliability[j], rep(parent_reliability[j], 2L))) * sd_score[ix]
    S[ix, ix] <- S[ix, ix] + E * tcrossprod(residual_sd)
  }
  gy <- drop(L %*% R %*% p[11:13])
  rbind(cbind(S, gy), c(gy, p[14]))
}

esa_covariance_jacobian <- function(p, step = 1e-5) {
  vapply(seq_along(p), function(j) {
    plus <- minus <- p
    plus[j] <- plus[j] + step
    minus[j] <- minus[j] - step
    upper <- esa_covariance(plus)
    lower <- esa_covariance(minus)
    ((upper - lower) / (2 * step))[lower.tri(upper, diag = TRUE)]
  }, numeric(28L))
}

esa_recurrence <- function(rho, h, n) {
  w <- c(sqrt(1 - h), sqrt(h))
  parent <- diag(2L)
  advance <- function(S) {
    cov_total <- drop(S %*% w)
    spouse <- rho * tcrossprod(cov_total) / sum(w * cov_total)
    list(spouse = spouse, child = .5 * (S + spouse) + .5 * diag(2L))
  }
  # An independent converged recurrence supplies equilibrium too.
  iterations <- if (is.infinite(n)) 500L else n - 1L
  if (iterations > 0L) for (j in seq_len(iterations)) parent <- advance(parent)$child
  next_state <- advance(parent)
  pc <- .5 * (parent + next_state$spouse)
  list(parent = parent, child = next_state$child, spouse = next_state$spouse,
    pc = pc, rho = next_state$spouse[1, 1] / parent[1, 1],
    eta = pc[1, 1] / sqrt(parent[1, 1] * next_state$child[1, 1]),
    total_parent = sum(w * drop(parent %*% w)))
}

esa_correct <- function(p, h) {
  rho <- p[1]
  eta <- p[2]
  growth <- ((1 + rho) / (2 * eta))^2
  vp <- 1 / (2 * growth - 1 - rho)
  vc <- vp * growth
  tp <- (vp - h) / (1 - h)
  rt <- rho * vp / (vp - h)
  tc <- .5 * (1 + rt) * tp + .5
  kgpc <- .5 * vp * (1 + rho)
  ktpc <- .5 * tp * (1 + rt)
  KG <- matrix(c(vc, kgpc, kgpc, kgpc, vp, rho * vp, kgpc, rho * vp, vp), 3L)
  KT <- matrix(c(tc, ktpc, ktpc, ktpc, tp, rt * tp, ktpc, rt * tp, tp), 3L)
  raw <- drop(solve(KT, KG %*% (p[11:13] / sqrt(c(vc, vp, vp))))) / sqrt(1 - h)
  c(raw * sqrt(c(tc, tp, tp)), rt, raw, rt * tp, growth, vp, tp, rt)
}

test_that("the 14 covariance coordinates are generically identified", {
  for (rho in c(-.2, 0, .3)) for (eta in if (rho > 0) c(.58, .61) else c(.48, .61)) {
    p <- c(rho, eta, .7, .85, log(c(.7, 1.4, .9, 1.8, .8, 1.2)), .3, .08, .04, 1)
    # Include zero current assortment after a non-random historical state;
    # it does not have the same degeneracy as rho=0, eta=.5 below.
    J <- esa_covariance_jacobian(p)
    singular <- svd(J)$d
    expect_equal(sum(singular > max(singular) * 1e-7), 14L)
    expect_equal(J, esa_covariance_jacobian(p, 5e-6), tolerance = 3e-8)
    p[11:13] <- 0
    singular <- svd(esa_covariance_jacobian(p))$d
    expect_equal(sum(singular > max(singular) * 1e-7), 14L)
  }
})

test_that("random mating without outcome association loses one reliability direction", {
  p <- c(0, .5, .7, .85, rep(0, 6L), 0, 0, 0, 1)
  singular <- svd(esa_covariance_jacobian(p))$d
  expect_equal(sum(singular > max(singular) * 1e-7), 13L)
  expect_lt(tail(singular, 1L), 1e-9)
  shifted <- p
  shifted[3:4] <- p[3:4] * exp(c(.05, -.05))
  expect_equal(esa_covariance(p), esa_covariance(shifted), tolerance = 1e-12)
  # A nonzero outcome association identifies the ratio of the two loadings.
  p[11:13] <- c(.3, .08, .04)
  singular <- svd(esa_covariance_jacobian(p))$d
  expect_equal(sum(singular > max(singular) * 1e-7), 14L)
})

test_that("estimated-state reconstruction agrees with an independent matrix recurrence", {
  for (rho in c(-.6, .4)) for (h in c(0, .4, .8)) for (n in c(1, 3, 10, Inf)) {
    expected <- esa_recurrence(rho, h, n)
    got <- RAVEL:::.hidden_estimated_state(expected$rho, expected$eta, expected$eta, h)
    expect_false(is.null(got))
    parent <- matrix(c(got$vG, got$d * sqrt(got$vG * got$vH),
      got$d * sqrt(got$vG * got$vH), got$vH), 2L)
    spouse <- matrix(c(got$rho_G * got$vG, got$q * sqrt(got$vG * got$vH),
      got$q * sqrt(got$vG * got$vH), got$rho_H * got$vH), 2L)
    child <- .5 * (parent + spouse) + .5 * diag(2L)
    expect_equal(parent, expected$parent, tolerance = 1e-10)
    expect_equal(spouse, expected$spouse, tolerance = 1e-10)
    expect_equal(child, expected$child, tolerance = 1e-10)
    expect_equal(got$estimated_state$B, expected$child[1, 1] / expected$parent[1, 1], tolerance = 1e-10)
    expect_equal(got$estimated_state$vTp, expected$total_parent, tolerance = 1e-10)
    expect_equal(got$rho_T, rho, tolerance = 1e-10)
    expect_equal(got$estimated_state$effective_generations, n, tolerance = 1e-5)
  }
})

test_that("estimated-state uncertainty agrees with independent Gaussian Fisher information", {
  h <- .4
  pi_se <- .02
  N <- 5000L
  for (generation in c(3L, Inf)) {
  reference <- esa_recurrence(.4, h, generation)
  vp <- reference$parent[1, 1]
  p <- c(reference$rho, reference$eta, sqrt(vp / (vp + c(.4, 1.6))),
    rep(0, 6L), .3, .08, .04, 1)
  S <- esa_covariance(p)
  W <- solve(S)
  derivative <- lapply(seq_along(p), function(j) {
    plus <- minus <- p
    plus[j] <- plus[j] + 1e-5
    minus[j] <- minus[j] - 1e-5
    (esa_covariance(plus) - esa_covariance(minus)) / 2e-5
  })
  info <- outer(seq_along(p), seq_along(p), Vectorize(function(i, j)
    sum(diag(W %*% derivative[[i]] %*% W %*% derivative[[j]])) / 2))
  J <- vapply(seq_len(15L), function(j) {
    plus <- minus <- c(p, h)
    plus[j] <- plus[j] + 1e-6
    minus[j] <- minus[j] - 1e-6
    (esa_correct(plus[1:14], plus[15]) - esa_correct(minus[1:14], minus[15])) / 2e-6
  }, numeric(12L))
  V <- matrix(0, 15L, 15L)
  V[1:14, 1:14] <- solve(info) / N
  V[15, 15] <- pi_se^2
  expected_covariance <- J %*% V %*% t(J)
  x <- outer(seq_len(N), seq_len(7L), function(i, j) sin(i * j) + cos(i * (j + .5)))
  Q <- qr.Q(qr(scale(x, center = TRUE, scale = FALSE)))
  dat <- as.data.frame(sqrt(N) * Q %*% chol(S))
  names(dat) <- c(paste0("PGS1_", c("child", "mother", "father")),
    paste0("PGS2_", c("child", "mother", "father")), "Y")
  RAVEL:::.ensure_lavaan_runtime()
  fit <- ravel_fit(dat, h = h, fit_constraints = "estimated_state")
  expect_true(fit$converged)
  # The public path uses MLR sandwich SEs. The deterministic exact-moment
  # fixture is not Gaussian in fourth moments, so refit its population moments
  # with ordinary expected-information ML for this Gaussian Fisher comparison.
  dimnames(S) <- list(names(dat), names(dat))
  fit$fit <- lavaan::sem(fit$model, sample.cov = S, sample.nobs = N,
    sample.cov.rescale = FALSE, std.lv = TRUE, estimator = "ML",
    se = "standard", information = "expected", warn = FALSE)
  expect_true(RAVEL:::.ravel_constraint_fit_valid(fit$fit))
  pe <- lavaan::parameterEstimates(fit$fit)
  labels <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f")
  theta <- setNames(pe$est[match(labels, pe$label)], labels)
  fit$beta_vis <- theta[1:3]
  fit$latent_cov_vis <- theta[4:6]
  corrected <- ravel_correct_pi(fit, pi = h, pi_se = pi_se)
  state <- attr(corrected, "estimated_state")
  expect_equal(c(corrected$estimate, state$summary$estimate), esa_correct(p, h), tolerance = 3e-5)
  expect_equal(c(corrected$se, state$summary$se), sqrt(diag(expected_covariance)), tolerance = 2e-4)
  expect_equal(unname(attr(corrected, "covariance")), expected_covariance[1:8, 1:8], tolerance = 3e-5)
  expect_equal(unname(state$covariance), expected_covariance[9:12, 9:12], tolerance = 3e-5)
  }
})
