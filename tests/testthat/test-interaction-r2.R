# R2 and variance decompositions in genetic x covariate interaction models
# must describe the whole genetic predictor b'G + sum_j g_j'(G x_j), so they
# cannot depend on interaction_center.

r2_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else
    rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  expr
}

# Independent standard normal G (child), x and e with Y = G x + e:
# the genetic predictor explains exactly 1/2 of Var(Y) and the main effect at
# centre c is c, so a main-effect-only R2 would be c^2 / 2.
r2_pure_gxe <- function(N = 8000L, seed = 31L) r2_with_seed(seed, {
  rho <- .2; rel <- c(.5, .3)
  z <- MASS::mvrnorm(N, c(0, 0), matrix(c(1, rho, rho, 1), 2))
  Gm <- z[, 1]; Gf <- z[, 2]; Gc <- (Gm + Gf) / 2 + sqrt(.5 * (1 - rho)) * rnorm(N)
  P <- list()
  for (k in 1:2) {
    em <- rnorm(N); ef <- rnorm(N); ec <- (em + ef) / 2 + sqrt(.5) * rnorm(N)
    P[[paste0("PGS", k, "_mother")]] <- sqrt(rel[k]) * Gm + sqrt(1 - rel[k]) * em
    P[[paste0("PGS", k, "_father")]] <- sqrt(rel[k]) * Gf + sqrt(1 - rel[k]) * ef
    P[[paste0("PGS", k, "_child")]] <- sqrt(rel[k]) * Gc + sqrt(1 - rel[k]) * ec
  }
  x <- rnorm(N)
  data.frame(FID = seq_len(N), as.data.frame(P), Y = Gc * x + rnorm(N), x = x)
})

test_that("interaction-model R2 is the full genetic predictor and ignores the centre", {
  skip_on_cran()
  d <- r2_pure_gxe()
  fits <- lapply(c(0, 1, 2), function(cc) suppressWarnings(ravel_fit(d, interactions = "x",
    interaction_roles = "child", interaction_center = c(x = cc))))
  r2 <- vapply(fits, `[[`, numeric(1), "R2_model")
  expect_equal(r2, rep(r2[1], 3), tolerance = 1e-4)
  expect_lt(abs(r2[1] - .5), .03)
  for (f in fits) {
    expect_equal(sum(f$components_vis), f$s2_genetic_vis, tolerance = 1e-12)
    expect_true(all(c("interaction", "cross_main_interaction") %in% names(f$components_vis)))
    expect_equal(f$components_vis, fits[[1]]$components_vis, tolerance = 1e-3)
  }
  # pi = 0: the total-scale quantity uses the same predictor.
  r2t <- vapply(fits, `[[`, numeric(1), "R2_total_assumed")
  expect_equal(r2t, rep(r2t[1], 3), tolerance = 1e-4)
  expect_lt(abs(r2t[1] - .5), .03)
})

# Hidden variation (pi = 0.4), two rounds of assortment, both parental
# siblings, and two moderators acting on every role's total liability.
r2_sib_data <- function(N = 6000L, seed = 8L) {
  s <- ravel_simulate_trios(N = N, m = 300, m_hidden = 300, hidden_fraction = .4, rho_MF = .4,
    grandparent_rho_MF = .4, parental_sibs = "both", seed = seed)
  z <- function(x) as.numeric(scale(x))
  Tz <- cbind(z(s$G_child), z(s$G_mother), z(s$G_father), z(s$G_mat_sib), z(s$G_pat_sib))
  r2_with_seed(seed, {
    age <- stats::runif(N, 4, 16); x2 <- stats::rnorm(N, 3, 2)
    Pg <- drop(Tz %*% c(.30, .10, .05, .04, 0)) + drop((Tz * (age - 10)) %*% c(.02, .01, 0, .01, 0)) +
      drop((Tz * (x2 - 3)) %*% c(-.03, 0, .02, 0, .01))
    Y <- Pg + .04 * age - .1 * x2 + rnorm(N, sd = .8)
    list(data = data.frame(s[grep("^PGS|^FID", names(s))], Y = Y, age = age, x2 = x2),
         oracle = var(Pg) / var(Y))
  })
}

test_that("R2 is centre-invariant with parental siblings and several moderators", {
  skip_on_cran()
  sim <- r2_sib_data()
  fit_at <- function(center) suppressWarnings(ravel_fit(sim$data, pi = .4, am_n_generations = 2L,
    parental_sibs = "both", interactions = c("age", "x2"), interaction_center = center))
  a <- fit_at(c(age = 10, x2 = 3))
  b <- fit_at(c(age = 4, x2 = -1))
  expect_true(a$converged && b$converged)
  expect_equal(b$R2_model, a$R2_model, tolerance = 1e-4)
  expect_equal(b$components_vis, a$components_vis, tolerance = 1e-3)
  expect_equal(b$R2_total_assumed, a$R2_total_assumed, tolerance = 1e-4)
  expect_equal(b$components_total_assumed, a$components_total_assumed, tolerance = 1e-3)
  expect_equal(sum(a$components_total_assumed), a$s2_genetic_total_assumed, tolerance = 1e-12)
  # The total-liability R2 recovers the realised share of the generating predictor.
  expect_lt(abs(a$R2_total_assumed - sim$oracle), .03)
})
