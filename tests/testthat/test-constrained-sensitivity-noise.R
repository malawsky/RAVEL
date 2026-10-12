test_that("pi sensitivity retries larger steps without relaxing its accuracy guard", {
  labels <- c("bD","bm","bf","rhoMF","eta_m","eta_f")
  base <- setNames(c(.3,.08,.04,.2,.58,.58),labels)
  derivative <- setNames(c(-.002,.001,.001,.002,-.004,-.004),labels)
  h0 <- .4
  object <- list(constraint_spec=list(mode="n_gen",h=h0,n_generations=3L,depends_on_pi=TRUE),
    constraint_refit=list(h=h0),beta_vis=base[1:3],latent_cov_vis=base[4:6])
  calls <- 0L
  testthat::local_mocked_bindings(ravel_fit=function(h) {
    calls <<- calls+1L
    theta <- base+derivative*(h-h0)
    # Bounded optimizer noise at the finest step creates a spurious derivative.
    if(abs(abs(h-h0)-.0005)<1e-12) theta[1] <- theta[1]+sign(h-h0)*1e-6
    list(converged=TRUE,beta_vis=theta[1:3],latent_cov_vis=theta[4:6])
  },.package="RAVEL")
  got <- RAVEL:::.ravel_constraint_core_h_sensitivity(object)
  expect_equal(got,derivative,tolerance=1e-10)
  expect_equal(calls,6L)
})

test_that("pi sensitivity continues to reject inconsistent derivatives", {
  labels <- c("bD","bm","bf","rhoMF","eta_m","eta_f")
  base <- setNames(c(.3,.08,.04,.2,.58,.58),labels)
  h0 <- .4
  object <- list(constraint_spec=list(mode="n_gen",h=h0,n_generations=3L,depends_on_pi=TRUE),
    constraint_refit=list(h=h0),beta_vis=base[1:3],latent_cov_vis=base[4:6])
  testthat::local_mocked_bindings(ravel_fit=function(h) {
    theta <- base
    theta[1] <- theta[1]+sign(h-h0)*1e-4
    list(converged=TRUE,beta_vis=theta[1:3],latent_cov_vis=theta[4:6])
  },.package="RAVEL")
  expect_error(RAVEL:::.ravel_constraint_core_h_sensitivity(object),"numerically unstable")
})

test_that("forward pi sensitivity stays inside the domain near zero", {
  labels <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f")
  base <- setNames(c(.3, .08, .04, .2, .58, .58), labels)
  derivative <- setNames(c(-.002, .001, .001, .002, -.004, -.004), labels)
  h0 <- .0001
  object <- list(constraint_spec = list(mode = "n_gen", h = h0,
    n_generations = 3L, depends_on_pi = TRUE), constraint_refit = list(h = h0),
    beta_vis = base[1:3], latent_cov_vis = base[4:6])
  called <- numeric()
  testthat::local_mocked_bindings(ravel_fit = function(h) {
    called <<- c(called, h)
    theta <- base + derivative * (h - h0)
    list(converged = TRUE, beta_vis = theta[1:3], latent_cov_vis = theta[4:6])
  }, .package = "RAVEL")
  got <- RAVEL:::.ravel_constraint_core_h_sensitivity(object)
  expect_equal(got, derivative, tolerance = 1e-10)
  expect_length(called, 2L)
  expect_true(all(called > h0 & called < 1))
})

test_that("upper-boundary pi sensitivity caps do not repeat identical retry steps", {
  labels <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f")
  base <- setNames(c(.3, .08, .04, .2, .58, .58), labels)
  h0 <- .9999
  object <- list(constraint_spec = list(mode = "n_gen", h = h0,
    n_generations = 3L, depends_on_pi = TRUE), constraint_refit = list(h = h0),
    beta_vis = base[1:3], latent_cov_vis = base[4:6])
  called <- numeric()
  testthat::local_mocked_bindings(ravel_fit = function(h) {
    called <<- c(called, h)
    theta <- base
    theta[1] <- theta[1] + sign(h - h0) * 1e-4
    list(converged = TRUE, beta_vis = theta[1:3], latent_cov_vis = theta[4:6])
  }, .package = "RAVEL")
  expect_error(RAVEL:::.ravel_constraint_core_h_sensitivity(object), "numerically unstable")
  expect_length(called, 4L)
  expect_length(unique(called), 4L)
  expect_true(all(called > 0 & called < 1))
})
