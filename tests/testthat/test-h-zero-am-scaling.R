am_variance <- function(rho, generation) {
  value <- 1
  if (generation > 0L) {
    for (i in seq_len(generation)) {
      value <- ((1 + rho) / 2) * value + 1 / 2
    }
  }
  value
}

scale_estimates <- function(beta_vis, eta_m, eta_f, rho, mode, generations = 1L) {
  RAVEL:::.compute_scale_estimates(
    beta_vis = beta_vis,
    cov_cm = eta_m,
    cov_cf = eta_f,
    cov_mf = rho,
    h = 0,
    am_correction = mode,
    am_n_generations = generations
  )
}

expect_h_zero_map <- function(result, beta_vis, variances, visible_covariances) {
  expected_map <- diag(sqrt(variances))

  expect_equal(result$var_total_assumed, variances, tolerance = 1e-10)
  expect_equal(
    result$beta_total_assumed,
    beta_vis / sqrt(variances),
    tolerance = 1e-10
  )
  expect_equal(result$beta_total_assumed_std, beta_vis, tolerance = 1e-10)
  expect_equal(result$latent_corr_total_assumed, visible_covariances, tolerance = 1e-10)
  expect_equal(unname(result$hidden_correction_details$M), expected_map, tolerance = 1e-10)
  expect_equal(
    unname(result$hidden_correction_details$Sigma_xh),
    matrix(0, 3, 3),
    tolerance = 1e-10
  )
}

test_that("h = 0 recovers the former identity under random mating", {
  beta_vis <- c(bD = 0.30, bm = 0.10, bf = -0.05)
  visible_covariances <- c(rhoMF = 0, eta_m = 0.5, eta_f = 0.5)
  result <- scale_estimates(beta_vis, 0.5, 0.5, 0, "n_gen", 1L)

  expect_h_zero_map(
    result,
    beta_vis,
    c(child = 1, mother = 1, father = 1),
    visible_covariances
  )
  expect_equal(result$beta_total_assumed, beta_vis, tolerance = 1e-10)
  expect_equal(result$latent_cov_total_assumed, visible_covariances, tolerance = 1e-10)
})

test_that("h = 0 retains raw scaling after one mating round", {
  rho <- 0.4
  var_child <- 1 + rho / 2
  eta <- (1 + rho) / (2 * sqrt(var_child))
  beta_vis <- c(bD = 0.30, bm = 0.10, bf = -0.05)
  result <- scale_estimates(beta_vis, eta, eta, rho, "n_gen", 1L)

  expect_h_zero_map(
    result,
    beta_vis,
    c(child = var_child, mother = 1, father = 1),
    c(rhoMF = rho, eta_m = eta, eta_f = eta)
  )
  expect_equal(result$hidden_correction_details$rho_T_corr_assumed, rho, tolerance = 1e-10)
  expect_equal(result$hidden_correction_details$eta_total_assumed,
               c(eta_m = (1 + rho) / 2, eta_f = (1 + rho) / 2),
               tolerance = 1e-10)
})

test_that("h = 0 retains raw scaling after finite repeated assortment", {
  rho <- 0.4
  generations <- 5L
  var_parent <- am_variance(rho, generations - 1L)
  var_child <- am_variance(rho, generations)
  eta_raw <- 0.5 * (1 + rho) * var_parent
  eta <- eta_raw / sqrt(var_child * var_parent)
  beta_vis <- c(bD = 0.30, bm = 0.10, bf = -0.05)
  result <- scale_estimates(beta_vis, eta, eta, rho, "n_gen", generations)

  expect_h_zero_map(
    result,
    beta_vis,
    c(child = var_child, mother = var_parent, father = var_parent),
    c(rhoMF = rho, eta_m = eta, eta_f = eta)
  )
  expect_equal(
    result$latent_cov_total_assumed,
    c(rhoMF = rho * var_parent, eta_m = eta_raw, eta_f = eta_raw),
    tolerance = 1e-10
  )
})

test_that("h = 0 retains raw scaling at equilibrium", {
  rho <- 0.4
  var_equilibrium <- 1 / (1 - rho)
  eta <- (1 + rho) / 2
  beta_vis <- c(bD = 0.30, bm = 0.10, bf = -0.05)
  result <- scale_estimates(beta_vis, eta, eta, rho, "equilibrium")

  expect_h_zero_map(
    result,
    beta_vis,
    setNames(rep(var_equilibrium, 3), c("child", "mother", "father")),
    c(rhoMF = rho, eta_m = eta, eta_f = eta)
  )
  expect_equal(
    result$latent_cov_total_assumed,
    c(
      rhoMF = rho * var_equilibrium,
      eta_m = eta * var_equilibrium,
      eta_f = eta * var_equilibrium
    ),
    tolerance = 1e-10
  )
})

test_that("h = 0 correction agrees with the analytical expected-value solver", {
  for (mode in c("n_gen", "equilibrium")) {
    generations <- if (mode == "n_gen") 5L else 1L
    expected <- ravel_expected_values(
      beta_total = c(bD = 0.31, bm = 0.12, bf = -0.04),
      beta_total_scale = "raw",
      h = 0,
      rho_T_corr = 0.4,
      tm = 0.5,
      tf = 0.5,
      am_correction = mode,
      am_n_generations = generations
    )
    visible <- expected$latent_cov_vis_expected
    result <- scale_estimates(
      expected$beta_vis_expected,
      visible["eta_m"],
      visible["eta_f"],
      visible["rhoMF"],
      mode,
      generations
    )

    expect_equal(result$var_total_assumed, expected$var_total_expected, tolerance = 1e-10)
    expect_equal(result$beta_total_assumed, expected$beta_total_recovered, tolerance = 1e-10)
    expect_equal(result$beta_total_assumed_std,
                 expected$beta_total_recovered_std, tolerance = 1e-10)
    expect_equal(result$latent_cov_total_assumed,
                 expected$latent_cov_total_expected, tolerance = 1e-10)
    expect_equal(result$latent_corr_total_assumed,
                 expected$latent_corr_total_expected, tolerance = 1e-10)
  }
})

test_that("the pre-patch h > 0 correction fixture is unchanged", {
  result <- RAVEL:::.hidden_correction_trio_core(
    beta_vis = c(bD = 0.31, bm = 0.12, bf = -0.04),
    cov_cm = 0.58,
    cov_cf = 0.55,
    cov_mf = 0.25,
    h = 0.4,
    am_correction = "n_gen",
    n_generations = 5L
  )

  expect_equal(
    result$beta_total_assumed,
    c(bD = 0.351852954380001, bm = 0.100031277732837, bf = -0.0900526093223272),
    tolerance = 1e-12
  )
  expect_equal(
    result$var_total_assumed,
    c(child = 1.45018683476103, mother = 1.44999499933893, father = 1.44999499933893),
    tolerance = 1e-12
  )
  expect_equal(result$rho_T_corr_assumed, 0.364942925142842, tolerance = 1e-12)
  expect_equal(result$rho_total_assumed, 0.529165416501298, tolerance = 1e-12)
  expect_equal(
    result$eta_total_assumed,
    c(eta_m = 0.921494859910108, eta_f = 0.884318850137703),
    tolerance = 1e-12
  )
  expect_equal(
    unname(result$Sigma_xh),
    matrix(c(
      0.160503552546864, 0.176940978756283, 0.178211598550692,
      0.176481575024480, 0.180082423003538, 0.211765442978727,
      0.177748895823488, 0.211765442978727, 0.180082423003538
    ), 3, 3),
    tolerance = 1e-12
  )
  expect_equal(
    unname(result$M),
    matrix(c(
      0.881049870808248, 0.0976658424536943, 0.0988375506003557,
      -4.12782103275799e-17, 0.966097033738216, 0.122194525965315,
      -2.79142733045595e-17, 0.122194525965315, 0.966097033738216
    ), 3, 3),
    tolerance = 1e-12
  )
})

test_that("named uncertainty blocks are populated at h = 0", {
  dat <- ravel_simulate_trios(
    N = 1000,
    m = 30,
    h2_direct = 0.20,
    h2_mat = 0.05,
    h2_pat = 0,
    rho_MF = 0.30,
    PGS1_reliability = 0.60,
    PGS2_reliability = 0.40,
    hidden_fraction = 0,
    seed = 1901
  )
  result <- ravel_fit(
    dat,
    h = 0,
    am_correction = "n_gen",
    am_n_generations = 1L,
    n_boot = 12,
    seed = 1902
  )

  expect_true(result$converged)
  expect_false(is.null(result$lavaan_inference$total_raw))
  expect_false(is.null(result$lavaan_inference$total_std))
  expect_equal(result$lavaan_inference$total_std$estimate,
               result$lavaan_inference$visible$estimate, tolerance = 1e-8)
  expect_equal(result$lavaan_inference$total_std$cov,
               result$lavaan_inference$visible$cov, tolerance = 1e-6)
  expect_gte(result$bootstrap$n_successful, 10)
  expect_false(is.null(result$bootstrap$total_raw))
  expect_false(is.null(result$bootstrap$total_std))
  expect_equal(result$bootstrap$total_std$estimate,
               result$bootstrap$visible$estimate, tolerance = 1e-8)
})
