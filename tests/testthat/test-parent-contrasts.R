pc_cache <- new.env(parent = emptyenv())

pc_data <- function() {
  if (!exists("dat", pc_cache, inherits = FALSE))
    pc_cache$dat <- ravel_simulate_trios(N = 12000, m = 200, h2_direct = .09, h2_mat = .0225,
      h2_pat = .0025, rho_MF = .3, seed = 12)
  pc_cache$dat
}

pc_rows <- function(fit, ...) {
  df <- as.data.frame(fit, include_vanilla = FALSE, parent_contrasts = TRUE, ...)
  df[df$parameter %in% c("beta_mid", "delta_M_mid", "delta_P_mid"), ]
}

test_that("parental contrasts are opt-in", {
  fit <- suppressWarnings(ravel_fit(pc_data()))
  expect_false(any(c("beta_mid", "delta_M_mid") %in% as.data.frame(fit)$parameter))
  expect_false(any(grepl("Parental midpoint", capture.output(print_simple(fit)))))
  expect_output(print_simple(fit, parent_contrasts = TRUE), "Parental midpoint and deviations")
  expect_output(print_ravel_results(fit, parent_contrasts = TRUE), "maternal deviation")
})

test_that("visible contrasts use the joint covariance and match the likelihood-ratio test", {
  d <- pc_data()
  fit <- suppressWarnings(ravel_fit(d))
  rows <- pc_rows(fit)
  expect_identical(rows$scale, rep("visible", 3))
  b <- fit$beta_vis
  V <- lavaan::vcov(fit$fit)[c("bm", "bf"), c("bm", "bf")]
  w <- c(1, -1)
  expect_equal(rows$estimate, c(b[["bm"]] + b[["bf"]], b[["bm"]] - b[["bf"]], b[["bf"]] - b[["bm"]]),
    tolerance = 1e-10)
  expect_equal(rows$se[2], sqrt(drop(t(w) %*% V %*% w)), tolerance = 1e-8)
  expect_equal(rows$se[1], sqrt(drop(t(c(1, 1)) %*% V %*% c(1, 1))), tolerance = 1e-8)
  expect_equal(rows$p_value[2], rows$p_value[3])
  # Same test as constraining bm == bf in the model (ML, normal-theory).
  f1 <- lavaan::sem(fit$model, data = d, std.lv = TRUE, fixed.x = FALSE, estimator = "ML", missing = "fiml")
  f0 <- lavaan::sem(paste(fit$model, "bm == bf", sep = "\n"), data = d, std.lv = TRUE, fixed.x = FALSE,
    estimator = "ML", missing = "fiml")
  lr <- lavaan::lavTestLRT(f1, f0)[2, "Chisq diff"]
  pe <- lavaan::parameterEstimates(f1)
  Vml <- lavaan::vcov(f1)[c("bm", "bf"), c("bm", "bf")]
  z2 <- (pe$est[pe$label == "bm"][1] - pe$est[pe$label == "bf"][1])^2 / drop(t(w) %*% Vml %*% w)
  expect_equal(lr, z2, tolerance = .1)
})

test_that("corrected contrasts use the corrected covariance and external pi", {
  d <- pc_data()
  fit <- suppressWarnings(ravel_fit(d, pi = .3))
  rows <- pc_rows(fit, include_raw = TRUE)
  expect_setequal(unique(rows$scale), c("visible", "corrected_std", "corrected_raw"))
  C <- fit$lavaan_inference$total_std$cov[c("bm", "bf"), c("bm", "bf")]
  std <- rows[rows$scale == "corrected_std", ]
  expect_equal(std$estimate[2], fit$beta_total_assumed_std[["bm"]] - fit$beta_total_assumed_std[["bf"]],
    tolerance = 1e-10)
  expect_equal(std$se[2], sqrt(C[1, 1] + C[2, 2] - 2 * C[1, 2]), tolerance = 1e-8)
  withpi <- suppressWarnings(ravel_fit(d, pi = .3, pi_se = .05))
  std_pi <- pc_rows(withpi)
  std_pi <- std_pi[std_pi$scale == "corrected_std", ]
  expect_equal(std_pi$estimate, std$estimate, tolerance = 1e-10)
  expect_true(all(std_pi$se >= std$se - 1e-10))
  expect_true(all(grepl("pi", std_pi$uncertainty_source)))
})

test_that("bootstrap contrasts come from the joint draws", {
  fit <- suppressWarnings(ravel_fit(pc_data()[1:4000, ], pi = .3, n_boot = 20, seed = 1))
  rows <- pc_rows(fit, uncertainty = "bootstrap")
  draws <- fit$bootstrap$core_draws
  ok <- stats::complete.cases(draws)
  expect_equal(rows$se[rows$scale == "visible"][2], stats::sd(draws[ok, "bm"] - draws[ok, "bf"]), tolerance = 1e-10)
  td <- fit$bootstrap$total_std_draws
  ok <- stats::complete.cases(td)
  expect_equal(rows$se[rows$scale == "corrected_std"][1], stats::sd(td[ok, "bm"] + td[ok, "bf"]), tolerance = 1e-10)
  expect_true(all(grepl("bootstrap", rows$uncertainty_source)))
})
