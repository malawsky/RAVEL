# Independent equilibrium Gaussian generator; no RAVEL simulation/correction
# helpers are used to construct its measurement or phenotype covariance.
score_diagnostic_data <- function(n = 1500L, seed = 9102L) {
  set.seed(seed)
  role <- c("child", "mother", "father")
  G <- MASS::mvrnorm(n, rep(0, 3), matrix(c(1, .65, .65, .65, 1, .3, .65, .3, 1), 3))
  K <- matrix(c(1, .5, .5, .5, 1, 0, .5, 0, 1), 3)
  d <- data.frame(FID = seq_len(n), M_ID = paste0("m", seq_len(n)), IID = paste0("i", seq_len(n)))
  for (j in 1:2) {
    X <- G + MASS::mvrnorm(n, rep(0, 3), c(.4, 1.3)[j] * K)
    d[paste0("PGS", j, "_", role)] <- X
  }
  d$sex <- rbinom(n, 1, .5)
  d$Y <- drop(G %*% c(.32, .13, .08)) + .2 * (d$sex - .5) + rnorm(n)
  d
}

score_diagnostic_cache <- new.env(parent = emptyenv())
score_diagnostic_example <- function() {
  if (!exists("data", score_diagnostic_cache, inherits = FALSE))
    score_diagnostic_cache$data <- score_diagnostic_data()
  score_diagnostic_cache$data
}

test_that("unadjusted refit runs despite covariates and preserves genetic restrictions", {
  d <- score_diagnostic_example()
  expect_warning(z <- ravel_score_diagnostic(d, fit_constraints = "equilibrium", pi = .4,
    covariates = "sex"), "Covariates were not used to compute the one-factor test.", fixed = TRUE)
  expect_identical(z$status, "ok", info = z$reason)
  expect_identical(z$auxiliary$covariates_omitted, "sex")
  expect_false(z$auxiliary$primary_model_test)
  expect_length(z$auxiliary_fit$fixed_effects, 0L)
  expect_identical(z$auxiliary_fit$constraint_spec$mode, "equilibrium")
  expect_equal(z$auxiliary_fit$h, .4)
  expect_equal(z$N, nrow(d))
  expect_length(unique(z$selection$clusters_used), nrow(d))
  expect_true(is.finite(z$p_value))
  expect_equal(z$df, 3)
  expect_output(print(z), "auxiliary unadjusted independent cross-sectional")
  score_diagnostic_cache$unadjusted <- z
})

test_that("primary fits supply settings and custom data maps without mutation", {
  d <- score_diagnostic_example()
  names(d)[names(d) == "Y"] <- "trait"
  names(d)[names(d) == "PGS1_child"] <- "score_child_A"
  primary <- suppressWarnings(ravel_fit(d, col_outcome = "trait", col_pgs1_child = "score_child_A",
    fixed_effects = "sex", fit_constraints = "equilibrium", pi = .3))
  expect_true(primary$converged)
  before <- serialize(primary, NULL)
  z <- suppressWarnings(ravel_score_diagnostic(d, primary_fit = primary))
  expect_identical(z$status, "ok", info = z$reason)
  expect_identical(serialize(primary, NULL), before)
  expect_equal(z$auxiliary$column_map[["Y"]], "trait")
  expect_equal(z$auxiliary$column_map[["PGS1_child"]], "score_child_A")
  expect_equal(z$auxiliary_fit$h, .3)
  expect_identical(z$auxiliary_fit$constraint_spec$mode, "equilibrium")
  expect_identical(z$auxiliary$covariates_omitted, "sex")
})

test_that("FID selection is reproducible under row order and preserves the RNG", {
  d <- score_diagnostic_example()
  d$FID <- (seq_len(nrow(d)) - 1L) %/% 2L
  d$M_ID <- paste0("m", d$FID)  # two individuals per maternal family
  d <- d[rep(seq_len(nrow(d)), each = 2L), ]
  d$age <- rep(c(8, 10), nrow(d) / 2L)
  d$Y <- d$Y + .05 * (d$age - 9)
  set.seed(214)
  rng <- .Random.seed
  z <- suppressWarnings(ravel_score_diagnostic(d, fit_constraints = "equilibrium", pi = .4,
    interactions = c("age", "sex"), age = "age", target_age = 9.8, selection_seed = 71L))
  expect_identical(.Random.seed, rng)
  expect_identical(z$status, "ok", info = z$reason)
  expect_equal(z$N, 750)
  expect_true(all(d$age[z$selection$rows_used] == 10))
  expect_length(unique(z$selection$clusters_used), 750)
  expect_null(z$auxiliary_fit$interactions)
  expect_setequal(z$auxiliary$interactions_omitted, c("age", "sex"))
  q <- d[rev(seq_len(nrow(d))), ]
  zz <- suppressWarnings(ravel_score_diagnostic(q, fit_constraints = "equilibrium", pi = .4,
    age = "age", target_age = 9.8, selection_seed = 71L))
  expect_identical(z$selection$individuals_used, zz$selection$individuals_used)
  expect_equal(z$statistic, zz$statistic, tolerance = 1e-7)
})

test_that("selection precedes completeness filtering and missing selected rows are not replaced", {
  d <- score_diagnostic_example()
  d$FID <- (seq_len(nrow(d)) - 1L) %/% 2L
  base <- suppressWarnings(ravel_score_diagnostic(d, fit_constraints = "equilibrium", pi = .4))
  expect_identical(base$status, "ok", info = base$reason)
  chosen <- base$selection$rows_selected[1:5]
  d$Y[chosen] <- NA_real_
  z <- suppressWarnings(ravel_score_diagnostic(d, fit_constraints = "equilibrium", pi = .4))
  expect_identical(z$selection$rows_selected, base$selection$rows_selected)
  expect_equal(z$selection$N_incomplete_selected, 5)
  expect_setequal(z$selection$rows_incomplete, chosen)
  expect_equal(z$N, base$N - 5L)
})

test_that("invalid family structure and ambiguous repeated records have explicit statuses", {
  d <- score_diagnostic_example()[1:80, ]
  d$M_ID[2] <- d$M_ID[1]
  expect_identical(ravel_score_diagnostic(d)$status, "invalid_design")
  d <- score_diagnostic_example()[1:80, ]
  d$FID <- rep(seq_len(40), each = 2)
  d$IID <- NULL
  expect_identical(ravel_score_diagnostic(d)$status, "needs_subject_id")
  d <- score_diagnostic_example()[1:80, ]
  d <- d[rep(seq_len(nrow(d)), each = 2), ]
  expect_identical(ravel_score_diagnostic(d)$status, "needs_occasion")
  d$age <- 10
  expect_identical(ravel_score_diagnostic(d, age = "age")$status, "ambiguous_occasion")
  d$visit <- rep(c(1, 2), 80)
  expect_warning(z <- ravel_score_diagnostic(d, occasion = "visit", occasion_value = 3),
    "0 complete observations", fixed = TRUE)
  expect_identical(z$status, "insufficient_data")
  expect_equal(z$selection$N_selected, 0)
  expect_equal(length(z$selection$clusters_without_selected_visit), 80)
})

test_that("the helper never substitutes a constrained model for none", {
  z <- suppressWarnings(ravel_score_diagnostic(score_diagnostic_example(), fit_constraints = "none", pi = .4))
  expect_identical(z$status, "unsupported_free_measurement", info = z$reason)
  expect_identical(z$auxiliary_fit$constraint_spec$mode, "none")
  expect_true(is.na(z$p_value))
})

test_that("public binary diagnostics report nominal p-values with calibration warnings", {
  skip_if_not(exists(".ravel_test_binary_score_signal", envir = asNamespace("RAVEL"), inherits = FALSE))
  d <- score_diagnostic_example()
  d$Y <- as.numeric(d$Y > 0)
  warnings <- character()
  z <- withCallingHandlers(ravel_score_diagnostic(d, outcome = "binary_liability",
    fit_constraints = "equilibrium", pi = .4, covariates = "sex"), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  expect_identical(z$auxiliary_fit$outcome, "binary_liability")
  expect_identical(z$status, "ok", info = z$reason)
  expect_identical(z$calibration_status, "uncalibrated")
  expect_equal(z$p_value, z$diagnostics$exploratory$nominal_p_value)
  expect_false(z$diagnostics$exploratory$calibrated)
  expect_true(is.finite(z$diagnostics$exploratory$nominal_p_value))
  expect_match(z$diagnostics$exploratory$label, "excess rejections", fixed = TRUE)
  expect_equal(z$statistic, z$diagnostics$exploratory$statistic)
  expect_equal(z$df, z$diagnostics$exploratory$df)
  expect_true(length(z$fits) > 0L)
  expect_true(any(grepl("8.0%", warnings, fixed = TRUE)))
  expect_true(any(grepl("Covariates were not used", warnings, fixed = TRUE)))
  expect_true(any(grepl("calibration warning", z$auxiliary$warnings, fixed = TRUE)))
  expect_match(z$test, "probit", ignore.case = TRUE)
  printed <- capture.output(print(z))
  expect_true(any(grepl("calibration warning", printed, fixed = TRUE)))
  expect_true(any(grepl("; p:", printed, fixed = TRUE)))
})
