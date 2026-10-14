.vpg_fixture <- function(binary = FALSE) {
  set.seed(481)
  n <- 800
  d <- data.frame(adjustment = rnorm(n))
  for (k in 1:2) {
    m <- rnorm(n); f <- .4 * m + rnorm(n)
    x <- list(child = (m + f) / 2 + rnorm(n), mother = m, father = f,
      mat_sib = .6 * m + rnorm(n), pat_sib = .6 * f + rnorm(n))
    for (role in names(x)) d[[paste0("PGS", k, "_", role)]] <- x[[role]]
  }
  d$Y <- .3 * d$PGS1_child + .4 * d$PGS1_mother - .1 * d$PGS1_father +
    .2 * d$PGS1_mat_sib + .1 * d$adjustment + rnorm(n)
  if (binary) d$Y <- as.numeric(d$Y > 0)
  d$PGS1_mat_sib[1:30] <- NA_real_
  d$PGS2_mat_sib[1:55] <- NA_real_
  d$PGS1_pat_sib[31:70] <- NA_real_
  d$PGS2_pat_sib[56:100] <- NA_real_
  d
}

.vpg_check_refit <- function(tab, fit, mother, other, parameter, level = .9) {
  d <- stats::model.frame(fit)
  d$midpoint <- (d[[mother]] + d[[other]]) / 2
  d$deviation <- (d[[mother]] - d[[other]]) / 2
  rhs <- c(setdiff(names(d)[-1], c(mother, other, "midpoint", "deviation")), "midpoint", "deviation")
  formula <- stats::reformulate(rhs, response = names(d)[1])
  ref <- if (inherits(fit, "glm")) stats::glm(formula, d, family = binomial("probit")) else stats::lm(formula, d)
  expected <- summary(ref)$coefficients["deviation", ]
  row <- tab[tab$parameter == parameter, ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$estimate, unname(expected[1]), tolerance = 1e-9)
  expect_equal(row$se, unname(expected[2]), tolerance = 1e-8)
  expect_equal(row$p_value, unname(expected[4]), tolerance = 1e-8)
  critical <- if (inherits(fit, "glm")) qnorm(1 - (1 - level) / 2) else qt(1 - (1 - level) / 2, ref$df.residual)
  expect_equal(c(row$ci_low, row$ci_high), unname(expected[1] + c(-1, 1) * critical * expected[2]), tolerance = 1e-8)
}

test_that("vanilla contrasts match reparameterized regressions for either sibling side", {
  for (binary in c(FALSE, TRUE)) for (mode in c("none", "maternal", "paternal", "both")) {
    d <- .vpg_fixture(binary)
    out <- .fit_vanilla_pgs_regression(d, if (binary) "binary" else "continuous", "adjustment", .9, mode)
    for (version in c("PGS1", "PGS2")) {
      x <- out[[version]]
      expect_true(x$converged)
      .vpg_check_refit(x$parent_contrasts, x$fit, ".pgs_mother", ".pgs_father", "delta_M_mid")
      expect_equal(x$parent_contrasts$estimate[2], -x$parent_contrasts$estimate[3])
      expect_equal(x$parent_contrasts$se[2], x$parent_contrasts$se[3])
      expected <- c("beta_child", "beta_mother", "beta_father")
      if (mode %in% c("maternal", "both")) {
        expected <- c(expected, "beta_MS", "delta_M")
        .vpg_check_refit(x$table, x$fit, ".pgs_mother", ".pgs_mat_sib", "delta_M")
      }
      if (mode %in% c("paternal", "both")) {
        expected <- c(expected, "beta_FS", "delta_F")
        .vpg_check_refit(x$table, x$fit, ".pgs_father", ".pgs_pat_sib", "delta_F")
      }
      expect_setequal(x$table$parameter, expected)
      expect_equal(x$N, nrow(model.frame(x$fit)))
    }
    if (mode == "both") {
      expect_equal(out$PGS1$N, 730)
      expect_equal(out$PGS2$N, 700)
    }
  }
})

test_that("non-estimable vanilla contrasts are unavailable rather than falsely precise", {
  d <- .vpg_fixture()
  d$PGS1_father <- d$PGS1_mother
  out <- .fit_vanilla_pgs_regression(d, "continuous", character(), .95)
  expect_true(all(is.na(out$PGS1$parent_contrasts$se)))
  d$PGS1_mat_sib <- NA_real_
  failed <- .fit_vanilla_pgs_regression(d, "continuous", character(), .95, "maternal")
  expect_false(failed$PGS1$converged)
  expect_true(all(is.na(failed$PGS1$table$se)))
  expect_true("delta_M" %in% failed$PGS1$table$parameter)
})

test_that("report vanilla contrasts preserve units and use actual complete-case N", {
  for (binary in c(FALSE, TRUE)) {
    d <- .vpg_fixture(binary)
    cols <- grep("^PGS", names(d), value = TRUE)
    scaling <- data.frame(model_column = cols, child_reference_sd = 1, child_reference_mean = 0,
      scaled_variance = vapply(d[cols], var, numeric(1), na.rm = TRUE),
      scaled_mean = vapply(d[cols], mean, numeric(1), na.rm = TRUE))
    for (version in 1:2) {
      out <- .ravel_report_vanilla(d, version, "adjustment", scaling, 1, 0,
        if (binary) "binary_liability" else "continuous", .9, "Y", "both", TRUE)
      cols_one <- paste0("PGS", version, "_", c("child", "mother", "father", "mat_sib", "pat_sib"))
      form <- reformulate(c(cols_one, "adjustment"), response = "Y")
      fit <- if (binary) glm(form, d, family = binomial("probit")) else lm(form, d)
      tab <- subset(out$coefficients, scale == "raw_PGS_units")
      .vpg_check_refit(tab, fit, cols_one[2], cols_one[3], "delta_M_mid")
      .vpg_check_refit(tab, fit, cols_one[2], cols_one[4], "delta_M")
      expect_true(all(out$coefficients$N == nobs(fit)))
      expect_equal(out$fit_statistics[["N"]], nobs(fit))
      own <- d
      own[cols_one] <- lapply(own[cols_one], function(x) as.numeric(scale(x)))
      ref <- if (binary) glm(form, own, family = binomial("probit")) else lm(form, own)
      tab <- subset(out$coefficients, scale == "own_role_PGS_SD" &
        outcome_scale != "model_implied_marginal_probit_liability_SD")
      .vpg_check_refit(tab, ref, cols_one[2], cols_one[3], "delta_M_mid")
      .vpg_check_refit(tab, ref, cols_one[3], cols_one[5], "delta_F")
    }
  }
})

.vpg_sem_cache <- new.env(parent = emptyenv())
.vpg_sem_data <- function() {
  if (!exists("data", .vpg_sem_cache, inherits = FALSE)) {
    d <- ravel_simulate_trios(N = 3000, m = 80, parental_sibs = "both", rho_MF = .2, seed = 48)
    d$PGS1_mat_sib[1:100] <- NA_real_
    d$PGS2_mat_sib[1:160] <- NA_real_
    .vpg_sem_cache$data <- d
  }
  .vpg_sem_cache$data
}

test_that("live tables and printers expose both contrast families and model-specific N", {
  fit <- suppressWarnings(ravel_fit(.vpg_sem_data(), parental_sibs = "both"))
  expect_true(fit$converged)
  tab <- as.data.frame(fit, parent_contrasts = TRUE)
  for (reg in c("RAVEL", "PGS1", "PGS2")) {
    rows <- tab[tab$regression == reg, ]
    expect_true(all(c("delta_M", "delta_F", "delta_M_mid", "delta_P_mid", "beta_mid") %in% rows$parameter))
    expect_true(all(rows$N == c(RAVEL = 3000, PGS1 = 2900, PGS2 = 2840)[[reg]]))
  }
  expect_false("delta_M_mid" %in% as.data.frame(fit)$parameter)
  expect_output(print_simple(fit, parent_contrasts = TRUE), "PGS1 delta_M_mid")
  expect_output(print_ravel_results(fit), "PGS2 delta_M")
  plot <- ravel_plot_report(tab, scale = "vanilla")
  expect_true(all(c("beta_D", "beta_M", "beta_P", "delta_M", "delta_F", "delta_M_mid") %in% plot$data$parameter))
  old <- fit
  old$vanilla_regression$parent_contrasts <- NULL
  expect_equal(as.data.frame(old, parent_contrasts = TRUE), tab)
})

test_that("dataframe report exports and bootstraps both vanilla contrast families", {
  d <- .vpg_sem_data()
  cols <- grep("^PGS", names(d), value = TRUE)
  out <- ravel_dataframe_report(d, pgs_columns = setNames(cols, cols), parental_sibs = "both",
    residual_models = "fixed", parent_contrasts = TRUE, make_plot = FALSE,
    n_boot = 10L, seed = 94)
  expect_true(all(out$model_status$status == "ok"))
  for (reg in c("RAVEL", "PGS1", "PGS2")) {
    tab <- subset(out$summary, regression == reg)
    expect_true(all(c("delta_M", "delta_F", "delta_M_mid", "beta_mid") %in% tab$parameter))
  }
  for (reg in c("PGS1", "PGS2")) {
    tab <- out$vanilla[[reg]]$coefficients
    draws <- out$bootstrap$vanilla_draws[[reg]]
    key <- .ravel_report_row_key(tab)
    for (sc in unique(tab$scale)) {
      rows <- which(tab$scale == sc & tab$parameter == "delta_M")
      parent <- which(tab$scale == sc & tab$parameter == "beta_M")
      sibling <- which(tab$scale == sc & tab$parameter == "beta_MS")
      expect_equal(draws[, key[rows]], draws[, key[parent]] - draws[, key[sibling]], tolerance = 1e-10)
    }
    expect_true(all(is.finite(out$vanilla[[reg]]$coefficients_bootstrap$se)))
  }
  plot <- ravel_plot_report(out)
  expect_true(all(c("beta_D", "beta_M", "beta_P", "delta_M", "delta_F", "delta_M_mid") %in% plot$data$parameter))
  expect_true(out$packet$settings$parent_contrasts)
})
