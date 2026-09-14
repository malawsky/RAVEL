test_that("covariate model syntax frees all genetic and covariate covariances", {
  for (covariates in list("PC1", c("PC1", "age", "sex"))) {
    model <- RAVEL:::.build_model(fixed_effects = covariates)
    pt <- lavaan::lavaanify(model, model.type = "sem", std.lv = TRUE,
                           fixed.x = FALSE, auto = TRUE)
    for (g in c("Gc", "Gm", "Gf")) {
      for (w in covariates) {
        hit <- pt$op == "~~" & ((pt$lhs == g & pt$rhs == w) | (pt$lhs == w & pt$rhs == g))
        expect_equal(sum(hit & pt$free > 0), 1L, info = paste(g, "~~", w))
      }
    }
    if (length(covariates) > 1L) {
      pairs <- combn(covariates, 2)
      for (j in seq_len(ncol(pairs))) {
        hit <- pt$op == "~~" & ((pt$lhs == pairs[1, j] & pt$rhs == pairs[2, j]) |
                                 (pt$rhs == pairs[1, j] & pt$lhs == pairs[2, j]))
        expect_equal(sum(hit & pt$free > 0), 1L)
      }
    }
    expect_equal(sum(pt$op == "~" & pt$lhs == "Y" & pt$rhs %in% covariates), length(covariates))
  }
})

test_that("NULL and empty covariate lists produce the same model", {
  expect_identical(RAVEL:::.build_model(fixed_effects = NULL),
                   RAVEL:::.build_model(fixed_effects = character()))
})

test_that("RAVEL and independent reference models recover known coefficients on the correct scale", {
  for (outcome in c("continuous", "binary_liability")) {
    for (correlated in c(FALSE, TRUE)) {
      fixture <- cov_test_fixture(outcome, correlated, n = 16000L)
      ref <- cov_test_fit(fixture, "reference")
      result <- cov_test_fit(fixture)
      cov_test_expect_reference(result, ref)
      if (cov_test_converged(ref)) {
        b <- lavaan::coef(ref$fit)[names(fixture$beta)]
        se <- sqrt(diag(lavaan::vcov(ref$fit))[names(b)])
        expect_true(all(is.finite(se) & se > 0))
        expect_lt(max(abs((b - fixture$beta) / se)), 5)
        expect_lt(max(abs(b - fixture$beta)), .06)
        expect_true(isTRUE(lavaan::lavInspect(ref$fit, "post.check")))
        if (outcome == "binary_liability") {
          expect_equal(lavaan::lavInspect(ref$fit, "cov.ov")["Y", "Y"], 1, tolerance = 1e-8)
          expect_equal(unname(result$beta_vis), unname(b), tolerance = 5e-5)
          expect_equal(unname(result$lavaan_inference$visible$se[seq_len(3L)]),
                       unname(se), tolerance = 5e-5)
        }
      }
    }
  }
})

test_that("the historically unstable binary case matches independent moments with simple starts", {
  fixture <- cov_test_fixture(n = 4000L)
  result <- cov_test_fit(fixture)
  if (cov_test_converged(result)) {
    reference <- cov_test_reference_fit(fixture, start = "simple")
    cov_test_expect_reference(result, list(fit = reference,
      converged = isTRUE(lavaan::lavInspect(reference, "converged"))))
    expect_true(all(is.finite(result$lavaan_inference$visible$se)))
  }
})

test_that("a covariate fit cannot return estimates when both numerical attempts are inadmissible", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  calls <- 0L
  testthat::local_mocked_bindings(
    sem = function(...) { calls <<- calls + 1L; structure(list(), class = "invalid_test_fit") },
    lavInspect = function(object, what, ...) switch(what, converged = TRUE, post.check = FALSE),
    .package = "lavaan"
  )
  fixture <- cov_test_fixture()
  expect_warning(result <- RAVEL::ravel_fit(fixture$dat, outcome = "binary_liability",
    fixed_effects = fixture$covariates, fixed_crossperson_corr = cov_test_fixed_corr),
    "admissible converged fit")
  expect_identical(calls, 2L)
  expect_false(result$converged)
  expect_null(result$beta_vis)
})

for (outcome in c("continuous", "binary_liability")) {
  local({
    type <- outcome
    test_that(paste(type, "without covariates preserves the original fit and SEs"), {
      fixture <- cov_test_fixture(type, with_covariates = FALSE)
      result <- cov_test_fit(fixture)
      reference <- cov_test_fit(fixture, "reference")
      cov_test_expect_reference(result, reference)
      if (isTRUE(result$converged)) {
        expect_equal(result$R2_model, cov_test_genetic_fraction(reference$fit), tolerance = 1e-5)
        expect_equal(result$beta_total_assumed_std, result$beta_vis, tolerance = 1e-8)
        expect_equal(result$lavaan_inference$total_std$cov,
                     result$lavaan_inference$visible$cov, tolerance = 1e-6)
      }
    })

    for (correlated in c(FALSE, TRUE)) {
      local({
        correlated_case <- correlated
        label <- if (correlated_case) "correlated" else "independent"
        test_that(paste(type, "with", label, "covariates matches the reference SEM"), {
          fixture <- cov_test_fixture(type, correlated_case)
          result <- cov_test_fit(fixture)
          reference <- cov_test_fit(fixture, "reference")
          cov_test_expect_reference(result, reference)
          if (isTRUE(result$converged)) {
            options <- lavaan::lavInspect(result$fit, "options")
            expect_false(options$fixed.x)
            expect_false(options$conditional.x)
            expect_identical(result$fixed_effects, fixture$covariates)
            if (type == "binary_liability") {
              expect_identical(lavaan::lavNames(result$fit, "ov.ord"), "Y")
              expect_equal(options$parameterization, "delta")
            }
          }
        })
      })
    }

    test_that(paste(type, "genetic variance fraction includes covariates in its denominator"), {
      fixture <- cov_test_fixture(type)
      result <- cov_test_fit(fixture)
      reference <- cov_test_fit(fixture, "reference")
      if (cov_test_converged(result)) {
        expect_equal(result$R2_model, cov_test_genetic_fraction(reference$fit), tolerance = 1e-5)
        if (type == "binary_liability") {
          expect_equal(result$R2_total_assumed, result$R2_model, tolerance = 1e-8)
        } else {
          expect_equal(result$R2_total_assumed,
                       result$s2_genetic_total_assumed / var(fixture$dat$Y), tolerance = 1e-8)
        }
      }
    })

    test_that(paste(type, "visible SEs, confidence intervals and p-values use the joint fit"), {
      fixture <- cov_test_fixture(type)
      result <- cov_test_fit(fixture)
      if (cov_test_converged(result)) {
        labels <- c("bD", "bm", "bf", "rhoMF")
        V <- lavaan::vcov(result$fit)[labels, labels]
        b <- lavaan::coef(result$fit)[labels]
        se <- sqrt(diag(V))
        expect_equal(result$lavaan_inference$visible$cov, V, tolerance = 1e-6)
        tbl <- RAVEL::ravel_to_dataframe(result, uncertainty = "lavaan", include_vanilla = FALSE)
        expect_equal(nrow(tbl), 4L)
        expect_equal(tbl$estimate, unname(b), tolerance = 1e-8)
        expect_equal(tbl$se, unname(se), tolerance = 1e-6)
        expect_equal(tbl$ci_low, unname(b - qnorm(.975) * se), tolerance = 1e-6)
        expect_equal(tbl$ci_high, unname(b + qnorm(.975) * se), tolerance = 1e-6)
        expect_equal(tbl$p_value, unname(2 * pnorm(-abs(b / se))), tolerance = 1e-6)
        expect_equal(result$lavaan_inference$total_std$estimate,
                     result$lavaan_inference$visible$estimate, tolerance = 1e-8)
        expect_equal(result$lavaan_inference$total_std$cov, V, tolerance = 1e-6)
      }
    })

    test_that(paste(type, "PGS comparison remains an ordinary covariate-adjusted regression"), {
      fixture <- cov_test_fixture(type)
      result <- cov_test_fit(fixture)
      if (cov_test_converged(result)) {
        for (k in 1:2) {
          d <- fixture$dat
          cols <- paste0("PGS", k, "_", c("child", "mother", "father"))
          d[cols] <- lapply(d[cols], function(x) as.numeric(scale(x)))
          if (type == "continuous") d$Y <- as.numeric(scale(d$Y))
          formula <- reformulate(c(cols, fixture$covariates), response = "Y")
          fit <- if (type == "continuous") lm(formula, d) else glm(formula, d, family = binomial("probit"))
          expected <- summary(fit)$coefficients[cols, ]
          got <- result$vanilla_regression[[paste0("PGS", k)]]
          expect_equal(unname(got$coefficients), unname(expected[, 1]), tolerance = 1e-10)
          expect_equal(unname(got$se), unname(expected[, 2]), tolerance = 1e-10)
          expect_equal(unname(got$p_value), unname(expected[, 4]), tolerance = 1e-10)
        }
      }
    })
  })
}

test_that("binary fitting retains the fixed shared-weight residual correlations", {
  result <- cov_test_fit(cov_test_fixture())
  if (cov_test_converged(result)) {
    expect_equal(result$crossperson_resid_corr, cov_test_fixed_corr, tolerance = 1e-6)
    pt <- lavaan::parTable(result$fit)
    expect_equal(sum(pt$op == "=~" & pt$free > 0), 6L)
    expect_equal(sum(pt$op == "=="), 6L)
    covariates <- c("PC1", "age", "sex")
    for (g in c("Gc", "Gm", "Gf")) {
      for (w in covariates) {
        hit <- pt$op == "~~" & ((pt$lhs == g & pt$rhs == w) | (pt$rhs == g & pt$lhs == w))
        expect_equal(sum(hit & pt$free > 0), 1L)
      }
    }
    expect_false(any(pt$op == "~~" & grepl("^PGS", pt$lhs) & pt$rhs %in% covariates & pt$free > 0))
  }
})

test_that("binary covariates support one predictor and explicit categorical dummy columns", {
  for (covariates in list("sex", c("PC1", "group_B", "group_C"))) {
    fixture <- cov_test_fixture()
    fixture$id <- paste0(fixture$id, paste(covariates, collapse = "_"))
    group <- rep(c("A", "B", "C"), length.out = nrow(fixture$dat))
    fixture$dat$group_B <- as.numeric(group == "B")
    fixture$dat$group_C <- as.numeric(group == "C")
    fixture$covariates <- covariates
    cov_test_expect_reference(cov_test_fit(fixture), cov_test_fit(fixture, "reference"))
  }
})

test_that("binary covariate order and positive PGS rescaling do not change genetic coefficients", {
  fixture <- cov_test_fixture()
  result <- cov_test_fit(fixture)
  altered <- fixture
  altered$id <- paste0(fixture$id, "_reordered_rescaled")
  altered$covariates <- rev(fixture$covariates)
  for (k in 1:2) {
    cols <- paste0("PGS", k, "_", c("child", "mother", "father"))
    altered$dat[cols] <- lapply(altered$dat[cols], function(x) x * c(2, .7)[k])
  }
  original_data <- altered$dat
  got <- cov_test_fit(altered)
  if (cov_test_converged(result) && cov_test_converged(got)) {
    expect_true(isTRUE(lavaan::lavInspect(got$fit, "post.check")))
    expect_equal(got$beta_vis, result$beta_vis, tolerance = 1e-4)
    expect_equal(got$lavaan_inference$visible$cov,
                 result$lavaan_inference$visible$cov, tolerance = 1e-4)
    expect_equal(got$R2_model, result$R2_model, tolerance = 1e-4)
    expect_identical(altered$dat, original_data)
  }
})

test_that("binary missingness uses complete cases including the covariates", {
  fixture <- cov_test_fixture(n = 16000L)
  fixture$id <- paste0(fixture$id, "_missing")
  fixture$dat$PC1[1:10] <- NA_real_
  fixture$dat$PGS1_father[11:20] <- NA_real_
  fixture$dat$Y[21:30] <- NA_real_
  result <- cov_test_fit(fixture)
  complete <- fixture
  complete$id <- paste0(fixture$id, "_complete")
  complete$dat <- complete$dat[complete.cases(complete$dat), ]
  reference <- cov_test_fit(complete, "reference")
  cov_test_expect_reference(result, reference)
  if (isTRUE(result$converged)) {
    expect_equal(as.numeric(lavaan::lavInspect(result$fit, "nobs")), nrow(complete$dat))
  }
})

test_that("reversing the binary outcome reverses coefficients but not their SEs or genetic fraction", {
  fixture <- cov_test_fixture()
  result <- cov_test_fit(fixture)
  reversed <- fixture
  reversed$id <- paste0(fixture$id, "_outcome_reversed")
  reversed$dat$Y <- 1 - fixture$dat$Y
  got <- cov_test_fit(reversed)
  if (cov_test_converged(result) && cov_test_converged(got)) {
    expect_equal(got$beta_vis, -result$beta_vis, tolerance = 1e-4)
    expect_equal(got$latent_cov_vis, result$latent_cov_vis, tolerance = 1e-4)
    expect_equal(got$lavaan_inference$visible$se,
                 result$lavaan_inference$visible$se, tolerance = 1e-4)
    expect_equal(got$R2_model, result$R2_model, tolerance = 1e-4)
  }
})

test_that("binary joint fitting preserves the alternative residual and parental-path options", {
  fixture <- cov_test_fixture()
  for (crossperson in c("free", "equal_pc")) {
    result <- cov_test_fit(fixture, crossperson = crossperson, fixed_corr = FALSE)
    reference <- cov_test_fit(fixture, "reference", crossperson = crossperson, fixed_corr = FALSE)
    cov_test_expect_reference(result, reference)
  }
  result <- cov_test_fit(fixture, fix_parent_paths = TRUE)
  reference <- cov_test_fit(fixture, "reference", fix_parent_paths = TRUE)
  cov_test_expect_reference(result, reference)
  if (isTRUE(result$converged)) expect_equal(unname(result$beta_vis[c("bm", "bf")]), c(0, 0))
})

test_that("binary bootstrap refits keep the covariates and joint model", {
  fixture <- cov_test_fixture()
  result <- cov_test_fit(fixture, n_boot = 12L)
  if (cov_test_converged(result)) {
    expect_false(is.null(result$bootstrap))
    if (!is.null(result$bootstrap)) {
      expect_gte(result$bootstrap$n_successful, 10L)
      set.seed(99173L)
      draws <- replicate(12L, sample(seq_len(nrow(fixture$dat)), nrow(fixture$dat), replace = TRUE), simplify = FALSE)
      b <- matrix(NA_real_, 12, 4, dimnames = list(NULL, c("bD", "bm", "bf", "rhoMF")))
      for (i in seq_along(draws)) {
        replicate_fixture <- fixture
        replicate_fixture$dat <- fixture$dat[draws[[i]], ]
        fit <- cov_test_reference_fit(replicate_fixture)
        V <- lavaan::vcov(fit)
        usable <- isTRUE(lavaan::lavInspect(fit, "converged")) &&
          isTRUE(lavaan::lavInspect(fit, "post.check")) && is.matrix(V) &&
          all(is.finite(V)) && all(diag(V)[colnames(b)] > 0)
        if (usable) b[i, ] <- lavaan::coef(fit)[colnames(b)]
      }
      expect_equal(result$bootstrap$n_successful, sum(complete.cases(b)))
      expect_equal(result$bootstrap$visible$se, apply(b, 2, sd, na.rm = TRUE), tolerance = 1e-4)
      ci <- apply(b, 2, quantile, probs = c(.025, .975), na.rm = TRUE, names = FALSE)
      expect_equal(unname(result$bootstrap$visible$ci), unname(ci), tolerance = 1e-4)
      expect_equal(result$bootstrap$total_std$se, result$bootstrap$visible$se, tolerance = 1e-6)
      tbl <- RAVEL::ravel_to_dataframe(result, uncertainty = "bootstrap", include_vanilla = FALSE)
      expect_equal(tbl$se, unname(result$bootstrap$visible$se), tolerance = 1e-8)
    }
  }
})

test_that("hidden-correction inference remains finite with no covariates", {
  fixture <- cov_test_fixture(with_covariates = FALSE)
  for (regime in c("n_gen", "equilibrium")) {
    result <- cov_test_fit(fixture, h = .4, regime = regime)
    if (cov_test_converged(result)) {
      for (block in c("visible", "total_raw", "total_std")) {
        inference <- result$lavaan_inference[[block]]
        expect_true(all(is.finite(inference$estimate)))
        expect_true(all(is.finite(inference$se) & inference$se > 0))
      }
      expected_fraction <- result$s2_genetic_total_assumed /
        lavaan::lavInspect(result$fit, "cov.ov")["Y", "Y"]
      expect_equal(result$R2_total_assumed, expected_fraction, tolerance = 1e-8)
    }
  }
})
