test_that("bootstrap extraction uses original fitted points for all three scales", {
  keys <- c("bD", "bm", "bf", "rhoMF")
  block <- list(se = setNames(rep(.02, 4L), keys),
    ci = matrix(rep(c(-.1, .5), 4), 2L, dimnames = list(c("lower", "upper"), keys)))
  result <- structure(list(converged = TRUE, h = .4,
    beta_vis = c(bD = .3, bm = .1, bf = .02), rho_vis = .1,
    beta_total_assumed = c(bD = .31, bm = .09, bf = .01), cov_total_assumed = .2,
    beta_total_assumed_std = c(bD = .4, bm = .12, bf = .015), rho_total_assumed = .15,
    bootstrap = list(source = "bootstrap", n_successful = 12L,
      visible = block, total_raw = block, total_std = block), outcome_variable = "trait"),
    class = c("ravel_results", "list"))
  before <- result
  got <- ravel_to_dataframe(result, uncertainty = "bootstrap", include_raw = TRUE,
                            include_vanilla = FALSE)
  expect_equal(got$estimate, c(.3, .1, .02, .1, .4, .12, .015, .15, .31, .09, .01, .2))
  expect_equal(got$se, rep(.02, 12L))
  expect_equal(got$ci_low, rep(-.1, 12L))
  expect_equal(got$ci_high, rep(.5, 12L))
  expect_equal(got$p_value, 2 * pnorm(-abs(got$estimate / got$se)))
  expect_true(all(got$uncertainty_source == "bootstrap (12 reps)"))
  expect_true(all(got$p_value_method == "normal_approx_from_bootstrap_se"))
  expect_true(all(got$outcome_variable == "trait"))
  expect_equal(ravel_to_dataframe(result, uncertainty = "auto", include_raw = TRUE,
                                  include_vanilla = FALSE), got)
  expect_identical(result, before)
})
