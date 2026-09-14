test_that("bundled example has documented raw inputs and reproducible metadata", {
  env <- new.env()
  data("ravel_example", package = "RAVEL", envir = env)
  df <- env$ravel_example
  expect_s3_class(df, "data.frame")
  expect_equal(dim(df), c(4000L, 18L))
  expect_false(anyNA(df))
  expect_equal(anyDuplicated(df$FID), 0L)
  expect_setequal(unique(df$Y_binary), c(0L, 1L))
  expect_true(mean(df$Y_binary) > .07 && mean(df$Y_binary) < .13)
  expect_identical(levels(df$sex), c("female", "male"))
  expect_identical(levels(df$site), c("A", "B", "C"))
  scores <- grep("^pgs", names(df), value = TRUE)
  expect_length(scores, 6L)
  expect_true(all(vapply(df[scores], is.numeric, logical(1))))
  expect_true(all(vapply(df[scores], var, numeric(1)) > 1.5))
  expect_equal(attr(df, "simulation")$hidden_fraction, .4)
  expect_equal(attr(df, "simulation")$rho_total, .4)
  expect_equal(attr(df, "simulation")$seed, 20260913L)
})

test_that("installed example fits both outcomes and retains correction inputs", {
  env <- new.env()
  source(system.file("examples", "dataframe_workflow.R", package = "RAVEL"),
         local = env)
  for (out in list(env$out, env$out_binary)) {
    expect_true(all(out$model_status$status == "ok"))
    expect_true(all(is.finite(out$summary$estimate)))
    expect_true(all(is.finite(out$summary$se)))
    expect_setequal(unique(out$summary$regression), c("RAVEL", "PGS1", "PGS2"))
    expect_equal(nrow(out$covariate_coding), 10L)
    expect_equal(unname(diag(out$pgs_covariance)), out$pgs_variances$variance)
    for (packet in out$packets) {
      expect_equal(dim(packet$vcov_core), c(6L, 6L))
      expect_true(nzchar(packet$settings$model_syntax))
      expect_true(nrow(packet$parameter_estimates) > 0L)
    }
    for (version in 1:2) {
      child <- paste0("pgs", version, "_child")
      for (role in c("child", "mother", "father")) {
        name <- paste0("PGS", version, "_", role)
        raw_name <- paste0("pgs", version, "_", role)
        actual <- out$pgs_variances$variance[out$pgs_variances$model_column == name]
        expect_equal(actual, var(env$df[[raw_name]]) / var(env$df[[child]]))
      }
    }
    corrected <- ravel_correct_grid(out, hidden_fraction = c(.4, .5),
      n_generations = c(1, 15, Inf), covariate_mapping_confirmed = TRUE)
    expect_equal(nrow(corrected), 96L)
    expect_true(all(is.finite(corrected$estimate)))
    expect_true(all(is.finite(corrected$se)))
    saved <- tempfile(fileext = ".rds")
    saveRDS(out, saved)
    restored <- ravel_correct_grid(saved, hidden_fraction = .4, n_generations = 1,
                                  covariate_mapping_confirmed = TRUE)
    direct <- ravel_correct_grid(out, hidden_fraction = .4, n_generations = 1,
                                covariate_mapping_confirmed = TRUE)
    expect_identical(restored, direct)
    unlink(saved)
  }
})
