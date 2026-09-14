test_that("report inference options and output paths are validated before fitting", {
  expect_error(ravel_dataframe_report(data.frame(), n_boot = 1), "at least 10")
  expect_error(ravel_dataframe_report(data.frame(), uncertainty = "bootstrap"), "requires")
  expect_error(ravel_dataframe_report(data.frame(), n_boot = NA_real_), "n_boot")
  expect_error(ravel_dataframe_report(data.frame(), seed = -1), "seed")
  expect_error(ravel_dataframe_report(data.frame(), report_name = "../escape", out_dir = tempdir()), "safe folder")
  expect_error(ravel_dataframe_report(data.frame(), make_plot = NA), "make_plot")
})

test_that("bootstrap SEs and intervals use draws without changing point estimates", {
  draws <- cbind(seq(.1, .3, length.out = 20), seq(-.1, .1, length.out = 20))
  tab <- data.frame(estimate = c(.19, .01), se = 1, ci_low = 0, ci_high = 1,
                    p_value = .5, inference = "model_based")
  result <- RAVEL:::.ravel_report_boot_table(tab, draws, .95, 16)
  expect_identical(result$estimate, tab$estimate)
  expect_equal(result$se, apply(draws, 2, sd), ignore_attr = TRUE)
  expect_equal(result$ci_low, apply(draws, 2, quantile, .025), ignore_attr = TRUE)
  expect_equal(result$ci_high, apply(draws, 2, quantile, .975), ignore_attr = TRUE)
  expect_equal(result$p_value, 2 * pnorm(-abs(tab$estimate / result$se)))
  expect_equal(result$n_boot_success, c(20L, 20L))
  draws[1:5, 2] <- NA
  failed <- RAVEL:::.ravel_report_boot_table(tab, draws, .95, 16)
  expect_equal(failed$uncertainty_status, c("ok", "insufficient_bootstrap"))
  expect_true(all(is.na(failed[2, c("se", "ci_low", "ci_high", "p_value")])))
  expect_identical(failed$estimate, tab$estimate)
})

test_that("detailed visible packet labels follow the chosen inference", {
  visible <- data.frame(parameter = c("beta_D", "beta_M"), estimate = c(.2, .1),
    se = .1, uncertainty_source = "lavaan/delta-method",
    p_value_method = "normal_approx_from_lavaan_se", inference = "lavaan_delta_wald")
  selected <- data.frame(parameter = c("beta_M", "beta_D"), se = c(.03, .02),
    inference = "bootstrap_percentile_CI_normal_p", n_boot_success = 20L)
  result <- RAVEL:::.ravel_report_visible_inference(visible, selected, "bootstrap")
  expect_identical(result$estimate, visible$estimate)
  expect_equal(result$se, c(.02, .03))
  expect_true(all(result$uncertainty_source == "bootstrap"))
  expect_true(all(result$p_value_method == "normal_approx_from_bootstrap_se"))
  expect_true(all(result$n_boot_success == 20L))
  expect_identical(RAVEL:::.ravel_report_visible_inference(visible, visible, "lavaan"), visible)
})

test_that("whole-trio resampling is paired, reproducible and retains six parameters", {
  roles <- c("child", "mother", "father")
  dat <- as.data.frame(setNames(rep(list(seq_len(100)), 6),
    unlist(lapply(1:2, function(k) paste0("PGS", k, "_", roles)))))
  dat$Y <- dat$PGS1_child * 2
  tab <- data.frame(regression = c("RAVEL", "RAVEL", "PGS1", "PGS2"),
    residual_model = c("fixed", "free", "not_applicable", "not_applicable"),
    parameter = "beta_D", scale = "test", outcome_scale = "test", estimate = 50.5,
    se = 1, ci_low = 0, ci_high = 100, p_value = .5, inference = "model_based")
  vanilla <- list(PGS1 = list(coefficients = tab[3, ]), PGS2 = list(coefficients = tab[4, ]))
  template <- list(summary = tab, vanilla = vanilla,
    packets = list(fixed = list(), free = list()))
  fake <- function(df, ...) {
    expect_equal(df$Y, df$PGS1_child * 2)
    expect_equal(df$PGS2_father, df$PGS1_child)
    expect_gt(nrow(df), length(unique(df$PGS1_child)))
    value <- mean(df$PGS1_child)
    answer <- template
    answer$summary$estimate <- value
    for (nm in names(vanilla)) answer$vanilla[[nm]]$coefficients$estimate <- value
    for (nm in names(answer$packets)) answer$packets[[nm]]$theta <-
      setNames(value + seq_len(6), RAVEL:::.report_core_labels)
    answer$model_status <- data.frame(residual_model = c("fixed", "free"), status = "ok", error = "")
    answer
  }
  runner <- RAVEL:::.ravel_report_bootstrap
  env <- new.env(parent = environment(runner))
  env$ravel_dataframe_report <- fake
  environment(runner) <- env
  set.seed(100)
  before <- .Random.seed
  run <- function() runner(template, dat, 10, 314, .95, "continuous", "joint", character(), c("fixed", "free"))
  first <- run()
  expect_identical(.Random.seed, before)
  second <- run()
  expect_identical(first, second)
  expect_identical(first$bootstrap$core_draws$fixed, first$bootstrap$core_draws$free)
  expect_identical(colnames(first$bootstrap$core_draws$fixed), RAVEL:::.report_core_labels)
  expect_equal(dim(first$bootstrap$core_draws$fixed), c(10L, 6L))
  expect_identical(first$summary_bootstrap$estimate, template$summary$estimate)
  expect_true(all(first$bootstrap$status$status == "ok"))
})

test_that("plots retain chosen intervals and distinguish binary outcome scales", {
  tab <- expand.grid(parameter = c("beta_D", "beta_M", "beta_P"),
    regression = c("PGS1", "PGS2", "RAVEL"), stringsAsFactors = FALSE)
  tab$outcome_variable <- "Y"
  tab$residual_model <- ifelse(tab$regression == "RAVEL", "fixed", "not_applicable")
  tab$estimate <- .2
  tab$ci_low <- .1
  tab$ci_high <- .3
  out <- list(summary = tab, settings = list(outcome_type = "binary_liability", ci_level = .95,
                                            uncertainty = "bootstrap"))
  plot <- ravel_plot_report(out)
  expect_s3_class(plot, "ravel_report_plot")
  expect_true(plot$binary)
  expect_equal(plot$data$ci_low, out$summary$ci_low)
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path, width = 14, height = 8)
  tryCatch(expect_invisible(print(plot)), finally = grDevices::dev.off())
  expect_gt(file.info(path)$size, 1000)
  unlink(path)
  directory <- tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  RAVEL:::.ravel_report_save_plot(plot, directory)
  paths <- file.path(directory, c("diagnostic_estimates.pdf", "diagnostic_estimates.png"))
  expect_true(all(file.exists(paths)))
  expect_true(all(file.info(paths)$size > 1000))
})

test_that("later corrections use joint bootstrap draws and reject insufficient draws", {
  theta <- setNames(c(.3, .1, .02, .1, .55, .55), RAVEL:::.report_core_labels)
  draws <- matrix(rep(theta, each = 20), 20, 6, dimnames = list(NULL, names(theta)))
  draws[, 1] <- draws[, 1] + seq(-.02, .02, length.out = 20)
  draws[, 4] <- draws[, 4] + seq(-.005, .005, length.out = 20)
  packet <- list(schema = "ravel_visible_transfer_v1", diagnostics = list(admissible = TRUE),
    settings = list(uncertainty = "bootstrap", residual_model = "fixed"),
    outcome = list(name = "Y"),
    software = list(transform_fingerprint = RAVEL:::.report_transform_fingerprint()),
    theta = theta, vcov_core = matrix(0, 6, 6, dimnames = list(names(theta), names(theta))),
    bootstrap = list(core_draws = draws, minimum_success = 16))
  corrected <- ravel_correct_packet(packet, h = 0, am_n_generations = 1)
  expect_equal(corrected$estimate[1:3], unname(theta[1:3]))
  expect_equal(corrected$se[1:3], unname(apply(draws[, 1:3], 2, sd)))
  expect_true(all(corrected$n_boot_success == 20))
  packet$bootstrap$core_draws[1:5, ] <- NA
  insufficient <- ravel_correct_packet(packet, h = 0, am_n_generations = 1)
  expect_true(all(is.na(insufficient$se)))
  expect_true(all(insufficient$uncertainty_status == "insufficient_bootstrap"))
  packet$bootstrap <- NULL
  expect_error(ravel_correct_packet(packet, h = 0), "no joint")
  packet$settings$requires_covariate_mapping_review <- TRUE
  expect_error(ravel_correct_packet(packet, h = 0, uncertainty = "lavaan"), "conditional")
})
test_that("binary covariate reports distinguish validated coefficients from baseline fit indices", {
  fixture <- cov_test_fixture(n = 4000L)
  scores <- grep("^PGS", names(fixture$dat), value = TRUE)
  out <- ravel_dataframe_report(fixture$dat, covariates = fixture$covariates,
    pgs_columns = setNames(scores, scores), outcome_type = "binary_liability",
    make_plot = FALSE)
  expect_true(all(out$model_status$status == "ok"))
  for (packet in out$packets) {
    expect_identical(packet$settings$binary_covariate_moments,
                     "conditional_probit_joint_marginal_liability")
    expect_match(packet$diagnostics$fit_measures_note, "independence baseline")
    names <- names(packet$unvalidated_baseline_fit_measures)
    expect_gt(length(names), 0L)
    expect_true(all(is.na(packet$fit_measures[names])))
    expect_true(is.finite(packet$fit_measures["chisq.scaled"]))
    expect_true(is.finite(packet$fit_measures["df.scaled"]))
    expect_true(all(is.finite(diag(packet$vcov_core))))
  }
})
