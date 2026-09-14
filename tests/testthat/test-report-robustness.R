robustness_packet <- function() {
  theta <- setNames(c(.3, .1, .02, .1, .55, .55), RAVEL:::.report_core_labels)
  draws <- sweep(outer(seq(-.01, .01, length.out = 20), c(1, .3, -.2, .1, .15, .2)),
                 2, theta, "+")
  colnames(draws) <- names(theta)
  list(schema = "ravel_visible_transfer_v1", diagnostics = list(admissible = TRUE),
    settings = list(uncertainty = "bootstrap", residual_model = "fixed"),
    outcome = list(name = "Y"), theta = theta,
    vcov_core = matrix(0, 6, 6, dimnames = list(names(theta), names(theta))),
    software = list(transform_fingerprint = RAVEL:::.report_transform_fingerprint()),
    bootstrap = list(core_draws = draws, minimum_success = 16L))
}

test_that("packet bootstrap transformations align named columns without changing inputs", {
  packet <- robustness_packet()
  shuffled <- packet
  shuffled$bootstrap$core_draws <- shuffled$bootstrap$core_draws[, c(3, 6, 1, 4, 2, 5)]
  before <- shuffled
  for (h in c(0, .4)) for (regime in c("n_gen", "equilibrium")) {
    expect_equal(ravel_correct_packet(shuffled, h, regime, am_n_generations = 1),
                 ravel_correct_packet(packet, h, regime, am_n_generations = 1))
  }
  expect_identical(shuffled, before)
  shuffled$bootstrap$core_draws[1:5, ] <- NA_real_
  got <- ravel_correct_packet(shuffled, .4, am_n_generations = 1)
  expect_true(all(got$n_boot_success == 15L))
  expect_true(all(got$uncertainty_status == "insufficient_bootstrap"))
  expect_true(all(is.na(got$se)))
})

test_that("packet bootstrap matrices require exactly six unique numeric named columns", {
  packet <- robustness_packet()
  x <- packet$bootstrap$core_draws
  duplicate <- x
  colnames(duplicate)[2] <- colnames(duplicate)[1]
  missing_name <- x
  colnames(missing_name)[1] <- NA_character_
  empty_name <- x
  colnames(empty_name)[1] <- ""
  bad <- list(unname(x), as.data.frame(x), x[, 1:5], cbind(x, extra = 0),
    duplicate, missing_name, empty_name, matrix(as.character(x), nrow(x),
      dimnames = dimnames(x)), x > 0, as.numeric(x), x + 1i, x[FALSE, , drop = FALSE],
    array(x, dim = c(10L, 6L, 2L)))
  for (draws in bad) {
    packet$bootstrap$core_draws <- draws
    expect_error(ravel_correct_packet(packet, .4), "numeric matrix.*unique column names")
  }
  packet$bootstrap$core_draws <- NULL
  expect_error(ravel_correct_packet(packet, .4), "no joint")
  packet$bootstrap$core_draws <- "not used"
  expect_s3_class(ravel_correct_packet(packet, .4, uncertainty = "lavaan"), "data.frame")
})

robustness_data <- function() {
  i <- seq_len(120)
  columns <- unlist(lapply(1:2, function(k) paste0("PGS", k, "_", c("child", "mother", "father"))))
  dat <- as.data.frame(setNames(lapply(seq_along(columns), function(j)
    sin(i * j / 7) + j * cos(i / 13)), columns))
  dat$Y <- sin(i / 3) + i / 120
  dat
}

robustness_report <- function(fail_at = NULL, attempts = new.env(parent = emptyenv())) {
  attempts$models <- character()
  report <- RAVEL::ravel_dataframe_report
  env <- new.env(parent = environment(report))
  env$.ravel_report_vanilla <- function(dat, version, fixed_effects, scaling, response_scale,
                                       response_center, outcome_type, ci_level, outcome_name) {
    model <- paste0("PGS", version)
    attempts$models <- c(attempts$models, model)
    if (identical(fail_at, model)) stop("Injected ", model, " failure")
    tab <- data.frame(outcome_variable = outcome_name, regression = model,
      scale = "own_role_PGS_SD", outcome_scale = "prepared_phenotype_SD",
      parameter = c("beta_D", "beta_M", "beta_P"), estimate = c(.3, .1, .02),
      se = .01, ci_low = 0, ci_high = .4, p_value = .1, N = nrow(dat), inference = "lm_t")
    list(coefficients = tab, fit_statistics = c(N = nrow(dat)))
  }
  env$.ravel_report_sem <- function(dat, common_packet, specification, fixed_corr, wanted,
                                   vanilla_summary, vanilla_statistics, stats_table, condition) {
    attempts$models <- c(attempts$models, specification)
    tab <- data.frame(outcome_variable = common_packet$outcome$name, regression = "RAVEL",
      scale = "visible", outcome_scale = "prepared_phenotype_SD",
      parameter = c("beta_D", "beta_M", "beta_P", "rho_visible"),
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_, N = nrow(dat), inference = "unavailable",
      residual_model = specification, fit_status = "failed")
    list(summary = tab[wanted], fit_statistics = stats_table(c(N = nrow(dat)), "RAVEL", specification),
      packet = NULL, diagnostics = list(admissible = FALSE, converged = FALSE, error = "Unavailable test SEM"),
      residuals = data.frame(residual_model = specification, parameter = names(fixed_corr),
        covariance = NA_real_, correlation = NA_real_))
  }
  if (identical(fail_at, "preprocessing"))
    env$.ravel_report_covariates <- function(...) stop("Injected preprocessing failure")
  environment(report) <- env
  report
}

robustness_run <- function(report, ...) {
  dat <- robustness_data()
  report(dat, pgs_columns = setNames(names(dat)[1:6], names(dat)[1:6]), make_plot = FALSE, ...)
}

test_that("report errors retain the failed stage and previously completed model status", {
  for (stage in c("preprocessing", "PGS1", "PGS2")) {
    attempts <- new.env(parent = emptyenv())
    report <- robustness_report(stage, attempts)
    err <- tryCatch(robustness_run(report), error = identity)
    expect_s3_class(err, "ravel_report_error")
    expect_identical(err$report_stage, stage)
    expected <- c(fixed = "not_attempted", free = "not_attempted", PGS1 = "not_attempted", PGS2 = "not_attempted")
    if (stage != "preprocessing") expected[stage] <- "failed"
    if (stage == "PGS2") expected["PGS1"] <- "ok"
    expect_identical(err$report_model_status, expected)
    expect_identical(attempts$models, switch(stage, preprocessing = character(), PGS1 = "PGS1",
                                            PGS2 = c("PGS1", "PGS2")))
    expect_match(conditionMessage(err), paste("Injected", stage, "failure"))
  }
})

test_that("aborted bootstrap reports retain truthful attempts but discard every estimate", {
  template <- robustness_run(robustness_report())
  dat <- robustness_data()
  for (stage in c("preprocessing", "PGS1", "PGS2", "unclassified")) {
    report <- if (stage == "unclassified") function(...) stop("Injected preprocessing failure") else
      robustness_report(stage)
    runner <- RAVEL:::.ravel_report_bootstrap
    env <- new.env(parent = environment(runner))
    env$ravel_dataframe_report <- report
    environment(runner) <- env
    expect_warning(out <- runner(template, dat, 10L, 314L, .95, "continuous", "joint",
                                 character(), c("fixed", "free")), "insufficient usable")
    expected <- c(fixed = "not_attempted", free = "not_attempted", PGS1 = "not_attempted", PGS2 = "not_attempted")
    if (stage %in% c("PGS1", "PGS2")) expected[stage] <- "failed"
    if (stage == "PGS2") expected["PGS1"] <- "ok"
    status <- out$bootstrap$status
    expect_identical(status$status, unname(expected[status$model]))
    expect_true(all(status$error_stage == if (stage == "unclassified") "unknown" else stage))
    expect_false(any(status$draw_retained))
    expect_true(all(is.na(out$bootstrap$summary_draws)))
    expect_true(all(vapply(out$bootstrap$core_draws, function(x) all(is.na(x)), logical(1))))
    expect_true(all(vapply(out$bootstrap$vanilla_draws, function(x) all(is.na(x)), logical(1))))
    expect_identical(out$summary_bootstrap$estimate, template$summary$estimate)
    expect_true(all(out$summary_bootstrap$uncertainty_status == "insufficient_bootstrap"))
    expect_true(all(is.na(out$summary_bootstrap$se)))
  }
})

test_that("nested preclassified report errors propagate without losing stage metadata", {
  original <- tryCatch(robustness_run(robustness_report("PGS2")), error = identity)
  expect_s3_class(original, "ravel_report_error")
  report <- robustness_report()
  environment(report)$.ravel_report_covariates <- function(...) signalCondition(original)
  propagated <- tryCatch(robustness_run(report), error = identity)
  expect_identical(propagated, original)
  expect_identical(propagated$report_stage, "PGS2")
  expect_identical(propagated$report_model_status[c("PGS1", "PGS2")],
                   c(PGS1 = "ok", PGS2 = "failed"))
})

test_that("export manifests list each existing file once and survive RDS round trips", {
  directory <- tempfile()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  report <- robustness_report()
  out <- robustness_run(report, out_dir = directory)
  expect_identical(anyDuplicated(unname(out$files)), 0L)
  expect_identical(anyDuplicated(names(out$files)), 0L)
  expect_true(all(file.exists(out$files)))
  expect_identical(unname(out$files), list.files(normalizePath(directory), full.names = TRUE))
  expect_equal(sum(names(out$files) == "ravel_report.rds"), 1L)
  restored <- readRDS(out$files[["ravel_report.rds"]])
  expect_identical(restored, out)
  expect_error(robustness_run(report, out_dir = directory), "nothing will be overwritten")
})

test_that("plot export failure preserves estimates and a correct saved manifest", {
  directory <- tempfile()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  report <- robustness_report()
  environment(report)$.ravel_report_save_plot <- function(...) stop("Injected plot export failure")
  dat <- robustness_data()
  expect_warning(out <- report(dat, pgs_columns = setNames(names(dat)[1:6], names(dat)[1:6]),
    out_dir = directory, make_plot = TRUE), "Estimates were saved, but plot export failed")
  expect_identical(anyDuplicated(unname(out$files)), 0L)
  expect_true(all(file.exists(out$files)))
  expect_identical(readRDS(out$files[["ravel_report.rds"]]), out)
})

test_that("completed bootstrap reports distinguish retained draws from unavailable models", {
  report <- robustness_report()
  template <- robustness_run(report)
  runner <- RAVEL:::.ravel_report_bootstrap
  env <- new.env(parent = environment(runner))
  env$ravel_dataframe_report <- report
  environment(runner) <- env
  expect_warning(out <- runner(template, robustness_data(), 10L, 314L, .95, "continuous", "joint",
    character(), c("fixed", "free")), "insufficient usable")
  status <- out$bootstrap$status
  pgs <- status$model %in% c("PGS1", "PGS2")
  expect_true(all(status$status[pgs] == "ok"))
  expect_true(all(status$draw_retained[pgs]))
  expect_true(all(status$error_stage[pgs] == ""))
  expect_true(all(status$status[!pgs] == "failed"))
  expect_false(any(status$draw_retained[!pgs]))
  expect_identical(status$error_stage[!pgs], status$model[!pgs])
  expect_true(all(vapply(out$bootstrap$vanilla_draws, function(x) all(is.finite(x)), logical(1))))
})
