.sib_report_cache <- new.env(parent = emptyenv())

.sib_report_data <- function() {
  if (!exists("data", .sib_report_cache, inherits = FALSE)) {
    dat <- ravel_simulate_trios(N = 6000, m = 80, parental_sibs = "both",
      rho_MF = .2, lineage_effect_mat = .1, seed = 48)
    dat[seq_len(700), c("PGS1_mat_sib", "PGS2_mat_sib")] <- NA_real_
    dat[701:1200, c("PGS1_pat_sib", "PGS2_pat_sib")] <- NA_real_
    .sib_report_cache$data <- dat
  }
  .sib_report_cache$data
}

.sib_report <- function(mode = "both", constraints = "none", bootstrap = FALSE) {
  key <- paste(mode, constraints, bootstrap, sep = "_")
  if (!exists(key, .sib_report_cache, inherits = FALSE)) {
    dat <- .sib_report_data()
    columns <- c("PGS1_child", "PGS1_mother", "PGS1_father", "PGS2_child", "PGS2_mother", "PGS2_father",
      .parental_sib_columns(.parental_sib_roles(mode)))
    .sib_report_cache[[key]] <- ravel_dataframe_report(dat, parental_sibs = mode,
      pgs_columns = setNames(columns, columns), residual_models = "fixed",
      fit_constraints = constraints, h = .3, am_n_generations = 1L,
      n_boot = if (bootstrap) 10L else 0L, seed = 91,
      make_plot = FALSE)
  }
  .sib_report_cache[[key]]
}

test_that("sibling reports retain incomplete relatives and transfer every estimand", {
  for (mode in c("maternal", "paternal", "both")) {
    out <- .sib_report(mode)
    expect_identical(out$model_status$status, "ok")
    roles <- .parental_sib_roles(mode)
    labels <- c(.report_core_labels, .sib_core_labels(roles))
    params <- unname(.parental_sib_coef_map(roles))
    expect_identical(names(out$packet$theta), labels)
    expect_identical(rownames(out$packet$vcov_core), labels)
    expect_identical(out$packet$settings$parental_sibs, mode)
    expect_equal(out$packet$settings$n_complete, 6000)
    expect_true(all(out$summary$N == 6000))
    expect_true(all(params %in% out$summary$parameter))
    expect_true(all(is.finite(out$summary$se)))
    expected_residuals <- .parental_sib_residual_benchmarks(roles)
    expect_equal(out$packet$parental_sibs$resid_corr[names(expected_residuals)], expected_residuals, tolerance = 1e-5)
    expect_true(all(is.finite(out$scaling$scaled_variance)))
    expect_equal(out$packet$parental_sibs$n_with_sib,
      c(mat_sib = 5300, pat_sib = 5500)[unname(c(ms = "mat_sib", fs = "pat_sib")[roles])])
    expect_true(all(unname(c(ms = "mat_sib", fs = "pat_sib")[roles]) %in% out$reliability$fixed$role))
    corrected <- ravel_correct_packet(out$packet, .3, am_n_generations = 1)
    for (scale in c("corrected_raw", "corrected_std")) {
      table <- corrected[corrected$scale == scale, ]
      expect_true(all(params %in% table$parameter))
      expect_true(all(is.finite(table$se)))
    }
    pi_zero <- ravel_correct_pi(out$packet, pi = .3, pi_se = 0, am_n_generations = 1)
    at <- match(paste(corrected$scale, corrected$parameter), paste(pi_zero$scale, pi_zero$parameter))
    expect_equal(pi_zero$estimate[at], corrected$estimate, tolerance = 1e-7)
    expect_equal(pi_zero$se[at], corrected$se, tolerance = 1e-7)
    pi_uncertain <- ravel_correct_pi(out$packet, pi = .3, pi_se = .04, am_n_generations = 1)
    expect_true(all(pi_uncertain$se[at] >= corrected$se - 1e-8))
    expect_output(print_simple(pi_uncertain), params[2])
    expect_output(print_ravel_results(pi_uncertain), params[1])
  }
})

test_that("packet corrections match refitting on the report's prepared phenotype and scores", {
  dat <- .sib_report_data()
  out <- .sib_report("both")
  for (i in seq_len(nrow(out$scaling))) {
    row <- out$scaling[i, ]
    dat[[row$model_column]] <- (dat[[row$source_column]] - row$child_reference_mean) / row$child_reference_sd
  }
  # Reproduce the recorded preprocessing exactly; machine-precision changes in
  # a separately standardized outcome can change nonlinear optimizer stopping.
  dat$Y <- (dat$Y - out$packet$outcome$preparation_center) / out$packet$outcome$preparation_sd
  fit <- suppressWarnings(ravel_fit(dat, parental_sibs = "both", pi = .3,
    fixed_crossperson_corr = c(e1_mf = 0, e1_mc = .5, e1_fc = .5, e2_mf = 0, e2_mc = .5, e2_fc = .5,
      .parental_sib_residual_benchmarks(c("ms", "fs")))))
  live <- ravel_correct_pi(fit, pi = .3, pi_se = .04)
  packet <- ravel_correct_pi(out$packet, pi = .3, pi_se = .04, am_n_generations = 1)
  expect_equal(packet$estimate, live$estimate, tolerance = 1e-6)
  expect_equal(packet$se, live$se, tolerance = 1e-6)
  expect_equal(packet$ci_low, live$ci_low, tolerance = 1e-6)
})

test_that("default sibling mappings, custom mappings and binary family selection are explicit", {
  dat <- .sib_report_data()
  for (nm in grep("^PGS[12]_", names(dat), value = TRUE)) names(dat)[names(dat) == nm] <- tolower(nm)
  out <- ravel_dataframe_report(dat, parental_sibs = "maternal", residual_models = "fixed", make_plot = FALSE)
  expect_equal(out$summary$estimate, .sib_report("maternal")$summary$estimate)
  expect_error(ravel_dataframe_report(dat, parental_sibs = "maternal",
    pgs_columns = c(PGS1_child = "pgs1_child"), make_plot = FALSE), "expected score columns")
  dat$Y <- as.numeric(dat$Y > stats::median(dat$Y))
  expect_warning(binary <- ravel_dataframe_report(dat, parental_sibs = "both",
    outcome_type = "binary_liability", residual_models = "fixed", make_plot = FALSE),
    "exclude families with missing requested sibling scores")
  expect_identical(binary$model_status$status, "ok")
  expect_true(all(binary$summary$N == 4800))
  row <- binary$descriptives$exclusions[binary$descriptives$exclusions$reason == "missing_or_nonfinite_sibling_PGS", ]
  expect_equal(row$n_excluded_sequential, 1200)
})

test_that("sibling correction grids and plots preserve derived contrasts including failed cells", {
  out <- .sib_report("both")
  corrected <- ravel_correct_grid(out, c(.2, .3), n_generations = c(1L, 2L, Inf), pi_se = .04)
  expect_equal(nrow(corrected), 2 * 3 * 16)
  expect_true(all(corrected$correction_status == "ok"))
  plot <- ravel_plot_report(corrected)
  expect_setequal(unique(plot$data$parameter), c("beta_D", "beta_M", "beta_P", "beta_MS", "delta_M", "beta_FS", "delta_F"))
  expect_true(all(plot$data$scale == "corrected_std"))
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path, width = 14, height = 14)
  tryCatch(print(plot), finally = grDevices::dev.off())
  expect_gt(file.info(path)$size, 0)
  unlink(path)
  bad <- out
  bad$packets$fixed <- NULL
  bad$packets <- list(fixed = NULL)
  expect_warning(failed <- ravel_correct_grid(bad, .3, 1), "Some corrections")
  expect_equal(nrow(failed), 16)
  expect_true(all(c("delta_M", "delta_F") %in% failed$parameter))
  expect_true(all(failed$correction_status == "unavailable_model"))
})

test_that("whole-family bootstrap packets retain sibling covariance and missingness", {
  out <- .sib_report("maternal", bootstrap = TRUE)
  expect_identical(colnames(out$packet$bootstrap$core_draws), names(out$packet$theta))
  expect_equal(ncol(out$packet$bootstrap$core_draws), 10)
  corrected <- ravel_correct_packet(out$packet, .3, am_n_generations = 1)
  expect_true(all(corrected$n_boot_success == 10))
  expect_true(all(is.finite(corrected$se)))
  uncertain <- ravel_correct_pi(out$packet, pi = .3, pi_se = .04, am_n_generations = 1)
  expect_true(all(is.finite(uncertain$se)))
  expect_true(all(c("beta_MS", "delta_M") %in% uncertain$parameter))
  zero <- ravel_correct_packet(out$packet, 0, am_n_generations = 1)
  zero <- zero[zero$scale == "corrected_std" & zero$parameter != "rho_corrected", ]
  visible <- out$summary[out$summary$regression == "RAVEL", ]
  visible <- visible[match(zero$parameter, visible$parameter), ]
  expect_equal(zero$estimate, visible$estimate, tolerance = 1e-8)
  expect_equal(zero$se, visible$se, tolerance = 1e-8)
  expect_equal(zero$ci_low, visible$ci_low, tolerance = 1e-8)
  expect_equal(zero$ci_high, visible$ci_high, tolerance = 1e-8)
})

test_that("sibling report export retains full covariance and coefficient tables", {
  dat <- .sib_report_data()
  columns <- c("PGS1_child", "PGS1_mother", "PGS1_father", "PGS2_child", "PGS2_mother", "PGS2_father",
    "PGS1_mat_sib", "PGS2_mat_sib")
  root <- tempfile()
  on.exit(unlink(root, recursive = TRUE))
  out <- ravel_dataframe_report(dat, parental_sibs = "maternal", pgs_columns = setNames(columns, columns),
    residual_models = "fixed", out_dir = root, make_plot = FALSE)
  saved <- readRDS(file.path(root, "ravel_report.rds"))
  tab <- utils::read.csv(file.path(root, "summary.csv"))
  expect_true(all(c("beta_MS", "delta_M") %in% tab$parameter))
  expect_identical(saved$packet$theta, out$packet$theta)
  expect_equal(ncol(saved$packet$vcov_core), 10)
})

test_that("constrained sibling packets retain history and match live pi propagation", {
  dat <- ravel_simulate_trios(N = 8000, m = 120, m_hidden = 120, hidden_fraction = .4,
    h2_direct = .09, h2_mat = .0225, h2_pat = .0025, rho_MF = .3, grandparent_rho_MF = .3,
    parental_sibs = "maternal", lineage_effect_mat = .10, seed = 31)
  columns <- c("PGS1_child", "PGS1_mother", "PGS1_father", "PGS2_child", "PGS2_mother", "PGS2_father",
    "PGS1_mat_sib", "PGS2_mat_sib")
  for (mode in c("n_gen", "equilibrium", "estimated_state")) {
    out <- ravel_dataframe_report(dat, parental_sibs = "maternal", pgs_columns = setNames(columns, columns),
      fit_constraints = mode, h = .4, am_n_generations = 2, make_plot = FALSE)
    expect_identical(out$model_status$status, "ok", info = mode)
    expect_identical(out$packet$settings$constraint_spec$mode, mode)
    corrected <- ravel_correct_packet(out$packet, .4)
    expect_true(all(is.finite(corrected$se)), info = mode)
    prepared <- dat
    for (i in seq_len(nrow(out$scaling))) {
      row <- out$scaling[i, ]
      prepared[[row$model_column]] <- (prepared[[row$source_column]] - row$child_reference_mean) / row$child_reference_sd
    }
    prepared$Y <- (dat$Y - out$packet$outcome$preparation_center) / out$packet$outcome$preparation_sd
    live <- suppressWarnings(ravel_fit(prepared, parental_sibs = "maternal", pi = .4, pi_se = .02,
      fit_constraints = mode, am_n_generations = 2))
    sensitivity <- if (mode == "n_gen") attr(live$pi_correction, "pi_uncertainty")$core_pi_sensitivity else NULL
    if (mode == "n_gen") expect_error(ravel_correct_pi(out$packet, pi = .4, pi_se = .02), "sensitivity")
    pi_packet <- ravel_correct_pi(out$packet, pi = .4, pi_se = .02, core_pi_sensitivity = sensitivity)
    expect_equal(pi_packet$estimate, live$pi_correction$estimate, tolerance = 1e-6, info = mode)
    expect_equal(pi_packet$se, live$pi_correction$se, tolerance = 1e-6, info = mode)
    grid <- ravel_correct_grid(out, hidden_fraction = .4, pi_se = .02, core_pi_sensitivity = sensitivity)
    expect_true(all(grid$correction_status == "ok"), info = mode)
    expect_equal(grid$se, pi_packet$se, tolerance = 1e-8, info = mode)
  }
})


test_that("generation-one reports accept exactly fixed sibling correlation parameters", {
  out <- .sib_report("maternal", constraints = "n_gen")
  expect_identical(out$model_status$status, "ok")
  expect_equal(out$packet$theta[["kap_m"]], .5, tolerance = 1e-7)
  expect_lt(abs(out$packet$vcov_core["kap_m", "kap_m"]), 1e-10)
  corrected <- ravel_correct_packet(out$packet, .3)
  expect_true(all(is.finite(corrected$se)))
})
