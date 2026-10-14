# Integration checks for downstream consumers of the opt-in inference features.
feature_reporting_cache <- new.env(parent = emptyenv())

feature_reporting_data <- function() {
  env <- new.env()
  data("ravel_example", package = "RAVEL", envir = env)
  env$ravel_example
}

feature_reporting_fit <- function() {
  if (!exists("fit", feature_reporting_cache, inherits = FALSE)) {
    dat <- feature_reporting_data()
    names(dat) <- sub("^pgs", "PGS", names(dat))
    feature_reporting_cache$fit <- ravel_fit(dat, h = .4, fit_constraints = "estimated_state")
  }
  feature_reporting_cache$fit
}

feature_reporting_report <- function(binary = FALSE, covariates = FALSE, bootstrap = FALSE) {
  key <- paste("report", binary, covariates, bootstrap, sep = "_")
  if (!exists(key, feature_reporting_cache, inherits = FALSE)) {
    feature_reporting_cache[[key]] <- ravel_dataframe_report(feature_reporting_data(),
      outcome = if (binary) "Y_binary" else "Y",
      outcome_type = if (binary) "binary_liability" else "continuous",
      covariates = if (covariates) c("age", "sex") else character(),
      fit_constraints = "estimated_state", h = .4, am_n_generations = 99L,
      n_boot = if (bootstrap) 10L else 0L, seed = 9184L, make_plot = FALSE)
  }
  feature_reporting_cache[[key]]
}

test_that("fit tables and zero-pi-uncertainty corrections retain the same inference", {
  fit <- feature_reporting_fit()
  expect_true(fit$converged)
  tab <- ravel_to_dataframe(fit, include_raw = TRUE, include_vanilla = FALSE)
  expect_equal(as.data.frame(fit, include_raw = TRUE, include_vanilla = FALSE), tab)
  expect_equal(tab$estimate[tab$scale == "visible"], unname(c(fit$beta_vis, fit$rho_vis)))
  pi <- ravel_correct_pi(fit, pi = .4, pi_se = 0)
  normalized <- tab
  normalized$parameter[normalized$parameter == "beta_F"] <- "beta_P"
  key <- function(x) paste(x$scale, x$parameter, sep = ":")
  index <- match(key(pi), key(normalized))
  expect_false(anyNA(index))
  expect_equal(normalized$estimate[index], pi$estimate, tolerance = 1e-8)
  expect_equal(normalized$se[index], pi$se, tolerance = 1e-5)
  expect_equal(normalized$ci_low[index], pi$ci_low, tolerance = 1e-5)
  expect_equal(normalized$ci_high[index], pi$ci_high, tolerance = 1e-5)
})

test_that("live corrected effects retain the input outcome units", {
  fit <- feature_reporting_fit()
  original <- ravel_correct_pi(fit, .4, pi_se = .02)
  dat <- feature_reporting_data()
  names(dat) <- sub("^pgs", "PGS", names(dat))
  dat$Y <- 2 * dat$Y
  scaled_fit <- ravel_fit(dat, h = .4, fit_constraints = "estimated_state")
  expect_true(scaled_fit$converged)
  scaled <- ravel_correct_pi(scaled_fit, .4, pi_se = .02)
  expect_true(all(original$outcome_scale == "input_phenotype_units"))
  expect_true(all(scaled$outcome_scale == "input_phenotype_units"))
  effects <- grepl("^beta_", original$parameter)
  expect_equal(scaled$estimate[effects], 2 * original$estimate[effects], tolerance = 1e-4)
  expect_equal(scaled$se[effects], 2 * original$se[effects], tolerance = 2e-4)
  expect_equal(scaled$estimate[!effects], original$estimate[!effects], tolerance = 1e-4)
})

test_that("aggregate correction round trips preserve pi uncertainty and state diagnostics", {
  report <- feature_reporting_report()
  expect_identical(report$model_status$status, "ok")
  expect_true(report$packet$diagnostics$estimated_state_identified)
  fixed <- ravel_correct_packet(report$packet, h = .4)
  certain <- ravel_correct_pi(report$packet, pi = .4, pi_se = 0)
  uncertain <- ravel_correct_pi(report$packet, pi = .4, pi_se = .02)
  expect_equal(fixed$estimate, certain$estimate, tolerance = 1e-8)
  expect_equal(fixed$se, certain$se, tolerance = 1e-6)
  expect_equal(uncertain$estimate, certain$estimate)
  expect_true(any(uncertain$se > certain$se + 1e-6))
  expect_equal(uncertain$pi_se, rep(.02, 8L))
  expect_identical(unique(uncertain$pi_uncertainty_source), "se")
  expect_true(all(is.na(uncertain$am_n_generations)))
  expect_false(attr(uncertain, "estimated_state")$history_imposed)
  expect_equal(nrow(attr(uncertain, "estimated_state")$summary), 4L)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(report$packet, path)
  expect_equal(ravel_correct_pi(path, pi = .4, pi_se = .02), uncertain)
  saveRDS(uncertain, path)
  expect_identical(readRDS(path), uncertain)
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv(uncertain, csv, row.names = FALSE)
  restored <- utils::read.csv(csv)
  expect_equal(restored$estimate, uncertain$estimate, tolerance = 1e-12)
  expect_equal(restored$se, uncertain$se, tolerance = 1e-12)
  expect_identical(restored$am_correction, uncertain$am_correction)
  expect_equal(restored$pi_se, uncertain$pi_se)
})

test_that("joint continuous and binary reports keep mapping review and correction metadata", {
  for (binary in c(FALSE, TRUE)) {
    report <- feature_reporting_report(binary = binary, covariates = TRUE)
    expect_identical(report$model_status$status, "ok", info = paste("binary", binary))
    packet <- report$packet
    expect_true(packet$settings$requires_covariate_mapping_review)
    expect_identical(packet$settings$constraint_spec$mode, "estimated_state")
    expect_true(packet$diagnostics$estimated_state_identified)
    expect_no_error(ravel_correct_packet(packet, h = .4))
    corrected <- ravel_correct_pi(packet, pi = .4, pi_se = .02)
    expect_equal(ravel_correct_pi(packet, pi = .4, pi_se = .02, covariate_mapping_confirmed = TRUE), corrected)
    expect_true(all(is.finite(corrected$estimate)))
    expect_true(all(is.finite(corrected$se)))
    expect_identical(unique(corrected$outcome_variable), if (binary) "Y_binary" else "Y")
    expect_identical(packet$outcome$type, if (binary) "binary_liability" else "continuous")
    expect_true(all(grepl(if (binary) "SEM_probit_liability" else "prepared_phenotype_SD",
      report$summary$outcome_scale[report$summary$regression == "RAVEL"])))
  }
})

test_that("bootstrap reporting carries joint covariance into external-pi uncertainty", {
  report <- feature_reporting_report(bootstrap = TRUE)
  expect_identical(report$model_status$status, "ok")
  packet <- report$packet
  expect_identical(packet$settings$uncertainty, "bootstrap")
  corrected <- ravel_correct_pi(packet, pi = .4, pi_se = .02)
  expect_true(all(corrected$inference == "bootstrap_covariance_pi_delta_wald"))
  expect_true(all(corrected$n_boot_success >= 10L))
  expect_true(all(is.finite(corrected$se)))
  expect_equal(attr(corrected, "pi_uncertainty")$raw_scale_se, .02)
  expect_true(all(is.finite(attr(corrected, "estimated_state")$summary$se)))
  expect_equal(report$summary, report$summary_bootstrap)
  expect_true(all(report$summary$inference == "bootstrap_percentile_CI_normal_p"))
})

test_that("fitted tables and print methods disclose estimated-state assumptions", {
  fit <- feature_reporting_fit()
  tab <- ravel_to_dataframe(fit, include_raw = TRUE, include_vanilla = FALSE)
  expect_identical(attr(tab, "constraint_spec")$mode, "estimated_state")
  expect_identical(attr(tab, "settings")$am_correction, "estimated_state")
  expect_true(is.na(attr(tab, "settings")$am_n_generations))
  expect_equal(attr(tab, "settings")$h, .4)
  expect_false(attr(tab, "estimated_state")$history_imposed)
  expect_true(is.finite(attr(tab, "estimated_state")$growth))
  for (printer in list(print_simple, print_ravel_results)) {
    text <- paste(capture.output(printer(fit)), collapse = "\n")
    expect_match(text, "estimated_state", fixed = TRUE)
    expect_match(text, "[Gg]rowth")
  }
  report <- feature_reporting_report()
  expect_true(is.na(report$settings$am_n_generations))
  expect_true(is.na(report$packet$settings$am_n_generations))
  correction <- ravel_correct_pi(report$packet, .4, pi_se = .02)
  for (printer in list(print_simple, print_ravel_results)) {
    text <- paste(capture.output(printer(correction)), collapse = "\n")
    expect_match(text, "estimated_state", fixed = TRUE)
    expect_match(text, "0.02", fixed = TRUE)
    expect_match(text, "[Pp]i")
  }
})

test_that("corrected tables plot all three roles on one explicit scale", {
  report <- feature_reporting_report()
  corrected <- ravel_correct_pi(report$packet, pi = .4, pi_se = .02, ci_level = .9)
  plot <- ravel_plot_report(corrected)
  expect_s3_class(plot, "ravel_report_plot")
  expect_setequal(plot$data$parameter, c("beta_D", "beta_M", "beta_P"))
  expect_equal(nrow(plot$data), 3L)
  expect_true(all(plot$data$scale == "corrected_std"))
  expect_false(any(grepl("Visible", plot$data$model, fixed = TRUE)))
  expect_equal(plot$ci_level, .9)
  expect_equal(plot$data$pi_se, rep(.02, 3L))
  want <- corrected[corrected$scale == "corrected_std" & corrected$parameter != "rho_corrected", ]
  ix <- match(plot$data$parameter, want$parameter)
  expect_equal(plot$data$estimate, want$estimate[ix])
  expect_equal(plot$data$ci_low, want$ci_low[ix])
  expect_equal(plot$data$ci_high, want$ci_high[ix])
  raw <- ravel_plot_report(corrected, scale = "corrected_raw")
  expect_true(all(raw$data$scale == "corrected_raw"))
  expect_equal(nrow(raw$data), 3L)
  expect_error(ravel_plot_report(corrected, scale = "unavailable"), "scale|Scale")
  legacy <- ravel_plot_report(report)
  expect_equal(nrow(legacy$data), 9L)
  expect_true(all(legacy$data$scale %in% c("visible", "own_role_PGS_SD")))
  # Legacy beta_F remains accepted, while plotting normalizes the role name.
  from_fit <- ravel_plot_report(ravel_to_dataframe(feature_reporting_fit(),
    include_raw = TRUE, include_vanilla = FALSE), scale = "corrected_std")
  expect_setequal(from_fit$data$parameter, c("beta_D", "beta_M", "beta_P"))
  expect_equal(nrow(from_fit$data), 3L)
  binary <- ravel_correct_pi(feature_reporting_report(binary = TRUE, covariates = TRUE)$packet,
    .4, pi_se = .02, covariate_mapping_confirmed = TRUE)
  binary_plot <- ravel_plot_report(binary)
  expect_true(binary_plot$binary)
  pdf <- tempfile(fileext = ".pdf")
  on.exit(unlink(pdf))
  grDevices::pdf(pdf, width = 14, height = 8)
  tryCatch(expect_invisible(print(binary_plot)), finally = grDevices::dev.off())
  expect_gt(file.info(pdf)$size, 1000)
})

test_that("report export creates usable graphics and round-trips the entire manifest", {
  directory <- tempfile()
  on.exit(unlink(directory, recursive = TRUE))
  report <- ravel_dataframe_report(feature_reporting_data(), fit_constraints = "estimated_state",
    h = .4, out_dir = directory, make_plot = TRUE)
  expect_identical(report$model_status$status, "ok")
  expect_true(all(c("diagnostic_estimates.pdf", "diagnostic_estimates.png", "ravel_report.rds") %in% names(report$files)))
  expect_true(all(file.exists(report$files)))
  expect_equal(anyDuplicated(unname(report$files)), 0L)
  restored <- readRDS(report$files[["ravel_report.rds"]])
  expect_identical(restored, report)
  expect_true(restored$packet$diagnostics$estimated_state_identified)
  expect_identical(restored$settings$fit_constraints, "estimated_state")
  png <- report$files[["diagnostic_estimates.png"]]
  expect_gt(file.info(png)$size, 1000)
  expect_identical(readBin(png, "raw", 8L), as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)))
  pdf <- report$files[["diagnostic_estimates.pdf"]]
  expect_identical(rawToChar(readBin(pdf, "raw", 4L)), "%PDF")
})

test_that("nonidentification remains visible in report diagnostics", {
  R <- matrix(c(1, .5, .5, .5, 1, 0, .5, 0, 1), 3L)
  product <- .7 * .85
  S <- rbind(cbind(R, product * R), cbind(product * R, R))
  S <- rbind(cbind(S, 0), c(rep(0, 6L), 1))
  N <- 600L
  x <- outer(seq_len(N), seq_len(7L), function(i, j) sin(i * j) + cos(i * (j + .5)))
  Q <- qr.Q(qr(scale(x, center = TRUE, scale = FALSE)))
  dat <- as.data.frame(sqrt(N) * Q %*% chol(S))
  names(dat) <- c(paste0("pgs1_", c("child", "mother", "father")),
    paste0("pgs2_", c("child", "mother", "father")), "Y")
  report <- ravel_dataframe_report(dat, h = .4, fit_constraints = "estimated_state", make_plot = FALSE)
  expect_identical(report$model_status$status, "failed")
  expect_null(report$packet)
  expect_match(tolower(report$model_status$error), "identif|uncertainty|singular|covariance")
  expect_identical(report$diagnostics$fixed$error, report$model_status$error)
})

test_that("zero hidden fraction keeps founder-scale raw bootstrap aliases", {
  dat <- feature_reporting_data()
  names(dat) <- sub("^pgs", "PGS", names(dat))
  fit <- ravel_fit(dat, h = 0, fit_constraints = "estimated_state", n_boot = 10L, seed = 9184L)
  expect_true(fit$converged)
  expect_gte(fit$bootstrap$n_successful, 10L)
  boot <- fit$bootstrap
  expect_gt(max(abs(fit$beta_total_assumed - fit$beta_vis)), 1e-4)
  for (role in c(D = "bD", M = "bm", F = "bf")) {
    prefix <- c(bD = "bD", bm = "bM", bf = "bF")[[role]]
    expect_equal(boot[[paste0(prefix, "_se")]], boot$total_raw$se[role])
    expect_equal(boot[[paste0(prefix, "_ci")]], boot$total_raw$ci[, role])
    expect_equal(boot[[paste0("visible_", prefix, "_se")]], boot$visible$se[role])
    expect_equal(boot[[paste0("total_std_", prefix, "_se")]], boot$total_std$se[role])
  }
  expect_equal(boot$cov_Tmf_se, boot$total_raw$se["rhoMF"])
  expect_equal(boot$cov_Tmf_ci, boot$total_raw$ci[, "rhoMF"])
  expect_equal(boot$corr_Tmf_se, boot$total_std$se["rhoMF"])
  tab <- ravel_to_dataframe(fit, uncertainty = "bootstrap", include_raw = TRUE, include_vanilla = FALSE)
  raw <- tab[tab$scale == "corrected_raw", ]
  expect_equal(nrow(raw), 4L)
  expect_equal(raw$se, unname(boot$total_raw$se))
  expect_equal(raw$ci_low, unname(boot$total_raw$ci["lower", ]))
  expect_equal(raw$ci_high, unname(boot$total_raw$ci["upper", ]))
})
