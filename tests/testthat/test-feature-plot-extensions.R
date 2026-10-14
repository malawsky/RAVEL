# Effect plots must not silently discard coefficient families returned by tables.
test_that("plots retain interaction and parental contrast coefficients", {
  parameters <- c("beta_D", "beta_M", "beta_F", "rho_MF", "beta_MS", "delta_M",
    "beta_mid", "delta_M_mid", "delta_P_mid", "beta_D:age", "beta_F:age",
    "beta_MS:age", "delta_M:age", "beta_D:age2")
  tab <- data.frame(parameter = parameters, scale = "corrected_std",
    estimate = seq_along(parameters) / 100, se = .01,
    ci_low = seq_along(parameters) / 100 - .02,
    ci_high = seq_along(parameters) / 100 + .02,
    uncertainty_source = "lavaan cluster-robust + pi (delta-Wald)")
  attr(tab, "settings") <- list(cluster = "FID", ci_level = .95)
  p <- ravel_plot_report(tab)
  expected <- sub("^beta_F(?=:|$)", "beta_P", setdiff(parameters, "rho_MF"), perl = TRUE)
  expect_setequal(p$data$parameter, expected)
  expect_identical(p$uncertainty, "pi_delta_wald")
  expect_match(p$caption, "external pi uncertainty")
  labels <- .ravel_plot_parameter_labels(p$data$parameter)
  expect_identical(unname(labels["beta_P:age"]), "Father x age")
  expect_identical(unname(labels["delta_M:age"]), "Mother-specific contrast x age")
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path, width = 14, height = 20)
  tryCatch(expect_no_error(print(p)), finally = grDevices::dev.off())
  expect_gt(file.info(path)$size, 0)
  unlink(path)

  # An interaction-only selection is a legitimate effect plot.
  ix <- grepl(":", tab$parameter, fixed = TRUE)
  p_ix <- ravel_plot_report(tab[ix, ])
  expect_equal(nrow(p_ix$data), sum(ix))
})

test_that("plot captions preserve clustered inference", {
  tab <- data.frame(parameter = "beta_D:age", scale = "visible", estimate = .02,
    ci_low = .01, ci_high = .03, uncertainty_source = "lavaan cluster-robust")
  attr(tab, "settings") <- list(cluster = "FID", ci_level = .95)
  expect_match(ravel_plot_report(tab)$caption, "cluster-robust")
  tab$uncertainty_source <- "bootstrap (500 reps)"
  expect_match(ravel_plot_report(tab)$caption, "cluster bootstrap")
  expect_false(grepl("whole-trio", ravel_plot_report(tab)$caption))
})


test_that("main-effect-only plots have no empty interaction panel", {
  tab <- data.frame(parameter = c("beta_D", "beta_M", "beta_F"),
    scale = "visible", estimate = c(.3, .1, .05), ci_low = c(.2, 0, -.05),
    ci_high = c(.4, .2, .15), uncertainty_source = "lavaan/delta-method")
  p <- ravel_plot_report(tab)
  expect_identical(names(.ravel_plot_parameter_labels(p$data$parameter)),
    c("beta_D", "beta_M", "beta_P"))
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path, width = 14, height = 5)
  tryCatch(expect_no_error(print(p)), finally = grDevices::dev.off())
  unlink(path)
})
