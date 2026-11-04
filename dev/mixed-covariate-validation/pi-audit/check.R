# Independent wider-step external-pi sensitivity checks with correlated PCs
# and a genetic-by-sex effect. Run from the package root.
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-mixed-validation.R")
source("tests/testthat/helper-mixed-covariates.R")
destination <- "dev/mixed-covariate-validation/pi-audit"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
hashes <- tools::md5sum(c("dev/mixed-covariate-validation/pi-audit/check.R",
  "tests/testthat/helper-mixed-covariates.R", list.files("R", pattern = "mixed", full.names = TRUE)))
results <- attempts <- list()
for (i in seq_along(modes <- c("estimated_state", "equilibrium", "n_gen"))) {
  d <- mixed_covariate_data(N_families = 500L, seed = 72021L + i,
    contract = mixed_covariate_contract(n_pc = 1L, sex_pc_association = .5),
    mode = modes[i], n_generations = 3L, family_variance = .5,
    sex_beta = c(child = .12, mother = -.06, father = .02))
  args <- list(dat = d, pi = .4, pi_se = .03, fit_constraints = modes[i],
    am_n_generations = 3L, fixed_effects = c("PC1", "sex", "age"),
    interactions = "sex", interaction_center = c(sex = 0), diagnostic = FALSE)
  f <- do.call(ravel_mixed_fit, args)
  sides <- lapply(c(-.005, .005), function(delta) {
    side <- args; side$pi <- .4 + delta; side$pi_se <- NULL; side$se <- FALSE
    side$start <- f$parameters; side$control <- list(pgtol = 1e-9, factr = 10)
    do.call(ravel_mixed_fit, side)
  })
  moments_unchanged <- all(vapply(sides, function(side)
    identical(side$settings$covariate$eta, f$settings$covariate$eta), logical(1)))
  derivative <- (sides[[2L]]$table$estimate - sides[[1L]]$table$estimate) / .01
  default <- if (is.null(f$pi_sensitivity)) rep(NA_real_, length(derivative)) else
    unlist(f$pi_sensitivity, use.names = FALSE)
  z <- data.frame(mode = modes[i], parameter = f$table$parameter, scale = f$table$scale,
    default_derivative = default, independent_derivative = derivative,
    discrepancy_in_reported_SE = .03 * abs(default - derivative) / f$table$se)
  D <- grep("^C_D_", names(f$parameters))
  attempts[[i]] <- data.frame(mode = modes[i], main_converged = f$converged,
    main_status = f$inference_status, minus_converged = sides[[1L]]$converged,
    plus_converged = sides[[2L]]$converged, covariate_moments_identical = moments_unchanged,
    D_parameter_count_refitted = length(D),
    maximum_D_change = max(abs(sides[[2L]]$parameters[D] - sides[[1L]]$parameters[D])))
  results[[i]] <- z
  utils::write.csv(do.call(rbind, attempts), file.path(destination, "attempts.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, results), file.path(destination, "derivatives.csv"), row.names = FALSE)
  cat(modes[i], f$inference_status, "Q unchanged", moments_unchanged, "\n")
}
writeLines(sub("[[:space:]]+$", "", c(capture.output(sessionInfo()), "", capture.output(hashes))),
           file.path(destination, "session.txt"))
print(do.call(rbind, attempts), row.names = FALSE)
