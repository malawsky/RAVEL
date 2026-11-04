# Sparse-sibling covariance regression checks and one independent fitted example.
# Run from package root: Rscript dev/check-sparse-sibling-covariates.R
# All data are synthetic; the single fitted example does not establish coverage.
pkgload::load_all(".", quiet = TRUE)
RAVEL:::.ensure_lavaan_runtime()
source("tests/testthat/helper-mixed-validation.R")
source("tests/testthat/helper-mixed-covariates.R")
arguments <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  x <- arguments[startsWith(arguments, paste0("--", name, "="))]
  if (length(x)) sub(paste0("^--", name, "="), "", tail(x, 1L)) else default
}
destination <- option("output", "dev/mixed-validation/sparse-sibling-covariates")
do_fit <- identical(tolower(option("fit", "true")), "true")
fit_scenario <- match.arg(option("fit-scenario", "ordinary"), c("ordinary", "shared", "shared_long"))
moment_replicates <- as.integer(option("moment-replicates", "0"))
stopifnot(moment_replicates >= 0L)
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
source_files <- c("dev/check-sparse-sibling-covariates.R",
  "tests/testthat/helper-mixed-covariates.R", "tests/testthat/helper-mixed-validation.R",
  list.files("R", pattern = "mixed", full.names = TRUE))
source_hashes <- tools::md5sum(source_files)
writeLines(sub("[[:space:]]+$", "", c(capture.output(sessionInfo()),
  paste("Arguments:", paste(arguments, collapse = " ")), "Source hashes at launch:",
  capture.output(source_hashes))), file.path(destination, "session.txt"))
checks <- list()
check <- function(ok, label, value = NA_real_) {
  checks[[length(checks) + 1L]] <<- data.frame(check = label,
    status = if (isTRUE(ok)) "pass" else "FAIL", value = value)
  write.csv(do.call(rbind, checks), file.path(destination, "checks.csv"), row.names = FALSE)
  if (!isTRUE(ok)) stop(label)
}
moments <- RAVEL:::.ravel_mixed_covariate_moments
synthetic <- function(nf, np, q, rho, seed, shared = integer()) {
  set.seed(seed)
  sizes <- c(rep(1L, nf - np), rep(2L, np))
  family <- rep(seq_len(nf), sizes)
  between <- matrix(rnorm(nf * q), nf, q)
  C <- sqrt(rho) * between[family, , drop = FALSE] +
    sqrt(1 - rho) * matrix(rnorm(length(family) * q), length(family), q)
  if (length(shared)) C[, shared] <- between[family, shared, drop = FALSE]
  C <- scale(C)
  lapply(split(seq_len(nrow(C)), family), function(i) C[i, , drop = FALSE])
}
old_moments <- function(values) {
  n <- sum(vapply(values, nrow, integer(1)))
  pairs <- sum(vapply(values, function(x) nrow(x) * (nrow(x) - 1L), numeric(1)))
  mean <- Reduce(`+`, lapply(values, colSums)) / n
  self <- Reduce(`+`, lapply(values, crossprod)) / n - tcrossprod(mean)
  sibling <- Reduce(`+`, lapply(values, function(x)
    tcrossprod(colSums(x)) - crossprod(x))) / pairs - tcrossprod(mean)
  list(Qself = self, Qsib = sibling)
}
minimum <- function(x) min(eigen((x + t(x)) / 2, symmetric = TRUE, only.values = TRUE)$values)

# High-dimensional analogue of the user's 95% singleton / 5% two-child design.
values <- synthetic(10000L, 500L, 20L, .9, 712L)
old <- old_moments(values)
check(minimum(old$Qself - old$Qsib) < 0, "old_empirical_subtraction_reproduces_failure",
  minimum(old$Qself - old$Qsib))
new <- moments(values)
check(minimum(new$Qself - new$Qsib) > 0, "twenty_PC_within_covariance_positive_definite",
  minimum(new$Qself - new$Qsib))
check(minimum(new$Qself + new$Qsib) > 0, "twenty_PC_pair_mean_covariance_positive_definite",
  minimum(new$Qself + new$Qsib))
check(all(is.finite(new$moment_scores)), "twenty_PC_first_stage_influence_finite")
check(max(abs(colSums(new$moment_scores))) < 1e-4,
  "twenty_PC_first_stage_scores_centered", max(abs(colSums(new$moment_scores))))

# Perturb a singleton and a sibling family. The estimator's reported influence
# must include the re-estimation of both covariance blocks and the common mean.
small <- synthetic(800L, 40L, 3L, .8, 173L)
base <- moments(small)
weights <- rep(1, length(small))
size <- vapply(small, nrow, integer(1))
rows <- c(which(size == 1L)[3L], which(size == 2L)[3L])
influence_results <- lapply(rows, function(i) {
  h <- 1e-3
  plus <- minus <- weights; plus[i] <- plus[i] + h; minus[i] <- minus[i] - h
  derivative <- (moments(small, weights = plus, influence = FALSE)$eta -
    moments(small, weights = minus, influence = FALSE)$eta) / (2 * h)
  expected <- base$moment_scores[i, ] / length(small)
  error <- max(abs(derivative - expected))
  check(error < 3e-5, paste0("family_weight_influence_nchild_", size[i]), error)
  data.frame(family = i, children = size[i], parameter = names(base$eta),
    finite_difference = unname(derivative), influence = unname(expected), row.names = NULL)
})
write.csv(do.call(rbind, influence_results), file.path(destination, "influence.csv"), row.names = FALSE)

# A genuinely shared covariate remains in the model instead of acquiring an
# artificial within-family variance or being deleted from the outcome design.
shared_values <- synthetic(800L, 40L, 3L, .8, 174L, shared = 2L)
shared <- moments(shared_values)
check(identical(as.integer(shared$shared), 2L), "exact_family_shared_covariate_identified")
check(max(abs((shared$Qself - shared$Qsib)[2L, ])) < 1e-10,
  "exact_family_shared_covariate_has_zero_within_variance")
check(all(is.finite(shared$moment_scores)), "shared_covariate_influence_finite")

if (moment_replicates > 0L) {
  contract <- mixed_covariate_contract(n_pc = 2L, pc_between = .95,
    pc_within = .015, baseline_pc_association = .2)
  oracle <- mixed_covariate_family(n_children = 2L, contract = contract)
  target <- c(PC1_variance = oracle$Q_self[1L, 1L],
    PC1_sibling_covariance = oracle$Q_sib[1L, 1L],
    PC1_within_variance = oracle$Q_self[1L, 1L] - oracle$Q_sib[1L, 1L],
    sex_variance = oracle$Q_self[3L, 3L], age_variance = oracle$Q_self[4L, 4L])
  attempts <- estimates <- list()
  for (replicate in seq_len(moment_replicates)) {
    seed <- 812000L + replicate
    d <- mixed_covariate_data(N_families = 2000L,
      sib_sizes = c(rep(1L, 1900L), rep(2L, 100L)), seed = seed, contract = contract)
    C <- as.matrix(d[c("PC1", "PC2", "sex", "age")])
    values <- lapply(split(seq_len(nrow(d)), d$M_ID), function(i) C[i, , drop = FALSE])
    f <- tryCatch(moments(values), error = identity)
    error <- inherits(f, "error")
    attempts[[replicate]] <- data.frame(replicate = replicate, seed = seed, error = error,
      old_within_min_eigenvalue = minimum(with(old_moments(values), Qself - Qsib)),
      note = if (error) conditionMessage(f) else "")
    point <- se <- rep(NA_real_, length(target))
    if (!error) {
      h <- f$moment_scores; q <- ncol(C)
      diagonal <- match(paste(seq_len(q), seq_len(q)),
        paste(f$triangle[, 1L], f$triangle[, 2L]))
      self_h <- h[, f$index$self[diagonal], drop = FALSE] -
        sweep(h[, f$index$mean, drop = FALSE], 2L, 2 * f$mean, "*")
      sib_h <- h[, f$index$sibling[diagonal], drop = FALSE] -
        sweep(h[, f$index$mean, drop = FALSE], 2L, 2 * f$mean, "*")
      H <- cbind(self_h[, 1L], sib_h[, 1L], self_h[, 1L] - sib_h[, 1L],
        self_h[, 3L], self_h[, 4L])
      point <- c(f$Qself[1L, 1L], f$Qsib[1L, 1L], f$Qself[1L, 1L] - f$Qsib[1L, 1L],
        f$Qself[3L, 3L], f$Qself[4L, 4L])
      se <- sqrt(colSums(H^2)) / length(values)
    }
    estimates[[replicate]] <- data.frame(replicate = replicate, parameter = names(target),
      truth = unname(target), estimate = point, se = se,
      covered_95 = is.finite(se) & abs(point - target) <= qnorm(.975) * se, row.names = NULL)
    write.csv(do.call(rbind, attempts), file.path(destination, "moment-attempts.csv"), row.names = FALSE)
    write.csv(do.call(rbind, estimates), file.path(destination, "moment-estimates.csv"), row.names = FALSE)
  }
  results <- do.call(rbind, estimates)
  summary <- do.call(rbind, lapply(split(results, results$parameter), function(z) {
    usable <- is.finite(z$estimate) & is.finite(z$se)
    data.frame(parameter = z$parameter[1L], attempts = nrow(z), regular = sum(usable),
      mean_bias = mean(z$estimate[usable] - z$truth[usable]),
      empirical_sd = sd(z$estimate[usable]), mean_se = mean(z$se[usable]),
      coverage_all_attempts = mean(z$covered_95),
      coverage_regular = mean(z$covered_95[usable]))
  }))
  write.csv(summary, file.path(destination, "moment-summary.csv"), row.names = FALSE)
  print(summary, row.names = FALSE)
}

if (do_fit) {
  long <- fit_scenario == "shared_long"
  nf <- if (long) 2000L else 4000L
  np <- as.integer(.05 * nf)
  fit_seed <- switch(fit_scenario, ordinary = 712L, shared = 713L, shared_long = 714L)
  d <- mixed_covariate_data(N_families = nf,
    sib_sizes = c(rep(1L, nf - np), rep(2L, np)), seed = fit_seed,
    contract = mixed_covariate_contract(n_pc = 2L, pc_between = .95,
      pc_within = .015, baseline_pc_association = .2),
    visits = if (long) 3L else 1L, family_variance = if (long) 0 else .5,
    individual_variance = if (long) .3 else 0,
    sex_beta = if (long) c(child = .12, mother = -.06, father = .02) else
      c(child = 0, mother = 0, father = 0))
  truth <- attr(d, "truth")
  if (fit_scenario != "ordinary") {
    latent <- attr(d, "latent")
    set.seed(713001L)
    family <- match(latent$M_ID, unique(latent$M_ID))
    PC <- .15 * latent$G_mother * truth$population$visible_raw_sd[2L] +
      .10 * latent$G_father * truth$population$visible_raw_sd[3L] +
      rnorm(length(unique(family)))[family]
    new_pc <- PC[match(d$IID, latent$IID)]
    d$Y <- d$Y + truth$fixed_effects["PC1"] * (new_pc - d$PC1)
    d$PC1 <- new_pc
    stopifnot(all(is.finite(d$PC1)), all(is.finite(d$Y)))
  }
  # Exercise factor IDs with unused levels, including the user's original setup.
  d$M_ID <- factor(d$M_ID, levels = c(unique(d$M_ID), "unused_family"))
  d$FID <- factor(d$FID)
  C <- scale(as.matrix(d[truth$covariates]))
  family_values <- lapply(split(seq_len(nrow(d)), droplevels(d$M_ID)),
    function(i) C[i, , drop = FALSE])
  previous <- old_moments(family_values)
  if (fit_scenario == "ordinary") check(minimum(previous$Qself - previous$Qsib) < 0,
    "actual_fit_fixture_failed_old_covariance", minimum(previous$Qself - previous$Qsib))
  notes <- character()
  start_time <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(ravel_mixed_fit(d, pi = truth$population$pi,
    fixed_effects = truth$covariates, covariate_model = "joint", family_random = !long,
    individual_random = long, interactions = if (long) "sex" else NULL,
    interaction_center = if (long) c(sex = 0) else NULL,
    diagnostic = FALSE, se = TRUE), warning = function(w) {
      notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity)
  elapsed <- proc.time()[["elapsed"]] - start_time
  error <- inherits(f, "error")
  attempt <- data.frame(scenario = fit_scenario, seed = fit_seed, N_families = nf, N_rows = nrow(d),
    singleton_families = nf - np, sibling_families = np, error = error,
    converged = !error && isTRUE(f$converged),
    inference_status = if (error) "error" else f$inference_status,
    elapsed_seconds = elapsed,
    note = paste(c(if (error) conditionMessage(f) else f$messages, notes), collapse = "; "))
  write.csv(attempt, file.path(destination, "attempts.csv"), row.names = FALSE)
  if (!error) {
    labels <- c("beta_D", "beta_M", "beta_F")
    expected <- truth$coefficients$corrected_std[c(labels, if (long) paste0(labels, ":sex"))]
    estimates <- f$table[f$table$scale == "corrected_std" & f$table$parameter %in% names(expected), ]
    estimates$truth <- unname(expected[estimates$parameter])
    estimates$covered_95 <- estimates$ci_low <= estimates$truth & estimates$ci_high >= estimates$truth
    write.csv(estimates, file.path(destination, "estimates.csv"), row.names = FALSE)
    write.csv(f$variance_components, file.path(destination, "variance-components.csv"), row.names = FALSE)
    write.csv(data.frame(parameter = names(f$parameters), estimate = f$parameters),
      file.path(destination, "parameters.csv"), row.names = FALSE)
  }
  check(!error, "sparse_sibling_fit_completes")
  check(isTRUE(f$converged), "sparse_sibling_fit_converged")
  check(identical(f$inference_status, "ok"), "sparse_sibling_fit_regular_inference")
  check(all(is.finite(estimates$se) & estimates$se > 0), "sparse_sibling_effect_SEs_finite")
  check(!is.null(f$uncertainty$nuisance_vcov) &&
      all(is.finite(f$uncertainty$nuisance_vcov)), "first_stage_uncertainty_propagated_to_FID_SEs")
  if (fit_scenario != "ordinary") {
    check(identical(as.integer(f$settings$covariate$shared), 1L),
      "full_fit_keeps_exact_shared_covariate")
    check("PC1" %in% f$settings$covariates, "shared_covariate_kept_in_outcome_adjustment")
  }
  if (long) {
    check(isTRUE(f$settings$individual_random) && !f$settings$family_random,
      "longitudinal_individual_effect_only")
    check(all(paste0(labels, ":sex") %in% estimates$parameter),
      "longitudinal_genetic_sex_interactions_reported")
  }
}
print(do.call(rbind, checks), row.names = FALSE)
