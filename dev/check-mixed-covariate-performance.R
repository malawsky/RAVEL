# Reproducible runtime/memory-shape check for the child-PC conditional extension.
# Run from the repository root. The default matches approximately 33k outcome
# rows / 30k maternal families, with 20 PCs, sex and age as fixed effects.
args <- commandArgs(trailingOnly = TRUE)
argument <- function(name, default) {
  value <- args[startsWith(args, paste0("--", name, "="))]
  if (length(value)) sub(paste0("^--", name, "="), "", value[1L]) else default
}
n_family <- as.integer(argument("families", "30000"))
n_pc <- as.integer(argument("pcs", "20"))
do_fit <- !"--no-fit" %in% args
do_se <- "--se" %in% args
output <- argument("output", "dev/mixed-validation/covariate-performance")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source_files <- c(list.files("R", pattern = "\\.R$", full.names = TRUE),
  "tests/testthat/helper-mixed-validation.R", "tests/testthat/helper-mixed-covariates.R",
  "dev/check-mixed-covariate-performance.R")
utils::write.csv(data.frame(file = source_files, md5 = unname(tools::md5sum(source_files))),
  file.path(output, "source-md5.csv"), row.names = FALSE)
for (file in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(file)
source("tests/testthat/helper-mixed-validation.R")
source("tests/testthat/helper-mixed-covariates.R")
timings <- data.frame(stage = character(), elapsed_seconds = numeric(),
  stringsAsFactors = FALSE)
record <- function(stage, elapsed) {
  timings <<- rbind(timings, data.frame(stage = stage, elapsed_seconds = elapsed))
  utils::write.csv(timings, file.path(output, "timings.csv"), row.names = FALSE)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stage, elapsed, "seconds\n")
  flush.console()
}
started <- proc.time()[["elapsed"]]
sizes <- c(rep(1L, n_family - ceiling(n_family / 10)), rep(2L, ceiling(n_family / 10)))
d <- mixed_covariate_data(N_families = n_family, sib_sizes = sizes, seed = 724109L,
  contract = mixed_covariate_contract(n_pc = n_pc))
fixed <- c(paste0("PC", seq_len(n_pc)), "sex", "age")
truth <- attr(d, "truth")
utils::write.csv(data.frame(scale = "corrected_raw", parameter = c("beta_D", "beta_M", "beta_F"),
  estimate = unname(truth$beta_raw)), file.path(output, "truth.csv"), row.names = FALSE)
record("generate", proc.time()[["elapsed"]] - started)
mapping <- c(child1 = "PGS1_child", child2 = "PGS2_child", mother1 = "PGS1_mother",
  mother2 = "PGS2_mother", father1 = "PGS1_father", father2 = "PGS2_father")
started <- proc.time()[["elapsed"]]
p <- .ravel_mixed_prepare(d, "continuous", "Y", mapping, "M_ID", "IID", "FID",
  fixed, NULL, NULL, TRUE, FALSE, NULL, "estimated_state", 1L, .4, TRUE,
  covariate_model = "joint")
record("prepare_genetic_records", proc.time()[["elapsed"]] - started)
started <- proc.time()[["elapsed"]]
p <- .ravel_mixed_covariate_prepare(p)
record("prepare_covariate_moments", proc.time()[["elapsed"]] - started)
started <- proc.time()[["elapsed"]]
p <- .ravel_mixed_covariate_compress(p)
record("compress_covariance_patterns", proc.time()[["elapsed"]] - started)
shape <- data.frame(outcomes = nrow(d), families = n_family, pcs = n_pc,
  covariate_features = p$spec$covariate$q, covariance_patterns = length(p$groups),
  covariance_parameters = length(p$start), mean_parameters = length(p$spec$mean_names),
  moment_parameters = length(p$spec$covariate$eta),
  prepared_size_bytes = as.numeric(object.size(p)),
  dense_design_cross_present = any(vapply(p$groups, function(x) !is.null(x$design_cross), logical(1))))
utils::write.csv(shape, file.path(output, "shape.csv"), row.names = FALSE)
ev <- .ravel_mixed_covariate_evaluator(p, p$spec)
started <- proc.time()[["elapsed"]]
value <- ev(p$start)
stopifnot(value$valid)
record("one_likelihood", proc.time()[["elapsed"]] - started)
started <- proc.time()[["elapsed"]]
gradient <- ev(p$start, derivatives = TRUE)
stopifnot(gradient$valid, all(is.finite(gradient$gradient)))
record("one_likelihood_and_gradient", proc.time()[["elapsed"]] - started)
rm(ev, value, gradient, p); invisible(gc())
if (do_fit) {
  started <- proc.time()[["elapsed"]]
  fit <- ravel_mixed_fit(d, fixed_effects = fixed, fit_constraints = "estimated_state",
    pi = .4, covariate_model = "joint", family_random = TRUE,
    se = do_se, diagnostic = FALSE, keep_data = FALSE)
  record(if (do_se) "public_fit_with_SE" else "public_fit_without_SE", proc.time()[["elapsed"]] - started)
  utils::write.csv(data.frame(converged = fit$converged,
    optimizer_converged = fit$optimizer_converged,
    inference_status = fit$inference_status, optimizer_code = fit$optimizer$convergence,
    objective_evaluations = unname(fit$optimizer$counts[["function"]]),
    gradient_evaluations = unname(fit$optimizer$counts[["gradient"]]),
    max_gradient_per_family = max(abs(fit$gradient)) / fit$n_families,
    rho_estimate = fit$state$rho, parent_variance_estimate = fit$state$parent,
    result_size_bytes = as.numeric(object.size(fit))), file.path(output, "fit.csv"), row.names = FALSE)
  utils::write.csv(fit$table, file.path(output, "coefficients.csv"), row.names = FALSE)
  writeLines(fit$messages, file.path(output, "messages.txt"))
  stopifnot(fit$converged, !do_se || identical(fit$inference_status, "ok"))
}
session <- capture.output(sessionInfo())
writeLines(sub("[[:blank:]]+$", "", session), file.path(output, "session.txt"))
cat("Performance check completed.\n")
