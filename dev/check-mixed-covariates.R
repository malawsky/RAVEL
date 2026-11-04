# Independent recovery checks for the optional correlated-child-covariate model.
# Rscript dev/check-mixed-covariates.R --families=1500 --se=true
# Add --resume=true to continue every pending seed without replacing failures.
# Explicit retries use a fresh output directory, --replicate-ids=5,7 and
# --warm-start-from=path/to/original/parameters.csv.
# All attempted fits remain in the output, including boundaries and failures.
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-mixed-validation.R")
source("tests/testthat/helper-mixed-covariates.R")
arguments <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  x <- arguments[startsWith(arguments, paste0("--", name, "="))]
  if (length(x)) sub(paste0("^--", name, "="), "", tail(x, 1L)) else default
}
families <- as.integer(option("families", "1000"))
replicates <- as.integer(option("replicates", "1"))
start_replicate <- as.integer(option("start-replicate", "1"))
resume <- identical(tolower(option("resume", "false")), "true")
with_se <- identical(tolower(option("se", "true")), "true")
destination <- option("output", file.path("dev", "mixed-covariate-validation", "smoke"))
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
scenarios <- list(
  continuous_pc = list(data = list(), fit = list()),
  continuous_equilibrium = list(data = list(mode = "equilibrium"),
    fit = list(fit_constraints = "equilibrium")),
  continuous_n_gen = list(data = list(mode = "n_gen", n_generations = 3L),
    fit = list(fit_constraints = "n_gen", am_n_generations = 3L)),
  long_age_sex = list(data = list(visits = 3L, individual_variance = .2,
    age_beta = c(child = .025, mother = -.008, father = 0),
    sex_beta = c(child = .12, mother = -.06, father = .02)),
    fit = list(individual_random = TRUE, interactions = c("age", "sex"),
      interaction_center = c(age = 10, sex = 0))),
  long_random_slope = list(data = list(visits = 4L, individual_variance = .2,
    slope_variance = .003, intercept_slope_covariance = .012,
    age_beta = c(child = .025, mother = -.008, father = 0),
    sex_beta = c(child = .12, mother = -.06, father = .02)),
    fit = list(individual_random = TRUE, random_slope = "age", interactions = c("age", "sex"),
      interaction_center = c(age = 10, sex = 0))),
  binary_pc = list(data = list(outcome = "binary_liability", residual_variance = 1),
    fit = list(outcome = "binary_liability")),
  external_pi = list(data = list(), fit = list(pi_se = .03), external_pi_sd = .03),
  external_pi_equilibrium = list(data = list(mode = "equilibrium"),
    fit = list(fit_constraints = "equilibrium", pi_se = .03), external_pi_sd = .03),
  external_pi_n_gen = list(data = list(mode = "n_gen", n_generations = 3L),
    fit = list(fit_constraints = "n_gen", am_n_generations = 3L, pi_se = .03), external_pi_sd = .03),
  extended_fid = list(data = list(families_per_FID = 3L, extended_fraction = .6), fit = list()),
  twenty_pc = list(data = list(contract = mixed_covariate_contract(n_pc = 20L,
    sex_pc_association = .5)), fit = list()))
selected <- strsplit(option("scenarios", "continuous_pc"), ",", fixed = TRUE)[[1L]]
replicate_ids <- seq.int(start_replicate, replicates)
if (nzchar(option("replicate-ids", "")))
  replicate_ids <- as.integer(strsplit(option("replicate-ids", ""), ",", fixed = TRUE)[[1L]])
stopifnot(families >= 50L, replicates >= 1L, start_replicate >= 1L,
  start_replicate <= replicates, !anyNA(replicate_ids), !anyDuplicated(replicate_ids),
  all(replicate_ids >= 1L & replicate_ids <= replicates), all(selected %in% names(scenarios)))
source_hashes <- tools::md5sum(c("dev/check-mixed-covariates.R",
  "tests/testthat/helper-mixed-covariates.R", "tests/testthat/helper-mixed-validation.R",
  list.files("R", pattern = "mixed", full.names = TRUE)))
attempts <- estimates <- states <- parameters <- list()
table_names <- c("attempts", "estimates", "states", "parameters")
checkpoint <- file.path(destination, ".checkpoint.rds")
if (!resume && (file.exists(file.path(destination, "attempts.csv")) || file.exists(checkpoint)))
  stop("Output already contains attempted fits. Use --resume=true or a fresh output directory.")
if (resume) {
  previous <- if (file.exists(checkpoint)) readRDS(checkpoint) else
    stats::setNames(lapply(table_names, function(name) {
      path <- file.path(destination, paste0(name, ".csv"))
      if (file.exists(path)) utils::read.csv(path, stringsAsFactors = FALSE) else data.frame()
    }), table_names)
  for (name in table_names) if (nrow(previous[[name]])) {
    if (name == "attempts" && !"source_run_id" %in% names(previous[[name]]))
      previous[[name]]$source_run_id <- "historical_pre_natural_polish"
    assign(name, list(previous[[name]]))
  }
  if (length(attempts)) {
    existing <- attempts[[1L]]
    stopifnot(all(existing$N_families[existing$scenario %in% selected] == families),
      all(existing$se_requested[existing$scenario %in% selected] == with_se))
    for (name in c("estimates", "states"))
      stopifnot(setequal(paste(existing$scenario, existing$replicate),
        paste(previous[[name]]$scenario, previous[[name]]$replicate)))
  }
}
completed <- if (length(attempts)) with(do.call(rbind, attempts), paste(scenario, replicate)) else character()
warm_path <- option("warm-start-from", "")
warm_parameters <- if (nzchar(warm_path)) utils::read.csv(warm_path, stringsAsFactors = FALSE) else NULL
source_run_id <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
run_directory <- file.path(destination, "source-runs", source_run_id)
dir.create(run_directory, recursive = TRUE, showWarnings = FALSE)
provenance <- c(paste("source_run_id:", source_run_id), paste("resume:", resume),
  paste("arguments:", paste(arguments, collapse = " ")),
  paste("completed records retained:", length(completed)),
  capture.output(sessionInfo()), "", "Source hashes:", capture.output(source_hashes))
writeLines(sub("[[:space:]]+$", "", provenance), file.path(run_directory, "session.txt"))
writeLines("running", file.path(run_directory, "status.txt"))
write_progress <- function() {
  tables <- stats::setNames(lapply(table_names, function(name) {
    value <- get(name, inherits = TRUE)
    if (length(value)) do.call(rbind, value) else data.frame()
  }), table_names)
  # One atomic checkpoint protects against interruption between CSV exports.
  # It contains only result tables, not data or fitted objects, and is removed
  # after a clean completion once all CSVs have been exported.
  temporary <- file.path(destination, ".checkpoint-pending.rds")
  saveRDS(tables, temporary)
  if (!file.rename(temporary, checkpoint)) stop("Cannot replace the simulation checkpoint.")
  for (name in table_names) if (nrow(tables[[name]])) {
    temporary <- file.path(destination, paste0(".", name, "-pending.csv"))
    utils::write.csv(tables[[name]], temporary, row.names = FALSE)
    if (!file.rename(temporary, file.path(destination, paste0(name, ".csv"))))
      stop("Cannot replace simulation output table: ", name)
  }
}
write_progress()
for (scenario_name in selected) for (replicate in replicate_ids) {
  if (paste(scenario_name, replicate) %in% completed) next
  scenario <- scenarios[[scenario_name]]
  seed <- 481000L + match(scenario_name, names(scenarios)) * 10000L + replicate
  data_arguments <- modifyList(list(N_families = families, seed = seed,
    contract = mixed_covariate_contract(sex_pc_association = .5)), scenario$data)
  d <- do.call(mixed_covariate_data, data_arguments)
  truth <- attr(d, "truth")
  args <- modifyList(list(dat = d, pi = truth$population$pi,
    fixed_effects = truth$covariates, covariate_model = "joint",
    diagnostic = FALSE, se = with_se), scenario$fit)
  if (!is.null(warm_parameters)) {
    w <- warm_parameters[warm_parameters$scenario == scenario_name & warm_parameters$replicate == replicate, ]
    if (!nrow(w)) stop("No saved start parameters for ", scenario_name, " replicate ", replicate)
    args$start <- stats::setNames(w$estimate, w$parameter)
  }
  if (!is.null(scenario$external_pi_sd)) {
    set.seed(seed + 1L)
    args$pi <- stats::rnorm(1L, truth$population$pi, scenario$external_pi_sd)
    stopifnot(args$pi > 0, args$pi < 1)
  }
  messages <- character()
  start_time <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(do.call(ravel_mixed_fit, args), warning = function(w) {
    messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = identity)
  elapsed <- proc.time()[["elapsed"]] - start_time
  error <- inherits(f, "error")
  attempts[[length(attempts) + 1L]] <- data.frame(scenario = scenario_name,
    replicate = replicate, seed = seed, N_families = families, N_rows = nrow(d),
    pi_supplied = args$pi, se_requested = with_se, error = error,
    converged = !error && f$converged,
    optimizer_converged = !error && f$optimizer_converged,
    inference_status = if (error) "error" else f$inference_status,
    elapsed_seconds = elapsed,
    note = paste(c(if (error) conditionMessage(f) else f$messages, messages), collapse = "; "),
    source_run_id = source_run_id)
  wanted <- c("beta_D", "beta_M", "beta_F", unlist(lapply(args$interactions,
    function(z) paste0(c("beta_D", "beta_M", "beta_F"), ":", z))))
  expected <- do.call(rbind, lapply(names(truth$coefficients), function(scale)
    data.frame(scale = scale, parameter = wanted,
      truth = unname(truth$coefficients[[scale]][wanted]))))
  expected$estimate <- expected$se <- NA_real_
  target_state <- c(rho = truth$population$rho_total,
    parent = truth$population$parent_variance, child = truth$population$Sigma_TT["child", "child"])
  state_estimate <- rep(NA_real_, length(target_state))
  if (!error) {
    idx <- match(paste(expected$scale, expected$parameter), paste(f$table$scale, f$table$parameter))
    expected$estimate <- f$table$estimate[idx]; expected$se <- f$table$se[idx]
    state_estimate <- unlist(f$state[names(target_state)], use.names = FALSE)
    parameters[[length(parameters) + 1L]] <- data.frame(scenario = scenario_name,
      replicate = replicate, parameter = names(f$parameters), estimate = unname(f$parameters))
  }
  states[[length(states) + 1L]] <- data.frame(scenario = scenario_name, replicate = replicate,
    parameter = names(target_state), truth = target_state, estimate = state_estimate)
  expected$finite_point <- (!error && isTRUE(f$converged)) & is.finite(expected$estimate)
  expected$finite_inference <- expected$finite_point & is.finite(expected$se) & expected$se > 0
  expected$covered <- expected$finite_inference & abs(expected$estimate - expected$truth) <= qnorm(.975) * expected$se
  estimates[[length(estimates) + 1L]] <- cbind(scenario = scenario_name,
    replicate = replicate, expected)
  write_progress()
  completed <- c(completed, paste(scenario_name, replicate))
  cat(scenario_name, replicate, "/", replicates,
      if (error) "error" else f$inference_status, sprintf("%.2fs", elapsed), "\n")
}
writeLines("completed", file.path(run_directory, "status.txt"))
writeLines(sub("[[:space:]]+$", "", provenance), file.path(destination, "session.txt"))
unlink(checkpoint)
print(do.call(rbind, attempts), row.names = FALSE)
cat("All attempted fits retained in", normalizePath(destination), "\n")
