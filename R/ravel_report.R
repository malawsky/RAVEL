.report_core_labels <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f")

.report_packet_roles <- function(packet) {
  .parental_sib_roles(packet$settings$parental_sibs %||% "none")
}

.report_packet_labels <- function(packet) {
  c(.report_core_labels, .sib_core_labels(.report_packet_roles(packet)))
}

.report_parameters <- function(roles, scale) {
  c("beta_D", "beta_M", "beta_P",
    if (scale == "visible") "rho_visible" else if (scale == "corrected_std") "rho_corrected" else "cov_corrected",
    unname(.parental_sib_coef_map(roles)))
}

.report_transform_fingerprint <- function() {
  ns <- asNamespace("RAVEL")
  funs <- sort(ls(ns, all.names = TRUE))
  funs <- funs[grepl("^\\.(hidden|compute_scale|standardize_total|latent_corr|finite_difference|build_lavaan|ravel_estimated_state|sib_|parental_sib_(coef_map|roles)|report_packet|report_parameters)", funs)]
  path <- tempfile()
  on.exit(unlink(path))
  writeLines(unlist(lapply(funs, function(nm) c(nm, deparse(get(nm, ns))))), path)
  unname(tools::md5sum(path))
}

#' @export
ravel_correct_packet <- function(packet, h, am_correction = c("n_gen", "equilibrium", "estimated_state"),
                                  am_n_generations = 15L, ci_level = .95,
                                  covariate_mapping_confirmed = FALSE,
                                  uncertainty = c("auto", "lavaan", "bootstrap")) {
  if (is.character(packet)) packet <- readRDS(packet)
  history <- packet$settings$constraint_spec
  roles <- .report_packet_roles(packet)
  labels <- .report_packet_labels(packet)
  if (!is.null(history$mode) && history$mode %in% c("equilibrium", "n_gen", "estimated_state")) {
    if (missing(am_correction)) am_correction <- history$mode
    if (missing(am_n_generations) && history$mode == "n_gen")
      am_n_generations <- history$n_generations
  }
  am_correction <- match.arg(am_correction)
  if (am_correction == "estimated_state" && !identical(history$mode, "estimated_state"))
    stop("Estimated-state correction requires a packet from fit_constraints = 'estimated_state'.")
  if (am_correction == "estimated_state" && !isTRUE(packet$diagnostics$estimated_state_identified))
    stop("Estimated-state correction requires a packet with a validated identification check; regenerate the report.")
  if (am_correction == "estimated_state" &&
      abs(packet$theta["eta_m"] - packet$theta["eta_f"]) > 1e-5)
    stop("Estimated-state correction requires equal maternal and paternal transmission correlations.")
  uncertainty <- match.arg(uncertainty)
  .ravel_check_constraint_correction(packet$settings$constraint_spec, h,
    am_correction, am_n_generations)
  if (uncertainty == "auto") uncertainty <- if (identical(packet$settings$uncertainty, "bootstrap"))
    "bootstrap" else "lavaan"
  if (uncertainty == "bootstrap" && is.null(packet$bootstrap$core_draws))
    stop("This packet has no joint core-parameter bootstrap draws")
  if (uncertainty == "bootstrap") {
    core_draws <- packet$bootstrap$core_draws
    if (!is.matrix(core_draws) || !is.numeric(core_draws) || is.complex(core_draws) ||
        nrow(core_draws) == 0L ||
        is.null(colnames(core_draws)) || anyNA(colnames(core_draws)) ||
        anyDuplicated(colnames(core_draws)) ||
        !setequal(colnames(core_draws), labels))
      stop("Bootstrap core_draws must be a numeric matrix with at least one row and exactly these unique column names: ",
           paste(labels, collapse = ", "))
    core_draws <- core_draws[, labels, drop = FALSE]
  }
  stopifnot(identical(packet$schema, "ravel_visible_transfer_v1"),
            isTRUE(packet$diagnostics$admissible), h >= 0, h < 1,
            ci_level > 0, ci_level < 1)
  if (!identical(packet$software$transform_fingerprint, .report_transform_fingerprint())) {
    stop("The installed correction code differs from the collaborator's; align package versions before proceeding")
  }
  transform <- utils::getFromNamespace(".compute_scale_estimates", "RAVEL")
  jacobian <- utils::getFromNamespace(".finite_difference_jacobian", "RAVEL")
  theta <- packet$theta[labels]
  if (!is.numeric(packet$theta) || is.complex(packet$theta) ||
      is.null(names(packet$theta)) || anyDuplicated(names(packet$theta)) ||
      !setequal(names(packet$theta), labels) || any(!is.finite(theta)))
    stop("Packet core estimates must include exactly the parameters required by its parental_sibs setting")
  V <- .ravel_pi_validate_covariance(packet$vcov_core, labels, "Core covariance")
  evaluate <- function(x, block) {
    x <- setNames(as.numeric(x), labels)
    result <- transform(beta_vis = x[c("bD", "bm", "bf")], cov_cm = x["eta_m"],
      cov_cf = x["eta_f"], cov_mf = x["rhoMF"], h = h,
      am_correction = am_correction, am_n_generations = am_n_generations,
      sib = .sib_input_from_theta(x, roles))
    if (block == "corrected_std") {
      c(result$beta_total_assumed_std, rhoMF = result$rho_total_assumed, result$sib_total$std)
    } else {
      c(result$beta_total_assumed, rhoMF = result$cov_total_assumed, result$sib_total$raw)
    }
  }
  out <- lapply(c("corrected_std", "corrected_raw"), function(block) {
    fun <- function(x) evaluate(x, block)
    est <- fun(theta)
    if (any(!is.finite(est))) stop("This hidden-fraction/AM assumption is incompatible with the fitted moments")
    if (uncertainty == "bootstrap") {
      draws <- t(apply(core_draws, 1, function(x)
        tryCatch(if (all(is.finite(x))) fun(x) else rep(NA_real_, length(est)),
                 error = function(e) rep(NA_real_, length(est)))))
      draws[!apply(is.finite(draws), 1, all), ] <- NA_real_
      interval <- .ravel_report_boot_interval(est, draws, ci_level, packet$bootstrap$minimum_success)
    } else {
      J <- jacobian(fun, theta)
      if (any(!is.finite(J))) stop("Correction derivatives unavailable near an admissibility boundary")
      propagated <- J %*% V %*% t(J)
      if (any(!is.finite(propagated)) || any(diag(propagated) < -1e-10)) stop("Invalid propagated covariance")
      se <- sqrt(pmax(diag(propagated), 0))
      z <- stats::qnorm(1 - (1 - ci_level) / 2)
      interval <- data.frame(se = unname(se), ci_low = unname(est - z * se),
        ci_high = unname(est + z * se), p_value = unname(2 * stats::pnorm(-abs(est / se))),
        n_boot_success = NA_integer_, uncertainty_status = "ok")
    }
    cbind(data.frame(outcome_variable = packet$outcome$name,
      outcome_scale = if (identical(packet$outcome$type, "binary_liability"))
        "SEM_probit_liability" else "prepared_phenotype_SD",
      residual_model = if (is.null(packet$settings$residual_model)) "unspecified" else packet$settings$residual_model,
      scale = block,
      parameter = .report_parameters(roles, block),
      estimate = unname(est),
      inference = if (uncertainty == "bootstrap") "bootstrap_percentile_CI_normal_p" else "lavaan_delta_wald",
      h = h, am_correction = am_correction,
      am_n_generations = if (am_correction == "n_gen") am_n_generations else NA_integer_), interval)
  })
  result <- do.call(rbind, out)
  attr(result, "ci_level") <- ci_level
  attr(result, "constraint_spec") <- history
  if (am_correction == "estimated_state") {
    state <- .hidden_estimated_state(theta["rhoMF"], theta["eta_m"], theta["eta_f"], h)
    attr(result, "estimated_state") <- state$estimated_state
  }
  result
}

.ravel_report_vanilla <- function(dat, version, fixed_effects, scaling, response_scale,
                                  response_center, outcome_type, ci_level, outcome_name) {
  roles <- c("child", "mother", "father")
  pgs <- paste0("PGS", version, "_", roles)
  formula <- stats::reformulate(c(pgs, fixed_effects), response = "Y")
  fit <- if (outcome_type == "continuous") stats::lm(formula, data = dat) else
    stats::glm(formula, data = dat, family = stats::binomial(link = "probit"))
  b <- stats::coef(fit)
  V <- stats::vcov(fit)
  if (any(!is.finite(b)) || any(!is.finite(V)) ||
      (outcome_type != "continuous" && (!isTRUE(fit$converged) || isTRUE(fit$boundary)))) {
    stop("Vanilla PGS", version, " regression did not return a usable fit")
  }
  make_table <- function(est, covariance, predictor_scale, response_unit, normal = FALSE) {
    se <- sqrt(diag(covariance))
    statistic <- est / se
    critical <- if (outcome_type == "continuous" && !normal)
      stats::qt(1 - (1 - ci_level) / 2, fit$df.residual) else stats::qnorm(1 - (1 - ci_level) / 2)
    p <- if (outcome_type == "continuous" && !normal)
      2 * stats::pt(-abs(statistic), fit$df.residual) else 2 * stats::pnorm(-abs(statistic))
    parameter <- names(est)
    parameter[match(pgs, names(est))] <- c("beta_D", "beta_M", "beta_P")
    data.frame(outcome_variable = outcome_name, regression = paste0("PGS", version),
      scale = predictor_scale, outcome_scale = response_unit, parameter = parameter,
      estimate = unname(est), se = unname(se), ci_low = unname(est - critical * se),
      ci_high = unname(est + critical * se), p_value = unname(p), N = nrow(dat),
      inference = if (outcome_type == "continuous" && !normal) "lm_t" else
        if (normal) "glm_delta_wald_conditional_on_sample_moments" else "glm_probit_wald",
      stringsAsFactors = FALSE)
  }
  response_unit <- if (outcome_type == "continuous") "prepared_phenotype_SD" else "probit_residual_SD"
  common <- make_table(b, V, "child_reference_PGS_SD", response_unit)
  S <- stats::cov(stats::model.matrix(fit))
  vcov_maps <- list(child_reference = V)
  beta_maps <- list(child_reference = b)
  all_tables <- list(common)
  s <- scaling[match(pgs, scaling$model_column), ]
  A <- diag(length(b))
  dimnames(A) <- list(names(b), names(b))
  A[,] <- A * response_scale
  A[cbind(match(pgs, names(b)), match(pgs, names(b)))] <- response_scale / s$child_reference_sd
  A["(Intercept)", pgs] <- -response_scale * s$child_reference_mean / s$child_reference_sd
  raw <- drop(A %*% b)
  names(raw) <- names(b)
  raw["(Intercept)"] <- raw["(Intercept)"] + response_center
  rawV <- A %*% V %*% t(A)
  all_tables[[length(all_tables) + 1L]] <- make_table(raw, rawV, "raw_PGS_units",
    if (outcome_type == "continuous") "prepared_phenotype_raw_units" else "probit_residual_SD")
  beta_maps$raw <- raw
  vcov_maps$raw <- rawV
  A <- diag(length(b))
  dimnames(A) <- list(names(b), names(b))
  A[cbind(match(pgs, names(b)), match(pgs, names(b)))] <- sqrt(s$scaled_variance)
  A["(Intercept)", pgs] <- s$scaled_mean
  role_b <- drop(A %*% b)
  names(role_b) <- names(b)
  roleV <- A %*% V %*% t(A)
  all_tables[[length(all_tables) + 1L]] <- make_table(role_b, roleV, "own_role_PGS_SD", response_unit)
  beta_maps$own_role <- role_b
  vcov_maps$own_role <- roleV
  liability_variance <- NA_real_
  if (outcome_type == "binary_liability") {
    liability_variance <- 1 + drop(t(b) %*% S %*% b)
    divisor <- sqrt(liability_variance)
    J <- diag(length(b)) / divisor - outer(b, drop(S %*% b)) / divisor^3
    dimnames(J) <- list(names(b), names(b))
    stdV <- J %*% V %*% t(J)
    stdb <- b / divisor
    all_tables[[length(all_tables) + 1L]] <- make_table(stdb, stdV,
      "child_reference_PGS_SD", "model_implied_marginal_probit_liability_SD", normal = TRUE)
    beta_maps$marginal_liability <- stdb
    vcov_maps$marginal_liability <- stdV
    own_std <- drop(A %*% stdb)
    names(own_std) <- names(b)
    own_stdV <- A %*% stdV %*% t(A)
    all_tables[[length(all_tables) + 1L]] <- make_table(own_std, own_stdV,
      "own_role_PGS_SD", "model_implied_marginal_probit_liability_SD", normal = TRUE)
    beta_maps$own_role_marginal_liability <- own_std
    vcov_maps$own_role_marginal_liability <- own_stdV
  }
  measures <- c(N = stats::nobs(fit), df_residual = fit$df.residual,
                 log_likelihood = as.numeric(stats::logLik(fit)), AIC = stats::AIC(fit),
                 BIC = stats::BIC(fit), deviance = stats::deviance(fit))
  if (outcome_type == "continuous") {
    sm <- summary(fit)
    measures <- c(measures, R2 = sm$r.squared, adjusted_R2 = sm$adj.r.squared,
                  residual_variance = sm$sigma^2)
  } else {
    measures <- c(measures, null_deviance = fit$null.deviance,
                  deviance_fraction_reduction = 1 - fit$deviance / fit$null.deviance,
                  model_implied_liability_variance = liability_variance)
  }
  list(coefficients = do.call(rbind, all_tables), vcov = vcov_maps, estimates = beta_maps,
       design_covariance = S, fit_statistics = measures,
       formula = paste(deparse(formula), collapse = " "), converged = TRUE)
}

.ravel_report_covariates <- function(dat, covariates, categorical_covariates = character()) {
  if (anyDuplicated(categorical_covariates) || !all(categorical_covariates %in% covariates))
    stop("categorical_covariates must be a unique subset of covariates")
  automatic <- covariates[vapply(dat[covariates], function(x)
    is.factor(x) || is.character(x) || is.logical(x) ||
      (is.numeric(x) && setequal(unique(x), c(0, 1))), logical(1))]
  categorical <- union(automatic, categorical_covariates)
  encoded <- dat[covariates]
  contrasts <- counts <- list()
  for (nm in categorical) {
    value <- encoded[[nm]]
    encoded[[nm]] <- if (is.factor(value))
      factor(value, levels = levels(droplevels(value)), ordered = FALSE) else factor(value)
    if (nlevels(encoded[[nm]]) < 2L) stop("Categorical covariate has fewer than two observed levels: ", nm)
    contrasts[[nm]] <- stats::contr.treatment(levels(encoded[[nm]]), base = 1L)
    freq <- base::table(encoded[[nm]])
    counts[[nm]] <- data.frame(variable = nm, level = names(freq), n = as.integer(freq),
      proportion = as.numeric(freq) / nrow(dat), reference = seq_along(freq) == 1L,
      stringsAsFactors = FALSE)
  }
  C <- if (length(covariates)) stats::model.matrix(stats::reformulate(covariates), encoded,
    contrasts.arg = if (length(contrasts)) contrasts else NULL) else
      matrix(1, nrow(dat), 1L, dimnames = list(NULL, "(Intercept)"))
  if (any(!is.finite(C)) || qr(C)$rank != ncol(C)) stop("Covariate design must be finite and full rank")
  assignment <- attr(C, "assign")
  original <- if (ncol(C) > 1L) covariates[assignment[-1]] else character()
  map <- data.frame(model_column = if (ncol(C) > 1L) paste0("adjustment_", seq_len(ncol(C) - 1L)) else character(),
    original_term = colnames(C)[-1], original_variable = original,
    coding = ifelse(original %in% categorical, "treatment_dummy", "numeric"),
    reference_level = vapply(original, function(nm)
      if (nm %in% categorical) levels(encoded[[nm]])[1] else NA_character_, character(1)),
    stringsAsFactors = FALSE)
  empty <- data.frame(variable = character(), level = character(), n = integer(),
                      proportion = numeric(), reference = logical())
  list(design = C, map = map, categorical = categorical,
       counts = if (length(counts)) do.call(rbind, counts) else empty)
}

.ravel_report_numeric_descriptives <- function(dat, columns) {
  if (!length(columns)) return(data.frame(variable = character(), n = integer(), mean = numeric(),
    sd = numeric(), variance = numeric(), minimum = numeric(), q05 = numeric(), median = numeric(),
    q95 = numeric(), maximum = numeric()))
  do.call(rbind, lapply(columns, function(nm) {
    x <- dat[[nm]]
    x <- x[is.finite(x)]
    data.frame(variable = nm, n = length(x), mean = mean(x), sd = stats::sd(x), variance = stats::var(x),
      minimum = min(x), q05 = unname(stats::quantile(x, .05)), median = stats::median(x),
      q95 = unname(stats::quantile(x, .95)), maximum = max(x))
  }))
}

#' @export
ravel_dataframe_report <- function(df, outcome = "Y", covariates = character(),
                                    categorical_covariates = character(),
                                    outcome_type = c("continuous", "binary_liability"),
                                    continuous_adjustment = c("joint", "residualize_Y"),
                                    pgs_columns = c(PGS1_child = "pgs1_child", PGS1_mother = "pgs1_mother",
                                      PGS1_father = "pgs1_father", PGS2_child = "pgs2_child",
                                      PGS2_mother = "pgs2_mother", PGS2_father = "pgs2_father"),
                                    ci_level = .95, out_dir = NULL,
                                    residual_models = c("fixed", "free"),
                                    uncertainty = c("auto", "lavaan", "bootstrap"),
                                    n_boot = 0L, seed = NULL, make_plot = TRUE,
                                    report_name = NULL,
                                    fit_constraints = c("none", "equilibrium", "n_gen", "estimated_state"),
                                    h = 0, am_n_generations = 1L,
                                    parental_sibs = c("none", "maternal", "paternal", "both")) {
  fit_constraints <- match.arg(fit_constraints)
  parental_sibs <- match.arg(parental_sibs)
  sib_roles <- .parental_sib_roles(parental_sibs)
  if (fit_constraints != "none") {
    if (missing(residual_models)) residual_models <- "fixed"
    if (!identical(residual_models, "fixed"))
      stop("Constrained fitting requires residual_models = 'fixed'.")
    if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h < 0 || h >= 1)
      stop("h must be a single finite number in [0, 1).")
    if (fit_constraints == "n_gen" && (!is.numeric(am_n_generations) ||
        length(am_n_generations) != 1L || !is.finite(am_n_generations) ||
        am_n_generations < 1 || am_n_generations > .Machine$integer.max ||
        am_n_generations != floor(am_n_generations)))
      stop("am_n_generations must be a positive integer for constrained finite-generation fitting.")
  }
  report_stage <- "preprocessing"
  report_model_status <- setNames(rep("not_attempted", length(residual_models) + 2L),
                                  c(residual_models, "PGS1", "PGS2"))
  withCallingHandlers({
  stopifnot(requireNamespace("RAVEL", quietly = TRUE), requireNamespace("lavaan", quietly = TRUE))
  outcome_type <- match.arg(outcome_type)
  continuous_adjustment <- match.arg(continuous_adjustment)
  uncertainty <- match.arg(uncertainty)
  if (length(n_boot) != 1L || !is.numeric(n_boot) || !is.finite(n_boot) ||
      n_boot != as.integer(n_boot) || n_boot < 0 || (n_boot > 0 && n_boot < 10))
    stop("n_boot must be zero or an integer of at least 10; use 500 or more for final inference")
  if (uncertainty == "auto") uncertainty <- if (n_boot > 0) "bootstrap" else "lavaan"
  if (uncertainty == "bootstrap" && n_boot == 0) stop("Bootstrap uncertainty requires n_boot >= 10")
  if (!is.null(seed) && (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed))) stop("seed must be a nonnegative integer")
  if (length(make_plot) != 1L || !is.logical(make_plot) || is.na(make_plot)) stop("make_plot must be TRUE or FALSE")
  if (!is.null(report_name)) {
    if (is.null(out_dir) || length(report_name) != 1L || is.na(report_name) ||
        !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", report_name))
      stop("report_name requires out_dir and must be a single safe folder name")
    out_dir <- file.path(out_dir, report_name)
  }
  if (!is.character(residual_models) || !length(residual_models) || anyNA(residual_models) ||
      anyDuplicated(residual_models) || !all(residual_models %in% c("fixed", "free")))
    stop("residual_models must select 'fixed', 'free', or both without duplicates")
  trio_columns <- c("PGS1_child", "PGS1_mother", "PGS1_father", "PGS2_child", "PGS2_mother", "PGS2_father")
  sib_columns <- .parental_sib_columns(sib_roles)
  expected <- c(trio_columns, sib_columns)
  if (missing(pgs_columns) && length(sib_columns))
    pgs_columns <- c(pgs_columns, stats::setNames(tolower(sib_columns), sib_columns))
  if (!is.data.frame(df) || anyDuplicated(names(df))) stop("df must have unique column names")
  if (!setequal(names(pgs_columns), expected) || anyDuplicated(unname(pgs_columns)))
    stop("pgs_columns must uniquely map the expected score columns: ", paste(expected, collapse = ", "))
  pgs_columns <- pgs_columns[expected]
  required <- c(unname(pgs_columns), outcome, covariates)
  if (length(outcome) != 1L || anyDuplicated(required) || !all(required %in% names(df)))
    stop("Check the outcome, PGS column mappings and covariate names")
  if (!is.numeric(ci_level) || length(ci_level) != 1L || !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1)
    stop("ci_level must lie between zero and one")
  if (!is.null(out_dir) && (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir) ||
      !nzchar(out_dir) || file.exists(out_dir)))
    stop("Choose a new output directory; nothing will be overwritten")
  if (any(make.names(covariates) != covariates)) stop("Use syntactically valid covariate names")
  if (!all(vapply(df[c(unname(pgs_columns), outcome)], is.numeric, logical(1))))
    stop("Scores and outcome must be numeric; binary outcomes must be explicitly coded 0/1")
  missing_field <- lapply(required, function(nm) {
    value <- df[[nm]]
    is.na(value) | if (is.numeric(value)) !is.finite(value) else rep(FALSE, length(value))
  })
  names(missing_field) <- required
  missing_pgs <- Reduce(`|`, missing_field[unname(pgs_columns[trio_columns])])
  missing_y <- missing_field[[outcome]]
  missing_covariates <- if (length(covariates)) Reduce(`|`, missing_field[covariates]) else rep(FALSE, nrow(df))
  keep <- !missing_pgs & !missing_y & !missing_covariates
  exclusion_counts <- data.frame(reason = c("missing_or_nonfinite_PGS", "missing_or_nonfinite_outcome",
    "missing_or_nonfinite_covariate"),
    n_any = c(sum(missing_pgs), sum(missing_y), sum(missing_covariates)),
    n_excluded_sequential = c(sum(missing_pgs), sum(!missing_pgs & missing_y),
                              sum(!missing_pgs & !missing_y & missing_covariates)))
  if (outcome_type == "binary_liability" && length(sib_columns)) {
    missing_sib <- Reduce(`|`, missing_field[unname(pgs_columns[sib_columns])])
    excluded_sib <- sum(keep & missing_sib)
    exclusion_counts <- rbind(exclusion_counts, data.frame(reason = "missing_or_nonfinite_sibling_PGS",
      n_any = sum(missing_sib), n_excluded_sequential = excluded_sib))
    keep <- keep & !missing_sib
    if (excluded_sib) warning("Binary reports exclude families with missing requested sibling scores from all models.")
  }
  exclusion_counts$n_remaining <- nrow(df) - cumsum(exclusion_counts$n_excluded_sequential)
  missingness <- data.frame(variable = required, n_input = nrow(df),
    n_missing_or_nonfinite = vapply(missing_field, sum, numeric(1)), stringsAsFactors = FALSE)
  source <- droplevels(df[keep, , drop = FALSE])
  n <- nrow(source)
  if (length(sib_columns)) {
    for (nm in unname(pgs_columns[sib_columns])) source[[nm]][!is.finite(source[[nm]])] <- NA_real_
  }
  if (n < 100L) stop("Fewer than 100 complete trios")
  for (id in intersect(c("FID", "child_id", "mother_id", "father_id"), names(source))) {
    values <- source[[id]][!is.na(source[[id]])]
    if (anyDuplicated(values)) stop("Repeated ", id, ": this workflow assumes independent families")
  }
  if (outcome_type == "binary_liability" && !setequal(unique(source[[outcome]]), c(0, 1)))
    stop("Binary outcome must contain both zero and one")
  design <- .ravel_report_covariates(source, covariates, categorical_covariates)
  C <- design$design
  covariate_map <- design$map
  fe <- covariate_map$model_column
  dat <- data.frame(FID = seq_len(n), Y = source[[outcome]])
  if (length(fe)) for (j in seq_along(fe)) dat[[fe[j]]] <- C[, j + 1L]
  scaling <- list()
  for (k in 1:2) {
    ref <- source[[pgs_columns[[paste0("PGS", k, "_child")]]]]
    center <- mean(ref)
    divisor <- stats::sd(ref)
    if (!is.finite(divisor) || divisor <= 0) stop("Child PGS", k, " has no variance")
    for (role in c("child", "mother", "father", unname(c(ms = "mat_sib", fs = "pat_sib")[sib_roles]))) {
      nm <- paste0("PGS", k, "_", role)
      raw <- source[[pgs_columns[[nm]]]]
      if (!is.finite(stats::var(raw, na.rm = TRUE)) || stats::var(raw, na.rm = TRUE) <= 0) stop(nm, " has no variance")
      dat[[nm]] <- (raw - center) / divisor
      scaling[[nm]] <- data.frame(outcome_variable = outcome, model_column = nm,
        source_column = pgs_columns[[nm]], PGS = paste0("PGS", k), role = role,
        raw_mean = mean(raw, na.rm = TRUE), raw_sd = stats::sd(raw, na.rm = TRUE), raw_variance = stats::var(raw, na.rm = TRUE),
        child_reference_mean = center, child_reference_sd = divisor,
        scaled_mean = mean(dat[[nm]], na.rm = TRUE), scaled_variance = stats::var(dat[[nm]], na.rm = TRUE))
    }
  }
  scaling <- do.call(rbind, scaling)
  rownames(scaling) <- NULL
  scaling$vanilla_reference_mean <- scaling$raw_mean
  scaling$vanilla_reference_sd <- scaling$raw_sd
  raw_y_mean <- mean(dat$Y)
  raw_y_sd <- stats::sd(dat$Y)
  response_center <- 0
  response_scale <- 1
  y_adjustment <- NULL
  if (outcome_type == "continuous") {
    if (continuous_adjustment == "residualize_Y") {
      residual_fit <- stats::lm.fit(C, dat$Y)
      y_adjustment <- residual_fit$coefficients
      dat$Y <- residual_fit$residuals
      fe <- character()
    }
    response_center <- mean(dat$Y)
    response_scale <- stats::sd(dat$Y)
    if (!is.finite(response_scale) || response_scale <= 0) stop("Prepared outcome has no variance")
    dat$Y <- (dat$Y - response_center) / response_scale
  }
  fixed_corr <- c(e1_mf = 0, e1_mc = .5, e1_fc = .5, e2_mf = 0, e2_mc = .5, e2_fc = .5,
    .parental_sib_residual_benchmarks(sib_roles))
  if (length(fe)) {
    model <- utils::getFromNamespace(".build_model", "RAVEL")(
      fixed_effects = fe, fixed_crossperson_corr = fixed_corr, parental_sibs = parental_sibs)
    if (!all(vapply(c("Gc", "Gm", "Gf"), function(g)
      grepl(paste(g, "~~", fe[1]), model, fixed = TRUE), logical(1)))) {
      stop("This workflow requires the local covariate candidate; the installed RAVEL model lacks genetic-covariate covariances")
    }
  }
  wanted <- c("outcome_variable", "regression", "scale", "outcome_scale", "parameter",
              "estimate", "se", "ci_low", "ci_high", "p_value", "N", "inference",
              "residual_model", "fit_status")
  vanilla <- setNames(vector("list", 2L), c("PGS1", "PGS2"))
  for (k in 1:2) {
    report_stage <- paste0("PGS", k)
    report_model_status[report_stage] <- "running"
    vanilla[[k]] <- .ravel_report_vanilla(dat, k, fe, scaling,
      response_scale, response_center, outcome_type, ci_level, outcome)
    report_model_status[report_stage] <- "ok"
  }
  report_stage <- "report"
  vanilla_summary <- do.call(rbind, lapply(vanilla, function(x) {
    tab <- x$coefficients
    tab <- tab[tab$scale == "own_role_PGS_SD" &
      tab$outcome_scale != "model_implied_marginal_probit_liability_SD" &
      tab$parameter %in% c("beta_D", "beta_M", "beta_P"), , drop = FALSE]
    tab$residual_model <- "not_applicable"
    tab$fit_status <- "ok"
    tab[wanted]
  }))
  stats_table <- function(x, name, residual_model = "not_applicable")
    data.frame(outcome_variable = outcome, regression = name, residual_model = residual_model,
      statistic = names(x), value = unname(x), stringsAsFactors = FALSE)
  vanilla_statistics <- do.call(rbind, lapply(names(vanilla), function(nm)
    stats_table(vanilla[[nm]]$fit_statistics, nm)))
  covariate_map$mean <- if (nrow(covariate_map)) colMeans(C[, -1, drop = FALSE]) else numeric()
  covariate_map$sd <- if (nrow(covariate_map)) apply(C[, -1, drop = FALSE], 2, stats::sd) else numeric()
  sample_counts <- data.frame(statistic = c("n_input_rows", "n_complete_trios", "n_excluded_rows",
    "n_recorded_unique_families", "n_covariate_variables", "n_covariate_design_columns"),
    value = c(nrow(df), n, nrow(df) - n,
      if ("FID" %in% names(source)) length(unique(source$FID[!is.na(source$FID)])) else NA_real_,
      length(covariates), ncol(C) - 1L))
  if (outcome_type == "binary_liability") sample_counts <- rbind(sample_counts,
    data.frame(statistic = c("n_cases", "n_controls", "case_fraction"),
               value = c(sum(dat$Y == 1), sum(dat$Y == 0), mean(dat$Y))))
  numeric_columns <- c(unname(pgs_columns), outcome,
    setdiff(covariates[vapply(source[covariates], is.numeric, logical(1))], design$categorical))
  descriptives <- list(sample_counts = sample_counts, exclusions = exclusion_counts,
    missingness = missingness,
    numeric = .ravel_report_numeric_descriptives(source, numeric_columns),
    categorical = design$counts)
  moment_columns <- c(expected, "Y", covariate_map$model_column)
  input_moments <- list(mean = colMeans(dat[moment_columns], na.rm = TRUE),
    covariance = stats::cov(dat[moment_columns], use = "pairwise.complete.obs"))
  if (length(sib_columns)) {
    counts <- vapply(sib_roles, function(role) sum(stats::complete.cases(dat[.parental_sib_columns(role)])), integer(1))
    sample_counts <- rbind(sample_counts, data.frame(statistic = paste0("n_with_", unname(c(ms = "maternal_sib", fs = "paternal_sib")[sib_roles])), value = counts))
    descriptives$sample_counts <- sample_counts
  }
  pgs_covariance <- input_moments$covariance[expected, expected, drop = FALSE]
  pgs_variances <- data.frame(outcome_variable = outcome, model_column = expected,
    variance = unname(diag(pgs_covariance)), N = colSums(!is.na(dat[expected])))
  fingerprint_path <- tempfile()
  on.exit(unlink(fingerprint_path), add = TRUE)
  writeLines(c(deparse(utils::getFromNamespace("ravel_fit", "RAVEL")),
               deparse(utils::getFromNamespace(".build_outcome", "RAVEL"))), fingerprint_path)
  common_packet <- list(schema = "ravel_visible_transfer_v1",
    input_moments = input_moments, pgs_covariance = pgs_covariance, pgs_variances = pgs_variances,
    raw_input_moments = list(mean = colMeans(cbind(source[c(unname(pgs_columns), outcome)], C[, -1, drop = FALSE]), na.rm = TRUE),
      covariance = stats::cov(cbind(source[c(unname(pgs_columns), outcome)], C[, -1, drop = FALSE]), use = "pairwise.complete.obs")),
    pgs_summary = scaling, covariate_map = covariate_map, covariates = covariates,
    categorical_covariates = design$categorical, descriptives = descriptives,
    vanilla = vanilla,
    outcome = list(name = outcome, type = outcome_type,
      adjustment = if (outcome_type == "continuous") continuous_adjustment else "joint_probit",
      raw_mean = raw_y_mean, raw_sd = raw_y_sd, preparation_center = response_center,
      preparation_sd = response_scale, adjustment_coefficients = y_adjustment,
      prepared_variance = stats::var(dat$Y),
      counts_0_1 = if (outcome_type == "binary_liability") base::table(dat$Y) else NULL,
      sample_prevalence = if (outcome_type == "binary_liability") mean(dat$Y) else NULL),
    settings = list(n_initial = nrow(df), n_complete = n, excluded = nrow(df) - n,
      h = if (fit_constraints == "none") 0 else h, n_boot = 0, ci_level = ci_level,
      fit_constraints = fit_constraints, parental_sibs = parental_sibs,
      moment_missingness = if (length(sib_columns) && outcome_type == "continuous")
        "pairwise observed descriptive moments; SEM uses FIML" else "complete observed inputs",
      am_n_generations =
        if (fit_constraints %in% c("equilibrium", "estimated_state")) NA_integer_ else am_n_generations,
      fixed_effects = fe, requires_covariate_mapping_review = length(fe) > 0,
      independent_families_assumed = TRUE,
      pgs_scaling_sem = "same child reference mean and SD for all roles within each score version",
      pgs_scaling_vanilla_summary = "each role standardized using its own mean and SD",
      identifiers_checked = intersect(c("FID", "child_id", "mother_id", "father_id"), names(source)),
      fitting_engine = "RAVEL::ravel_fit"),
    software = list(RAVEL_version = as.character(utils::packageVersion("RAVEL")),
      lavaan_version = as.character(utils::packageVersion("lavaan")), R_version = R.version.string,
      transform_fingerprint = .report_transform_fingerprint(),
      sem_fingerprint = unname(tools::md5sum(fingerprint_path))))
  fitted <- setNames(vector("list", length(residual_models)), residual_models)
  for (specification in residual_models) {
    report_stage <- specification
    report_model_status[specification] <- "running"
    fitted[[specification]] <- .ravel_report_sem(dat, common_packet, specification, fixed_corr, wanted,
      vanilla_summary, vanilla_statistics, stats_table, kappa(C))
    report_model_status[specification] <- if (isTRUE(fitted[[specification]]$diagnostics$admissible)) "ok" else "failed"
  }
  report_stage <- "report"
  summary <- do.call(rbind, c(lapply(fitted, `[[`, "summary"), list(vanilla_summary)))
  statistics <- do.call(rbind, c(lapply(fitted, `[[`, "fit_statistics"), list(vanilla_statistics)))
  rownames(summary) <- rownames(statistics) <- NULL
  packets <- lapply(fitted, `[[`, "packet")
  diagnostics <- lapply(fitted, `[[`, "diagnostics")
  status <- do.call(rbind, lapply(residual_models, function(nm) {
    d <- diagnostics[[nm]]
    data.frame(residual_model = nm, status = if (isTRUE(d$admissible)) "ok" else "failed",
      converged = d$converged, admissible = d$admissible, error = d$error, stringsAsFactors = FALSE)
  }))
  residuals <- do.call(rbind, lapply(fitted, `[[`, "residuals"))
  out <- list(summary = summary, scaling = scaling, fit_statistics = statistics,
              pgs_covariance = pgs_covariance, pgs_variances = pgs_variances,
              descriptives = descriptives, covariate_coding = covariate_map,
              reliability = lapply(packets, `[[`, "reliability"), residual_parameters = residuals,
              diagnostics = diagnostics, model_status = status, vanilla = vanilla,
              packets = packets, packet = packets[["fixed"]])
  out$summary_model_based <- out$summary
  out["summary_bootstrap"] <- list(NULL)
  out["bootstrap"] <- list(NULL)
  out$settings <- list(uncertainty = uncertainty, n_boot = n_boot, seed = seed,
    ci_level = ci_level, outcome_type = outcome_type, fit_constraints = fit_constraints, parental_sibs = parental_sibs,
    h = if (fit_constraints == "none") 0 else h, am_n_generations =
      if (fit_constraints %in% c("equilibrium", "estimated_state")) NA_integer_ else am_n_generations)
  if (n_boot > 0) {
    boot_data <- source[unname(pgs_columns)]
    names(boot_data) <- expected
    boot_data$Y <- source[[outcome]]
    if (nrow(covariate_map)) for (j in seq_len(nrow(covariate_map)))
      boot_data[[covariate_map$model_column[j]]] <- C[, j + 1L]
    out <- .ravel_report_bootstrap(out, boot_data, n_boot, seed, ci_level,
      outcome_type, continuous_adjustment, covariate_map$model_column, residual_models)
    if (uncertainty == "bootstrap") out$summary <- out$summary_bootstrap
  }
  for (nm in residual_models) if (!is.null(out$packets[[nm]])) {
    out$packets[[nm]]$settings$uncertainty <- uncertainty
    out$packets[[nm]]$settings$n_boot <- n_boot
    out$packets[[nm]]$summary_model_based <- out$packets[[nm]]$summary
    selected <- out$summary$regression != "RAVEL" | out$summary$residual_model == nm
    out$packets[[nm]]$summary <- out$summary[selected, , drop = FALSE]
    out$packets[[nm]]$visible_table_model_based <- out$packets[[nm]]$visible_table
    selected_visible <- out$summary[out$summary$regression == "RAVEL" & out$summary$residual_model == nm, ]
    out$packets[[nm]]$visible_table <- .ravel_report_visible_inference(
      out$packets[[nm]]$visible_table, selected_visible, uncertainty)
    out$packets[[nm]]$vanilla <- out$vanilla
  }
  out["packet"] <- list(out$packets[["fixed"]])
  out["plot"] <- list(if (make_plot) ravel_plot_report(out) else NULL)
  out$files <- character()
  if (!is.null(out_dir)) {
    if (dir.exists(out_dir)) stop("Output directory already exists; nothing will be overwritten")
    if (!dir.create(out_dir, recursive = TRUE)) stop("Could not create output directory")
    out_dir <- normalizePath(out_dir, mustWork = TRUE)
    utils::write.csv(out$summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
    utils::write.csv(out$summary_model_based, file.path(out_dir, "summary_model_based.csv"), row.names = FALSE)
    if (!is.null(out$bootstrap)) {
      utils::write.csv(out$summary_bootstrap, file.path(out_dir, "summary_bootstrap.csv"), row.names = FALSE)
      utils::write.csv(out$bootstrap$status, file.path(out_dir, "bootstrap_status.csv"), row.names = FALSE)
    }
    utils::write.csv(scaling, file.path(out_dir, "pgs_scaling.csv"), row.names = FALSE)
    utils::write.csv(statistics, file.path(out_dir, "fit_statistics.csv"), row.names = FALSE)
    utils::write.csv(pgs_covariance, file.path(out_dir, "pgs_input_covariance.csv"), row.names = TRUE)
    utils::write.csv(pgs_variances, file.path(out_dir, "pgs_input_variances.csv"), row.names = FALSE)
    utils::write.csv(status, file.path(out_dir, "model_status.csv"), row.names = FALSE)
    utils::write.csv(residuals, file.path(out_dir, "residual_covariances_correlations.csv"), row.names = FALSE)
    for (nm in residual_models) if (!is.null(packets[[nm]])) {
      utils::write.csv(packets[[nm]]$vcov_core,
        file.path(out_dir, paste0("core_parameter_vcov_", nm, ".csv")), row.names = TRUE)
      utils::write.csv(packets[[nm]]$parameter_estimates,
        file.path(out_dir, paste0("SEM_parameters_model_based_", nm, ".csv")), row.names = FALSE)
      utils::write.csv(packets[[nm]]$reliability,
        file.path(out_dir, paste0("reliability_", nm, ".csv")), row.names = FALSE)
    }
    utils::write.csv(sample_counts, file.path(out_dir, "sample_counts.csv"), row.names = FALSE)
    utils::write.csv(exclusion_counts, file.path(out_dir, "exclusions.csv"), row.names = FALSE)
    utils::write.csv(missingness, file.path(out_dir, "missingness.csv"), row.names = FALSE)
    utils::write.csv(descriptives$numeric, file.path(out_dir, "numeric_descriptives.csv"), row.names = FALSE)
    utils::write.csv(descriptives$categorical, file.path(out_dir, "categorical_descriptives.csv"), row.names = FALSE)
    utils::write.csv(covariate_map, file.path(out_dir, "covariate_coding.csv"), row.names = FALSE)
    for (nm in names(out$vanilla)) {
      utils::write.csv(out$vanilla[[nm]]$coefficients,
        file.path(out_dir, paste0("vanilla_coefficients_", nm, ".csv")), row.names = FALSE)
      if (!is.null(out$vanilla[[nm]]$coefficients_bootstrap))
        utils::write.csv(out$vanilla[[nm]]$coefficients_bootstrap,
          file.path(out_dir, paste0("vanilla_coefficients_bootstrap_", nm, ".csv")), row.names = FALSE)
    }
    saveRDS(out, file.path(out_dir, "ravel_report.rds"))
    if (make_plot) tryCatch(.ravel_report_save_plot(out$plot, out_dir), error = function(e)
      warning("Estimates were saved, but plot export failed: ", conditionMessage(e)))
    paths <- list.files(out_dir, full.names = TRUE)
    out$files <- stats::setNames(paths, basename(paths))
    saveRDS(out, file.path(out_dir, "ravel_report.rds"))
  }
  out
  }, error = function(e) {
    if (inherits(e, "ravel_report_error")) stop(e)
    report_model_status[report_model_status == "running"] <- "failed"
    stop(structure(list(message = conditionMessage(e), call = conditionCall(e),
      report_stage = report_stage, report_model_status = report_model_status),
      class = unique(c("ravel_report_error", class(e)))))
  })
}

.ravel_report_row_key <- function(tab) {
  columns <- intersect(c("regression", "scale", "outcome_scale", "parameter", "residual_model"), names(tab))
  do.call(paste, c(tab[columns], sep = "|"))
}

.ravel_report_visible_inference <- function(visible, selected, uncertainty) {
  order <- match(visible$parameter, selected$parameter)
  for (column in intersect(c("se", "ci_low", "ci_high", "p_value", "inference",
                             "n_boot_success", "uncertainty_status"), names(selected)))
    visible[[column]] <- selected[[column]][order]
  if (uncertainty == "bootstrap") {
    visible$uncertainty_source <- "bootstrap"
    visible$p_value_method <- "normal_approx_from_bootstrap_se"
  }
  visible
}

.ravel_report_boot_interval <- function(estimate, draws, ci_level, minimum_success) {
  stopifnot(is.matrix(draws), ncol(draws) == length(estimate))
  result <- lapply(seq_along(estimate), function(j) {
    x <- draws[is.finite(draws[, j]), j]
    valid <- length(x) >= minimum_success && is.finite(estimate[j])
    se <- if (valid) stats::sd(x) else NA_real_
    bounds <- if (valid) stats::quantile(x, c((1 - ci_level) / 2, (1 + ci_level) / 2),
                                       names = FALSE) else c(NA_real_, NA_real_)
    data.frame(se = se, ci_low = bounds[1], ci_high = bounds[2],
      p_value = if (valid && se > 0) 2 * stats::pnorm(-abs(estimate[j] / se)) else NA_real_,
      n_boot_success = length(x), uncertainty_status = if (valid) "ok" else "insufficient_bootstrap")
  })
  do.call(rbind, result)
}

.ravel_report_boot_table <- function(tab, draws, ci_level, minimum_success) {
  interval <- .ravel_report_boot_interval(tab$estimate, draws, ci_level, minimum_success)
  for (nm in names(interval)) tab[[nm]] <- interval[[nm]]
  tab$inference <- "bootstrap_percentile_CI_normal_p"
  tab
}

.ravel_report_bootstrap <- function(out, dat, n_boot, seed, ci_level, outcome_type,
                                    continuous_adjustment, covariates, residual_models) {
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  replicate_seeds <- sample.int(.Machine$integer.max, n_boot)
  minimum_success <- max(10L, ceiling(.8 * n_boot))
  allocate <- function(keys) matrix(NA_real_, n_boot, length(keys),
    dimnames = list(as.character(seq_len(n_boot)), keys))
  summary_draws <- allocate(.ravel_report_row_key(out$summary))
  labels <- c(.report_core_labels, .sib_core_labels(.parental_sib_roles(out$settings$parental_sibs %||% "none")))
  core_draws <- stats::setNames(lapply(residual_models, function(x) allocate(labels)), residual_models)
  vanilla_draws <- lapply(out$vanilla, function(x) allocate(.ravel_report_row_key(x$coefficients)))
  status <- vector("list", n_boot)
  pgs_names <- names(dat)[grepl("^PGS[12]_", names(dat))]
  pgs_map <- stats::setNames(pgs_names, pgs_names)
  for (i in seq_len(n_boot)) {
    set.seed(replicate_seeds[i])
    sample <- dat[sample.int(nrow(dat), nrow(dat), replace = TRUE), , drop = FALSE]
    fitted <- tryCatch(ravel_dataframe_report(sample, outcome = "Y", covariates = covariates,
      pgs_columns = pgs_map, outcome_type = outcome_type, continuous_adjustment = continuous_adjustment,
      residual_models = residual_models, ci_level = ci_level, make_plot = FALSE,
      fit_constraints = out$settings$fit_constraints %||% "none",
      parental_sibs = out$settings$parental_sibs %||% "none",
      h = out$settings$h %||% 0, am_n_generations = out$settings$am_n_generations %||% 1L), error = identity)
    if (inherits(fitted, "error")) {
      models <- c(residual_models, "PGS1", "PGS2")
      fit_status <- setNames(rep("not_attempted", length(models)), models)
      error_stage <- "unknown"
      if (inherits(fitted, "ravel_report_error")) {
        fit_status[] <- fitted$report_model_status[models]
        error_stage <- fitted$report_stage
      }
      status[[i]] <- data.frame(replicate = i, model = c(residual_models, "PGS1", "PGS2"),
        status = unname(fit_status), error = conditionMessage(fitted),
        error_stage = error_stage, draw_retained = FALSE)
      next
    }
    summary_draws[i, ] <- fitted$summary$estimate[match(colnames(summary_draws), .ravel_report_row_key(fitted$summary))]
    for (nm in residual_models) if (!is.null(fitted$packets[[nm]]))
      core_draws[[nm]][i, ] <- fitted$packets[[nm]]$theta[labels]
    for (nm in names(vanilla_draws)) {
      tab <- fitted$vanilla[[nm]]$coefficients
      vanilla_draws[[nm]][i, ] <- tab$estimate[match(colnames(vanilla_draws[[nm]]), .ravel_report_row_key(tab))]
    }
    status[[i]] <- rbind(data.frame(replicate = i, model = fitted$model_status$residual_model,
      status = fitted$model_status$status, error = fitted$model_status$error,
      error_stage = ifelse(fitted$model_status$status == "ok", "", fitted$model_status$residual_model),
      draw_retained = fitted$model_status$status == "ok"),
      data.frame(replicate = i, model = c("PGS1", "PGS2"), status = "ok", error = "",
        error_stage = "", draw_retained = TRUE))
  }
  out$bootstrap <- list(n_requested = n_boot, minimum_success = minimum_success, seed = seed,
    replicate_seeds = replicate_seeds, status = do.call(rbind, status),
    summary_draws = summary_draws, core_draws = core_draws, vanilla_draws = vanilla_draws,
    resampling = "whole families with complete core trios, retaining sibling missingness; preprocessing and all models refitted per draw",
    interval = "percentile", p_value = "normal approximation using bootstrap SE; not a bootstrap null test")
  out$summary_bootstrap <- .ravel_report_boot_table(out$summary, summary_draws, ci_level, minimum_success)
  for (nm in names(out$vanilla)) out$vanilla[[nm]]$coefficients_bootstrap <-
    .ravel_report_boot_table(out$vanilla[[nm]]$coefficients, vanilla_draws[[nm]], ci_level, minimum_success)
  for (nm in residual_models) if (!is.null(out$packets[[nm]])) out$packets[[nm]]$bootstrap <-
    list(n_requested = n_boot, minimum_success = minimum_success, seed = seed,
      core_draws = core_draws[[nm]], status = out$bootstrap$status[out$bootstrap$status$model == nm, ])
  if (any(out$summary_bootstrap$uncertainty_status != "ok"))
    warning("Some estimates have insufficient usable bootstrap fits; see summary_bootstrap and bootstrap$status")
  out
}

#' @export
ravel_plot_report <- function(out, scale = NULL) {
  direct <- is.data.frame(out)
  tab <- if (direct) out else out$summary
  # Legacy hand-built report objects predate explicit predictor-scale columns.
  if (!direct && is.data.frame(tab) && !"scale" %in% names(tab))
    tab$scale <- if ("regression" %in% names(tab))
      ifelse(tab$regression == "RAVEL", "visible", "own_role_PGS_SD") else "visible"
  if (!is.data.frame(tab) || !all(c("parameter", "scale", "estimate", "ci_low", "ci_high") %in% names(tab)))
    stop("Supply a RAVEL report or a coefficient table with scale, parameter, estimate and interval columns.")
  settings <- if (direct) attr(out, "settings") else out$settings
  if (is.null(scale) && direct) scale <- if ("corrected_std" %in% tab$scale) "corrected_std" else "visible"
  if (!is.null(scale)) {
    if (!is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% tab$scale)
      stop("scale must name one available scale in the coefficient table.")
    tab <- tab[tab$scale == scale, , drop = FALSE]
  }
  tab$parameter[tab$parameter == "beta_F"] <- "beta_P"
  tab <- tab[tab$parameter %in% c("beta_D", "beta_M", "beta_P", "beta_MS", "beta_FS", "delta_M", "delta_F"), , drop = FALSE]
  if (!nrow(tab)) stop("No coefficient estimates to plot")
  if (!"regression" %in% names(tab)) tab$regression <- "RAVEL"
  if (!"residual_model" %in% names(tab)) tab$residual_model <- "unspecified"
  if (!"outcome_variable" %in% names(tab)) tab$outcome_variable <- "Y"
  if (length(unique(tab$outcome_variable)) != 1L)
    stop("Plot one outcome at a time; subset the table before plotting.")
  corrected <- tab$scale %in% c("corrected_std", "corrected_raw")
  if (!"am_correction" %in% names(tab) && !is.null(settings$am_correction))
    tab$am_correction <- settings$am_correction
  if (!any(c("am_n_generations", "n_generations") %in% names(tab)) && !is.null(settings$am_n_generations))
    tab$am_n_generations <- settings$am_n_generations
  if (!any(c("pi", "h", "hidden_fraction") %in% names(tab)) && !is.null(settings$h)) tab$h <- settings$h
  tab$model <- ifelse(tab$regression == "RAVEL",
    paste0(ifelse(corrected, "Corrected RAVEL", "Visible RAVEL"), " (", tab$residual_model, ")"),
    paste0(tab$regression, " trio regression"))
  # Keep each sensitivity scenario on its own row; never overlay distinct pi/n values.
  for (nm in intersect(c("am_correction", "am_n_generations", "n_generations",
                         "pi", "h", "hidden_fraction", "pi_se"), names(tab))) {
    show <- corrected & !is.na(tab[[nm]])
    if (any(show)) {
      label <- switch(nm, am_n_generations = "n", n_generations = "n",
                      hidden_fraction = "pi", h = "pi", pi_se = "SE(pi)", nm)
      value <- format(tab[[nm]][show], trim = TRUE)
      value <- if (nm == "am_correction") gsub("_", " ", value, fixed = TRUE) else paste0(label, "=", value)
      value[value == "n gen"] <- "finite generations"
      tab$model[show] <- paste0(tab$model[show], "; ", value)
    }
  }
  mode <- settings$fit_constraints
  if (!is.null(mode) && mode != "none") {
    show <- tab$regression == "RAVEL" & !corrected
    tab$model[show] <- paste0(tab$model[show], "; ", gsub("_", " ", mode, fixed = TRUE))
  }
  key <- paste(tab$model, tab$parameter, sep = "|")
  if (anyDuplicated(key)) stop("Multiple estimates share a plot row; select one outcome, scale and uncertainty source.")
  inference <- unique(tab$inference %||% tab$uncertainty_source)
  pi_delta <- any(grepl("pi_delta_wald", inference))
  uncertainty <- if (pi_delta) "pi_delta_wald" else settings$uncertainty %||%
    if (any(grepl("bootstrap", inference))) "bootstrap" else "lavaan"
  level <- attr(out, "ci_level") %||% settings$ci_level %||%
    attr(out, "pi_uncertainty")$correction_ci_level %||% unique(tab$ci_level)
  binary <- identical(settings$outcome_type, "binary_liability") ||
    any(grepl("probit|liability", tab$outcome_scale %||% ""))
  raw <- all(tab$scale == "corrected_raw")
  caption <- paste0(if (length(level) == 1L && is.finite(level)) paste0(round(100 * level), "% ") else "",
    "intervals: ", if (pi_delta) "joint delta-Wald, including external pi uncertainty" else
    if (uncertainty == "mixed") "mixed model-based and bootstrap percentile; see inference column" else
    if (uncertainty == "bootstrap") "whole-trio bootstrap (percentile)" else "model-based (lm/glm; lavaan delta)")
  structure(list(data = tab, outcome = unique(tab$outcome_variable),
    ci_level = level, uncertainty = uncertainty, inference = inference, caption = caption,
    scale = scale, raw = raw, binary = binary,
    pi_uncertainty = attr(out, "pi_uncertainty"),
    estimated_state = attr(out, "estimated_state")), class = "ravel_report_plot")
}

#' @export
print.ravel_report_plot <- function(x, ...) {
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old))
  groups <- if (x$binary) intersect(c("PGS", "RAVEL"),
    unique(ifelse(x$data$regression == "RAVEL", "RAVEL", "PGS"))) else "all"
  roles <- c(beta_D = "Offspring", beta_M = "Mother", beta_P = "Father",
    beta_MS = "Maternal sibling", delta_M = "Mother-specific contrast",
    beta_FS = "Paternal sibling", delta_F = "Father-specific contrast")
  roles <- roles[names(roles) %in% x$data$parameter]
  corrected <- any(x$data$scale %in% c("corrected_std", "corrected_raw"))
  label_width <- if (corrected) 26L else 29L
  columns <- min(3L, length(roles))
  rows_per_group <- ceiling(length(roles) / columns)
  graphics::par(mfrow = c(length(groups) * rows_per_group, columns), mar = c(3.8, if (corrected) 9.4 else 8.4, 2.7, 1.1),
                oma = c(if (x$binary) 3.2 else 2.3, 0, 2.5, 0), family = "sans", las = 1)
  colors <- c(PGS1 = "#a56b28", PGS2 = "#cf9f60", fixed = "#175a83", free = "#5195b6")
  for (group in groups) {
    tab <- x$data
    if (group != "all") tab <- tab[(tab$regression == "RAVEL") == (group == "RAVEL"), , drop = FALSE]
    models <- unique(tab$model)
    lim <- range(c(0, tab$estimate, tab$ci_low, tab$ci_high), finite = TRUE)
    padding <- max(diff(lim) * .08, .01)
    lim <- lim + c(-padding, padding)
    for (role in names(roles)) {
      d <- tab[tab$parameter == role, , drop = FALSE]
      yy <- length(models) + 1L - match(d$model, models)
      graphics::plot(NA_real_, NA_real_, xlim = lim, ylim = c(.4, length(models) + .6),
        xlab = if (isTRUE(x$raw)) "Effect per founder-scale total liability" else
          if (identical(x$scale, "corrected_std")) "Effect per current total-liability SD" else
          if (x$binary) if (group == "PGS") "Probit residual-SD units" else "SEM liability units" else
          "Standardized coefficient", ylab = "", axes = FALSE, main = roles[[role]])
      graphics::axis(1)
      labels <- vapply(models, function(label) paste(strwrap(label, width = label_width), collapse = "\n"), character(1))
      graphics::axis(2, at = seq_along(models), labels = rev(labels), tick = FALSE, cex.axis = .75)
      graphics::abline(v = 0, col = "#969696", lty = 3)
      if (!nrow(d) && group == "PGS") graphics::text(mean(lim), mean(seq_along(models)),
        "Not estimated by\ntrio regression", cex = .8, col = "#777777")
      col <- unname(colors[ifelse(d$regression == "RAVEL", d$residual_model, d$regression)])
      col[is.na(col)] <- "#175a83"
      graphics::segments(d$ci_low, yy, d$ci_high, yy, col = col, lwd = 2)
      graphics::points(d$estimate, yy, pch = ifelse(d$regression == "RAVEL", 16, 17), col = col, cex = 1.15)
      bad <- !is.finite(d$estimate)
      if (any(bad)) graphics::text(mean(lim), yy[bad], "Fit unavailable", cex = .8, col = "#777777")
      no_ci <- is.finite(d$estimate) & (!is.finite(d$ci_low) | !is.finite(d$ci_high))
      if (any(no_ci)) graphics::text(d$estimate[no_ci], yy[no_ci] + .22, "CI unavailable", cex = .65)
      graphics::box(bty = "l", col = "#bdbdbd")
    }
    if (length(roles) %% columns) for (i in seq_len(columns - length(roles) %% columns)) graphics::plot.new()
  }
  graphics::mtext(paste(x$outcome, collapse = ", "), side = 3, outer = TRUE, line = .7, font = 2, cex = 1.1)
  caption <- x$caption %||% paste0(round(100 * x$ci_level), "% intervals: ", x$uncertainty)
  graphics::mtext(caption, side = 1, outer = TRUE, line = .5, cex = .85)
  if (x$binary && length(groups) > 1L) graphics::mtext("Separate outcome scales: do not compare magnitudes between rows directly.",
    side = 1, outer = TRUE, line = 1.8, cex = .8)
  invisible(x)
}

.ravel_report_save_plot <- function(plot, out_dir) {
  label_width <- if (any(plot$data$scale %in% c("corrected_std", "corrected_raw"))) 26L else 29L
  lines <- sum(vapply(unique(plot$data$model), function(x) length(strwrap(x, width = label_width)), integer(1)))
  separate_outcomes <- plot$binary && any(plot$data$regression != "RAVEL") && any(plot$data$regression == "RAVEL")
  height <- max(if (separate_outcomes) 8 else 4.5, 2 + .28 * lines) *
    ceiling(length(unique(plot$data$parameter)) / 3L)
  path <- file.path(out_dir, "diagnostic_estimates.pdf")
  grDevices::pdf(path, width = 14, height = height)
  tryCatch(print(plot), finally = grDevices::dev.off())
  path <- file.path(out_dir, "diagnostic_estimates.png")
  # capabilities("cairo") describes the build, not whether its shared libraries
  # can load. Try native devices and verify the file rather than trusting it.
  types <- unique(c(if (isTRUE(capabilities("aqua"))) "quartz", getOption("bitmapType"),
    if (isTRUE(capabilities("cairo"))) "cairo", if (.Platform$OS.type == "windows") "windows"))
  failures <- character()
  for (type in types) {
    before <- grDevices::dev.cur()
    ok <- tryCatch(withCallingHandlers({
      grDevices::png(path, width = 14, height = height, units = "in", res = 160, type = type)
      print(plot)
      grDevices::dev.off()
      file.exists(path) && isTRUE(file.info(path)$size > 0)
    }, warning = function(w) stop(conditionMessage(w))), error = function(e) {
      failures <<- c(failures, paste0(type, ": ", conditionMessage(e)))
      FALSE
    }, finally = if (grDevices::dev.cur() != before) grDevices::dev.off())
    if (ok) return(invisible(NULL))
    if (file.exists(path)) unlink(path)
  }
  stop("No usable PNG device: ", paste(failures, collapse = "; "))
}

.ravel_report_sem <- function(dat, common_packet, specification, fixed_corr, wanted,
                              vanilla_summary, vanilla_statistics, stats_table, condition) {
  messages <- character()
  diagnostics <- list(converged = NA, admissible = FALSE, warnings = character(), error = "")
  outcome <- common_packet$outcome
  settings <- common_packet$settings
  chosen_corr <- if (specification == "fixed") fixed_corr else NULL
  tryCatch(withCallingHandlers({
    constraint_mode <- settings$fit_constraints %||% "none"
    correction_mode <- if (constraint_mode == "none") "n_gen" else constraint_mode
    fit <- RAVEL::ravel_fit(dat, h = settings$h %||% 0, outcome = outcome$type,
      fixed_effects = settings$fixed_effects, crossperson = "free",
      fixed_crossperson_corr = chosen_corr, n_boot = 0, ci_level = settings$ci_level,
      am_correction = correction_mode,
      am_n_generations = if (constraint_mode == "none") 1L else settings$am_n_generations,
      fit_constraints = constraint_mode, parental_sibs = settings$parental_sibs %||% "none")
    diagnostics$converged <- isTRUE(fit$converged)
    if (!diagnostics$converged) stop(fit$error %||% "RAVEL did not converge")
    V <- lavaan::vcov(fit$fit)
    roles <- .report_packet_roles(common_packet)
    labels <- .report_packet_labels(common_packet)
    if (!is.matrix(V) || !all(labels %in% rownames(V))) stop("Missing joint parameter covariance matrix")
    V_core <- V[labels, labels, drop = FALSE]
    theta <- c(fit$beta_vis[c("bD", "bm", "bf")], fit$latent_cov_vis[c("rhoMF", "eta_m", "eta_f")],
      .sib_theta(fit$parental_sibs, roles))[labels]
    postcheck <- suppressWarnings(isTRUE(lavaan::lavInspect(fit$fit, "post.check")))
    if (!postcheck || any(!is.finite(theta)) || any(!is.finite(V_core)) ||
        any(diag(V_core)[.report_core_labels] <= 0) || any(diag(V_core) < -1e-10) ||
        min(eigen((V_core + t(V_core)) / 2, symmetric = TRUE, only.values = TRUE)$values) < -1e-8)
      stop("RAVEL returned inadmissible estimates or unusable parameter uncertainty")
    if (specification == "fixed") {
      actual <- c(fit$crossperson_resid_corr, fit$parental_sibs$resid_corr)[names(fixed_corr)]
      diagnostics$maximum_fixed_correlation_error <- max(abs(actual - fixed_corr))
      if (any(!is.finite(actual)) || diagnostics$maximum_fixed_correlation_error > 1e-5)
        stop("Fitted residual correlations violate the requested fixed constraints; numerical review required")
    }
    visible <- RAVEL::ravel_to_dataframe(fit, uncertainty = "lavaan", include_vanilla = FALSE)
    visible <- visible[visible$scale == "visible", , drop = FALSE]
    visible$parameter[visible$parameter == "beta_F"] <- "beta_P"
    visible$outcome_variable <- outcome$name
    visible$outcome_scale <- if (outcome$type == "continuous") "prepared_phenotype_SD" else "SEM_probit_liability"
    visible$N <- nrow(dat)
    visible$inference <- "lavaan_delta_wald"
    visible$residual_model <- specification
    visible$fit_status <- "ok"
    pe <- lavaan::parameterEstimates(fit$fit, standardized = TRUE, ci = TRUE, level = settings$ci_level)
    measures <- tryCatch(lavaan::fitMeasures(fit$fit), error = function(e) {
      messages <<- c(messages, paste("Fit measures unavailable:", conditionMessage(e)))
      c(unavailable = NA_real_)
    })
    binary_covariates <- outcome$type == "binary_liability" && length(settings$fixed_effects) > 0L
    unvalidated_measures <- numeric()
    if (binary_covariates) {
      baseline_indices <- grepl("^(baseline\\.|cfi|tli|nnfi|rfi|nfi|pnfi|ifi|rni)", names(measures))
      unvalidated_measures <- measures[baseline_indices]
      measures[baseline_indices] <- NA_real_
      diagnostics$fit_measures_note <- paste(
        "Baseline-relative indices are unavailable: the default independence baseline",
        "is not validated for the conditional-probit moments and can contradict dummy coding.",
        "Adjusted chi-square and RMSEA are moment-matched approximations;",
        "coefficient-SE checks alone do not establish their null calibration.")
    }
    measures <- c(measures, N = nrow(dat), genetic_variance_fraction_visible = fit$R2_model)
    statistics <- stats_table(measures, "RAVEL", specification)
    options <- lavaan::lavInspect(fit$fit, "options")
    diagnostics$admissible <- TRUE
    if (constraint_mode == "estimated_state") diagnostics$estimated_state_identified <- TRUE
    diagnostics$warnings <- unique(messages)
    diagnostics$maximum_core_se <- max(sqrt(diag(V_core)))
    diagnostics$covariate_design_condition <- condition
    diagnostics$binary_covariate_validation <- if (binary_covariates)
      "conditional-probit moments; requires a linear, homoscedastic conditional liability model" else "not_applicable"
    packet <- common_packet
    packet$parental_sibs <- fit$parental_sibs
    packet$theta <- theta
    packet$vcov_core <- V_core
    packet$vcov_full <- V
    packet$parameter_table <- lavaan::parTable(fit$fit)
    packet$parameter_estimates <- pe
    packet$model_implied_covariance <- lavaan::lavInspect(fit$fit, "cov.all")
    packet$reliability <- fit$reliability
    packet$crossperson_resid <- fit$crossperson_resid
    packet$crossperson_resid_corr <- fit$crossperson_resid_corr
    packet$visible_table <- visible
    packet$summary <- rbind(visible[wanted], vanilla_summary)
    packet$fit_measures <- measures
    packet$unvalidated_baseline_fit_measures <- unvalidated_measures
    packet$fit_statistics <- rbind(statistics, vanilla_statistics)
    packet$diagnostics <- diagnostics
    packet$outcome$model_implied_variance <- lavaan::lavInspect(fit$fit, "cov.ov")["Y", "Y"]
    packet$settings$residual_model <- specification
    packet$settings$constraint_spec <- fit$constraint_spec
    packet$settings$binary_covariate_moments <- if (binary_covariates)
      "conditional_probit_joint_marginal_liability" else "not_applicable"
    packet$settings$fixed_crossperson_corr <- chosen_corr
    packet$settings$model_syntax <- fit$model
    packet$settings$lavaan_options <- options[c("estimator", "estimator.orig", "se", "parameterization",
                                               "fixed.x", "conditional.x", "std.lv")]
    residuals <- data.frame(outcome_variable = outcome$name, residual_model = specification,
      parameter = names(fixed_corr),
      covariance = unname(c(fit$crossperson_resid, fit$parental_sibs$resid_cov)[names(fixed_corr)]),
      correlation = unname(c(fit$crossperson_resid_corr, fit$parental_sibs$resid_corr)[names(fixed_corr)]))
    list(summary = visible[wanted], fit_statistics = statistics, residuals = residuals,
         diagnostics = diagnostics, packet = packet)
  }, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    diagnostics$admissible <- FALSE
    diagnostics$error <- conditionMessage(e)
    diagnostics$warnings <- unique(messages)
    tab <- data.frame(outcome_variable = outcome$name, regression = "RAVEL", scale = "visible",
      outcome_scale = if (outcome$type == "continuous") "prepared_phenotype_SD" else "SEM_probit_liability",
      parameter = .report_parameters(.report_packet_roles(common_packet), "visible"),
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
      N = nrow(dat), inference = "unavailable", residual_model = specification, fit_status = "failed")
    list(summary = tab[wanted], fit_statistics = stats_table(c(N = nrow(dat)), "RAVEL", specification),
      residuals = data.frame(outcome_variable = outcome$name, residual_model = specification,
        parameter = names(fixed_corr), covariance = NA_real_, correlation = NA_real_),
      diagnostics = diagnostics, packet = NULL)
  })
}
