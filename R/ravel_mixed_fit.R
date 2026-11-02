# Opt-in family likelihood. The established ravel_fit() path remains separate.

.ravel_mixed_prepare <- function(dat, outcome, col_outcome, pgs_columns,
    family_id, subject_id, cluster, fixed_effects, interactions, interaction_center,
    family_random, individual_random, random_slope, fit_constraints,
    am_n_generations, pi, covariate_mapping_confirmed, paternal_id = NULL) {
  if (!is.data.frame(dat)) stop("dat must be a data frame.")
  for (nm in c("family_id", "subject_id", "cluster")) {
    value <- get(nm)
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value))
      stop(nm, " must name one identifier column; use a distinct FID per family when extended families are independent.")
  }
  if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi) || pi < 0 || pi >= 1)
    stop("pi must be one number in [0, 1).")
  if (!fit_constraints %in% c("estimated_state", "n_gen", "equilibrium"))
    stop("The mixed fitter currently supports estimated_state, n_gen and equilibrium; the original fitter retains its other modes.")
  if (fit_constraints == "n_gen" && (length(am_n_generations) != 1L ||
      !is.finite(am_n_generations) || am_n_generations < 1 || am_n_generations != as.integer(am_n_generations)))
    stop("am_n_generations must be a positive integer.")
  for (flag in list(family_random, individual_random))
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) stop("Random-intercept options must be TRUE or FALSE.")
  if (!is.null(random_slope) && (!is.character(random_slope) || length(random_slope) != 1L))
    stop("random_slope must name one numeric time variable.")
  covariates <- unique(c(fixed_effects, interactions, random_slope))
  if (length(covariates) && !isTRUE(covariate_mapping_confirmed))
    stop("Set covariate_mapping_confirmed=TRUE only if the conditional genetic distribution is assumed independent of the supplied covariates/moderators. This assumption is not established by a small sample correlation.")
  if (anyDuplicated(interactions) || anyDuplicated(unname(pgs_columns)))
    stop("Interactions and PGS column mappings must be distinct.")
  needed <- unique(c(col_outcome, unname(pgs_columns), family_id, subject_id, cluster, covariates, paternal_id))
  absent <- setdiff(needed, names(dat))
  if (length(absent)) stop("Columns not found: ", paste(absent, collapse = ", "))
  if (any(vapply(dat[c(family_id, subject_id, cluster)], anyNA, logical(1))))
    stop("Individual, maternal-family and FID cluster identifiers must be observed.")
  if (any(vapply(dat[c(family_id, subject_id, cluster)], function(x) any(!nzchar(as.character(x))), logical(1))))
    stop("Identifiers must not be empty strings.")
  fam <- as.character(dat[[family_id]]); iid <- as.character(dat[[subject_id]])
  fid <- as.character(dat[[cluster]])
  mapping <- unique(data.frame(iid = iid, family = fam, fid = fid))
  if (anyDuplicated(mapping$iid)) stop("Each individual must belong to one maternal family and one FID.")
  family_mapping <- unique(mapping[c("family", "fid")])
  if (anyDuplicated(family_mapping$family)) stop("Each maternal family must lie within one FID; the two ID columns need not have identical values.")
  if (!is.null(paternal_id)) {
    fathers <- as.character(dat[[paternal_id]])
    if (anyNA(fathers) || any(!nzchar(fathers))) stop("Supplied paternal identifiers must be observed and nonempty.")
    paternal_mapping <- unique(data.frame(father = fathers, fid = fid))
    if (anyDuplicated(paternal_mapping$father))
      stop("Each father must lie within one FID, including when he has children in multiple maternal families.")
  }
  if (length(unique(fid)) < 2L) stop("At least two FID clusters are required.")
  numeric_columns <- unique(c(col_outcome, unname(pgs_columns), covariates))
  if (any(!vapply(dat[numeric_columns], is.numeric, logical(1)))) stop("Outcomes, PGSs and covariates must be numeric.")
  if (any(vapply(dat[numeric_columns], function(x) any(!is.na(x) & !is.finite(x)), logical(1))))
    stop("Numerical inputs must be finite or NA.")
  if (anyNA(dat[unname(pgs_columns)])) stop("This mixed likelihood currently requires complete PGS measurements; it does not impute missing scores.")
  binary <- outcome == "binary_liability"
  if (binary && !setequal(unique(stats::na.omit(dat[[col_outcome]])), c(0, 1)))
    stop("A binary outcome must contain both numeric 0 and 1.")
  counts <- table(iid)
  if (binary && any(counts > 1L)) stop("Binary longitudinal models are outside the supported scope.")
  if ((individual_random || !is.null(random_slope)) && !any(counts > 1L))
    stop("Individual random effects require repeated observations; a cross-sectional individual intercept is not separately identifiable.")
  children_per_family <- table(mapping$family)
  if (family_random && individual_random && !any(children_per_family > 1L))
    stop("Family and individual intercepts are indistinguishable when every family contains one individual.")
  if (family_random && !any(children_per_family > 1L) && !any(counts > 1L))
    stop("A cross-sectional family variance needs some families with multiple observed siblings.")
  # Shared genotypes are represented once, not as fresh measurements per visit.
  constant <- function(x) max(x) - min(x) <= 1e-10 * max(1, max(abs(x)))
  for (id in unique(iid)) for (nm in unname(pgs_columns))
    if (!constant(dat[iid == id, nm])) stop("PGSs must remain constant across visits for each individual: ", id)
  parent_columns <- unname(pgs_columns[grepl("mother|father", names(pgs_columns))])
  for (id in unique(fam)) {
    rows <- fam == id
    if (!is.null(paternal_id) && (anyNA(dat[[paternal_id]][rows]) ||
        length(unique(dat[[paternal_id]][rows])) != 1L))
      stop("The current genetic block requires full siblings; maternal family ", id,
           " contains missing or multiple paternal IDs. Apply an explicit full-sibling restriction before fitting.")
    if (any(!vapply(dat[rows, parent_columns, drop = FALSE], constant, logical(1))))
      stop("Parental PGSs differ within maternal family ", id,
           "; the current genetic block assumes full siblings with the same two parents.")
  }
  if (!is.null(paternal_id)) for (id in unique(fathers)) {
    columns <- unname(pgs_columns[grepl("father", names(pgs_columns))])
    if (any(!vapply(dat[fathers == id, columns, drop = FALSE], constant, logical(1))))
      stop("Paternal PGSs disagree for the same father across maternal families: ", id)
  }
  observed <- !is.na(dat[[col_outcome]])
  if (length(covariates)) observed <- observed & stats::complete.cases(dat[covariates])
  if (sum(observed) < 50L) stop("At least 50 observed outcomes with complete covariates are required.")
  observed_subjects <- unique(mapping[mapping$iid %in% iid[observed], c("iid", "family")])
  observed_counts <- table(iid[observed])
  observed_family_counts <- table(observed_subjects$family)
  if ((individual_random || !is.null(random_slope)) && !any(observed_counts > 1L))
    stop("Individual random effects require repeated observed outcomes after missing-data exclusions.")
  if (family_random && individual_random && !any(observed_family_counts > 1L))
    stop("The observed outcomes do not separate family and individual intercept variances.")
  if (family_random && !any(observed_family_counts > 1L) && !any(observed_counts > 1L))
    stop("The observed outcomes do not identify a family intercept variance.")
  if (!is.null(random_slope) && !any(vapply(split(dat[[random_slope]][observed], iid[observed]),
      function(x) length(unique(x)) > 1L, logical(1))))
    stop("An individual age slope requires within-individual time variation in observed outcomes.")
  centers <- scales <- setNames(numeric(length(covariates)), covariates)
  for (nm in covariates) {
    x <- dat[[nm]][observed]
    scales[nm] <- stats::sd(x)
    if (!is.finite(scales[nm]) || scales[nm] <= 0) stop("Covariate has no observed variation: ", nm)
    centers[nm] <- if (all(x %in% c(0, 1))) 0 else mean(x)
  }
  if (!is.null(interaction_center)) {
    if (!is.numeric(interaction_center) || is.null(names(interaction_center)) ||
        anyDuplicated(names(interaction_center)) || any(!is.finite(interaction_center)) ||
        !all(names(interaction_center) %in% covariates))
      stop("interaction_center must be a finite named numeric vector of supplied covariates.")
    centers[names(interaction_center)] <- interaction_center
  }
  adjusted <- dat
  for (nm in covariates) adjusted[[nm]] <- (dat[[nm]] - centers[nm]) / scales[nm]
  outcome_center <- if (binary) 0 else mean(dat[[col_outcome]][observed])
  outcome_scale <- if (binary) 1 else stats::sd(dat[[col_outcome]][observed])
  if (!is.finite(outcome_scale) || outcome_scale <= 0) stop("The outcome has no observed variation.")
  adjusted[[col_outcome]] <- (dat[[col_outcome]] - outcome_center) / outcome_scale
  score_scaling <- data.frame(column = unname(pgs_columns), center = 0, scale = 1)
  for (j in seq_along(pgs_columns)) {
    nm <- pgs_columns[j]
    unique_rows <- if (grepl("child", names(pgs_columns)[j])) !duplicated(iid) else !duplicated(fam)
    x <- dat[[nm]][unique_rows]
    score_scaling$center[j] <- mean(x); score_scaling$scale[j] <- stats::sd(x)
    if (!is.finite(score_scaling$scale[j]) || score_scaling$scale[j] <= 0)
      stop("PGS column has no variation: ", nm)
    adjusted[[nm]] <- (dat[[nm]] - score_scaling$center[j]) / score_scaling$scale[j]
  }
  beta_names <- unlist(lapply(c("main", interactions), function(nm)
    paste0("genetic_", c("child", "mother", "father"), "_", nm)))
  spec <- list(pi = pi, mode = fit_constraints, generations = am_n_generations,
    family_random = family_random, individual_random = individual_random,
    random_slope = random_slope, interactions = interactions, covariates = covariates,
    centers = centers, scales = scales, beta_names = beta_names,
    mean_names = c(paste0("score_mean", rep(1:2, each = 3), "_", rep(c("mother", "father", "child"), 2)),
                   "mean_intercept", if (length(covariates)) paste0("fixed_", covariates)),
    outcome_scale = outcome_scale, outcome_center = outcome_center, outcome = outcome)
  groups <- list(); keys <- new.env(parent = emptyenv())
  for (id in sort(unique(fam))) {
    rows <- which(fam == id); subjects <- sort(unique(iid[rows]))
    first <- rows[match(subjects, iid[rows])]
    yr <- rows[observed[rows]]
    yr <- yr[order(match(iid[yr], subjects), if (!is.null(random_slope)) adjusted[[random_slope]][yr] else yr)]
    child <- match(iid[yr], subjects); nchild <- length(subjects)
    p <- unlist(lapply(1:2, function(k) c(adjusted[first[1], pgs_columns[[paste0("mother", k)]]],
      adjusted[first[1], pgs_columns[[paste0("father", k)]]],
      adjusted[first, pgs_columns[[paste0("child", k)]]])), use.names = FALSE)
    X <- cbind(rep(1, length(yr)), if (length(covariates)) as.matrix(adjusted[yr, covariates, drop = FALSE]) else matrix(numeric(), length(yr), 0))
    moderators <- cbind(rep(1, length(yr)), if (length(interactions)) as.matrix(adjusted[yr, interactions, drop = FALSE]) else matrix(numeric(), length(yr), 0))
    slope <- if (is.null(random_slope)) numeric(length(yr)) else adjusted[[random_slope]][yr]
    score_index <- c(1, 2, rep(3, nchild), 4, 5, rep(6, nchild))
    D <- matrix(0, length(p) + length(yr), length(spec$mean_names))
    D[cbind(seq_along(p), score_index)] <- 1
    if (length(yr)) D[length(p) + seq_along(yr), 6 + seq_len(ncol(X))] <- X
    key <- paste(nchild, paste(child, collapse = ","),
      paste(format(c(D, moderators, slope), digits = 17, scientific = TRUE, trim = TRUE), collapse = ","), sep = "|")
    if (!exists(key, keys, inherits = FALSE)) {
      at <- length(groups) + 1L; assign(key, at, keys)
      groups[[at]] <- list(nchild = as.integer(nchild), child = child, moderators = moderators,
        slope = slope, D = D, n = 0L, observed = list(), families = character(), clusters = character())
    } else at <- get(key, keys)
    groups[[at]]$n <- groups[[at]]$n + 1L
    groups[[at]]$observed[[groups[[at]]$n]] <- c(p, adjusted[[col_outcome]][yr])
    groups[[at]]$families <- c(groups[[at]]$families, id)
    groups[[at]]$clusters <- c(groups[[at]]$clusters, unique(fid[rows]))
  }
  groups <- lapply(groups, function(z) {
    z$observed <- do.call(rbind, z$observed)
    z$sum <- colSums(z$observed); z$cross <- crossprod(z$observed)
    z
  })
  start <- c(state_rho = atanh(.25), if (fit_constraints == "estimated_state") c(state_variance = log(.7)))
  for (k in 1:2) for (role in c("mother", "father", "child"))
    start[paste0("loading", k, "_", role)] <- log(.65)
  start[c("noise1", "noise2")] <- log(1)
  start[beta_names] <- 0
  start[beta_names[1:3]] <- c(.25, .1, .1)
  start["residual_sd"] <- log(.7)
  if (family_random) start["family_sd"] <- log(.25)
  if (individual_random) start["individual_sd"] <- log(.3)
  if (!is.null(random_slope)) {
    start["slope_sd"] <- log(.1)
    if (individual_random) start["intercept_slope_cor"] <- 0
  }
  lower <- upper <- start
  lower[] <- -Inf; upper[] <- Inf
  logs <- grepl("loading|noise|_sd$|state_variance", names(start))
  lower[logs] <- -12; upper[logs] <- 8
  lower["state_rho"] <- atanh(-.98); upper["state_rho"] <- atanh(.98)
  if ("intercept_slope_cor" %in% names(start)) {
    lower["intercept_slope_cor"] <- atanh(-.995)
    upper["intercept_slope_cor"] <- atanh(.995)
  }
  list(groups = groups, spec = spec, start = start, lower = lower, upper = upper,
    data = dat, adjusted_data = adjusted, observed = observed, score_scaling = score_scaling,
    N = sum(observed), n_subjects = nrow(mapping), n_families = nrow(family_mapping),
    n_clusters = length(unique(fid)), mapping = mapping, col_outcome = col_outcome,
    pgs_columns = pgs_columns, family_id = family_id, subject_id = subject_id, cluster = cluster,
    paternal_id = paternal_id)
}

.ravel_mixed_effects <- function(par, spec, scale) {
  state <- .ravel_mixed_state(par, spec)
  roles <- c("child", "mother", "father")
  b <- matrix(par[spec$beta_names], nrow = 3L)
  multiplier <- c(1, 1 / spec$scales[spec$interactions]) * spec$outcome_scale
  b <- sweep(b, 2, multiplier, "*")
  if (scale == "corrected_raw") b <- b / sqrt(c(state$child, state$parent, state$parent))
  if (scale == "visible") {
    genetic <- .ravel_mixed_genetic(1L, par, spec)
    order <- c(3L, 1L, 2L)
    Traw <- genetic$T * outer(sqrt(c(state$parent, state$parent, state$child)),
                              sqrt(c(state$parent, state$parent, state$child)))
    K <- matrix(c(1, 0, .5, 0, 1, .5, .5, .5, 1), 3L)
    V <- (1 - spec$pi) * Traw + spec$pi * K
    GT <- sqrt(1 - spec$pi) * sweep(Traw, 2, sqrt(diag(Traw)), "/")
    GT <- sweep(GT, 1, sqrt(diag(V)), "/")
    G <- stats::cov2cor(V)
    b <- solve(G[order, order], GT[order, order]) %*% b
  }
  labels <- unlist(lapply(c("", spec$interactions), function(nm)
    paste0(c("beta_D", "beta_M", "beta_F"), if (nzchar(nm)) paste0(":", nm) else "")))
  stats::setNames(as.numeric(b), labels)
}

.ravel_mixed_variances <- function(par, spec) {
  ans <- c(residual = exp(2 * par[["residual_sd"]]))
  if (spec$family_random) ans["family_intercept"] <- exp(2 * par[["family_sd"]])
  if (spec$individual_random) ans["individual_intercept"] <- exp(2 * par[["individual_sd"]])
  if (!is.null(spec$random_slope)) {
    ans["individual_slope"] <- exp(2 * par[["slope_sd"]]) / spec$scales[[spec$random_slope]]^2
    if (spec$individual_random) ans["individual_intercept_slope_covariance"] <-
      exp(par[["individual_sd"]] + par[["slope_sd"]]) * tanh(par[["intercept_slope_cor"]]) /
        spec$scales[[spec$random_slope]]
  }
  ans * spec$outcome_scale^2
}

#' Fit an optional RAVEL family random-effects model
#'
#' A separate full-family likelihood for shared parents and repeated children.
#' Maternal families define the random intercept; FID defines robust inference.
#' Genetic measurement error is integrated jointly with the outcome model.
#' @export
ravel_mixed_fit <- function(dat, outcome = c("continuous", "binary_liability"),
    col_outcome = "Y", col_pgs1_child = "PGS1_child", col_pgs2_child = "PGS2_child",
    col_pgs1_mother = "PGS1_mother", col_pgs2_mother = "PGS2_mother",
    col_pgs1_father = "PGS1_father", col_pgs2_father = "PGS2_father",
    family_id = "M_ID", subject_id = "IID", cluster = "FID", paternal_id = NULL,
    fixed_effects = NULL, interactions = NULL, interaction_center = NULL,
    family_random = TRUE, individual_random = FALSE, random_slope = NULL,
    fit_constraints = c("estimated_state", "equilibrium", "n_gen"), am_n_generations = 1L,
    pi = 0, pi_se = NULL, covariate_mapping_confirmed = FALSE,
    diagnostic = TRUE, diagnostic_control = list(), control = list(), start = NULL,
    se = TRUE, keep_data = TRUE) {
  outcome <- match.arg(outcome); fit_constraints <- match.arg(fit_constraints)
  for (flag in list(diagnostic, se, keep_data))
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag))
      stop("diagnostic, se and keep_data must each be TRUE or FALSE.")
  if (!is.null(pi_se) && (!is.numeric(pi_se) || length(pi_se) != 1L ||
      !is.finite(pi_se) || pi_se < 0)) stop("pi_se must be a finite nonnegative number.")
  pgs_columns <- c(child1 = col_pgs1_child, child2 = col_pgs2_child,
    mother1 = col_pgs1_mother, mother2 = col_pgs2_mother,
    father1 = col_pgs1_father, father2 = col_pgs2_father)
  prepared <- .ravel_mixed_prepare(dat, outcome, col_outcome, pgs_columns,
    family_id, subject_id, cluster, fixed_effects, interactions, interaction_center,
    family_random, individual_random, random_slope, fit_constraints,
    am_n_generations, pi, covariate_mapping_confirmed, paternal_id)
  prepared$spec$pi_se <- pi_se
  if (outcome == "binary_liability") {
    if (!exists(".ravel_mixed_fit_binary", mode = "function"))
      stop("The binary mixed likelihood is not implemented in this prototype.")
    return(.ravel_mixed_fit_binary(prepared, pi_se = pi_se, diagnostic = diagnostic,
      diagnostic_control = diagnostic_control,
      control = control, start = start, se = se, keep_data = keep_data))
  }
  prepared$groups <- .ravel_mixed_compress_means(prepared$groups)
  spec <- prepared$spec
  par <- prepared$start
  if (!is.null(start)) {
    if (is.null(names(start)) || !all(names(start) %in% names(par)) || any(!is.finite(start)))
      stop("start must be a finite named vector of model parameters.")
    par[names(start)] <- start
  }
  evaluate <- .ravel_mixed_evaluator(prepared$groups, spec)
  initial <- evaluate(par)
  if (!isTRUE(initial$valid)) stop("Initial family likelihood is invalid: ", initial$reason)
  settings <- utils::modifyList(list(maxit = 500L, factr = 1e2, pgtol = 1e-6), control)
  fit <- stats::optim(par, function(x) evaluate(x)$value,
    function(x) evaluate(x, derivatives = TRUE)$gradient,
    method = "L-BFGS-B", lower = prepared$lower, upper = prepared$upper, control = settings)
  par <- fit$par
  final <- evaluate(par, derivatives = TRUE)
  if (!isTRUE(final$valid)) stop("Optimizer returned an invalid family likelihood: ", final$reason)
  certification <- .ravel_mixed_certification(par, final, prepared, spec)
  projected_gradient <- certification$projected_gradient
  gradient_ok <- certification$gradient_ok
  converged <- fit$convergence == 0L && isTRUE(final$valid) && is.finite(final$value) && gradient_ok
  messages <- character()
  if (!gradient_ok) messages <- c(messages, "The likelihood gradient, including its natural-variance scale, remains too large for a certified optimum.")
  boundary <- certification$boundary_parameters
  if (length(boundary)) messages <- c(messages, paste("Parameters near a boundary:",
    paste(boundary, collapse = ", "), "; ordinary Wald uncertainty is not certified."))
  uncertainty <- NULL
  if (converged && !length(boundary) && isTRUE(se)) {
    uncertainty <- tryCatch(.ravel_mixed_inference(par, prepared$groups, spec, evaluate),
      error = function(e) { messages <<- c(messages, conditionMessage(e)); NULL })
  }
  V <- if (is.null(uncertainty)) matrix(NA_real_, length(par), length(par),
    dimnames = list(names(par), names(par))) else uncertainty$vcov[names(par), names(par), drop = FALSE]
  maps <- lapply(c("visible", "corrected_std", "corrected_raw"), function(scale)
    .ravel_mixed_delta(function(x) .ravel_mixed_effects(x, spec, scale), par))
  names(maps) <- c("visible", "corrected_std", "corrected_raw")
  pi_derivatives <- NULL
  variance_pi_derivative <- NULL
  # Unlike the original post-fit correction, this joint likelihood depends on
  # pi through the outcome covariance in every assortment mode.
  if (converged && !length(boundary) && !is.null(pi_se) && pi_se > 0 && isTRUE(se)) {
    step <- min(1e-3, max(1e-5, pi_se / 10), (1 - pi) / 4)
    values <- if (pi > step) c(pi - step, pi + step) else c(pi, pi + step)
    refits <- lapply(values, function(value) {
      s <- spec; s$pi <- value
      ev <- .ravel_mixed_evaluator(prepared$groups, s)
      z <- stats::optim(par, function(x) ev(x)$value,
        function(x) ev(x, derivatives = TRUE)$gradient, method = "L-BFGS-B",
        lower = prepared$lower, upper = prepared$upper, control = settings)
      ez <- ev(z$par, derivatives = TRUE)
      check <- if (isTRUE(ez$valid)) .ravel_mixed_certification(z$par, ez, prepared, s) else NULL
      list(par = z$par, spec = s, valid = z$convergence == 0L && isTRUE(ez$valid) &&
        isTRUE(check$gradient_ok) && !length(check$boundary_parameters))
    })
    if (all(vapply(refits, `[[`, logical(1), "valid"))) {
      pi_derivatives <- lapply(names(maps), function(scale)
        (.ravel_mixed_effects(refits[[2]]$par, refits[[2]]$spec, scale) -
         .ravel_mixed_effects(refits[[1]]$par, refits[[1]]$spec, scale)) / diff(values))
      names(pi_derivatives) <- names(maps)
      variance_pi_derivative <-
        (.ravel_mixed_variances(refits[[2]]$par, refits[[2]]$spec) -
         .ravel_mixed_variances(refits[[1]]$par, refits[[1]]$spec)) / diff(values)
    } else messages <- c(messages, "A pi-sensitivity refit failed; externally corrected uncertainty is unavailable.")
  }
  coefficient_vcov <- list()
  tables <- lapply(names(maps), function(scale) {
    z <- maps[[scale]]; covariance <- z$jacobian %*% V %*% t(z$jacobian)
    if (!is.null(pi_se) && pi_se > 0) {
      if (is.null(pi_derivatives)) covariance[,] <- NA_real_ else
        covariance <- covariance + pi_se^2 * tcrossprod(pi_derivatives[[scale]])
    }
    error <- sqrt(pmax(diag(covariance), 0))
    coefficient_vcov[[scale]] <<- covariance
    data.frame(parameter = names(z$value), scale = scale, estimate = unname(z$value), se = error,
      ci_low = z$value - stats::qnorm(.975) * error, ci_high = z$value + stats::qnorm(.975) * error,
      p_value = 2 * stats::pnorm(-abs(z$value / error)),
      uncertainty_source = if (is.null(uncertainty)) "unavailable" else uncertainty$method,
      outcome_scale = "input outcome units", row.names = NULL)
  })
  variances <- .ravel_mixed_delta(function(x) .ravel_mixed_variances(x, spec), par)
  variance_vcov <- variances$jacobian %*% V %*% t(variances$jacobian)
  if (!is.null(pi_se) && pi_se > 0) {
    if (is.null(variance_pi_derivative)) variance_vcov[,] <- NA_real_ else
      variance_vcov <- variance_vcov + pi_se^2 * tcrossprod(variance_pi_derivative)
  }
  variance_table <- data.frame(component = names(variances$value), estimate = unname(variances$value),
    se = sqrt(pmax(diag(variance_vcov), 0)),
    uncertainty_conditions = if (!is.null(pi_se) && pi_se > 0)
      "includes independent external pi uncertainty when sensitivity refits succeed" else "conditional on supplied pi", row.names = NULL)
  result <- structure(list(converged = converged, optimizer_converged = fit$convergence == 0L,
    inference_status = if (!se) "not_requested" else if (length(boundary)) "boundary" else
      if (!converged) "not_converged" else if (is.null(uncertainty)) "unavailable" else
      if (!is.null(pi_se) && pi_se > 0 && is.null(pi_derivatives)) "pi_sensitivity_unavailable" else "ok",
    table = do.call(rbind, tables), coefficient_vcov = coefficient_vcov,
    variance_components = variance_table, variance_vcov = variance_vcov,
    parameters = par, means = stats::setNames(final$means, spec$mean_names),
    vcov = if (is.null(uncertainty)) NULL else uncertainty$vcov,
    model_vcov = if (is.null(uncertainty)) NULL else uncertainty$model_vcov,
    uncertainty = uncertainty, pi_sensitivity = pi_derivatives,
    variance_pi_sensitivity = variance_pi_derivative,
    N = prepared$N, n_subjects = prepared$n_subjects, n_families = prepared$n_families,
    n_clusters = prepared$n_clusters, fit_statistics = c(logLik = -final$value,
      n_parameters = length(par) + length(final$means),
      AIC = 2 * final$value + 2 * (length(par) + length(final$means))),
    state = .ravel_mixed_state(par, spec), settings = spec, optimizer = fit,
    gradient = final$gradient, projected_gradient = projected_gradient,
    natural_variance_gradient = certification$natural_variance_gradient,
    projected_natural_variance_gradient = certification$projected_natural_variance_gradient,
    boundary_parameters = boundary, messages = messages, score_scaling = prepared$score_scaling,
    family_id = family_id, subject_id = subject_id, cluster = cluster,
    diagnostics = list(full_sibling_assumption = TRUE,
      paternal_id_checked = !is.null(paternal_id), design_groups = length(prepared$groups),
      score_measurements_counted_once_per_person_within_maternal_family = TRUE,
      likelihood = "Gaussian maternal-family block likelihood; FID-clustered inference",
      visible_scale = "Marginal trio projection, not the conditional regression on every family member's scores"),
    prepared = if (keep_data) prepared else NULL), class = "ravel_mixed_results")
  if (diagnostic) result$one_factor <- .ravel_mixed_run_diagnostic(prepared, diagnostic_control)
  result
}

#' @export
as.data.frame.ravel_mixed_results <- function(x, ...) x$table

#' @export
print.ravel_mixed_results <- function(x, ...) {
  cat("RAVEL optional mixed family model\n")
  cat(x$N, "outcomes;", x$n_subjects, "individuals;", x$n_families,
      "maternal families;", x$n_clusters, "FID clusters\n")
  cat("Optimizer converged:", x$optimizer_converged, "; inference:", x$inference_status, "\n")
  print(x$table[x$table$scale == "corrected_std", ], row.names = FALSE)
  if (length(x$messages)) cat(paste(x$messages, collapse = "\n"), "\n")
  invisible(x)
}
