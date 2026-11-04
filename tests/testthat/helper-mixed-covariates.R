# Independent correlated-covariate fixtures. Genetics are drawn by the existing
# primitive validation generator; generation and population oracles use no
# fitter or correction helper. The final audit function checks fitted internals.

mixed_covariate_contract <- function(n_pc = 2L, pc_loading = NULL,
    pc_between = .2, pc_within = .6, pc_correlation = .25,
    sex_pc_association = 0, baseline_pc_association = .4,
    baseline_between = .3, baseline_within = .7, age_center = 10) {
  stopifnot(n_pc >= 1L, n_pc == as.integer(n_pc), abs(pc_correlation) < 1,
    pc_between >= 0, pc_within > 0, baseline_between >= 0, baseline_within >= 0)
  if (is.null(pc_loading)) {
    j <- seq_len(n_pc)
    pc_loading <- cbind(child = .20 * cos(j - 1), mother = .16 * sin(j),
                        father = -.12 * cos(j))
    pc_loading <- pc_loading / sqrt(max(1, n_pc / 2))
  }
  pc_loading <- as.matrix(pc_loading)
  stopifnot(identical(dim(pc_loading), c(as.integer(n_pc), 3L)),
            all(is.finite(pc_loading)))
  colnames(pc_loading) <- c("child", "mother", "father")
  rownames(pc_loading) <- paste0("PC", seq_len(n_pc))
  correlation <- pc_correlation^abs(outer(seq_len(n_pc), seq_len(n_pc), "-"))
  sex_loading <- rep(0, n_pc); sex_loading[1L] <- sex_pc_association
  baseline_loading <- rep(0, n_pc); baseline_loading[1L] <- baseline_pc_association
  list(n_pc = as.integer(n_pc), pc_names = rownames(pc_loading),
    pc_loading = pc_loading, pc_between = pc_between * correlation,
    pc_within = pc_within * correlation, sex_loading = sex_loading,
    baseline_loading = baseline_loading, baseline_between = baseline_between,
    baseline_within = baseline_within, age_center = age_center)
}

mixed_covariate_family <- function(population = mixed_validation_population(),
    n_children = 2L, contract = mixed_covariate_contract(),
    covariates = c(contract$pc_names, "sex", "age"), c_values = NULL) {
  stopifnot(n_children >= 1L, n_children == as.integer(n_children))
  n <- as.integer(n_children); r <- n + 2L; p <- contract$n_pc
  # Primitive innovations rather than the fitted compact T/K covariance formula.
  primitive <- matrix(0, 2L * r, 2L * r)
  primitive[1:4, 1:4] <- population$parental
  diag(primitive)[-(1:4)] <- .5
  inheritance <- diag(2L * r)
  for (i in seq_len(n)) {
    rows <- 2L * (i + 1L) + 1:2
    inheritance[rows, 1:2] <- inheritance[rows, 3:4] <- diag(.5, 2L)
  }
  visible <- inheritance[seq.int(1L, 2L * r, by = 2L), , drop = FALSE]
  total <- kronecker(diag(r), matrix(population$w, nrow = 1L)) %*% inheritance
  VV <- visible %*% primitive %*% t(visible)
  TT <- total %*% primitive %*% t(total)
  TV <- total %*% primitive %*% t(visible)
  pc_map <- matrix(0, n * p, r)
  for (i in seq_len(n)) {
    rows <- (i - 1L) * p + seq_len(p)
    pc_map[rows, 1L] <- contract$pc_loading[, "mother"]
    pc_map[rows, 2L] <- contract$pc_loading[, "father"]
    pc_map[rows, i + 2L] <- contract$pc_loading[, "child"]
  }
  Qpc <- pc_map %*% VV %*% t(pc_map) +
    kronecker(matrix(1, n, n), contract$pc_between) +
    kronecker(diag(n), contract$pc_within)
  Vpc <- VV %*% t(pc_map); Tpc <- TV %*% t(pc_map)
  all_names <- c(contract$pc_names, "sex", "age")
  q <- length(all_names)
  gaussian_map <- matrix(0, n * q, n * p)
  added <- matrix(0, n * q, n * q)
  for (i in seq_len(n)) {
    rows <- (i - 1L) * q + seq_len(q)
    pcs <- (i - 1L) * p + seq_len(p)
    gaussian_map[rows[seq_len(p)], pcs] <- diag(p)
    gaussian_map[rows[p + 1L], pcs] <- contract$sex_loading
    gaussian_map[rows[p + 2L], pcs] <- contract$baseline_loading
    added[rows[p + 1L], rows[p + 1L]] <- 1
    for (j in seq_len(n)) {
      column <- (j - 1L) * q + p + 2L
      added[rows[p + 2L], column] <- contract$baseline_between +
        if (i == j) contract$baseline_within else 0
    }
  }
  Qgaussian <- gaussian_map %*% Qpc %*% t(gaussian_map) + added
  Dgaussian <- Vpc %*% t(gaussian_map)
  Egaussian <- Tpc %*% t(gaussian_map)
  sex_rows <- (seq_len(n) - 1L) * q + p + 1L
  other_rows <- setdiff(seq_len(n * q), sex_rows)
  Q <- Qgaussian; D <- Dgaussian; E <- Egaussian
  # Sex is 1(Z>0) with Z Gaussian conditional on PCs. These are exact binary
  # moments, and V|PC,sex remains Gaussian because sex adds no information
  # after all PCs have been observed.
  for (j in sex_rows) {
    divisor <- sqrt(2 * pi * Qgaussian[j, j])
    D[, j] <- Dgaussian[, j] / divisor
    E[, j] <- Egaussian[, j] / divisor
    Q[other_rows, j] <- Qgaussian[other_rows, j] / divisor
    Q[j, other_rows] <- Q[other_rows, j]
  }
  for (i in sex_rows) for (j in sex_rows) {
    correlation <- Qgaussian[i, j] / sqrt(Qgaussian[i, i] * Qgaussian[j, j])
    Q[i, j] <- asin(max(-1, min(1, correlation))) / (2 * pi)
  }
  keep <- unlist(lapply(seq_len(n), function(i)
    (i - 1L) * q + match(covariates, all_names)))
  stopifnot(!anyNA(keep), !anyDuplicated(covariates))
  Q <- Q[keep, keep, drop = FALSE]
  D <- D[, keep, drop = FALSE]; E <- E[, keep, drop = FALSE]
  labels <- unlist(lapply(seq_len(n), function(i) paste0("child", i, ":", covariates)))
  dimnames(Q) <- list(labels, labels); colnames(D) <- colnames(E) <- labels
  projection <- TV %*% solve(VV)
  hidden <- TT - projection %*% VV %*% t(projection)
  Qinv <- solve(Q)
  conditional_V <- VV - D %*% Qinv %*% t(D)
  conditional_T <- TT - E %*% Qinv %*% t(E)
  conditional_TV <- TV - E %*% Qinv %*% t(D)
  mean_c <- c(stats::setNames(rep(0, p), contract$pc_names), sex = .5,
              age = contract$age_center)[covariates]
  mean_V <- mean_T <- rep(0, r)
  if (!is.null(c_values)) {
    if (is.matrix(c_values)) {
      stopifnot(identical(dim(c_values), c(n, length(covariates))))
      c_values <- as.vector(t(c_values))
    }
    centered <- as.numeric(c_values) - rep(mean_c, n)
    stopifnot(length(centered) == ncol(D))
    mean_V <- drop(D %*% Qinv %*% centered)
    mean_T <- drop(E %*% Qinv %*% centered)
  }
  list(VV = VV, TT = TT, TV = TV, Q = Q, D = D, E = E,
    Q_self = Q[seq_along(covariates), seq_along(covariates), drop = FALSE],
    Q_sib = if (n > 1L) Q[seq_along(covariates), length(covariates) +
      seq_along(covariates), drop = FALSE] else NULL,
    projection = projection, hidden = hidden,
    conditional_V = conditional_V, conditional_T = conditional_T,
    conditional_TV = conditional_TV, mean_V = mean_V, mean_T = mean_T,
    mean_C = mean_c, covariates = covariates, n_children = n)
}

mixed_covariate_observed <- function(family, child_index,
    age = rep(10, length(child_index)), sex = rep(0, length(child_index)),
    beta = c(child = .3, mother = .12, father = -.04),
    age_beta = c(child = 0, mother = 0, father = 0),
    sex_beta = c(child = 0, mother = 0, father = 0),
    reliability = c(.65, .4), score_scales = matrix(1, 2L, 3L),
    family_variance = .25, individual_variance = 0, slope_variance = 0,
    intercept_slope_covariance = 0, residual_variance = .75,
    age_center = 10, fixed_mean = rep(0, length(child_index))) {
  n <- family$n_children; r <- n + 2L; m <- length(child_index)
  roles <- c("mother", "father", rep("child", n))
  role_index <- match(roles, c("child", "mother", "father"))
  loading <- rbind(diag(score_scales[1L, role_index], r),
                   diag(score_scales[2L, role_index], r))
  innovation <- diag(r)
  innovation[-(1:2), 1:2] <- .5
  diag(innovation)[-(1:2)] <- sqrt(.5)
  noise <- matrix(0, 2L * r, 2L * r)
  for (k in 1:2) {
    rows <- (k - 1L) * r + seq_len(r)
    A <- score_scales[k, role_index] * innovation * sqrt((1 - reliability[k]) / reliability[k])
    noise[rows, rows] <- tcrossprod(A)
  }
  age0 <- age - age_center
  b <- outer(rep(1, m), beta) + outer(age0, age_beta) + outer(sex, sex_beta)
  B <- matrix(0, m, r)
  B[, 1L] <- b[, "mother"]; B[, 2L] <- b[, "father"]
  B[cbind(seq_len(m), child_index + 2L)] <- b[, "child"]
  same <- outer(child_index, child_index, "==")
  R <- diag(residual_variance, m) + family_variance + same * (
    individual_variance + slope_variance * outer(age0, age0) +
      intercept_slope_covariance * outer(age0, age0, "+"))
  PP <- loading %*% family$conditional_V %*% t(loading) + noise
  PY <- loading %*% t(family$conditional_TV) %*% t(B)
  YY <- B %*% family$conditional_T %*% t(B) + R
  list(covariance = rbind(cbind(PP, PY), cbind(t(PY), YY)),
    mean = c(drop(loading %*% family$mean_V), fixed_mean + drop(B %*% family$mean_T)),
    loading = loading, outcome_loading = B)
}

mixed_covariate_data <- function(N_families = 500L, sib_sizes = 2L,
    visits = 1L, seed = 37111L, contract = mixed_covariate_contract(),
    pi = .4, rho_total = .3, parent_variance = 1.2,
    mode = "estimated_state", n_generations = 3L,
    beta = c(child = .3, mother = .12, father = -.04),
    age_beta = c(child = 0, mother = 0, father = 0),
    sex_beta = c(child = 0, mother = 0, father = 0),
    fixed_pc = rep(.1, contract$n_pc), age_effect = .03, sex_effect = .12,
    age_pattern = c("scheduled", "random"),
    family_variance = .25, individual_variance = 0, slope_variance = 0,
    intercept_slope_covariance = 0, residual_variance = .75,
    families_per_FID = 1L, extended_fraction = 0,
    reliability = c(.65, .4), score_scales = NULL,
    outcome = c("continuous", "binary_liability"), binary_intercept = -1.4) {
  outcome <- match.arg(outcome); age_pattern <- match.arg(age_pattern)
  roles <- c("child", "mother", "father")
  for (x in list(beta, age_beta, sex_beta))
    stopifnot(length(x) == 3L, setequal(names(x), roles), all(is.finite(x)))
  beta <- beta[roles]; age_beta <- age_beta[roles]; sex_beta <- sex_beta[roles]
  if (outcome == "binary_liability" && (residual_variance != 1 || any(visits != 1L)))
    stop("Binary fixtures require one visit and conditional residual variance 1.")
  stopifnot(length(fixed_pc) == contract$n_pc)
  d <- mixed_validation_data(N_families = N_families, sib_sizes = sib_sizes,
    visits = visits, seed = seed, pi = pi, rho_total = rho_total,
    parent_variance = parent_variance, mode = mode, n_generations = n_generations,
    beta = beta, age_effect = 0, sex_effect = 0, age_pattern = "scheduled",
    family_variance = family_variance, individual_variance = individual_variance,
    slope_variance = slope_variance, intercept_slope_covariance = intercept_slope_covariance,
    residual_variance = residual_variance, families_per_FID = families_per_FID,
    extended_fraction = extended_fraction, reliability = reliability, score_scales = score_scales)
  truth <- attr(d, "truth"); latent <- attr(d, "latent")
  previous_seed <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit(if (is.null(previous_seed)) {
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", previous_seed, .GlobalEnv), add = TRUE)
  set.seed(seed + 1000000L)
  family <- match(latent$M_ID, unique(latent$M_ID))
  subject <- match(d$IID, latent$IID)
  V <- sweep(as.matrix(latent[paste0("G_", c("child", "mother", "father"))]),
             2L, truth$population$visible_raw_sd, "*")
  between <- if (any(contract$pc_between != 0))
    .mixed_validation_mvn(N_families, contract$pc_between) else matrix(0, N_families, contract$n_pc)
  pc <- V %*% t(contract$pc_loading) + between[family, , drop = FALSE] +
    .mixed_validation_mvn(nrow(latent), contract$pc_within)
  colnames(pc) <- contract$pc_names
  sex <- as.numeric(drop(pc %*% contract$sex_loading) + stats::rnorm(nrow(latent)) > 0)
  baseline <- contract$age_center + drop(pc %*% contract$baseline_loading) +
    sqrt(contract$baseline_between) * stats::rnorm(N_families)[family] +
    sqrt(contract$baseline_within) * stats::rnorm(nrow(latent))
  counts <- tabulate(subject, nbins = nrow(latent))
  offset <- unlist(lapply(counts, function(k) {
    z <- if (k == 1L) 0 else if (age_pattern == "scheduled") seq(-3, 3, length.out = k) else
      stats::runif(k, -3, 3)
    z - mean(z)
  }), use.names = FALSE)
  age <- baseline[subject] + offset
  age0 <- age - contract$age_center
  T <- as.matrix(latent[paste0("T_", c("child", "mother", "father"))])
  liability <- d$Y + latent$individual_slope[subject] *
    (age0 - (d$age - truth$age_center)) +
    age0 * drop(T[subject, , drop = FALSE] %*% age_beta) +
    sex[subject] * drop(T[subject, , drop = FALSE] %*% sex_beta) +
    drop(pc[subject, , drop = FALSE] %*% fixed_pc) + age_effect * age0 + sex_effect * sex[subject]
  if (outcome == "binary_liability") liability <- liability + binary_intercept
  d$Y <- if (outcome == "binary_liability") as.numeric(liability > 0) else liability
  d$age <- age; d$sex <- sex[subject]; d$baseline_age <- baseline[subject]
  for (j in seq_len(contract$n_pc)) d[[contract$pc_names[j]]] <- pc[subject, j]
  labels <- c("beta_D", "beta_M", "beta_F")
  for (scale in names(truth$coefficients)) {
    transform <- function(b) switch(scale, corrected_raw = b,
      corrected_std = b * truth$population$total_sd,
      visible = drop(truth$population$total_to_visible %*% b))
    truth$coefficients[[scale]] <- c(stats::setNames(transform(beta), labels),
      stats::setNames(transform(age_beta), paste0(labels, ":age")),
      stats::setNames(transform(sex_beta), paste0(labels, ":sex")))
  }
  truth$age_beta_raw <- age_beta; truth$sex_beta_raw <- sex_beta
  truth$fixed_effects <- c(intercept = if (outcome == "binary_liability") binary_intercept else 0,
    stats::setNames(fixed_pc, contract$pc_names), age = age_effect, sex = sex_effect)
  truth$outcome <- outcome; truth$age_center <- contract$age_center
  truth$covariate_contract <- contract; truth$covariate_model <- "marginal Schur; hidden projection residual independent of covariates"
  truth$covariates <- c(contract$pc_names, "sex", "age")
  truth$prevalence_reference <- truth$expected_prevalence <- NA_real_
  truth$omitted_covariate_independence <- NULL
  latent$sex <- sex; latent$baseline_age <- baseline
  latent[contract$pc_names] <- as.data.frame(pc)
  attr(d, "truth") <- truth; attr(d, "latent") <- latent
  attr(d, "liability") <- liability
  d
}

# Independent unprofiled Gaussian calculations and finite differences. This
# audits the evaluator/inference, while the primitive oracle above separately
# audits the covariance supplied to that evaluator.
mixed_covariate_engine_audit <- function(fit) {
  internal <- function(name) getFromNamespace(name, "RAVEL")
  p <- fit$prepared; s <- fit$settings
  stopifnot(!is.null(p), identical(s$outcome, "continuous"), !is.null(s$covariate),
            !is.null(fit$uncertainty))
  p0 <- internal(".ravel_mixed_prepare")(p$data, s$outcome, p$col_outcome,
    p$pgs_columns, p$family_id, p$subject_id, p$cluster, s$covariates,
    s$interactions, s$centers, s$family_random, s$individual_random,
    s$random_slope, s$mode, s$generations, s$pi, FALSE, p$paternal_id, "joint")
  p0$spec$person_mean_covariates <- s$person_mean_covariates
  p0 <- internal(".ravel_mixed_covariate_prepare")(p0)
  raw <- p0$groups
  names(raw) <- vapply(raw, function(g) g$families[1L], character(1))
  designs <- lapply(p$groups, function(g)
    do.call(rbind, lapply(raw[g$families], `[[`, "D")))
  k <- length(fit$parameters); parameter <- c(fit$parameters, fit$means)
  family_nll <- function(value, spec = s) {
    unlist(lapply(seq_along(p$groups), function(i) {
      g <- p$groups[[i]]
      joint <- internal(".ravel_mixed_covariate_joint")(value[seq_len(k)], g, spec)
      ch <- chol(joint$S)
      residual <- g$observed - joint$offset -
        matrix(designs[[i]] %*% value[-seq_len(k)], g$n, byrow = TRUE)
      result <- .5 * (2 * sum(log(diag(ch))) + ncol(residual) * log(2 * pi) +
        rowSums((residual %*% chol2inv(ch)) * residual))
      stats::setNames(result, g$families)
    }), use.names = TRUE)
  }
  difference <- function(fun, value, step = 2e-5) {
    out <- matrix(0, length(fun(value)), length(value))
    for (j in seq_along(value)) {
      a <- b <- value; a[j] <- a[j] + step; b[j] <- b[j] - step
      out[, j] <- (fun(a) - fun(b)) / (2 * step)
    }
    out
  }
  information <- stats::optimHess(parameter, function(x) sum(family_nll(x)),
    control = list(ndeps = rep(1e-4, length(parameter))))
  scores <- -difference(family_nll, parameter)
  at <- match(fit$uncertainty$family_ids, names(family_nll(parameter)))
  scores <- scores[at, , drop = FALSE]
  eta <- s$covariate$eta
  eta_nll <- function(value) {
    spec <- s; spec$covariate$eta <- value
    sum(family_nll(parameter, spec))
  }
  evaluate <- internal(".ravel_mixed_covariate_evaluator")(p, s)
  eta_gradient <- evaluate(fit$parameters, means = fit$means, moments = TRUE)$moment_gradient
  cross_score <- matrix(0, length(parameter), length(eta))
  for (j in seq_along(eta)) {
    a <- b <- s; step <- 2e-4
    a$covariate$eta[j] <- a$covariate$eta[j] + step
    b$covariate$eta[j] <- b$covariate$eta[j] - step
    cross_score[, j] <- -(colSums(difference(function(x) family_nll(x, a), parameter)) -
      colSums(difference(function(x) family_nll(x, b), parameter))) / (2 * step)
  }
  # Saved fits deliberately omit the large family-by-moment influence matrix.
  # Reconstruct that stage from the retained raw data for the independent audit.
  cv <- p0$spec$covariate; triangle <- cv$triangle
  stopifnot(isTRUE(all.equal(cv$eta, s$covariate$eta, tolerance = 1e-12)))
  values <- lapply(cv$family_ids, function(id)
    matrix(raw[[id]]$C, nrow = raw[[id]]$nchild, byrow = TRUE))
  counts <- vapply(values, nrow, integer(1)); pairs <- counts * (counts - 1L)
  first <- do.call(rbind, lapply(values, colSums))
  second <- do.call(rbind, lapply(values, function(x) crossprod(x)[triangle]))
  sibling <- do.call(rbind, lapply(values, function(x)
    (tcrossprod(colSums(x)) - crossprod(x))[triangle]))
  weighted_moments <- function(weights) {
    c(drop(crossprod(weights, first)) / sum(weights * counts),
      drop(crossprod(weights, second)) / sum(weights * counts),
      if (cv$has_siblings) drop(crossprod(weights, sibling)) / sum(weights * pairs))
  }
  n <- length(values)
  moment_influence <- n * t(difference(weighted_moments, rep(1, n), step = 1e-3))
  moment_error <- max(abs(moment_influence - cv$moment_scores))
  moment_influence <- moment_influence[match(fit$uncertainty$family_ids, cv$family_ids), , drop = FALSE]
  effective <- scores + moment_influence %*% t(cross_score / n)
  cluster_scores <- rowsum(effective, fit$uncertainty$cluster_ids, reorder = FALSE)
  covariance <- solve(information) %*% crossprod(cluster_scores) %*% solve(information)
  errors <- c(information_relative_error = max(abs(information - fit$uncertainty$information)) / max(abs(information)),
    family_score_max_error = max(abs(scores - fit$uncertainty$family_scores)),
    moment_gradient_max_error = max(abs(drop(difference(eta_nll, eta)) - eta_gradient)),
    moment_influence_max_error = moment_error,
    cross_derivative_relative_error = max(abs(cross_score - fit$uncertainty$nuisance_score_sensitivity)) /
      max(abs(cross_score)),
    stacked_covariance_relative_error = max(abs(covariance - fit$uncertainty$vcov)) / max(abs(covariance)))
  stopifnot(all(is.finite(errors)))
  errors
}
