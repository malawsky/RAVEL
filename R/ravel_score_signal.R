# Optional joint-likelihood diagnostic for score-specific outcome association.

#' Evaluate an expression while retaining warnings and errors
#' @noRd
.ravel_score_signal_capture <- function(expr) {
  warnings <- character(0)
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = if (inherits(value, "error")) NULL else value,
    error = if (inherits(value, "error")) conditionMessage(value) else NULL,
    warnings = unique(warnings))
}

#' Construct the observable alternative without changing score restrictions
#' @noRd
.ravel_score_signal_syntax <- function(model, scores, latents, type = c("H1", "score")) {
  type <- match.arg(type)
  lines <- strsplit(model, "\n", fixed = TRUE)[[1]]
  # Outcome paths are spread across the trio and sibling builders. Remove all
  # of them, including defined parent-minus-sibling contrasts. History
  # definitions are retained because nonlinear constraints refer to them.
  keep <- !grepl("^\\s*Y\\s*~", lines) &
    !(grepl(":=", lines, fixed = TRUE) &
        grepl("\\b(bD|bm|bf|bms|bfs)\\b", lines, perl = TRUE))
  lines <- lines[keep]
  if (type == "H1") lines <- c(lines,
    paste("Y ~", paste(paste0("score_q", seq_along(scores), "*", scores), collapse = " + ")),
    paste0("Y ~~ 0*", latents))
  paste(lines, collapse = "\n")
}

#' Stable natural-unit scales for constrained numerical diagnostics
#' @noRd
.ravel_score_signal_scales <- function(fit) {
  pt <- lavaan::parameterTable(fit)
  S <- lavaan::lavInspect(fit, "cov.ov")
  if (is.list(S)) S <- S[[1]]
  sd <- sqrt(diag(S))
  latent <- lavaan::lavNames(fit, "lv")
  sd <- c(sd, stats::setNames(rep(1, length(latent)), latent))
  scale <- rep(1, max(pt$free))
  for (i in which(pt$free > 0L)) {
    z <- switch(pt$op[i], "=~" = sd[[pt$rhs[i]]],
      "~" = sd[[pt$lhs[i]]] / sd[[pt$rhs[i]]],
      "~~" = sd[[pt$lhs[i]]] * sd[[pt$rhs[i]]],
      "~1" = sd[[pt$lhs[i]]], 1)
    scale[pt$free[i]] <- z
  }
  list(parameter = scale, observed = sd[colnames(S)])
}

#' Constrained moment-Jacobian rank and tangent diagnostics
#' @noRd
.ravel_score_signal_rank <- function(fit, tolerance = 1e-7, weak_tolerance = 1e-5) {
  J <- lavaan::lavInspect(fit, "delta")
  if (is.list(J)) J <- do.call(rbind, J)
  J <- as.matrix(J)
  scales <- .ravel_score_signal_scales(fit)
  p <- ncol(J)
  if (length(scales$parameter) != p || any(!is.finite(J)))
    stop("The fitted moment Jacobian has an unsupported parameter representation.")
  # lavaan's covariance-moment rows use lower-triangular column order; free
  # observed means precede these rows. Scale both by their natural units.
  sd <- scales$observed
  covariance_scale <- outer(sd, sd)[lower.tri(diag(length(sd)), diag = TRUE)]
  moment_scale <- c(if (isTRUE(fit@Options$meanstructure)) sd else numeric(0), covariance_scale)
  if (length(moment_scale) != nrow(J) || any(!is.finite(moment_scale)) || any(moment_scale <= 0))
    stop("The moment ordering or marginal scales are unsupported.")
  D <- sweep(sweep(J, 1, moment_scale, "/"), 2, scales$parameter, "*")
  A <- fit@Model@con.jac
  equal <- attr(A, "ceq.idx")
  inequality <- attr(A, "cin.idx")
  inactive <- attr(A, "inactive.idx")
  active <- setdiff(inequality, inactive)
  if (!length(A)) A <- matrix(0, 0, p)
  if (ncol(A) != p || any(!is.finite(A))) stop("The constraint Jacobian is unavailable.")
  E <- if (length(equal)) A[equal, , drop = FALSE] else matrix(0, 0, p)
  E <- sweep(E, 2, scales$parameter, "*")
  norms <- sqrt(rowSums(E^2))
  E <- E[norms > 0, , drop = FALSE]
  if (nrow(E)) E <- E / norms[norms > 0]
  thresholds <- c(tolerance / 10, tolerance, tolerance * 10)
  erank <- 0L
  constraint_ranks <- rep(0L, length(thresholds))
  N <- diag(p)
  constraint_singular <- numeric(0)
  if (nrow(E)) {
    z <- svd(E, nu = 0, nv = p)
    constraint_singular <- z$d
    constraint_ranks <- vapply(thresholds, function(tol) sum(z$d > tol * max(z$d)), integer(1))
    erank <- constraint_ranks[2]
    if (erank >= p) stop("The constraint tangent has zero dimension.")
    N <- z$v[, seq.int(erank + 1L, p), drop = FALSE]
  }
  singular <- svd(D %*% N, nu = 0, nv = 0)$d
  dimension <- ncol(N)
  ranks <- vapply(thresholds, function(tol) sum(singular > tol * max(singular)), integer(1))
  ratio <- if (length(singular) < dimension) 0 else min(singular) / max(singular)
  # Where lavaan retains a linear equality basis, verify its compatibility.
  # Nonlinear fits may leave that slot empty; the robust comparison explicitly
  # constructs the final equality tangent below.
  K <- fit@Model@eq.constraints.K
  tangent_match <- TRUE
  if (nrow(E)) {
    tangent_match <- !length(K) || (is.matrix(K) && nrow(K) == p && ncol(K) == dimension)
    if (tangent_match && length(K)) {
      KN <- sweep(K, 1, scales$parameter, "/")
      tangent_match <- max(abs(E %*% KN)) <= 1e-5 * max(1, max(abs(KN)))
    }
  }
  list(rank = ranks[2], tangent_dimension = dimension,
    covariance_rank = ranks[2] - if (isTRUE(fit@Options$meanstructure)) length(sd) else 0L,
    constraint_rank = erank, singular_values = singular,
    constraint_singular_values = constraint_singular,
    constraint_ranks_at_tolerances = stats::setNames(constraint_ranks, thresholds),
    minimum_relative_singular_value = ratio, ranks_at_tolerances = stats::setNames(ranks, thresholds),
    rank_tolerance = tolerance, weak_tolerance = weak_tolerance,
    identified = all(ranks == dimension) && all(constraint_ranks == erank) && ratio >= weak_tolerance,
    active_inequalities = active, tangent_match = tangent_match,
    tangent = N, parameter_scale = scales$parameter)
}

#' Inspect regularity of a fitted diagnostic model
#' @noRd
.ravel_score_signal_validate_fit <- function(fit) {
  rank <- .ravel_score_signal_rank(fit)
  gradient <- as.numeric(lavaan::lavInspect(fit, "gradient"))
  if (length(gradient) != length(rank$parameter_scale)) stop("The score gradient is unavailable.")
  projected <- as.vector(crossprod(rank$tangent, gradient * rank$parameter_scale))
  pt <- lavaan::parameterTable(fit)
  pe <- lavaan::parameterEstimates(fit)
  values <- new.env(parent = baseenv())
  for (i in which(nzchar(pe$label))) assign(pe$label[i], pe$est[i], envir = values)
  eq_residual <- numeric(0)
  inequality_margin <- numeric(0)
  for (i in which(pt$op %in% c("==", ">", "<"))) {
    lhs <- eval(parse(text = pt$lhs[i]), envir = values)
    rhs <- eval(parse(text = pt$rhs[i]), envir = values)
    divisor <- 1 + max(abs(c(lhs, rhs)))
    if (pt$op[i] == "==") eq_residual <- c(eq_residual, abs(lhs - rhs) / divisor)
    else inequality_margin <- c(inequality_margin,
      if (pt$op[i] == ">") (lhs - rhs) / divisor else (rhs - lhs) / divisor)
  }
  # Simple inequalities may be compiled into parameter-table bounds.
  for (i in which(pt$free > 0L)) {
    unit <- rank$parameter_scale[pt$free[i]]
    if ("lower" %in% names(pt) && is.finite(pt$lower[i]))
      inequality_margin <- c(inequality_margin, (pt$est[i] - pt$lower[i]) / unit)
    if ("upper" %in% names(pt) && is.finite(pt$upper[i]))
      inequality_margin <- c(inequality_margin, (pt$upper[i] - pt$est[i]) / unit)
  }
  pd <- function(x) {
    if (is.list(x)) x <- x[[1]]
    if (!is.matrix(x) || any(!is.finite(x)) || any(diag(x) <= 0)) return(NA_real_)
    min(eigen(stats::cov2cor((x + t(x)) / 2), symmetric = TRUE, only.values = TRUE)$values)
  }
  eigenvalues <- c(observed = pd(lavaan::lavInspect(fit, "cov.ov")),
    latent = pd(lavaan::lavInspect(fit, "cov.lv")),
    residual = pd(lavaan::lavInspect(fit, "theta")))
  observed_variance <- diag(lavaan::lavInspect(fit, "cov.ov"))
  residual_variance <- diag(lavaan::lavInspect(fit, "theta"))
  minimum_residual_fraction <- min(residual_variance / observed_variance[names(residual_variance)])
  loading <- pe[pe$op == "=~", , drop = FALSE]
  sd <- .ravel_score_signal_scales(fit)$observed
  minimum_loading <- min(abs(loading$est / sd[loading$rhs]))
  max_gradient <- max(abs(projected))
  max_equality <- max(c(0, eq_residual))
  margin <- min(c(Inf, inequality_margin))
  converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
  reasons <- c(if (!converged) "optimizer did not converge",
    if (!is.finite(max_gradient) || max_gradient > 1e-4) "projected gradient is too large",
    if (!is.finite(max_equality) || max_equality > 1e-6) "equality constraints are not satisfied",
    if (length(rank$active_inequalities) || margin <= 1e-7) "an inequality boundary is active or too close",
    if (any(!is.finite(eigenvalues)) || any(eigenvalues <= 1e-8)) "a covariance matrix is not strictly positive definite",
    if (!is.finite(minimum_residual_fraction) || minimum_residual_fraction <= 1e-8)
      "a residual variance is too close to zero",
    if (!is.finite(minimum_loading) || minimum_loading <= 1e-5) "a standardized loading is too close to zero",
    if (!rank$identified) "the constrained moment Jacobian is singular or weakly identified",
    if (!rank$tangent_match) "the robust-test tangent does not match the final constraints")
  list(valid = !length(reasons), reason = paste(reasons, collapse = "; "),
    converged = converged, projected_gradient_max = max_gradient,
    equality_residual_max = max_equality, inequality_margin_min = margin,
    covariance_min_eigenvalues = eigenvalues, minimum_standardized_loading = minimum_loading,
    minimum_residual_variance_fraction = minimum_residual_fraction,
    thresholds = list(projected_gradient = 1e-4, equality_residual = 1e-6,
      inequality_margin = 1e-7, covariance_eigenvalue = 1e-8, standardized_loading = 1e-5,
      residual_variance_fraction = 1e-8),
    rank = rank)
}

#' Fit one diagnostic model with explicit starts and preserved likelihood options
#' @noRd
.ravel_score_signal_fit <- function(model, dat, reference, estimator, start_values = NULL, do_fit = TRUE,
                                    reference_units = NULL) {
  args <- list(model = model, data = dat, std.lv = TRUE, fixed.x = FALSE,
    conditional.x = FALSE, estimator = estimator, missing = "fiml",
    meanstructure = isTRUE(reference@Options$meanstructure),
    likelihood = reference@Options$likelihood, check.gradient = TRUE, warn = TRUE,
    control = list(iter.max = 20000L, rel.tol = 1e-10))
  initial <- do.call(lavaan::sem, c(args, list(do.fit = FALSE)))
  pt <- lavaan::parameterTable(initial)
  old <- lavaan::parameterTable(reference)
  key <- function(x) paste(x$lhs, x$op, x$rhs, sep = "\r")
  pos <- match(key(pt), key(old))
  good <- !is.na(pos) & pt$op %in% c("=~", "~", "~~", "~1")
  pt$ustart[good] <- old$est[pos[good]]
  if (!is.null(reference_units)) {
    units <- c(reference_units, stats::setNames(rep(1, length(lavaan::lavNames(reference, "lv"))),
      lavaan::lavNames(reference, "lv")))
    for (i in which(good)) {
      divisor <- switch(pt$op[i], "=~" = units[[pt$rhs[i]]] / units[[pt$lhs[i]]],
        "~" = units[[pt$lhs[i]]] / units[[pt$rhs[i]]],
        "~~" = units[[pt$lhs[i]]] * units[[pt$rhs[i]]],
        "~1" = units[[pt$lhs[i]]], 1)
      pt$ustart[i] <- pt$ustart[i] / divisor
    }
  }
  if (!is.null(start_values)) {
    selected <- match(key(pt), names(start_values))
    good <- !is.na(selected)
    pt$ustart[good] <- start_values[selected[good]]
  }
  start <- as.numeric(initial@Fit@x)
  for (i in which(pt$free > 0L & is.finite(pt$ustart))) start[pt$free[i]] <- pt$ustart[i]
  args$model <- pt
  args$start <- start
  args$do.fit <- do_fit
  do.call(lavaan::sem, args)
}

#' Resolve numerical trial warnings only by a clean, equivalent refit
#' @noRd
.ravel_score_signal_fit_checked <- function(model, dat, reference, estimator,
                                           start_values = NULL, reference_units = NULL) {
  first <- .ravel_score_signal_capture(.ravel_score_signal_fit(model, dat, reference,
    estimator, start_values = start_values, reference_units = reference_units))
  trial_warnings <- c("NaNs produced", "NA/NaN function evaluation")
  if (is.null(first$value) || !length(first$warnings) ||
      !all(first$warnings %in% trial_warnings)) return(first)
  before <- .ravel_score_signal_capture(.ravel_score_signal_validate_fit(first$value))
  if (!isTRUE(before$value$valid) || length(before$warnings)) return(first)
  # A finite, regular optimum can follow an inadmissible line-search proposal.
  # Do not discard that warning: require an independently clean reoptimization
  # from the final solution with the same likelihood and all regularity checks.
  second <- .ravel_score_signal_capture(.ravel_score_signal_fit(model, dat,
    first$value, estimator))
  after <- if (!is.null(second$value))
    .ravel_score_signal_capture(.ravel_score_signal_validate_fit(second$value)) else NULL
  difference <- if (isTRUE(after$value$valid))
    as.numeric(lavaan::fitMeasures(second$value, "logl") - lavaan::fitMeasures(first$value, "logl")) else NA_real_
  resolved <- !length(second$warnings) && !length(after$warnings) &&
    isTRUE(after$value$valid) && is.finite(difference) && abs(difference) <= 1e-5
  record <- list(initial_warnings = first$warnings, retry_warnings = second$warnings,
    validation_warnings = after$warnings, validation_error = after$error,
    retry_error = second$error,
    resolved = resolved, log_likelihood_difference = difference,
    likelihood_tolerance = 1e-5, initial_validation = before$value,
    retry_validation = after$value)
  answer <- if (resolved) second else first
  answer$retry <- record
  answer
}

#' Rebuild the supported null in common diagnostic units
#' @noRd
.ravel_score_signal_model <- function(result, units = NULL) {
  fixed <- result[["fixed_crossperson"]]
  if (!is.null(units) && length(fixed)) {
    role_names <- c(m = "mother", f = "father", c = "child", ms = "mat_sib", fs = "pat_sib")
    for (label in names(fixed)) {
      version <- substr(label, 2L, 2L)
      pair <- substring(label, 4L)
      roles <- if (grepl("_", pair, fixed = TRUE)) strsplit(pair, "_", fixed = TRUE)[[1]] else
        strsplit(pair, "", fixed = TRUE)[[1]]
      columns <- paste0("PGS", version, "_", role_names[roles])
      if (length(columns) != 2L || any(!columns %in% names(units)))
        stop("An unrecognized fixed residual covariance prevents unit conversion.")
      fixed[[label]] <- fixed[[label]] / prod(units[columns])
    }
  }
  mode <- result$constraint_spec$mode %||% "none"
  parental_sibs <- result$parental_sibs$mode %||% "none"
  model <- .build_model(crossperson = result$crossperson_structure %||% "free",
    fixed_crossperson = fixed, fixed_crossperson_corr = result[["fixed_crossperson_corr"]],
    parental_sibs = parental_sibs, constrained = mode != "none")
  if (mode != "none") model <- paste(model, .ravel_constraint_syntax(mode,
    result$constraint_spec$h, result$constraint_spec$n_generations,
    .parental_sib_roles(parental_sibs)), sep = "\n")
  if (!is.null(units) && mode != "none") {
    role_names <- c(m = "mother", f = "father", c = "child", ms = "mat_sib", fs = "pat_sib")
    roles <- c("m", "f", "c", .parental_sib_roles(parental_sibs))
    for (j in 1:2) for (role in roles) {
      sd <- units[[paste0("PGS", j, "_", role_names[[role]])]]
      model <- gsub(sprintf("l%s%d > 0.00000001", role, j),
        sprintf("l%s%d > %s", role, j, .format_lavaan_numeric(1e-8 / sd)), model, fixed = TRUE)
      model <- gsub(sprintf("th%d_%s > 0.00000001", j, role),
        sprintf("th%d_%s > %s", j, role, .format_lavaan_numeric(1e-8 / sd^2)), model, fixed = TRUE)
    }
  }
  model
}

#' Ignore only comments, starting values and whitespace when checking model identity
#' @noRd
.ravel_score_signal_model_key <- function(model) {
  lines <- strsplit(model, "\n", fixed = TRUE)[[1]]
  model <- paste(sub("#.*$", "", lines), collapse = "")
  model <- gsub("start\\([^)]*\\)\\*", "", model)
  gsub("[[:space:]]", "", model)
}

#' Supply the final nonlinear equality tangent to the robust delta comparison
#' @noRd
.ravel_score_signal_lrt_fit <- function(fit, rank) {
  # Some supported versions only populate this basis for linear constraints.
  # lavTestLRT's delta method otherwise uses unconstrained nuisance directions.
  # Work on a copy, convert the scaled basis back to parameter units, and
  # orthonormalize it; the original scientific and diagnostic fits stay intact.
  if (rank$constraint_rank > 0L) {
    basis <- sweep(rank$tangent, 1, rank$parameter_scale, "*")
    basis <- qr.Q(qr(basis), complete = FALSE)
    if (ncol(basis) != rank$tangent_dimension || any(!is.finite(basis)))
      stop("The final equality tangent could not be constructed.")
    fit@Model@eq.constraints <- TRUE
    fit@Model@eq.constraints.K <- basis
  }
  fit
}

#' Tangent-projected Satorra (2000) statistic computed directly
#'
#' Independent of whether lavTestLRT honours the supplied eq.constraints.K
#' slot: projects the information onto the final equality-constraint tangent
#' and forms the scaled-and-shifted statistic. Used as a runtime self-check.
#' @noRd
.ravel_score_signal_reference_lrt <- function(f0, f1, df) {
  nullspace <- function(A) {
    z <- svd(A, nu = 0, nv = ncol(A))
    rank <- sum(z$d > 1e-8 * max(z$d))
    z$v[, seq.int(rank + 1L, ncol(A)), drop = FALSE]
  }
  fits <- list(f0, f1)
  J <- lapply(fits, function(f) as.matrix(lavaan::lavInspect(f, "delta")))
  tangent <- lapply(1:2, function(i) {
    A <- fits[[i]]@Model@con.jac
    ceq <- attr(A, "ceq.idx")
    if (!length(A) || !length(ceq)) return(diag(ncol(J[[i]])))
    nullspace(A[ceq, , drop = FALSE])
  })
  H <- MASS::ginv(J[[2]] %*% tangent[[2]]) %*% (J[[1]] %*% tangent[[1]])
  A <- t(nullspace(t(H))) %*% t(tangent[[2]])
  P <- lavaan::lavTech(f1, "information")
  N <- tangent[[2]]
  Pn <- N %*% solve(crossprod(N, P %*% N)) %*% t(N)
  Q <- Pn %*% t(A) %*% solve(A %*% Pn %*% t(A)) %*% A %*% Pn
  W <- lavaan::lavTech(f1, "WLS.V")[[1]]
  Gamma <- lavaan::lavTech(f1, "Gamma")[[1]]
  UG <- W %*% J[[2]] %*% Q %*% t(J[[2]]) %*% W %*% Gamma
  scale <- sqrt(df / sum(diag(UG %*% UG)))
  shift <- df - scale * sum(diag(UG))
  raw <- 2 * (lavaan::fitMeasures(f1, "logl") - lavaan::fitMeasures(f0, "logl"))
  c(statistic = unname(raw * scale + shift), scale = unname(scale), shift = unname(shift))
}

#' Test score-specific phenotype associations
#'
#' Compare the shared latent-factor outcome model with unrestricted outcome
#' regression on every observed score, preserving the score model and all its
#' fitting constraints. The alternative coefficients are diagnostic only.
#'
#' @param result A full result from \code{ravel_fit()}, retaining its lavaan fit.
#' @param data Optional original fitting data. Supplied model variables must
#'   match the retained fitting sample exactly. Usually unnecessary.
#' @return A \code{ravel_score_signal_test} list with status, reason, inference,
#'   model metadata, fitted models, and regularity diagnostics. Unsupported or
#'   nonregular comparisons return missing inference with an explicit reason.
#' @details The initial scope is continuous complete independent families with
#'   no covariates and an identified score measurement model. The default free
#'   residual model is unsupported because its score-only model has a scale
#'   invariance. The test cannot certify unbiased effects or genetic correlation
#'   one: outcome signals proportional to both loadings evade this restriction.
#'   Inference conditions on any external pi used in fitting; post-fit pi
#'   uncertainty propagation does not enter this likelihood comparison.
#'
#'   Both diagnostic models use common sample-SD units, with fixed residual
#'   covariances, numerical bounds and initial values transformed accordingly.
#'   Marginal variances remain estimated. This makes the robust delta projection
#'   invariant to positive changes of score units. Returned log likelihoods are
#'   restored to the original input units; the retained diagnostic fits and
#'   their slopes use the common units. They do not replace the scientific fit.
#'
#'   MLR uses the scaled-and-shifted Satorra (2000) difference test with the delta
#'   Jacobian. Copies of both fits receive the final equality-constraint tangent
#'   so nonlinear measurement constraints enter that comparison. Diagnostic
#'   ranks and conditioning use the same final constraint Jacobian. The
#'   \code{crossperson = "equal_pc"} configuration is currently unsupported.
#'   Numerical NaN trial warnings are resolved only when a second fit from the
#'   admissible optimum is warning-free, passes all checks and reproduces its
#'   likelihood. Both attempts and the warning history are recorded under
#'   \code{diagnostics$refit_retries}; other warnings block inference.
#' @export
ravel_test_score_signal <- function(result, data = NULL) {
  if (!is.list(result)) result <- list()
  out <- structure(list(status = "not_run", reason = "", statistic = NA_real_,
    df = NA_real_, p_value = NA_real_, raw_statistic = NA_real_,
    test = NA_character_, adjustment = NULL,
    log_likelihood = c(H0 = NA_real_, H1 = NA_real_),
    model_df = c(H0 = NA_real_, H1 = NA_real_), N = result$N %||% NA_integer_,
    roles = c("child", "mother", "father", c(ms = "maternal_sibling", fs = "paternal_sibling")[result$parental_sibs$roles]),
    measurement_constraints = list(mode = result$constraint_spec$mode %||% "none",
      crossperson = result$crossperson_structure, fixed_covariance = result[["fixed_crossperson"]],
      fixed_correlation = result[["fixed_crossperson_corr"]]),
    mode = result$constraint_spec$mode %||% "none",
    external_pi = result$constraint_spec$h %||% result$h,
    history = result$constraint_spec, software = list(R = as.character(getRversion()),
      lavaan = as.character(utils::packageVersion("lavaan"))),
    diagnostics = list(), fits = list()), class = "ravel_score_signal_test")
  finish <- function(status, reason) {
    out$status <- status
    out$reason <- reason
    out
  }
  if (!inherits(result$fit, "lavaan"))
    return(finish("needs_fit", "A full fitted result is required; effect tables or transfer packets cannot reconstruct a likelihood test."))
  .ensure_lavaan_runtime()
  original <- result$fit
  if (identical(result$converged, FALSE) ||
      !isTRUE(lavaan::lavInspect(original, "converged")))
    return(finish("invalid_H0", "The supplied scientific fit did not converge or failed its fitting checks."))
  opt <- original@Options
  metadata <- result$score_signal_input
  if (!identical(result$outcome, "continuous") || length(lavaan::lavNames(original, "ov.ord")))
    return(finish("unsupported_outcome", "Only continuous ML or MLR outcomes are supported."))
  if (length(result$fixed_effects)) return(finish("unsupported_covariates", "Covariate comparisons have not been validated."))
  if (original@Data@ngroups != 1L) return(finish("unsupported_multigroup", "Multigroup comparisons have not been validated."))
  if (length(original@Data@cluster) || length(original@Data@sampling.weights) ||
      isTRUE(metadata$dependent_families))
    return(finish("unsupported_design", "Weighted or clustered observations are outside the independent-family scope."))
  if (isTRUE(metadata$imputed) || isTRUE(result$imputed))
    return(finish("unsupported_imputed", "Imputed relatives are outside the validated scope."))
  estimator <- opt$estimator.orig %||% opt$estimator
  if (!estimator %in% c("ML", "MLR")) return(finish("unsupported_estimator", "Only ML and MLR likelihoods are supported."))
  if (!isTRUE(opt$std.lv) || isTRUE(opt$fixed.x) || isTRUE(opt$conditional.x))
    return(finish("unsupported_options", "The joint score model requires std.lv=TRUE, fixed.x=FALSE and conditional.x=FALSE."))
  pt <- lavaan::parameterTable(original)
  latents <- lavaan::lavNames(original, "lv")
  expected_latents <- c("Gc", "Gm", "Gf", c(ms = "Gms", fs = "Gfs")[result$parental_sibs$roles])
  if (!setequal(latents, expected_latents)) return(finish("unsupported_model", "Unexpected latent family roles."))
  paths <- pt[pt$lhs == "Y" & pt$op == "~", , drop = FALSE]
  if (!setequal(paths$rhs, latents) || any(paths$free == 0L))
    return(finish("unsupported_fixed_paths", "Every included latent outcome path must be freely estimated under the null."))
  if (!identical(result$crossperson_structure %||% "free", "free"))
    return(finish("unsupported_measurement", "Equal parent-child residual covariance configurations have not been validated for this diagnostic."))
  if (out$mode == "none" && !length(result[["fixed_crossperson"]]) && !length(result[["fixed_crossperson_corr"]]))
    return(finish("unsupported_free_measurement", "The free residual score model has an observational scale invariance; an identified alternative has not been implemented."))
  dat <- tryCatch(lavaan::lavInspect(original, "data"), error = function(e) NULL)
  if (is.list(dat) && length(dat) == 1L) dat <- dat[[1]]
  if (!is.matrix(dat) || !nrow(dat)) return(finish("needs_data", "The full observed fitting matrix is unavailable; refit with raw complete data."))
  if (any(!is.finite(dat)) || isTRUE(metadata$incomplete) ||
      nrow(dat) != as.integer(lavaan::lavInspect(original, "nobs")))
    return(finish("unsupported_missingness", "Incomplete or dropped observations are outside the validated scope; no rows are dropped by this diagnostic."))
  if (!is.null(data)) {
    if (!is.data.frame(data)) return(finish("data_mismatch", "data must be the original fitting data frame."))
    mapped <- data
    for (nm in names(metadata$column_map)) {
      source <- metadata$column_map[[nm]]
      if (source %in% names(mapped)) mapped[[nm]] <- mapped[[source]]
    }
    if (!all(colnames(dat) %in% names(mapped)) ||
        !isTRUE(all.equal(unname(as.matrix(mapped[, colnames(dat), drop = FALSE])), unname(dat),
          tolerance = 0, check.attributes = FALSE)))
      return(finish("data_mismatch", "Supplied model variables must exactly match the observations and units used in the scientific fit."))
    if ("FID" %in% names(data) && (anyNA(data$FID) || anyDuplicated(data$FID)))
      return(finish("unsupported_design", "Repeated or missing family identifiers are outside the independent-family scope."))
    if (isTRUE(attr(data, "imputed"))) return(finish("unsupported_imputed", "Imputed data are outside the validated scope."))
  }
  dat <- as.data.frame(dat)
  supported_model <- .ravel_score_signal_capture(.ravel_score_signal_model(result))
  if (is.null(supported_model$value) || !identical(.ravel_score_signal_model_key(result$model),
      .ravel_score_signal_model_key(supported_model$value)))
    return(finish("unsupported_model", "The saved null differs from the supported generated RAVEL model; its constraints cannot safely be rebuilt in common units."))
  units <- vapply(dat, stats::sd, numeric(1))
  if (any(!is.finite(units)) || any(units <= 0))
    return(finish("invalid_data", "Every model variable must have finite positive sample variance."))
  # The covariance-nesting delta projection uses a Euclidean moment metric.
  # Fit both models in the same canonical units so arbitrary score units cannot
  # change that projection. Fixed covariances and nuisance starts are converted
  # with the data; genetic/history/correlation restrictions are unchanged.
  dat[] <- lapply(names(dat), function(nm) dat[[nm]] / units[[nm]])
  diagnostic_model <- .ravel_score_signal_model(result, units)
  scores <- setdiff(names(dat), "Y")
  if (length(scores) != 2L * length(latents)) return(finish("unsupported_model", "Each included role must have exactly two observed scores."))
  out$roles <- c("child", "mother", "father", c(ms = "maternal_sibling", fs = "paternal_sibling")[result$parental_sibs$roles])
  out$N <- nrow(dat)
  out$diagnostics$conditional_on_pi <- isTRUE(result$constraint_spec$depends_on_pi)
  out$diagnostics$warnings <- list()
  out$diagnostics$estimator <- estimator
  out$diagnostics$assumptions <- "Observed, complete independent families; externally imputed inputs cannot always be detected."
  out$diagnostics$input_standard_deviations <- units
  out$diagnostics$likelihood_units <- "Reported log likelihoods are restored to original input units; retained diagnostic fits use sample-SD units."
  h0 <- .ravel_score_signal_fit_checked(diagnostic_model, dat, original,
    estimator, reference_units = units)
  out$diagnostics$warnings$H0 <- h0$warnings
  out$diagnostics$refit_retries$H0 <- h0$retry
  if (is.null(h0$value)) return(finish("fit_error", paste("Null refit failed:", h0$error)))
  out$fits$H0 <- h0$value
  h0check <- .ravel_score_signal_capture(.ravel_score_signal_validate_fit(h0$value))
  out$diagnostics$warnings$H0_validation <- h0check$warnings
  out$diagnostics$H0 <- h0check$value
  if (is.null(h0check$value) || !isTRUE(h0check$value$valid))
    return(finish("invalid_H0", h0check$error %||% h0check$value$reason))
  score_model <- .ravel_score_signal_syntax(diagnostic_model, scores, latents, "score")
  score <- .ravel_score_signal_fit_checked(score_model, dat[scores], h0$value, estimator)
  out$diagnostics$warnings$score <- score$warnings
  out$diagnostics$refit_retries$score <- score$retry
  if (is.null(score$value)) return(finish("fit_error", paste("Score-only refit failed:", score$error)))
  out$fits$score <- score$value
  scorecheck <- .ravel_score_signal_capture(.ravel_score_signal_validate_fit(score$value))
  out$diagnostics$warnings$score_validation <- scorecheck$warnings
  out$diagnostics$score <- scorecheck$value
  if (is.null(scorecheck$value) || !isTRUE(scorecheck$value$valid))
    return(finish("invalid_score_model", scorecheck$error %||% scorecheck$value$reason))
  S0 <- lavaan::lavInspect(h0$value, "cov.ov")
  mu0 <- if (isTRUE(opt$meanstructure)) lavaan::lavInspect(h0$value, "mean.ov") else
    stats::setNames(rep(0, ncol(S0)), colnames(S0))
  q <- as.vector(solve(S0[scores, scores], S0[scores, "Y"]))
  residual <- as.numeric(S0["Y", "Y"] - crossprod(q, S0[scores, "Y"]))
  starts <- stats::setNames(c(q, residual, mu0[["Y"]] - sum(q * mu0[scores])),
    c(paste("Y", "~", scores, sep = "\r"), "Y\r~~\rY", "Y\r~1\r"))
  h1_model <- .ravel_score_signal_syntax(diagnostic_model, scores, latents, "H1")
  nested <- .ravel_score_signal_capture(.ravel_score_signal_fit(h1_model, dat, h0$value,
    estimator, starts, do_fit = FALSE))
  out$diagnostics$warnings$nesting <- nested$warnings
  if (is.null(nested$value)) return(finish("invalid_nesting", nested$error))
  mapped_S <- lavaan::lavInspect(nested$value, "cov.ov")[rownames(S0), colnames(S0)]
  mapped_mu <- if (isTRUE(opt$meanstructure)) lavaan::lavInspect(nested$value, "mean.ov")[names(mu0)] else mu0
  nesting_error <- max(abs((mapped_S - S0) / outer(sqrt(diag(S0)), sqrt(diag(S0)))),
    abs((mapped_mu - mu0) / sqrt(diag(S0))))
  out$diagnostics$nesting <- list(maximum_scaled_moment_error = nesting_error,
    alternative_outcome_residual = residual, mapped_fit = nested$value)
  if (!is.finite(nesting_error) || nesting_error > 1e-7 || !is.finite(residual) || residual <= 0)
    return(finish("invalid_nesting", "The constructive null-to-alternative moment mapping failed."))
  # On complete observations H1 factors into the fitted score likelihood and
  # unrestricted Gaussian Y|scores regression. Their separate MLEs provide
  # exact joint starts, while the joint refit still estimates every nuisance
  # parameter and produces its robust likelihood ingredients.
  design <- cbind(intercept = 1, as.matrix(dat[scores]))
  conditional <- stats::lm.fit(design, dat$Y)
  if (conditional$rank != ncol(design))
    return(finish("invalid_data", "The unrestricted conditional score regression is rank deficient."))
  h1_starts <- stats::setNames(c(conditional$coefficients[-1],
    mean(conditional$residuals^2), conditional$coefficients[1]),
    c(paste("Y", "~", scores, sep = "\r"), "Y\r~~\rY", "Y\r~1\r"))
  h1 <- .ravel_score_signal_fit_checked(h1_model, dat, score$value,
    estimator, h1_starts)
  out$diagnostics$warnings$H1 <- h1$warnings
  out$diagnostics$refit_retries$H1 <- h1$retry
  if (is.null(h1$value)) return(finish("fit_error", paste("Alternative refit failed:", h1$error)))
  out$fits$H1 <- h1$value
  h1check <- .ravel_score_signal_capture(.ravel_score_signal_validate_fit(h1$value))
  out$diagnostics$warnings$H1_validation <- h1check$warnings
  out$diagnostics$H1 <- h1check$value
  if (is.null(h1check$value) || !isTRUE(h1check$value$valid))
    return(finish("invalid_H1", h1check$error %||% h1check$value$reason))
  alternative_pt <- lavaan::parameterTable(h1$value)
  yp <- alternative_pt[alternative_pt$lhs == "Y" & alternative_pt$op == "~", , drop = FALSE]
  yc <- alternative_pt[alternative_pt$op == "~~" & xor(alternative_pt$lhs == "Y", alternative_pt$rhs == "Y"), , drop = FALSE]
  if (!setequal(yp$rhs, scores) || any(yp$free == 0L) || any(yc$free != 0L | yc$est != 0))
    return(finish("invalid_alternative", "The alternative contains unintended outcome paths or residual covariances."))
  rank_difference <- h1check$value$rank$rank - h0check$value$rank$rank
  expected <- length(latents)
  h0_increment <- h0check$value$rank$rank - scorecheck$value$rank$rank
  h1_increment <- h1check$value$rank$rank - scorecheck$value$rank$rank
  mean_increment <- as.integer(isTRUE(opt$meanstructure))
  out$diagnostics$rank_difference <- rank_difference
  out$log_likelihood <- c(H0 = as.numeric(lavaan::fitMeasures(h0$value, "logl")), H1 = as.numeric(lavaan::fitMeasures(h1$value, "logl")))
  out$model_df <- c(H0 = lavaan::fitMeasures(h0$value, "df"), H1 = lavaan::fitMeasures(h1$value, "df"))
  names(out$model_df) <- c("H0", "H1")
  out$diagnostics$canonical_log_likelihood <- out$log_likelihood
  out$log_likelihood <- out$log_likelihood - out$N * sum(log(units))
  if (rank_difference != expected || h0_increment != expected + 1L + mean_increment ||
      h1_increment != 2L * expected + 1L + mean_increment || diff(rev(out$model_df)) != expected)
    return(finish("invalid_rank_difference", "The actual constrained dimensions disagree with the identified shared-factor comparison."))
  raw <- 2 * (out$log_likelihood[["H1"]] - out$log_likelihood[["H0"]])
  out$raw_statistic <- raw
  if (!is.finite(raw) || raw < -1e-6) return(finish("negative_likelihood_improvement", "The alternative likelihood is materially worse than the null; inference is unavailable."))
  # A harmless negative value below floating-point resolution is recorded in
  # raw_statistic; only the numerical round-off is removed for the ML tail.
  if (any(lengths(out$diagnostics$warnings) > 0L))
    return(finish("fit_warning", "A refit or nesting calculation emitted a warning; inspect diagnostics$warnings before using this comparison."))
  if (estimator == "ML") {
    out$statistic <- max(0, raw)
    out$df <- expected
    out$p_value <- stats::pchisq(out$statistic, expected, lower.tail = FALSE)
    out$test <- "Gaussian ML likelihood-ratio test"
  } else {
    available <- names(formals(lavaan::lavTestLRT))
    robust_H0 <- .ravel_score_signal_lrt_fit(h0$value, h0check$value$rank)
    robust_H1 <- .ravel_score_signal_lrt_fit(h1$value, h1check$value$rank)
    args <- list(robust_H0, robust_H1, method = "satorra.2000")
    out$diagnostics$robust_tangent <- "Final equality-constraint Jacobian, projected in parameter units"
    delta_name <- if ("a_method" %in% available) "a_method" else "A.method"
    shifted_name <- if ("scaled_shifted" %in% available) "scaled_shifted" else "scaled.shifted"
    args[[delta_name]] <- "delta"
    args[[shifted_name]] <- TRUE
    if ("check_nested" %in% available) args$check_nested <- TRUE
    comparison <- .ravel_score_signal_capture(do.call(lavaan::lavTestLRT, args))
    out$diagnostics$warnings$test <- comparison$warnings
    out$diagnostics$test_table <- comparison$value
    if (is.null(comparison$value) || length(comparison$warnings))
      return(finish("invalid_robust_test", comparison$error %||% paste(comparison$warnings, collapse = "; ")))
    tab <- comparison$value
    row <- which(is.finite(tab[["Df diff"]]))
    heading <- paste(attr(tab, "heading"), collapse = " ")
    scale <- attr(tab, "scale")[row]
    shift <- attr(tab, "shift")[row]
    if (length(row) != 1L || !grepl("satorra.2000", heading, fixed = TRUE) ||
        length(scale) != 1L || !is.finite(scale) || scale <= 0 ||
        length(shift) != 1L || !is.finite(shift) || tab[["Df diff"]][row] != expected ||
        !is.finite(tab[["Chisq diff"]][row]) || tab[["Chisq diff"]][row] < 0 ||
        !is.finite(tab[["Pr(>Chisq)"]][row]))
      return(finish("invalid_robust_test", "A finite scaled-and-shifted delta-Jacobian comparison was not returned; no unscaled fallback is reported."))
    # lavTestLRT is driven through internal slots; recompute the projected
    # statistic directly and refuse to report if the two disagree.
    reference <- tryCatch(.ravel_score_signal_reference_lrt(h0$value, h1$value, expected),
      error = function(e) NULL)
    out$diagnostics$reference_statistic <- reference
    if (is.null(reference) || !all(is.finite(reference)) ||
        abs(reference[["statistic"]] - tab[["Chisq diff"]][row]) >
          1e-6 * max(1, abs(reference[["statistic"]])))
      return(finish("invalid_robust_test", paste0("lavTestLRT's robust statistic did not match the directly ",
        "computed tangent-projected statistic (lavaan ", as.character(utils::packageVersion("lavaan")),
        "); no robust p value is reported.")))
    out$statistic <- tab[["Chisq diff"]][row]
    out$df <- tab[["Df diff"]][row]
    out$p_value <- tab[["Pr(>Chisq)"]][row]
    out$test <- "MLR robust difference test: satorra.2000, delta Jacobian"
    out$adjustment <- list(method = "satorra.2000", jacobian = "delta", scaled_shifted = TRUE,
      tangent = "final equality-constraint Jacobian",
      scale = scale, shift = shift, heading = heading)
  }
  finish("ok", "Score-phenotype compatibility test under the retained measurement assumptions.")
}

#' Print a score-specific association diagnostic
#' @param x A score-specific association test result.
#' @param ... Unused.
#' @export
print.ravel_score_signal_test <- function(x, ...) {
  cat("Score-specific phenotype association test\n")
  cat("Status:", x$status, "\n")
  if (identical(x$status, "ok")) {
    cat(sprintf("Statistic: %.4f; df: %g; p: %.4g; N: %d\n", x$statistic, x$df, x$p_value, x$N))
    cat(x$test, "\n")
    if (isTRUE(x$diagnostics$conditional_on_pi)) cat("Inference is conditional on pi used in fitting.\n")
  }
  cat(x$reason, "\n")
  invisible(x)
}
