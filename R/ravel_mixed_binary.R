# Optional cross-sectional probit family likelihood. Repeated binary outcomes
# remain outside this interface. The genotype and outcome likelihood are joint.

.ravel_mixed_binary_theta <- function(par, covariance_names) {
  theta <- c(par[covariance_names], residual_sd = 0)
  if ("family_variance" %in% names(theta)) {
    theta["family_sd"] <- if (theta[["family_variance"]] == 0) -Inf else
      .5 * log(theta[["family_variance"]])
    theta <- theta[names(theta) != "family_variance"]
  }
  theta
}

.ravel_mixed_normal_probability <- function(upper, covariance, mean = rep(0, length(upper))) {
  d <- length(upper)
  if (!d) return(1)
  if (d == 1L) return(stats::pnorm((upper - mean) / sqrt(covariance[1L, 1L])))
  algorithm <- if (d <= 3L) mvtnorm::TVPACK(abseps = 1e-9) else mvtnorm::Miwa(steps = 128L)
  as.numeric(mvtnorm::pmvnorm(lower = rep(-Inf, d), upper = upper, mean = mean,
                            sigma = covariance, algorithm = algorithm))
}

.ravel_mixed_probit_probability <- function(mu, covariance, y, derivatives = TRUE) {
  m <- length(y)
  if (!m) return(list(log_probability = 0, mean_gradient = numeric(), covariance_gradient = matrix(0, 0, 0)))
  sd <- sqrt(diag(covariance)); sign <- 2 * y - 1
  b <- sign * mu / sd
  C <- covariance / outer(sd, sd) * outer(sign, sign)
  pr <- .ravel_mixed_normal_probability(b, C)
  if (!is.finite(pr) || pr <= 0 || pr > 1 + 1e-8) stop("Probit orthant probability could not be evaluated accurately.")
  if (!derivatives) return(list(log_probability = log(min(pr, 1))))
  db <- numeric(m); dc <- matrix(0, m, m)
  for (i in seq_len(m)) {
    others <- setdiff(seq_len(m), i)
    conditional <- C[others, others, drop = FALSE] - tcrossprod(C[others, i])
    db[i] <- stats::dnorm(b[i]) * .ravel_mixed_normal_probability(b[others], conditional, C[others, i] * b[i])
  }
  if (m > 1L) for (i in seq_len(m - 1L)) for (j in seq.int(i + 1L, m)) {
    rho <- C[i, j]; determinant <- 1 - rho^2
    density <- exp(-(b[i]^2 - 2 * rho * b[i] * b[j] + b[j]^2) / (2 * determinant)) /
      (2 * pi * sqrt(determinant))
    ij <- c(i, j); others <- setdiff(seq_len(m), ij)
    if (length(others)) {
      mapping <- C[others, ij, drop = FALSE] %*% solve(C[ij, ij, drop = FALSE])
      conditional <- C[others, others, drop = FALSE] - mapping %*% C[ij, others, drop = FALSE]
      density <- density * .ravel_mixed_normal_probability(b[others], conditional, drop(mapping %*% b[ij]))
    }
    dc[i, j] <- dc[j, i] <- density
  }
  gradient <- dc * outer(sign / sd, sign / sd) / (2 * pr)
  diag(gradient) <- -.5 * (b * db + rowSums(C * dc)) / (diag(covariance) * pr)
  list(log_probability = log(min(pr, 1)), mean_gradient = sign / sd * db / pr,
       covariance_gradient = gradient)
}

.ravel_mixed_binary_evaluator <- function(groups, spec, covariance_names) {
  cache <- new.env(parent = emptyenv())
  evaluate <- function(par, derivatives = TRUE, keep = FALSE) {
    if (exists("par", cache, inherits = FALSE) && identical(par, cache$par) &&
        (!derivatives || !is.null(cache$value$gradient))) return(cache$value)
    invalid <- function(e) list(value = 1e50, valid = FALSE, reason = conditionMessage(e),
      gradient = if (derivatives) stats::setNames(rep(0, length(par)), names(par)) else NULL)
    answer <- tryCatch({
      theta <- .ravel_mixed_binary_theta(par, covariance_names)
      alpha <- par[spec$mean_names]
      value <- 0; all_scores <- list(); work <- vector("list", length(groups))
      for (g in seq_along(groups)) {
        z <- groups[[g]]; np <- 2L * (z$nchild + 2L); ny <- length(z$child)
        ip <- seq_len(np); yi <- np + seq_len(ny)
        S <- .ravel_mixed_covariance(theta, z, spec)
        P <- S[ip, ip, drop = FALSE]; PY <- S[ip, yi, drop = FALSE]
        Q <- chol2inv(chol(P)); H <- Q %*% PY
        R <- S[yi, yi, drop = FALSE] - crossprod(PY, H)
        Dp <- z$D[ip, , drop = FALSE]; Dy <- z$D[yi, , drop = FALSE]
        residual <- sweep(z$observed[, ip, drop = FALSE], 2, drop(Dp %*% alpha))
        U <- residual %*% Q
        mu <- sweep(residual %*% H, 2, drop(Dy %*% alpha), "+")
        score <- matrix(0, z$n, length(par), dimnames = list(NULL, names(par)))
        gradient_mean <- matrix(0, z$n, ny)
        gradient_cov <- vector("list", z$n)
        log_probability <- numeric(z$n)
        for (i in seq_len(z$n)) {
          probability <- .ravel_mixed_probit_probability(mu[i, ], R, z$observed[i, yi], derivatives)
          log_probability[i] <- probability$log_probability
          if (derivatives) {
            gradient_mean[i, ] <- probability$mean_gradient
            gradient_cov[[i]] <- probability$covariance_gradient
          }
        }
        value <- value + .5 * (z$n * (np * log(2 * pi) + 2 * sum(log(diag(chol(P))))) +
                              sum(residual * U)) - sum(log_probability)
        if (derivatives) {
          score[, spec$mean_names] <- U %*% Dp + gradient_mean %*% (Dy - crossprod(H, Dp))
          for (name in covariance_names) {
            if (name == "family_variance") {
              # The natural variance has a nonvanishing score at zero. A log
              # variance can falsely satisfy a gradient tolerance near zero.
              score[, name] <- vapply(gradient_cov, sum, numeric(1))
              next
            }
            step <- 1e-5 * max(1, abs(theta[name])); plus <- minus <- theta
            plus[name] <- plus[name] + step; minus[name] <- minus[name] - step
            sp <- .ravel_mixed_covariance(plus, z, spec)
            sm <- .ravel_mixed_covariance(minus, z, spec)
            hp <- solve(sp[ip, ip, drop = FALSE], sp[ip, yi, drop = FALSE])
            hm <- solve(sm[ip, ip, drop = FALSE], sm[ip, yi, drop = FALSE])
            dr <- (sp[yi, yi, drop = FALSE] - crossprod(sp[ip, yi, drop = FALSE], hp) -
                    sm[yi, yi, drop = FALSE] + crossprod(sm[ip, yi, drop = FALSE], hm)) / (2 * step)
            dp <- (sp[ip, ip, drop = FALSE] - sm[ip, ip, drop = FALSE]) / (2 * step)
            dm <- residual %*% ((hp - hm) / (2 * step))
            score[, name] <- .5 * (rowSums((U %*% dp) * U) - sum(Q * dp)) +
              rowSums(gradient_mean * dm) + vapply(gradient_cov, function(x) sum(x * dr), numeric(1))
          }
        }
        all_scores[[g]] <- score
        work[[g]] <- list(P = P, conditional_covariance = R, conditional_mean = mu)
      }
      scores <- do.call(rbind, all_scores)
      list(value = value, valid = is.finite(value), gradient = if (derivatives) -colSums(scores) else NULL,
           family_scores = scores, work = work)
    }, error = invalid)
    cache$par <- par; cache$value <- answer; answer
  }
  evaluate
}

.ravel_mixed_fit_binary <- function(prepared, pi_se, diagnostic, control, start, se, keep_data,
                                    diagnostic_control = list()) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) stop("Binary mixed fitting requires the optional mvtnorm package.")
  spec <- prepared$spec
  if (spec$individual_random || !is.null(spec$random_slope)) stop("Cross-sectional binary fitting supports only the maternal-family random intercept.")
  if (max(vapply(prepared$groups, function(g) length(g$child), integer(1))) > 20L)
    stop("Deterministic probit integration currently supports at most 20 observed siblings per maternal family.")
  covariance_names <- setdiff(names(prepared$start), c("residual_sd", "family_sd"))
  par <- prepared$start[covariance_names]
  if (spec$family_random) {
    covariance_names <- c(covariance_names, "family_variance")
    par["family_variance"] <- exp(2 * prepared$start[["family_sd"]])
  }
  par <- c(par, stats::setNames(rep(0, length(spec$mean_names)), spec$mean_names))
  par["mean_intercept"] <- stats::qnorm(mean(prepared$data[[prepared$col_outcome]], na.rm = TRUE))
  lower <- c(prepared$lower[covariance_names], stats::setNames(rep(-Inf, length(spec$mean_names)), spec$mean_names))
  upper <- c(prepared$upper[covariance_names], stats::setNames(rep(Inf, length(spec$mean_names)), spec$mean_names))
  names(lower) <- names(upper) <- names(par)
  if (spec$family_random) { lower["family_variance"] <- 0; upper["family_variance"] <- Inf }
  if (!is.null(start)) {
    if ("family_sd" %in% names(start) && spec$family_random) {
      start["family_variance"] <- exp(2 * start[["family_sd"]])
      start <- start[names(start) != "family_sd"]
    }
    if (is.null(names(start)) || any(!names(start) %in% names(par)) || any(!is.finite(start)))
      stop("start must be a finite named parameter vector.")
    par[names(start)] <- start
  }
  settings <- utils::modifyList(list(maxit = 500L, factr = 1e2, pgtol = 1e-6), control)
  evaluate <- .ravel_mixed_binary_evaluator(prepared$groups, spec, covariance_names)
  initial <- evaluate(par)
  if (!isTRUE(initial$valid)) stop("Invalid initial probit family likelihood: ", initial$reason)
  optimize <- function(ev, initial) stats::optim(initial, function(x) ev(x)$value,
    function(x) ev(x, derivatives = TRUE)$gradient, method = "L-BFGS-B",
    lower = lower, upper = upper, control = settings)
  fit <- optimize(evaluate, par); par <- fit$par; final <- evaluate(par)
  if (!isTRUE(final$valid)) stop("Invalid final probit family likelihood: ", final$reason)
  projected_gradient <- final$gradient
  at_lower <- par - lower < 1e-8
  at_upper <- upper - par < 1e-8
  projected_gradient[(at_lower & projected_gradient > 0) | (at_upper & projected_gradient < 0)] <- 0
  converged <- fit$convergence == 0L && max(abs(projected_gradient)) / prepared$n_families < 1e-5
  messages <- character(); uncertainty <- NULL
  boundary_parameters <- names(par)[at_lower | at_upper]
  boundary <- length(boundary_parameters) > 0L
  if (boundary) messages <- c(messages, paste("Parameters at a boundary:",
    paste(boundary_parameters, collapse = ", "), "; ordinary Wald uncertainty is not certified."))
  if (!converged) messages <- c(messages, "The binary likelihood did not pass optimizer and projected-gradient checks.")
  if (converged && isTRUE(se) && !boundary) uncertainty <- tryCatch({
    hstep <- rep(1e-4, length(par)); names(hstep) <- names(par)
    if (spec$family_random) hstep["family_variance"] <- min(hstep["family_variance"],
      par[["family_variance"]] / 4)
    H <- stats::optimHess(par, function(x) evaluate(x)$value, function(x) evaluate(x)$gradient,
                         control = list(ndeps = hstep))
    H <- (H + t(H)) / 2
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= max(ev) * 1e-9) stop("Binary model information is singular or indefinite.")
    bread <- solve(H)
    clusters <- unlist(lapply(prepared$groups, `[[`, "clusters"), use.names = FALSE)
    scores <- rowsum(final$family_scores, clusters, reorder = FALSE)
    V <- bread %*% crossprod(scores) %*% bread
    dimnames(V) <- dimnames(bread) <- list(names(par), names(par))
    list(vcov = (V + t(V)) / 2, model_vcov = bread, information = H,
         family_scores = final$family_scores, cluster_scores = scores,
         family_ids = unlist(lapply(prepared$groups, `[[`, "families"), use.names = FALSE),
         cluster_ids = clusters,
         method = "FID-clustered CR0 sandwich")
  }, error = function(e) { messages <<- c(messages, conditionMessage(e)); NULL })
  V <- if (is.null(uncertainty)) matrix(NA_real_, length(par), length(par), dimnames = list(names(par), names(par))) else uncertainty$vcov
  maps <- lapply(c("visible", "corrected_std", "corrected_raw"), function(scale)
    .ravel_mixed_delta(function(x) .ravel_mixed_effects(x, spec, scale), par))
  names(maps) <- c("visible", "corrected_std", "corrected_raw")
  varfun <- function(x) {
    # Express the identified variance directly so derivatives remain valid at 0.
    out <- c(residual = 1)
    if (spec$family_random) out["family_intercept"] <- x[["family_variance"]]
    out
  }
  variances <- .ravel_mixed_delta(varfun, par)
  pi_derivatives <- NULL; variance_pi <- NULL
  if (converged && !boundary && isTRUE(se) && !is.null(pi_se) && pi_se > 0) {
    step <- min(1e-3, max(1e-5, pi_se / 10), (1 - spec$pi) / 4)
    values <- if (spec$pi > step) c(spec$pi - step, spec$pi + step) else c(spec$pi, spec$pi + step)
    refits <- lapply(values, function(value) {
      s <- spec; s$pi <- value
      ev <- .ravel_mixed_binary_evaluator(prepared$groups, s, covariance_names)
      f <- optimize(ev, par)
      ef <- ev(f$par)
      list(par = f$par, spec = s, valid = f$convergence == 0L && isTRUE(ef$valid) &&
        max(abs(ef$gradient)) / prepared$n_families < 1e-5 &&
        !any(f$par - lower < 1e-8 | upper - f$par < 1e-8))
    })
    if (all(vapply(refits, `[[`, logical(1), "valid"))) {
      pi_derivatives <- lapply(names(maps), function(scale)
        (.ravel_mixed_effects(refits[[2]]$par, refits[[2]]$spec, scale) -
         .ravel_mixed_effects(refits[[1]]$par, refits[[1]]$spec, scale)) / diff(values))
      names(pi_derivatives) <- names(maps)
      variance_pi <- (varfun(refits[[2]]$par) - varfun(refits[[1]]$par)) / diff(values)
    } else messages <- c(messages, "A binary pi-sensitivity refit failed; pi-propagated uncertainty is unavailable.")
  }
  coefficient_vcov <- list()
  tables <- lapply(names(maps), function(scale) {
    z <- maps[[scale]]; C <- z$jacobian %*% V %*% t(z$jacobian)
    if (!is.null(pi_se) && pi_se > 0) {
      if (is.null(pi_derivatives)) C[,] <- NA_real_ else C <- C + pi_se^2 * tcrossprod(pi_derivatives[[scale]])
    }
    coefficient_vcov[[scale]] <<- C
    error <- sqrt(pmax(diag(C), 0))
    data.frame(parameter = names(z$value), scale = scale, estimate = unname(z$value), se = error,
      ci_low = z$value - stats::qnorm(.975) * error, ci_high = z$value + stats::qnorm(.975) * error,
      p_value = 2 * stats::pnorm(-abs(z$value / error)),
      uncertainty_source = if (is.null(uncertainty)) "unavailable" else uncertainty$method,
      outcome_scale = "conditional probit residual SD = 1", row.names = NULL)
  })
  variance_vcov <- variances$jacobian %*% V %*% t(variances$jacobian)
  if (!is.null(pi_se) && pi_se > 0) {
    if (is.null(variance_pi)) variance_vcov[,] <- NA_real_ else
      variance_vcov <- variance_vcov + pi_se^2 * tcrossprod(variance_pi)
  }
  variance_table <- data.frame(component = names(variances$value), estimate = unname(variances$value),
    se = sqrt(pmax(diag(variance_vcov), 0)), uncertainty_conditions = if (!is.null(pi_se) && pi_se > 0)
      "independent external pi uncertainty propagated" else "conditional on supplied pi", row.names = NULL)
  # The occasion residual is fixed by probit identification, not estimated.
  variance_table$se[variance_table$component == "residual"] <- 0
  result <- structure(list(converged = converged, optimizer_converged = fit$convergence == 0L,
    inference_status = if (!se) "not_requested" else if (boundary) "boundary" else
      if (!converged) "not_converged" else if (is.null(uncertainty)) "unavailable" else
      if (!is.null(pi_se) && pi_se > 0 && is.null(pi_derivatives)) "pi_sensitivity_unavailable" else "ok",
    table = do.call(rbind, tables), coefficient_vcov = coefficient_vcov,
    variance_components = variance_table, variance_vcov = variance_vcov,
    parameters = par[covariance_names], means = par[spec$mean_names],
    vcov = if (is.null(uncertainty)) NULL else uncertainty$vcov,
    model_vcov = if (is.null(uncertainty)) NULL else uncertainty$model_vcov,
    uncertainty = uncertainty, pi_sensitivity = pi_derivatives, variance_pi_sensitivity = variance_pi,
    N = prepared$N, n_subjects = prepared$n_subjects, n_families = prepared$n_families, n_clusters = prepared$n_clusters,
    fit_statistics = c(logLik = -final$value, n_parameters = length(par), AIC = 2 * final$value + 2 * length(par)),
    state = .ravel_mixed_state(par, spec), settings = spec, optimizer = fit,
    gradient = final$gradient, projected_gradient = projected_gradient,
    boundary_parameters = boundary_parameters, messages = messages, score_scaling = prepared$score_scaling,
    family_id = prepared$family_id, subject_id = prepared$subject_id, cluster = prepared$cluster,
    diagnostics = list(full_sibling_assumption = TRUE, paternal_id_checked = !is.null(prepared$paternal_id),
      likelihood = "Joint PGS density and conditional multivariate-probit family probability",
      outcome_scale = "Conditional occasion residual variance fixed to 1; not legacy unit marginal liability variance",
      score_measurements_counted_once_per_person_within_maternal_family = TRUE),
    prepared = if (keep_data) prepared else NULL), class = "ravel_mixed_results")
  if (diagnostic) result$one_factor <- .ravel_mixed_run_diagnostic(prepared, diagnostic_control)
  result
}
