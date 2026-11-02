# Optional family likelihood. This module does not change the original SEM.

# Fixed-effect designs can differ while covariance designs agree. Compress
# their sufficient statistics without forcing each covariate value into its own
# covariance group. The binary likelihood retains its original unmerged groups.
.ravel_mixed_compress_means <- function(groups) {
  index <- new.env(parent = emptyenv()); merged <- list()
  for (z in groups) {
    key <- paste(z$nchild, paste(z$child, collapse = ","),
      paste(format(c(z$moderators, z$slope), digits = 17, scientific = TRUE, trim = TRUE), collapse = ","), sep = "|")
    if (!exists(key, index, inherits = FALSE)) {
      at <- length(merged) + 1L; assign(key, at, index)
      out <- z
      out$observed <- list(z$observed)
      out$design_columns <- lapply(seq_len(ncol(z$D)), function(j)
        matrix(rep(z$D[, j], each = z$n), z$n))
      merged[[at]] <- out
    } else {
      at <- get(key, index); out <- merged[[at]]
      out$observed[[length(out$observed) + 1L]] <- z$observed
      out$design_columns <- lapply(seq_len(ncol(z$D)), function(j)
        rbind(out$design_columns[[j]], matrix(rep(z$D[, j], each = z$n), z$n)))
      out$n <- out$n + z$n
      out$families <- c(out$families, z$families); out$clusters <- c(out$clusters, z$clusters)
      merged[[at]] <- out
    }
  }
  lapply(merged, function(z) {
    z$observed <- do.call(rbind, z$observed)
    z$sum <- colSums(z$observed); z$cross <- crossprod(z$observed)
    q <- ncol(z$D)
    z$design_cross <- lapply(seq_len(q), function(j)
      lapply(seq_len(q), function(k) crossprod(z$design_columns[[j]], z$design_columns[[k]])))
    z$design_outcome_cross <- lapply(z$design_columns, function(x) crossprod(x, z$observed))
    z
  })
}

.ravel_mixed_mean_matrix <- function(group, means) {
  if (is.null(group$design_columns))
    return(matrix(rep(drop(group$D %*% means), each = group$n), nrow = group$n))
  Reduce(`+`, Map(function(x, b) x * b, group$design_columns, means))
}

.ravel_mixed_state <- function(par, spec) {
  rho <- tanh(par[["state_rho"]])
  vp <- switch(spec$mode,
    estimated_state = .5 + exp(par[["state_variance"]]),
    equilibrium = 1 / (1 - rho),
    n_gen = if (abs(rho) < 1e-8) 1 else
      (1 - rho * ((1 + rho) / 2)^(spec$generations - 1L)) / (1 - rho))
  vc <- .5 * vp * (1 + rho) + .5
  vg <- spec$pi + (1 - spec$pi) * vp
  vgc <- spec$pi + (1 - spec$pi) * vc
  list(rho = rho, parent = vp, child = vc, visible_parent = vg,
       visible_child = vgc, rho_visible = (1 - spec$pi) * rho * vp / vg,
       eta = (1 - spec$pi) * .5 * vp * (1 + rho) / sqrt(vg * vgc) +
         spec$pi * .5 / sqrt(vg * vgc))
}

.ravel_mixed_genetic <- function(nchild, par, spec) {
  state <- .ravel_mixed_state(par, spec)
  n <- nchild + 2L
  # Order: mother, father, then each distinct child. The neutral process also
  # supplies the within-GWAS measurement-error kinship used by constrained RAVEL.
  K <- matrix(.5, n, n); diag(K) <- 1; K[1, 2] <- K[2, 1] <- 0
  T <- matrix(.5 * state$parent * (1 + state$rho), n, n)
  T[1:2, 1:2] <- state$parent * matrix(c(1, state$rho, state$rho, 1), 2)
  diag(T)[-(1:2)] <- state$child
  V <- (1 - spec$pi) * T + spec$pi * K
  sdT <- sqrt(diag(T))
  roles <- c("mother", "father", rep("child", nchild))
  a <- lapply(1:2, function(k) exp(par[paste0("loading", k, "_", roles)]))
  noise <- exp(par[paste0("noise", 1:2)])
  P <- matrix(0, 2L * n, 2L * n)
  for (k in 1:2) for (j in 1:2) {
    rows <- (k - 1L) * n + seq_len(n)
    cols <- (j - 1L) * n + seq_len(n)
    P[rows, cols] <- outer(a[[k]], a[[j]]) *
      (V + if (k == j) noise[k]^2 * K else 0)
  }
  PT <- rbind(a[[1]] * T, a[[2]] * T) * sqrt(1 - spec$pi)
  PT <- sweep(PT, 2, sdT, "/")
  list(P = P, PT = PT, T = T / outer(sdT, sdT), state = state)
}

.ravel_mixed_covariance <- function(par, group, spec, genetic = NULL) {
  if (is.null(genetic)) genetic <- .ravel_mixed_genetic(group$nchild, par, spec)
  n <- group$nchild + 2L
  m <- length(group$child)
  B <- matrix(0, m, n)
  beta <- matrix(par[spec$beta_names], nrow = 3L,
                 dimnames = list(c("child", "mother", "father"), NULL))
  effects <- group$moderators %*% t(beta)
  B[, 1] <- effects[, "mother"]
  B[, 2] <- effects[, "father"]
  B[cbind(seq_len(m), group$child + 2L)] <- effects[, "child"]
  R <- diag(exp(2 * par[["residual_sd"]]), m)
  if (spec$family_random) R <- R + exp(2 * par[["family_sd"]])
  same <- outer(group$child, group$child, "==")
  if (spec$individual_random) R <- R + same * exp(2 * par[["individual_sd"]])
  if (!is.null(spec$random_slope)) {
    slope <- exp(par[["slope_sd"]])
    R <- R + same * outer(group$slope, group$slope) * slope^2
    if (spec$individual_random) {
      covariance <- exp(par[["individual_sd"]]) * slope * tanh(par[["intercept_slope_cor"]])
      R <- R + same * outer(group$slope, group$slope, "+") * covariance
    }
  }
  PY <- genetic$PT %*% t(B)
  YY <- B %*% genetic$T %*% t(B) + R
  rbind(cbind(genetic$P, PY), cbind(t(PY), YY))
}

.ravel_mixed_evaluator <- function(groups, spec) {
  cache <- new.env(parent = emptyenv())
  mean_dimension <- ncol(groups[[1]]$D)
  evaluate <- function(par, derivatives = FALSE, keep = FALSE) {
    if (exists("par", cache, inherits = FALSE) && identical(par, cache$par) &&
        (!derivatives || !is.null(cache$value$gradient))) return(cache$value)
    invalid <- function(reason) list(value = 1e50, valid = FALSE, reason = reason,
                              gradient = if (derivatives) rep(0, length(par)) else NULL)
    answer <- tryCatch({
      genetic <- lapply(unique(vapply(groups, `[[`, integer(1), "nchild")),
        function(n) .ravel_mixed_genetic(n, par, spec))
      names(genetic) <- unique(vapply(groups, `[[`, integer(1), "nchild"))
      A <- matrix(0, mean_dimension, mean_dimension); b <- numeric(mean_dimension)
      work <- vector("list", length(groups)); constant <- 0
      for (g in seq_along(groups)) {
        z <- groups[[g]]
        S <- .ravel_mixed_covariance(par, z, spec, genetic[[as.character(z$nchild)]])
        C <- chol(S); Q <- chol2inv(C)
        QD <- Q %*% z$D
        if (is.null(z$design_columns)) {
          A <- A + z$n * crossprod(z$D, QD)
          b <- b + drop(crossprod(QD, z$sum))
        } else {
          for (j in seq_len(mean_dimension)) {
            b[j] <- b[j] + sum(Q * z$design_outcome_cross[[j]])
            for (k in seq_len(mean_dimension))
              A[j, k] <- A[j, k] + sum(Q * z$design_cross[[j]][[k]])
          }
        }
        logdet <- 2 * sum(log(diag(C)))
        constant <- constant + z$n * (logdet + nrow(S) * log(2 * pi)) + sum(Q * z$cross)
        work[[g]] <- list(S = S, Q = Q, QD = QD, logdet = logdet)
      }
      means <- drop(solve(A, b))
      for (g in seq_along(groups)) {
        residual <- groups[[g]]$observed - .ravel_mixed_mean_matrix(groups[[g]], means)
        work[[g]]$residual <- residual
        work[[g]]$scatter <- crossprod(residual)
      }
      ans <- list(value = .5 * (constant - sum(b * means)), valid = TRUE,
                  means = means, mean_information = A, work = work)
      if (derivatives) {
        grad <- numeric(length(par)); names(grad) <- names(par)
        # Envelope theorem: derivative of the profiled means does not enter the
        # total likelihood gradient. Covariance derivatives use symmetric steps.
        for (j in seq_along(par)) {
          step <- 1e-5 * max(1, abs(par[j]))
          plus <- minus <- par; plus[j] <- par[j] + step; minus[j] <- par[j] - step
          gp <- lapply(names(genetic), function(n) .ravel_mixed_genetic(as.integer(n), plus, spec))
          gm <- lapply(names(genetic), function(n) .ravel_mixed_genetic(as.integer(n), minus, spec))
          names(gp) <- names(gm) <- names(genetic)
          for (g in seq_along(groups)) {
            z <- groups[[g]]; w <- work[[g]]
            derivative <- (.ravel_mixed_covariance(plus, z, spec, gp[[as.character(z$nchild)]]) -
              .ravel_mixed_covariance(minus, z, spec, gm[[as.character(z$nchild)]])) / (2 * step)
            grad[j] <- grad[j] + .5 * sum((z$n * w$Q - w$Q %*% w$scatter %*% w$Q) * derivative)
          }
        }
        ans$gradient <- grad
      }
      ans
    }, error = function(e) invalid(conditionMessage(e)))
    cache$par <- par; cache$value <- answer
    answer
  }
  evaluate
}

# Scores in variance units remain informative when a log-SD score tends to
# zero simply because its variance is small. Compute these analytically rather
# than dividing a finite-difference log-SD score by an almost-zero variance.
.ravel_mixed_variance_scores <- function(par, groups, spec, evaluated) {
  parameters <- intersect(c("residual_sd", "family_sd", "individual_sd", "slope_sd"), names(par))
  score <- stats::setNames(numeric(length(parameters)), parameters)
  for (g in seq_along(groups)) {
    z <- groups[[g]]; w <- evaluated$work[[g]]
    m <- length(z$child)
    if (!m) next
    yi <- 2L * (z$nchild + 2L) + seq_len(m)
    W <- (z$n * w$Q - w$Q %*% w$scatter %*% w$Q)[yi, yi, drop = FALSE] / 2
    score["residual_sd"] <- score["residual_sd"] + sum(diag(W))
    if (spec$family_random) score["family_sd"] <- score["family_sd"] + sum(W)
    same <- outer(z$child, z$child, "==")
    if (spec$individual_random) {
      derivative <- same * 1
      if (!is.null(spec$random_slope)) derivative <- derivative + same *
        outer(z$slope, z$slope, "+") * tanh(par[["intercept_slope_cor"]]) *
        exp(par[["slope_sd"]] - par[["individual_sd"]]) / 2
      score["individual_sd"] <- score["individual_sd"] + sum(W * derivative)
    }
    if (!is.null(spec$random_slope)) {
      derivative <- same * outer(z$slope, z$slope)
      if (spec$individual_random) derivative <- derivative + same *
        outer(z$slope, z$slope, "+") * tanh(par[["intercept_slope_cor"]]) *
        exp(par[["individual_sd"]] - par[["slope_sd"]]) / 2
      score["slope_sd"] <- score["slope_sd"] + sum(W * derivative)
    }
  }
  score
}

.ravel_mixed_certification <- function(par, evaluated, prepared, spec) {
  at_lower <- par - prepared$lower < 1e-6
  at_upper <- prepared$upper - par < 1e-6
  projected <- evaluated$gradient
  projected[(at_lower & projected > 0) | (at_upper & projected < 0)] <- 0
  natural <- .ravel_mixed_variance_scores(par, prepared$groups, spec, evaluated)
  random <- intersect(c("family_sd", "individual_sd", "slope_sd"), names(natural))
  near_zero <- random[exp(2 * par[random]) < 1e-6]
  natural_projected <- natural
  # A near-zero positive score points out of the variance parameter space.
  # Treat it as a boundary approximation, and never issue ordinary Wald SEs.
  natural_projected[intersect(near_zero, names(natural)[natural > 0])] <- 0
  boundary <- unique(c(near_zero, names(par)[at_lower | at_upper]))
  list(projected_gradient = projected, natural_variance_gradient = natural,
    projected_natural_variance_gradient = natural_projected, boundary_parameters = boundary,
    gradient_ok = all(is.finite(c(projected, natural_projected))) &&
      max(abs(c(projected, natural_projected))) / prepared$n_families < 1e-5)
}

.ravel_mixed_inference <- function(par, groups, spec, evaluate) {
  base <- evaluate(par, derivatives = TRUE)
  p <- length(par); q <- length(base$means)
  Hprofile <- stats::optimHess(par, function(x) evaluate(x)$value,
    function(x) evaluate(x, derivatives = TRUE)$gradient)
  Hprofile <- (Hprofile + t(Hprofile)) / 2
  Jmean <- matrix(0, q, p)
  family_scores <- matrix(0, sum(vapply(groups, `[[`, integer(1), "n")), p + q)
  clusters <- character(nrow(family_scores)); families <- character(nrow(family_scores))
  offset <- 0L
  for (g in seq_along(groups)) {
    z <- groups[[g]]; w <- base$work[[g]]
    rows <- offset + seq_len(z$n); offset <- offset + z$n
    residual <- w$residual
    if (is.null(z$design_columns)) family_scores[rows, p + seq_len(q)] <- residual %*% w$QD else {
      U <- residual %*% w$Q
      for (j in seq_len(q)) family_scores[rows, p + j] <- rowSums(U * z$design_columns[[j]])
    }
    clusters[rows] <- z$clusters; families[rows] <- z$families
  }
  for (j in seq_along(par)) {
    step <- 1e-5 * max(1, abs(par[j]))
    plus <- minus <- par; plus[j] <- par[j] + step; minus[j] <- par[j] - step
    ep <- evaluate(plus); em <- evaluate(minus)
    if (!isTRUE(ep$valid) || !isTRUE(em$valid)) stop("Cannot differentiate the fitted likelihood.")
    Jmean[, j] <- (ep$means - em$means) / (2 * step)
    offset <- 0L
    for (g in seq_along(groups)) {
      z <- groups[[g]]; w <- base$work[[g]]
      rows <- offset + seq_len(z$n); offset <- offset + z$n
      derivative <- (ep$work[[g]]$S - em$work[[g]]$S) / (2 * step)
      residual <- w$residual
      U <- residual %*% w$Q
      family_scores[rows, j] <- .5 * (rowSums((U %*% derivative) * U) - sum(w$Q * derivative))
    }
  }
  A <- base$mean_information
  H <- rbind(cbind(Hprofile + crossprod(Jmean, A %*% Jmean), -crossprod(Jmean, A)),
             cbind(-A %*% Jmean, A))
  H <- (H + t(H)) / 2
  eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= max(eigenvalues) * 1e-9)
    stop("The observed information is not positive definite; regular inference is unavailable.")
  bread <- solve(H)
  cluster_scores <- rowsum(family_scores, clusters, reorder = FALSE)
  G <- nrow(cluster_scores)
  if (G < 2L) stop("At least two FID clusters are needed for robust inference.")
  # CR0 is the unadjusted sandwich. Report its convention rather than silently
  # applying a regression degrees-of-freedom correction to a latent likelihood.
  robust <- bread %*% crossprod(cluster_scores) %*% bread
  robust <- (robust + t(robust)) / 2
  labels <- c(names(par), spec$mean_names)
  dimnames(robust) <- dimnames(bread) <- dimnames(H) <- list(labels, labels)
  colnames(family_scores) <- colnames(cluster_scores) <- labels
  list(vcov = robust, model_vcov = bread, information = H, family_scores = family_scores,
       cluster_scores = cluster_scores, family_ids = families, cluster_ids = clusters,
       mean_jacobian = Jmean, profile_information = Hprofile, method = "FID-clustered CR0 sandwich")
}

.ravel_mixed_delta <- function(fun, par) {
  value <- fun(par)
  J <- matrix(0, length(value), length(par), dimnames = list(names(value), names(par)))
  for (j in seq_along(par)) {
    step <- 1e-5 * max(1, abs(par[j])); plus <- minus <- par
    plus[j] <- par[j] + step; minus[j] <- par[j] - step
    J[, j] <- (fun(plus) - fun(minus)) / (2 * step)
  }
  list(value = value, jacobian = J)
}
