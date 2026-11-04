# Conditional family likelihood for genetically associated covariates.
# This extension leaves the original independent-covariate evaluator unchanged.

# Store one compact fixed-effect block per observed response position. The old
# representation allocated a d by d matrix for every pair of mean coefficients;
# with continuously varying designs that scales as families * d^2 * q^2.
.ravel_mixed_covariate_compress <- function(prepared) {
  groups <- prepared$groups
  keys <- vapply(groups, function(z) paste(z$nchild, paste(z$child, collapse = ","),
    paste(format(c(z$moderators, z$slope), digits = 17, scientific = TRUE,
                 trim = TRUE), collapse = ","), sep = "|"), character(1))
  memberships <- split(seq_along(groups), factor(keys, levels = unique(keys)))
  prepared$groups <- unname(lapply(memberships, function(indices) {
    chunks <- groups[indices]
    z <- chunks[[1L]]
    z$n <- as.integer(sum(vapply(chunks, `[[`, integer(1), "n")))
    z$observed <- do.call(rbind, lapply(chunks, `[[`, "observed"))
    z$families <- unlist(lapply(chunks, `[[`, "families"), use.names = FALSE)
    z$clusters <- unlist(lapply(chunks, `[[`, "clusters"), use.names = FALSE)
    if (!is.null(z$C)) z$C <- do.call(rbind, lapply(chunks, `[[`, "C"))
    z$sum <- colSums(z$observed)
    z$cross <- crossprod(z$observed)
    z$mean_blocks <- lapply(seq_len(nrow(z$D)), function(row) {
      active <- which(Reduce(`|`, lapply(chunks, function(x) x$D[row, ] != 0)))
      values <- if (length(active)) do.call(rbind, lapply(chunks, function(x)
        matrix(rep(x$D[row, active], each = x$n), x$n, length(active)))) else
          matrix(numeric(), z$n, 0L)
      list(columns = active, values = values)
    })
    # Crossproducts are useful only if they are smaller than the raw designs.
    # In particular, never materialize q^2 products for singleton patterns.
    widths <- vapply(z$mean_blocks, function(x) length(x$columns), integer(1))
    if (sum(widths)^2 <= z$n * sum(widths)) {
      z$mean_cross <- lapply(z$mean_blocks, function(a)
        lapply(z$mean_blocks, function(b) crossprod(a$values, b$values)))
    } else z$mean_cross <- NULL
    if (z$n == 1L) z$mean_blocks <- NULL
    z$design_columns <- z$design_cross <- z$design_outcome_cross <- NULL
    z
  }))
  prepared
}

.ravel_mixed_covariate_mean <- function(group, means) {
  if (is.null(group$mean_blocks))
    return(matrix(rep(drop(group$D %*% means), each = group$n), group$n))
  out <- matrix(0, group$n, length(group$mean_blocks))
  for (i in seq_along(group$mean_blocks)) {
    b <- group$mean_blocks[[i]]
    if (length(b$columns)) out[, i] <- drop(b$values %*% means[b$columns])
  }
  out
}

.ravel_mixed_covariate_mean_information <- function(group, precision, q) {
  if (is.null(group$mean_blocks)) return(group$n * crossprod(group$D, precision %*% group$D))
  out <- matrix(0, q, q)
  blocks <- group$mean_blocks
  for (i in seq_along(blocks)) {
    a <- blocks[[i]]
    if (!length(a$columns)) next
    for (j in seq_along(blocks)) {
      b <- blocks[[j]]
      if (!length(b$columns) || precision[i, j] == 0) next
      cp <- if (is.null(group$mean_cross)) crossprod(a$values, b$values) else group$mean_cross[[i]][[j]]
      out[a$columns, b$columns] <- out[a$columns, b$columns, drop = FALSE] + precision[i, j] * cp
    }
  }
  out
}

.ravel_mixed_covariate_mean_scores <- function(group, precision_residual, q) {
  if (is.null(group$mean_blocks)) return(precision_residual %*% group$D)
  out <- matrix(0, group$n, q)
  for (i in seq_along(group$mean_blocks)) {
    b <- group$mean_blocks[[i]]
    if (length(b$columns)) out[, b$columns] <- out[, b$columns, drop = FALSE] +
      b$values * precision_residual[, i]
  }
  out
}

.ravel_mixed_covariate_kernel <- function(joint) {
  if (!is.null(joint$H)) return(joint$H)
  rbind(joint$blocks$L, joint$B %*% joint$blocks$R) %*%
    joint$blocks$D %*% joint$blocks$Qinv
}

# The first-stage moments are packed as means and raw second moments. Keep their
# exact packing in the covariate module; differentiating this small deterministic
# transformation does not require refitting either stage.
.ravel_mixed_covariate_eta_jacobian <- function(covariate) {
  if (exists(".ravel_mixed_covariate_unpack_jacobian", mode = "function"))
    return(.ravel_mixed_covariate_unpack_jacobian(covariate))
  physical <- function(eta) {
    x <- .ravel_mixed_covariate_unpack(covariate, eta)
    c(x$mean, as.vector(x$Qself), as.vector(x$Qsib))
  }
  eta <- covariate$eta
  J <- matrix(0, length(physical(eta)), length(eta))
  for (j in seq_along(eta)) {
    h <- 1e-5 * max(1, abs(eta[j])); a <- b <- eta
    a[j] <- a[j] + h; b[j] <- b[j] - h
    J[, j] <- (physical(a) - physical(b)) / (2 * h)
  }
  colnames(J) <- names(eta)
  J
}

.ravel_mixed_covariate_evaluator <- function(prepared, spec) {
  groups <- prepared$groups
  q <- length(spec$mean_names)
  counts <- unique(vapply(groups, `[[`, integer(1), "nchild"))
  cache <- new.env(parent = emptyenv())
  eta_jacobian <- NULL
  blocks_at <- function(par) {
    out <- lapply(counts, function(n) .ravel_mixed_covariate_blocks(par, spec, n))
    names(out) <- as.character(counts)
    out
  }
  joint_at <- function(par, group, blocks, include_offset = TRUE) {
    ans <- .ravel_mixed_covariate_joint(par, group, spec,
      blocks = blocks[[as.character(group$nchild)]], include_offset = include_offset)
    ans$kernel <- .ravel_mixed_covariate_kernel(ans)
    ans
  }
  evaluate <- function(par, derivatives = FALSE, keep = FALSE, means = NULL,
                       moments = FALSE) {
    if (exists("par", cache, inherits = FALSE) && identical(par, cache$par) &&
        identical(means, cache$means) &&
        (!derivatives || !is.null(cache$value$gradient)) &&
        (!moments || !is.null(cache$value$moment_gradient))) return(cache$value)
    invalid <- function(reason) list(value = 1e50, valid = FALSE, reason = reason,
      gradient = if (derivatives) rep(0, length(par)) else NULL)
    answer <- tryCatch({
      blocks <- blocks_at(par)
      A <- matrix(0, q, q); b <- numeric(q); constant <- 0
      work <- vector("list", length(groups))
      for (g in seq_along(groups)) {
        z <- groups[[g]]; j <- joint_at(par, z, blocks)
        ch <- chol(j$S); Q <- chol2inv(ch)
        centered <- z$C - rep(j$blocks$mean, z$nchild)[col(z$C)]
        adjusted <- z$observed - j$offset
        logdet <- 2 * sum(log(diag(ch)))
        A <- A + .ravel_mixed_covariate_mean_information(z, Q, q)
        b <- b + colSums(.ravel_mixed_covariate_mean_scores(z, adjusted %*% Q, q))
        constant <- constant + z$n * (logdet + nrow(Q) * log(2 * pi)) +
          sum((adjusted %*% Q) * adjusted)
        work[[g]] <- list(S = j$S, Q = Q, QD = Q %*% z$D,
          logdet = logdet, offset = j$offset, kernel = j$kernel,
          centered_C = centered, blocks = j$blocks, B = j$B, adjusted = adjusted)
      }
      fitted_means <- if (is.null(means)) drop(solve(A, b)) else means
      if (length(fitted_means) != q || any(!is.finite(fitted_means)))
        stop("The conditional mean coefficients are invalid.")
      for (g in seq_along(groups)) {
        w <- work[[g]]
        w$residual <- w$adjusted - .ravel_mixed_covariate_mean(groups[[g]], fitted_means)
        w$scatter <- crossprod(w$residual)
        w$precision_residual <- w$residual %*% w$Q
        w$mean_cross <- crossprod(w$precision_residual, w$centered_C)
        w$covariance_gradient <- .5 * (groups[[g]]$n * w$Q - w$Q %*% w$scatter %*% w$Q)
        w$adjusted <- NULL
        work[[g]] <- w
      }
      ans <- list(value = .5 * (constant - 2 * sum(b * fitted_means) +
        sum(fitted_means * drop(A %*% fitted_means))), valid = TRUE,
        means = fitted_means, mean_information = A, work = work)
      if (derivatives) {
        gradient <- setNames(numeric(length(par)), names(par))
        for (k in seq_along(par)) {
          h <- 1e-5 * max(1, abs(par[k])); plus <- minus <- par
          plus[k] <- plus[k] + h; minus[k] <- minus[k] - h
          bp <- blocks_at(plus); bm <- blocks_at(minus)
          for (g in seq_along(groups)) {
            jp <- joint_at(plus, groups[[g]], bp, FALSE)
            jm <- joint_at(minus, groups[[g]], bm, FALSE)
            dS <- (jp$S - jm$S) / (2 * h)
            dM <- (jp$kernel - jm$kernel) / (2 * h)
            gradient[k] <- gradient[k] + sum(work[[g]]$covariance_gradient * dS) -
              sum(work[[g]]$mean_cross * dM)
          }
        }
        ans$gradient <- gradient
      }
      if (moments) {
        pC <- spec[["covariate"]]$q
        gmean <- numeric(pC); gself <- gsib <- matrix(0, pC, pC)
        sym <- function(x) (x + t(x)) / 2
        for (g in seq_along(groups)) {
          z <- groups[[g]]; w <- work[[g]]
          J <- w$kernel
          ZU <- w$blocks$Qinv %*% t(w$mean_cross)
          Js <- matrix(0, nrow(J), pC); ZUs <- matrix(0, pC, nrow(J))
          diagonal <- matrix(0, pC, pC)
          for (i in seq_len(z$nchild)) {
            at <- (i - 1L) * pC + seq_len(pC)
            Ji <- J[, at, drop = FALSE]; ZUi <- ZU[at, , drop = FALSE]
            diagonal <- diagonal + crossprod(Ji, w$covariance_gradient %*% Ji) + sym(ZUi %*% Ji)
            Js <- Js + Ji; ZUs <- ZUs + ZUi
          }
          gself <- gself + diagonal
          gsib <- gsib + crossprod(Js, w$covariance_gradient %*% Js) + sym(ZUs %*% Js) - diagonal
          gmean <- gmean + drop(crossprod(Js, colSums(w$precision_residual)))
        }
        if (is.null(eta_jacobian)) eta_jacobian <<- .ravel_mixed_covariate_eta_jacobian(spec[["covariate"]])
        physical <- c(gmean, as.vector(gself), as.vector(gsib))
        ans$moment_gradient <- drop(crossprod(eta_jacobian, physical))
        names(ans$moment_gradient) <- names(spec[["covariate"]]$eta)
      }
      ans
    }, error = function(e) invalid(conditionMessage(e)))
    cache$par <- par; cache$means <- means; cache$value <- answer
    answer
  }
  evaluate
}

# Combine likelihood scores and the first-stage empirical covariate moments.
# The moment module supplies h_f such that eta_hat - eta = sum(h_f)/Nfamily.
# This retains uncertainty and correlation from estimating both stages on the
# same families before aggregation to the user-supplied extended-family FID.
.ravel_mixed_covariate_sandwich <- function(par, prepared, spec, information,
    family_scores, clusters, score_total_function = NULL,
    nuisance_gradient_function = NULL, family_ids = NULL) {
  if (is.null(family_ids)) family_ids <- unlist(lapply(prepared$groups, `[[`, "families"), use.names = FALSE)
  H <- (information + t(information)) / 2
  eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= max(eigenvalues) * 1e-9)
    stop("The observed information is not positive definite; regular inference is unavailable.")
  bread <- solve(H)
  cv <- spec[["covariate"]]
  at <- match(family_ids, cv$family_ids)
  if (anyNA(at) || anyDuplicated(family_ids) || length(at) != nrow(family_scores))
    stop("The first-stage covariate moments do not match the likelihood families.")
  h <- cv$moment_scores[at, , drop = FALSE]
  n <- length(at)
  eta <- cv$eta
  cross_score <- matrix(0, length(par), length(eta),
    dimnames = list(names(par), names(eta)))
  if (!is.null(nuisance_gradient_function)) {
    # By equality of mixed partials, differentiate the analytic negative
    # likelihood moment gradient in each likelihood parameter. This avoids
    # hundreds of full score-vector evaluations when 20+ PCs are supplied.
    for (j in seq_along(par)) {
      step <- 1e-5 * max(1, abs(par[j])); plus <- minus <- par
      plus[j] <- plus[j] + step; minus[j] <- minus[j] - step
      a <- nuisance_gradient_function(plus, spec)
      b <- nuisance_gradient_function(minus, spec)
      if (length(a) != length(eta) || length(b) != length(eta) || any(!is.finite(c(a, b))))
        stop("Cannot differentiate the first-stage moment adjustment.")
      cross_score[j, ] <- -(a - b) / (2 * step)
    }
  } else {
    if (is.null(score_total_function)) stop("A likelihood score function is required for the moment adjustment.")
    for (j in seq_along(eta)) {
      step <- 1e-5 * max(1, abs(eta[j])); plus <- minus <- spec
      plus[["covariate"]]$eta[j] <- eta[j] + step
      minus[["covariate"]]$eta[j] <- eta[j] - step
      a <- score_total_function(par, plus)
      b <- score_total_function(par, minus)
      if (length(a) != length(par) || length(b) != length(par) || any(!is.finite(c(a, b))))
        stop("Cannot differentiate the first-stage moment adjustment.")
      cross_score[, j] <- (a - b) / (2 * step)
    }
  }
  effective_scores <- family_scores + h %*% t(cross_score / n)
  if (anyDuplicated(clusters)) {
    cluster_scores <- rowsum(effective_scores, clusters, reorder = FALSE)
    cluster_moments <- rowsum(h, clusters, reorder = FALSE)
  } else {
    cluster_scores <- effective_scores
    cluster_moments <- h
  }
  if (nrow(cluster_scores) < 2L) stop("At least two FID clusters are needed for robust inference.")
  robust <- bread %*% crossprod(cluster_scores) %*% bread
  nuisance_vcov <- crossprod(cluster_moments) / n^2
  param_nuisance_cov <- bread %*% crossprod(cluster_scores, cluster_moments) / n
  labels <- names(par)
  dimnames(robust) <- dimnames(bread) <- dimnames(H) <- list(labels, labels)
  dimnames(nuisance_vcov) <- list(names(eta), names(eta))
  dimnames(param_nuisance_cov) <- list(labels, names(eta))
  colnames(family_scores) <- colnames(effective_scores) <- colnames(cluster_scores) <- labels
  list(vcov = (robust + t(robust)) / 2, model_vcov = bread, information = H,
    family_scores = family_scores, cluster_scores = cluster_scores,
    family_ids = family_ids, cluster_ids = clusters,
    nuisance_vcov = nuisance_vcov, param_nuisance_cov = param_nuisance_cov,
    nuisance_parameters = eta, nuisance_score_sensitivity = cross_score,
    first_stage_uncertainty = "Estimated covariate means and coherent self/sibling covariance moments included through stacked FID influence",
    method = "FID-clustered CR0 sandwich with estimated covariate moments")
}

.ravel_mixed_covariate_inference <- function(par, prepared, spec, evaluate) {
  groups <- prepared$groups
  base <- evaluate(par, derivatives = TRUE)
  p <- length(par); q <- length(base$means)
  Hprofile <- stats::optimHess(par, function(x) evaluate(x)$value,
    function(x) evaluate(x, derivatives = TRUE)$gradient)
  Hprofile <- (Hprofile + t(Hprofile)) / 2
  Jmean <- matrix(0, q, p)
  family_scores <- matrix(0, sum(vapply(groups, `[[`, integer(1), "n")), p + q)
  clusters <- families <- character(nrow(family_scores))
  offset <- 0L
  for (g in seq_along(groups)) {
    z <- groups[[g]]; w <- base$work[[g]]
    rows <- offset + seq_len(z$n); offset <- offset + z$n
    family_scores[rows, p + seq_len(q)] <-
      .ravel_mixed_covariate_mean_scores(z, w$precision_residual, q)
    clusters[rows] <- z$clusters; families[rows] <- z$families
  }
  for (j in seq_along(par)) {
    step <- 1e-5 * max(1, abs(par[j])); plus <- minus <- par
    plus[j] <- plus[j] + step; minus[j] <- minus[j] - step
    ep <- evaluate(plus); em <- evaluate(minus)
    if (!isTRUE(ep$valid) || !isTRUE(em$valid)) stop("Cannot differentiate the fitted conditional likelihood.")
    Jmean[, j] <- (ep$means - em$means) / (2 * step)
    offset <- 0L
    for (g in seq_along(groups)) {
      z <- groups[[g]]; w <- base$work[[g]]
      rows <- offset + seq_len(z$n); offset <- offset + z$n
      dS <- (ep$work[[g]]$S - em$work[[g]]$S) / (2 * step)
      dM <- (ep$work[[g]]$kernel - em$work[[g]]$kernel) / (2 * step)
      U <- w$precision_residual
      doffset <- w$centered_C %*% t(dM)
      family_scores[rows, j] <- .5 * (rowSums((U %*% dS) * U) - sum(w$Q * dS)) +
        rowSums(U * doffset)
    }
  }
  A <- base$mean_information
  H <- rbind(cbind(Hprofile + crossprod(Jmean, A %*% Jmean), -crossprod(Jmean, A)),
             cbind(-A %*% Jmean, A))
  parameters <- c(par, stats::setNames(base$means, spec$mean_names))
  # Both functions hold the nuisance mean coefficients fixed when varying the
  # other parameters. Reprofiling here would give the wrong joint cross block.
  total_score <- function(x, s) {
    ev <- .ravel_mixed_covariate_evaluator(prepared, s)
    value <- ev(x[seq_len(p)], derivatives = TRUE, means = x[p + seq_len(q)])
    if (!isTRUE(value$valid)) stop(value$reason)
    mean_score <- Reduce(`+`, lapply(seq_along(groups), function(g)
      colSums(.ravel_mixed_covariate_mean_scores(groups[[g]], value$work[[g]]$precision_residual, q))))
    c(-value$gradient, mean_score)
  }
  moment_gradient <- function(x, s) {
    # The same spec is used in this derivative route, so reuse the evaluator's
    # cached moment-coordinate Jacobian instead of constructing it per call.
    value <- evaluate(x[seq_len(p)], means = x[p + seq_len(q)], moments = TRUE)
    if (!isTRUE(value$valid)) stop(value$reason)
    value$moment_gradient
  }
  ans <- .ravel_mixed_covariate_sandwich(parameters, prepared, spec, H,
    family_scores, clusters, total_score, moment_gradient, families)
  ans$mean_jacobian <- Jmean
  ans$profile_information <- Hprofile
  ans
}
