# Coherent exchangeable covariate moments. A Gaussian working likelihood is
# used only for first-stage covariance estimation; family influence scores,
# rather than Gaussian model SEs, enter the final FID sandwich.
.ravel_mixed_covariate_moments <- function(family_values, weights = NULL,
    influence = TRUE) {
  sizes <- vapply(family_values, nrow, integer(1))
  F <- length(sizes); q <- ncol(family_values[[1L]])
  if (is.null(weights)) weights <- rep(1, F)
  if (length(weights) != F || any(!is.finite(weights)) || any(weights <= 0))
    stop("Covariate family weights must be finite and positive.")
  mass <- sum(weights); N <- sum(weights * sizes); maxn <- max(sizes)
  triangle <- which(lower.tri(matrix(0, q, q), diag = TRUE), arr.ind = TRUE)
  k <- nrow(triangle)
  index <- list(mean = seq_len(q), self = q + seq_len(k),
    sibling = if (maxn > 1L) q + k + seq_len(k) else integer())
  averages <- matrix(vapply(family_values, colMeans, numeric(q)), F, q, byrow = TRUE)
  mu <- colSums(averages * (weights * sizes)) / N
  centered <- sweep(averages, 2L, mu, "-")
  within <- lapply(family_values, function(x) crossprod(sweep(x, 2L, colMeans(x), "-")))
  Sw <- Reduce(`+`, Map(function(x, w) x * w, within, weights))
  self <- (Sw + crossprod(centered * sqrt(weights * sizes))) / N
  shared <- if (maxn > 1L) which(diag(Sw) <= 1e-20 * N * pmax(1, diag(self))) else integer()
  varying <- setdiff(seq_len(q), shared); v <- length(varying)
  rank_check <- function(x, label) {
    if (!length(x)) return(invisible(NULL))
    e <- eigen((x + t(x))/2, symmetric = TRUE, only.values = TRUE)$values
    if (min(e) <= 1e-10 * max(1, max(e)))
      stop("Covariate ", label, " lacks independent variation. Check redundant covariates or too few informative sibling families.")
  }
  rank_check(self, "person-level covariance")
  mean_influence <- sweep(averages, 2L, mu, "-") * (sizes / (N / mass))
  # No sibling covariance is estimated without observed siblings.
  if (maxn == 1L) {
    Qself <- self; Qsib <- matrix(0, q, q)
    eta <- c(mu, (Qself + tcrossprod(mu))[triangle])
    h <- if (influence) cbind(mean_influence, matrix(vapply(family_values, function(x)
      (crossprod(x) - (Qself + tcrossprod(mu)))[triangle], numeric(k)), F, k, byrow = TRUE)) else NULL
    diagnostics <- list(estimator = "Unique-child empirical covariance", converged = TRUE)
  } else {
    W <- matrix(0, q, q)
    if (v) {
      W[varying, varying] <- Sw[varying, varying, drop = FALSE] / sum(weights * (sizes - 1L))
      rank_check(W[varying, varying, drop = FALSE], "within-sibling contrast covariance")
    }
    # SPD W on its observed support and SPD M at the largest family size imply
    # valid covariance at every smaller family size, including singletons.
    iw <- which(lower.tri(matrix(0, v, v), diag = TRUE), arr.ind = TRUE)
    nw <- nrow(iw)
    pack <- function(x, ix) {
      if (!nrow(ix)) return(numeric())
      L <- t(chol(x)); z <- L[ix]; d <- ix[,1L] == ix[,2L]
      z[d] <- log(z[d]); z
    }
    unpack <- function(z, ix, dim) {
      L <- matrix(0, dim, dim)
      if (nrow(ix)) {
        d <- ix[,1L] == ix[,2L]; z[d] <- exp(z[d]); L[ix] <- z
      }
      list(L = L, S = tcrossprod(L))
    }
    decode <- function(phi) {
      a <- unpack(phi[seq_len(nw)], iw, v)
      b <- unpack(phi[nw + seq_len(k)], triangle, q)
      A <- matrix(0, q, q)
      if (v) A[varying, varying] <- a$S
      list(W = A, M = b$S, Lw = a$L, Lm = b$L)
    }
    chain <- function(G, L, ix) {
      if (!nrow(ix)) return(numeric())
      out <- (2 * G %*% L)[ix]; d <- ix[,1L] == ix[,2L]
      out[d] <- out[d] * diag(L); out
    }
    groups <- split(seq_len(F), sizes)
    summaries <- lapply(groups, function(at) list(n = sizes[at[1L]],
      count = sum(weights[at]), at = at,
      within = Reduce(`+`, Map(function(x, w) w*x, within[at], weights[at]))))
    evaluate <- function(phi, mean = mu, per_family = FALSE) {
      s <- decode(phi)
      GW <- matrix(0, q, q); GM <- GW; value <- 0
      wi <- matrix(0, q, q)
      if (v) wi[varying, varying] <- chol2inv(chol(s$W[varying, varying, drop=FALSE]))
      ldW <- if (v) 2 * sum(log(diag(s$Lw))) else 0
      grad_family <- if (per_family) matrix(0, F, nw + k) else NULL
      for (g in summaries) {
        a <- g$n/maxn; Mn <- (1-a)*s$W + a*s$M
        ch <- chol(Mn); mi <- chol2inv(ch)
        d <- sweep(averages[g$at,,drop=FALSE], 2L, mean, "-")
        Tn <- crossprod(d * sqrt(weights[g$at] * g$n))
        Hn <- .5 * (g$count*mi - mi %*% Tn %*% mi)
        HW <- .5 * (g$count*(g$n-1L)*wi - wi %*% g$within %*% wi)
        GW <- GW + HW + (1-a)*Hn; GM <- GM + a*Hn
        value <- value + .5 * (g$count*((g$n-1L)*ldW + 2*sum(log(diag(ch)))) +
          sum(wi*g$within) + sum(mi*Tn))
        if (per_family) for (j in seq_along(g$at)) {
          i <- g$at[j]; z <- drop(mi %*% d[j,])
          hi <- .5 * (mi - g$n*tcrossprod(z))
          hw <- .5 * ((g$n-1L)*wi - wi %*% within[[i]] %*% wi) + (1-a)*hi
          grad_family[i,] <- c(chain(hw[varying,varying,drop=FALSE],s$Lw,iw),
            chain(a*hi,s$Lm,triangle))
        }
      }
      gradient <- c(chain(GW[varying,varying,drop=FALSE],s$Lw,iw), chain(GM,s$Lm,triangle))/mass
      list(value = value/mass, gradient = gradient, family_gradient = grad_family)
    }
    if (length(groups) == 1L) {
      # Balanced designs have the exact unconstrained likelihood solution.
      M <- crossprod(centered * sqrt(weights * maxn))/mass
      rank_check(M, "family-mean covariance")
      phi <- c(pack(W[varying,varying,drop=FALSE], iw), pack(M, triangle))
      status <- 0L
    } else {
      start <- c(pack(W[varying,varying,drop=FALSE], iw), pack(maxn*self, triangle))
      fit <- stats::nlminb(start, function(x) evaluate(x)$value,
        function(x) evaluate(x)$gradient,
        control = list(iter.max = 1000L, eval.max = 1500L, rel.tol = 1e-12, x.tol = 1e-10))
      phi <- fit$par; status <- fit$convergence
      # A small Newton polish makes influence comparisons reproducible even
      # when the optimizer stops on relative objective change.
      for (iteration in seq_len(3L)) {
        ev <- evaluate(phi)
        if (max(abs(ev$gradient)) < 1e-10) break
        H <- stats::optimHess(phi, function(x) evaluate(x)$value, function(x) evaluate(x)$gradient,
          control = list(ndeps = rep(1e-5, length(phi))))
        step <- tryCatch(solve(H, ev$gradient), error = function(e) NULL)
        if (is.null(step)) break
        accepted <- FALSE
        for (fraction in 2^-(0:12)) {
          candidate <- phi - fraction*step
          ce <- tryCatch(evaluate(candidate), error = function(e) NULL)
          if (!is.null(ce) && ce$value <= ev$value + 1e-13) {
            phi <- candidate; accepted <- TRUE; break
          }
        }
        if (!accepted) break
      }
    }
    ev <- evaluate(phi, per_family = influence)
    if (any(!is.finite(ev$gradient)) || max(abs(ev$gradient)) > 1e-6)
      stop("The coherent covariate covariance fit did not converge; regular inference is unavailable.")
    s <- decode(phi)
    Qself <- (s$M + (maxn-1L)*s$W)/maxn
    Qsib <- (s$M - s$W)/maxn
    eta <- c(mu, (Qself+tcrossprod(mu))[triangle], (Qsib+tcrossprod(mu))[triangle])
    h <- NULL
    if (influence) {
      H <- stats::optimHess(phi, function(x) evaluate(x)$value, function(x) evaluate(x)$gradient,
          control = list(ndeps = rep(1e-5, length(phi))))
      H <- (H+t(H))/2
      he <- eigen(H, symmetric=TRUE, only.values=TRUE)$values
      if (min(he) <= 1e-9*max(he))
        stop("The coherent covariate covariance fit has singular information; more independent sibling variation is needed.")
      K <- matrix(0, length(phi), q)
      for (j in seq_len(q)) {
        plus <- minus <- mu; step <- 1e-5*max(1,abs(mu[j]))
        plus[j] <- plus[j]+step; minus[j] <- minus[j]-step
        K[,j] <- (evaluate(phi,plus)$gradient-evaluate(phi,minus)$gradient)/(2*step)
      }
      hp <- -(ev$family_gradient + mean_influence %*% t(K)) %*% solve(H)
      # Exact Cholesky-to-covariance Jacobian; no refitting in the SEM sandwich.
      J <- matrix(0, 2*k, length(phi))
      for (j in seq_along(phi)) {
        dw <- dm <- matrix(0,q,q)
        if (j <= nw) {
          dL <- matrix(0,v,v); row <- iw[j,1L]; col <- iw[j,2L]
          dL[row,col] <- if (row == col) s$Lw[row,col] else 1
          dw[varying,varying] <- tcrossprod(dL,s$Lw)+tcrossprod(s$Lw,dL)
        } else {
          dL <- matrix(0,q,q); row <- triangle[j-nw,1L]; col <- triangle[j-nw,2L]
          dL[row,col] <- if (row == col) s$Lm[row,col] else 1
          dm <- tcrossprod(dL,s$Lm)+tcrossprod(s$Lm,dL)
        }
        J[,j] <- c(((dm+(maxn-1L)*dw)/maxn)[triangle], ((dm-dw)/maxn)[triangle])
      }
      raw_mean <- matrix(vapply(seq_len(F), function(i)
        (tcrossprod(mu,mean_influence[i,])+tcrossprod(mean_influence[i,],mu))[triangle], numeric(k)), F, k, byrow = TRUE)
      h <- cbind(mean_influence, hp %*% t(J) + cbind(raw_mean,raw_mean))
      h <- sweep(h, 2L, colSums(h*weights)/mass, "-")
    }
    diagnostics <- list(estimator = "Coherent exchangeable Gaussian working covariance likelihood",
      converged = TRUE, optimizer_status = status, max_abs_gradient = max(abs(ev$gradient)),
      max_family_size = maxn, shared_columns = shared,
      within_min_eigenvalue = if (v) min(eigen(s$W[varying,varying,drop=FALSE],symmetric=TRUE)$values) else NA_real_,
      mean_min_eigenvalue = min(eigen(s$M,symmetric=TRUE)$values))
  }
  names(eta) <- c(paste0("C_mean_", seq_len(q)),
    paste0("C_self_",triangle[,1L],"_",triangle[,2L]),
    if (maxn > 1L) paste0("C_sibling_",triangle[,1L],"_",triangle[,2L]))
  if (!is.null(h)) dimnames(h) <- list(names(family_values), names(eta))
  list(mean = mu, Qself = Qself, Qsib = Qsib, shared = shared,
    eta = eta, triangle = triangle, index = index, moment_scores = h,
    diagnostics = diagnostics)
}
