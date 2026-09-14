.ravel_binary_moments <- function(data, model) {
  variables <- lavaan::lavNames(lavaan::lavaanify(model), "ov")
  predictors <- setdiff(variables, "Y")
  d <- data[variables]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (!all(vapply(d, is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(d))) || !setequal(unique(d$Y), c(0, 1))) {
    stop("Binary covariate fitting requires finite numeric predictors and both Y = 0 and Y = 1.")
  }
  X <- as.matrix(d[predictors])
  n <- nrow(X)
  mu <- colMeans(X)
  xc <- sweep(X, 2, mu)
  S <- crossprod(xc) / n
  sx <- sqrt(diag(S))
  if (n <= ncol(X) + 1L || any(sx <= 0)) {
    stop("Binary covariate fitting requires nonconstant predictors and more rows than parameters.")
  }
  Z <- cbind(1, sweep(xc, 2, sx, "/"))
  probit <- stats::glm.fit(Z, d$Y, family = stats::binomial(link = "probit"))
  if (!isTRUE(probit$converged) || isTRUE(probit$boundary) ||
      probit$rank != ncol(Z) || any(!is.finite(probit$coefficients))) {
    stop("Conditional probit moment fit did not converge to a full-rank, nonboundary solution.")
  }
  p <- probit$fitted.values
  density <- stats::dnorm(probit$linear.predictors)
  if (any(p <= 0 | p >= 1) || any(density == 0)) {
    stop("Conditional probit moment fit has separated or numerically extreme predictions.")
  }
  bread <- crossprod(Z, Z * (density^2 / (p * (1 - p)))) / n
  score <- Z * ((d$Y - p) * density / (p * (1 - p)))
  IFb <- score %*% solve(bread)
  IFdelta <- sweep(IFb[, -1, drop = FALSE], 2, sx, "/")
  delta <- probit$coefficients[-1] / sx
  Sdelta <- drop(S %*% delta)
  projection <- drop(xc %*% delta)
  variance <- 1 + sum(delta * Sdelta)
  sd_y <- sqrt(variance)
  covariance_y <- Sdelta / sd_y
  threshold <- -probit$coefficients[1] / sd_y

  # Joint influence functions retain uncertainty in the probit and X moments.
  IFsd <- (projection^2 - sum(delta * Sdelta) + 2 * drop(IFdelta %*% Sdelta)) / (2 * sd_y)
  IFcov <- (sweep(xc * projection, 2, Sdelta) + IFdelta %*% S) / sd_y -
    outer(IFsd, covariance_y / sd_y)
  IFthreshold <- -(IFb[, 1] + projection) / sd_y -
    threshold * IFsd / sd_y

  covariance <- matrix(0, length(variables), length(variables),
                       dimnames = list(variables, variables))
  covariance[predictors, predictors] <- S
  covariance["Y", "Y"] <- 1
  covariance[predictors, "Y"] <- covariance["Y", predictors] <- covariance_y
  means <- setNames(numeric(length(variables)), variables)
  means[predictors] <- mu

  # lavaan orders means/thresholds, numeric variances, then lower-triangle covariances.
  location_names <- ifelse(variables == "Y", "Y|t1", paste0(variables, "~1"))
  pairs <- which(lower.tri(covariance), arr.ind = TRUE)
  labels <- c(location_names, paste0(predictors, "~~", predictors),
              paste0(variables[pairs[, 2]], "~~", variables[pairs[, 1]]))
  IF <- matrix(0, n, length(labels), dimnames = list(NULL, labels))
  working <- numeric(length(labels))
  for (j in seq_along(labels)) {
    label <- labels[j]
    if (label == "Y|t1") {
      IF[, j] <- IFthreshold
      working[j] <- stats::pnorm(threshold) * stats::pnorm(-threshold) /
        stats::dnorm(threshold)^2
    } else if (grepl("~~", label, fixed = TRUE)) {
      pair <- strsplit(label, "~~", fixed = TRUE)[[1]]
      if ("Y" %in% pair) {
        IF[, j] <- IFcov[, match(setdiff(pair, "Y"), predictors)]
      } else {
        IF[, j] <- xc[, pair[1]] * xc[, pair[2]] - S[pair[1], pair[2]]
      }
      working[j] <- covariance[pair[1], pair[1]] * covariance[pair[2], pair[2]] +
        covariance[pair[1], pair[2]]^2
    } else {
      name <- sub("~1$", "", label)
      IF[, j] <- xc[, name]
      working[j] <- S[name, name]
    }
  }
  if (any(!is.finite(IF)) || any(!is.finite(working) | working <= 0)) {
    stop("Conditional probit moment covariance is not finite.")
  }
  IF <- sweep(IF, 2, colMeans(IF))
  th <- -means
  th["Y"] <- threshold
  attr(th, "th.idx") <- setNames(ifelse(variables == "Y", match("Y", variables), 0L),
                                location_names)
  # Positive working weights; robust SEs use the full, possibly singular NACOV.
  list(sample.cov = covariance, sample.mean = means, sample.th = th,
       sample.nobs = n, NACOV = crossprod(IF) / n,
       WLS.V = diag(1 / working), sample.cov.rescale = FALSE)
}
