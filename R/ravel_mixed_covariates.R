# Joint child-covariate extension for the optional mixed likelihood.
# All score measurements remain unique-person measurements. Covariates belong
# to children; conditioning on the whole sibling vector also informs parents.

.ravel_mixed_covariate_unpack <- function(covariate, eta = covariate$eta) {
  if (!is.numeric(eta) || length(eta) != length(covariate$eta) || any(!is.finite(eta)))
    stop("Covariate moments must be a finite vector matching the fitted moment specification.")
  q <- covariate$q; ix <- covariate$triangle
  symmetric <- function(value) {
    out <- matrix(0, q, q); out[ix] <- value
    out[cbind(ix[, 2L], ix[, 1L])] <- value
    out
  }
  mean <- unname(eta[covariate$index$mean])
  self <- symmetric(eta[covariate$index$self]) - tcrossprod(mean)
  sibling <- if (covariate$has_siblings)
    symmetric(eta[covariate$index$sibling]) - tcrossprod(mean) else matrix(0, q, q)
  list(mean = mean, Qself = self, Qsib = sibling, shared = covariate$shared)
}

# Derivative rows are mean, column-major vec(Qself), column-major vec(Qsib).
# The sibling rows remain present, with zero derivatives, when no siblings
# identify a sibling-covariance block.
.ravel_mixed_covariate_unpack_jacobian <- function(covariate, eta = covariate$eta) {
  value <- .ravel_mixed_covariate_unpack(covariate, eta)
  q <- covariate$q; ix <- covariate$triangle
  J <- matrix(0, q + 2L * q^2L, length(eta),
    dimnames = list(NULL, names(covariate$eta)))
  J[cbind(seq_len(q), covariate$index$mean)] <- 1
  for (j in seq_len(q)) {
    unit <- numeric(q); unit[j] <- 1
    derivative <- -(tcrossprod(unit, value$mean) + tcrossprod(value$mean, unit))
    J[q + seq_len(q^2L), covariate$index$mean[j]] <- c(derivative)
    if (covariate$has_siblings)
      J[q + q^2L + seq_len(q^2L), covariate$index$mean[j]] <- c(derivative)
  }
  for (j in seq_len(nrow(ix))) {
    derivative <- matrix(0, q, q)
    derivative[ix[j, 1L], ix[j, 2L]] <- derivative[ix[j, 2L], ix[j, 1L]] <- 1
    J[q + seq_len(q^2L), covariate$index$self[j]] <- c(derivative)
    if (covariate$has_siblings)
      J[q + q^2L + seq_len(q^2L), covariate$index$sibling[j]] <- c(derivative)
  }
  J
}

.ravel_mixed_covariate_inverse <- function(moments, nchild) {
  q <- length(moments$mean)
  inverse <- function(x, label) {
    x <- (x + t(x)) / 2
    e <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
    if (!all(is.finite(e)) || min(e) <= 1e-10 * max(1, max(e)))
      stop("The estimated covariate ", label, " covariance is singular or not positive definite. ",
        "Inspect covariate variation and the first-stage covariance fit; no nearest-positive-definite replacement is made.")
    chol2inv(chol(x))
  }
  mean_inverse <- inverse(moments$Qself + (nchild - 1L) * moments$Qsib,
    paste0("family-mean block for ", nchild, " children"))
  within_inverse <- mean_inverse
  if (nchild > 1L) {
    # Exactly shared columns have no sibling-contrast coordinates. Preserve
    # their family-mean information rather than inventing contrast variance.
    varying <- setdiff(seq_len(q), moments$shared)
    within_inverse <- matrix(0, q, q)
    if (length(varying)) within_inverse[varying, varying] <- inverse(
      (moments$Qself - moments$Qsib)[varying, varying, drop = FALSE],
      "within-sibling block")
  }
  Qinv <- kronecker(diag(nchild), within_inverse) +
    kronecker(matrix(1 / nchild, nchild, nchild), mean_inverse - within_inverse)
  list(Qinv = Qinv, Qwithin_inv = within_inverse, Qmean_inv = mean_inverse)
}

.ravel_mixed_covariate_prepare <- function(prepared) {
  columns <- prepared$spec$covariates
  if (!length(columns)) { prepared$spec$covariate <- NULL; return(prepared) }
  dat <- prepared$data
  iid <- as.character(dat[[prepared$subject_id]])
  fam <- as.character(dat[[prepared$family_id]])
  ids <- unique(iid)
  person_rows <- split(seq_len(nrow(dat)), factor(iid, levels = ids))
  first <- vapply(person_rows, `[`, integer(1), 1L)
  allowed_mean <- prepared$spec$person_mean_covariates
  if (is.null(allowed_mean)) allowed_mean <- intersect(c("age", prepared$spec$random_slope), columns)
  if (!is.character(allowed_mean) || anyNA(allowed_mean) || anyDuplicated(allowed_mean) ||
      any(!allowed_mean %in% columns))
    stop("person_mean_covariates must contain distinct supplied covariate column names.")
  C <- matrix(NA_real_, length(ids), length(columns), dimnames = list(ids, columns))
  features <- data.frame(column = columns, aggregation = "constant_within_person",
    center = NA_real_, scale = NA_real_, status = "included", missing_repeat_values = 0L,
    stringsAsFactors = FALSE)
  for (j in seq_along(columns)) {
    name <- columns[j]
    features$missing_repeat_values[j] <- sum(is.na(dat[[name]]))
    for (i in seq_along(person_rows)) {
      x <- dat[[name]][person_rows[[i]]]
      observed <- x[!is.na(x)]
      if (!length(observed) || any(!is.finite(observed)))
        stop("Child covariate '", name, "' needs a finite person-level value for individual '", ids[i],
          "'. Missing child covariate features are not silently dropped or imputed.")
      changes <- diff(range(observed)) > 1e-10 * max(1, max(abs(observed)))
      if (changes && !name %in% allowed_mean)
        stop("Child covariate '", name, "' changes across visits for individual '", ids[i],
          "'. PCs and sex must be stable; explicitly name time-varying columns in person_mean_covariates to condition scores on their person mean.")
      C[i, j] <- if (name %in% allowed_mean) mean(observed) else observed[1L]
    }
    if (name %in% allowed_mean) features$aggregation[j] <- "person_mean_of_supplied_visits"
    features$center[j] <- mean(C[, j]); features$scale[j] <- stats::sd(C[, j])
    if (!is.finite(features$scale[j]) || features$scale[j] <=
        1e-12 * max(1, max(abs(C[, j])))) features$status[j] <- "constant_between_persons"
  }
  keep <- which(features$status == "included")
  if (length(keep)) {
    C <- sweep(sweep(C[, keep, drop = FALSE], 2L, features$center[keep], "-"),
      2L, features$scale[keep], "/")
    rank <- qr(C, tol = 1e-10, LAPACK = FALSE)
    independent <- sort(rank$pivot[seq_len(rank$rank)])
    redundant <- setdiff(seq_along(keep), independent)
    if (length(redundant)) features$status[keep[redundant]] <- "linearly_redundant_person_feature"
    keep <- keep[independent]; C <- C[, independent, drop = FALSE]
  }
  prepared$covariate_features <- features
  if (!length(keep)) { prepared$spec$covariate <- NULL; return(prepared) }
  q <- ncol(C)
  families <- sort(unique(fam))
  person_family <- fam[first]
  person_groups <- split(seq_along(ids), factor(person_family, levels = families))
  family_values <- lapply(person_groups, function(rows) C[rows[order(ids[rows])], , drop = FALSE])
  sizes <- vapply(family_values, nrow, integer(1))
  has_siblings <- any(sizes > 1L)
  moments <- .ravel_mixed_covariate_moments(family_values)
  features$shared_within_family <- FALSE
  features$shared_within_family[keep[moments$shared]] <- TRUE
  prepared$covariate_features <- features
  eta <- moments$eta; triangle <- moments$triangle; index <- moments$index
  moment_scores <- moments$moment_scores
  roles <- c("mother", "father", "child_self", if (has_siblings) "child_sibling")
  parameter_names <- stats::setNames(lapply(roles, function(role)
    paste0("C_D_", role, "_", seq_len(q))), roles)
  if (has_siblings && length(moments$shared))
    parameter_names$child_sibling[moments$shared] <-
      parameter_names$child_self[moments$shared]
  covariate <- list(q = q, columns = features$column[keep], features = features,
    eta = eta, triangle = triangle, index = index, has_siblings = has_siblings,
    shared = moments$shared, estimation = moments$diagnostics,
    moment_scores = moment_scores, family_ids = families, n_children = sizes,
    parameter_names = parameter_names, person_mean_covariates = allowed_mean,
    assumption = "Hidden residual after marginal projection onto visible genetics is independent of child covariates; PGS measurement errors are independent of them.",
    covariance_model = "Exchangeable child self/sibling covariances; marginal genetic AM and pi parameters retained",
    moment_weighting = "Unique-child mean and coherent exchangeable working covariance likelihood; estimated first-stage family influence propagated jointly",
    person_mean_scope = "Person means use all supplied visit rows; no informative-attendance or cohort-selection model is implied")
  moments <- .ravel_mixed_covariate_unpack(covariate)
  for (n in unique(sizes)) .ravel_mixed_covariate_inverse(moments, n)
  group_families <- lapply(prepared$groups, `[[`, "families")
  family_index <- match(unlist(group_families, use.names = FALSE), families)
  group_ends <- cumsum(lengths(group_families))
  group_starts <- c(1L, utils::head(group_ends, -1L) + 1L)
  if (anyNA(family_index)) stop("Covariate and genetic family identities do not agree.")
  for (g in seq_along(prepared$groups)) {
    z <- prepared$groups[[g]]
    values <- family_values[family_index[seq.int(group_starts[g], group_ends[g])]]
    if (any(vapply(values, nrow, integer(1)) != z$nchild))
      stop("Covariate and genetic child ordering do not agree.")
    prepared$groups[[g]]$C <- do.call(rbind, lapply(values, function(value) c(t(value))))
  }
  prepared$spec$covariate <- covariate
  new_parameters <- unique(unlist(parameter_names, use.names = FALSE))
  prepared$start <- c(prepared$start, stats::setNames(rep(0, length(new_parameters)), new_parameters))
  prepared$lower <- c(prepared$lower, stats::setNames(rep(-Inf, length(new_parameters)), new_parameters))
  prepared$upper <- c(prepared$upper, stats::setNames(rep(Inf, length(new_parameters)), new_parameters))
  prepared
}

.ravel_mixed_covariate_blocks <- function(par, spec, nchild, genetic = NULL, eta = NULL) {
  covariate <- spec[["covariate"]]
  if (is.null(covariate)) return(NULL)
  if (is.null(eta)) eta <- covariate$eta
  moments <- .ravel_mixed_covariate_unpack(covariate, eta)
  inverse <- .ravel_mixed_covariate_inverse(moments, nchild)
  q <- covariate$q; n <- nchild + 2L
  if (is.null(genetic)) genetic <- .ravel_mixed_genetic(nchild, par, spec)
  values <- lapply(covariate$parameter_names, function(names) unname(par[names]))
  if (any(!is.finite(unlist(values)))) stop("Covariate genetic cross-covariance parameters must be finite.")
  mother <- values$mother; father <- values$father; own <- values$child_self
  sibling <- if (covariate$has_siblings) values$child_sibling else rep(0, q)
  D <- matrix(0, n, nchild * q)
  for (j in seq_len(nchild)) {
    columns <- (j - 1L) * q + seq_len(q)
    D[1L, columns] <- mother; D[2L, columns] <- father
    for (i in seq_len(nchild)) D[i + 2L, columns] <- if (i == j) own else sibling
  }
  roles <- c("mother", "father", rep("child", nchild))
  a1 <- exp(par[paste0("loading1_", roles)])
  a2 <- exp(par[paste0("loading2_", roles)])
  L <- rbind(diag(a1, n), diag(a2, n))
  A <- genetic$P[seq_len(n), n + seq_len(n), drop = FALSE] / outer(a1, a2)
  A <- (A + t(A)) / 2
  TV <- t(sweep(genetic$PT[seq_len(n), , drop = FALSE], 1L, a1, "/"))
  R <- t(solve(A, t(TV)))
  D_Qinv <- D %*% inverse$Qinv
  Schur <- D_Qinv %*% t(D); Schur <- (Schur + t(Schur)) / 2
  conditional_visible <- A - Schur
  if (min(eigen(conditional_visible, symmetric = TRUE, only.values = TRUE)$values) <= 1e-10)
    stop("Genetic-covariate cross-covariance leaves a non-positive conditional visible-genetic covariance.")
  conditional <- list(P = genetic$P - L %*% Schur %*% t(L),
    PT = genetic$PT - L %*% Schur %*% t(R),
    T = genetic$T - R %*% Schur %*% t(R), state = genetic$state)
  c(list(genetic = conditional, marginal_genetic = genetic, mean = moments$mean,
    R = R, L = L, D = D, D_Qinv = D_Qinv, Schur = Schur,
    marginal_visible = A, conditional_visible = conditional_visible,
    nchild = nchild, q = q, Qself = moments$Qself, Qsib = moments$Qsib), inverse)
}

.ravel_mixed_covariate_means <- function(blocks, C) {
  if (!is.matrix(C) || ncol(C) != blocks$nchild * blocks$q)
    stop("C must contain contiguous child-covariate blocks, one family per row.")
  centered <- sweep(C, 2L, rep(blocks$mean, blocks$nchild), "-")
  V <- centered %*% t(blocks$D_Qinv)
  list(V = V, P = V %*% t(blocks$L), T = V %*% t(blocks$R))
}

.ravel_mixed_covariate_joint <- function(par, group, spec, eta = NULL, blocks = NULL,
                                        include_offset = TRUE) {
  if (is.null(spec[["covariate"]])) return(list(S = .ravel_mixed_covariance(par, group, spec),
    offset = if (include_offset) matrix(0, group$n, ncol(group$observed)) else NULL,
    blocks = NULL, B = NULL, H = NULL))
  if (is.null(blocks)) blocks <- .ravel_mixed_covariate_blocks(par, spec, group$nchild, eta = eta)
  n <- group$nchild + 2L; m <- length(group$child)
  beta <- matrix(par[spec$beta_names], nrow = 3L)
  effects <- group$moderators %*% t(beta)
  B <- matrix(0, m, n)
  B[, 1L] <- effects[, 2L]; B[, 2L] <- effects[, 3L]
  B[cbind(seq_len(m), group$child + 2L)] <- effects[, 1L]
  F <- rbind(blocks$L, B %*% blocks$R)
  H <- F %*% blocks$D_Qinv
  offset <- if (include_offset) {
    centered <- sweep(group$C, 2L, rep(blocks$mean, group$nchild), "-")
    centered %*% t(H)
  } else NULL
  list(S = .ravel_mixed_covariance(par, group, spec, genetic = blocks$genetic),
    offset = offset, blocks = blocks, B = B, H = H)
}
