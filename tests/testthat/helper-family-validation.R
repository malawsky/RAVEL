# Independent generative reference for family-subset regressions. None of the
# construction or oracle calculations calls RAVEL's correction machinery.
family_validation_population <- function(
    pi = .4, rho = .3, parent_variance = 1.15,
    beta = c(child = .12, mother = .25, father = -.08,
             mat_sib = .03, pat_sib = -.02),
    lineage = c(maternal = .2, paternal = -.1),
    noise = c(.4, 1.6), residual_variance = .7) {
  roles <- c("child", "mother", "father", "mat_sib", "pat_sib")
  stopifnot(pi >= 0, pi < 1, abs(rho) < 1, parent_variance > .5,
            length(noise) == 2L, all(noise > 0), residual_variance > 0,
            length(lineage) == 2L)
  if (is.null(names(beta))) names(beta) <- roles
  stopifnot(setequal(names(beta), roles))
  beta <- beta[roles]
  w <- c(sqrt(1 - pi), sqrt(pi))
  S <- diag(2) + (parent_variance - 1) * tcrossprod(w)
  A <- S - diag(2) / 2
  K <- rho * parent_variance * tcrossprod(w)
  H <- A %*% solve(S)
  Q <- A - H %*% A

  # Primitive variables: parents; innovations in their respective parental
  # midpoints conditional on the parent; and three segregation innovations.
  P <- matrix(0, 25, 25)
  P[1:2, 1:2] <- P[3:4, 3:4] <- S
  P[1:2, 3:4] <- P[3:4, 1:2] <- K
  P[5:6, 5:6] <- P[7:8, 7:8] <- Q
  P[9:10, 9:10] <- P[11:12, 11:12] <- P[13:14, 13:14] <- diag(2) / 2
  Z <- matrix(0, 14, 25)
  Z[1:2, 1:2] <- Z[1:2, 3:4] <- diag(2) / 2
  Z[1:2, 13:14] <- diag(2)
  Z[3:4, 1:2] <- Z[5:6, 3:4] <- diag(2)
  Z[11:12, 1:2] <- Z[13:14, 3:4] <- H
  Z[11:12, 5:6] <- Z[13:14, 7:8] <- diag(2)
  Z[7:8, ] <- Z[11:12, ]
  Z[7:8, 9:10] <- diag(2)
  Z[9:10, ] <- Z[13:14, ]
  Z[9:10, 11:12] <- diag(2)

  # The two GWAS error processes are independent of one another and of all
  # true liabilities. Within each process, shared weights induce pedigree
  # relatedness, including 1/4 between the child and a parental sibling.
  kin <- matrix(c(1, .5, .5, .25, .25,
                  .5, 1, 0, .5, 0,
                  .5, 0, 1, 0, .5,
                  .25, .5, 0, 1, 0,
                  .25, 0, .5, 0, 1), 5, 5)
  P[15:19, 15:19] <- noise[1] * kin
  P[20:24, 20:24] <- noise[2] * kin
  P[25, 25] <- residual_variance

  Tmap <- kronecker(diag(7), matrix(w, 1)) %*% Z
  Vmap <- Z[seq(1, 9, 2), , drop = FALSE]
  G1 <- G2 <- Vmap
  G1[, 15:19] <- diag(5)
  G2[, 20:24] <- diag(5)
  ymap <- drop(c(beta, lineage) %*% Tmap)
  ymap[25] <- 1
  observed_map <- rbind(G1, G2, ymap)
  rownames(observed_map) <- c(paste0("PGS1_", roles),
                             paste0("PGS2_", roles), "Y")
  truth_map <- rbind(Tmap, Vmap)
  rownames(truth_map) <- c(paste0("T_", roles), "T_mat_lineage",
                          "T_pat_lineage", paste0("V_", roles))
  TT <- Tmap[1:5, , drop = FALSE] %*% P %*% t(Tmap[1:5, , drop = FALSE])
  GG <- Vmap %*% P %*% t(Vmap)
  dimnames(TT) <- dimnames(GG) <- list(roles, roles)
  cyT <- setNames(drop(Tmap[1:5, , drop = FALSE] %*% P %*% ymap), roles)
  cyG <- setNames(drop(Vmap %*% P %*% ymap), roles)
  list(roles = roles, primitive_cov = P, observed_map = observed_map,
       truth_map = truth_map, Sigma_TT = TT, Sigma_GG = GG,
       cov_TY = cyT, cov_GY = cyG, var_Y = drop(ymap %*% P %*% ymap),
       observed_cov = observed_map %*% P %*% t(observed_map),
       pi = pi, rho = rho, parent_variance = parent_variance,
       beta = beta, lineage = lineage, noise = noise)
}

family_validation_data <- function(population, N = 2000L, seed = 419L,
                                   exact = TRUE) {
  stopifnot(N > nrow(population$primitive_cov))
  set.seed(seed)
  z <- matrix(stats::rnorm(N * nrow(population$primitive_cov)), N)
  if (exact) {
    # Orthogonalize all primitive innovations together; exact covariance
    # fixtures are used for algebraic recovery, never for coverage claims.
    # lavaan's likelihood uses the N-denominator covariance. Match that
    # convention so exact-fit coefficients need no finite-N scale adjustment.
    z <- sqrt(N) * qr.Q(qr(cbind(1, z)))[, -1, drop = FALSE]
  }
  innovations <- z %*% chol(population$primitive_cov)
  observed <- innovations %*% t(population$observed_map)
  truth <- innovations %*% t(population$truth_map)
  dat <- as.data.frame(cbind(observed, truth))
  dat$FID <- seq_len(N)
  dat
}

family_validation_oracle <- function(population, outcome_roles) {
  S <- population$Sigma_TT[outcome_roles, outcome_roles, drop = FALSE]
  beta <- setNames(drop(solve(S, population$cov_TY[outcome_roles])), outcome_roles)
  std <- beta * sqrt(diag(S))
  contrasts <- c()
  for (pair in list(c("mother", "mat_sib"), c("father", "pat_sib"))) {
    if (all(pair %in% outcome_roles)) {
      nm <- if (pair[1] == "mother") "delta_m" else "delta_f"
      contrasts[nm] <- beta[pair[1]] - beta[pair[2]]
    }
  }
  list(beta_raw = beta, beta_std = std, contrasts_raw = contrasts,
       variance_explained = drop(t(beta) %*% S %*% beta),
       R2 = drop(t(beta) %*% S %*% beta) / population$var_Y)
}

family_validation_subsets <- function(focal = "mother") {
  optional <- if (focal == "mother") c("mat_sib", "father", "pat_sib") else
    c("pat_sib", "mother", "mat_sib")
  lapply(0:7, function(mask) c(focal, optional[
    as.logical(intToBits(mask)[seq_along(optional)])]))
}
