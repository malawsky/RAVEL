# Run from the repository root, or supply the destination .rda path.
args <- commandArgs(trailingOnly = TRUE)
destination <- if (length(args)) args[1] else "data/ravel_example.rda"

local({
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(20260913)
  n <- 4000L
  hidden_fraction <- .4
  rho_total <- .4
  beta <- c(beta_D = .45, beta_M = .12, beta_P = 0)
  reliability <- c(PGS1 = .5, PGS2 = .3)

  # One generation of the same Gaussian trio-moment model used in validation.
  parent_child <- (1 + rho_total) / 2
  ST <- matrix(c(1 + rho_total / 2, parent_child, parent_child,
                 parent_child, 1, rho_total, parent_child, rho_total, 1), 3L)
  SU <- matrix(c(1, .5, .5, .5, 1, 0, .5, 0, 1), 3L)
  SGraw <- (1 - hidden_fraction) * ST + hidden_fraction * SU
  sdG <- sqrt(diag(SGraw))
  SG <- SGraw / outer(sdG, sdG)
  SGT <- sqrt(1 - hidden_fraction) * sweep(ST, 1L, sdG, "/")
  T <- matrix(rnorm(n * 3L), n) %*% chol(ST)
  U <- matrix(rnorm(n * 3L), n) %*% chol(SU)
  Graw <- sqrt(1 - hidden_fraction) * T + sqrt(hidden_fraction) * U
  G <- sweep(Graw, 2L, sdG, "/")

  A <- rbind(c(.30, .10, .05), c(.05, .35, .05), c(.05, .10, .35),
             c(.20, .15, 0), c(.20, 0, .15), c(0, .25, .15), c(.15, .15, .15))
  noise_variance <- 1 - diag(A %*% SG %*% t(A))
  C <- G %*% t(A) + sweep(matrix(rnorm(n * 7L), n), 2L, sqrt(noise_variance), "*")
  SC <- A %*% SG %*% t(A) + diag(noise_variance)
  SGC <- SG %*% t(A)
  STC <- t(SGT) %*% t(A)
  sex <- rbinom(n, 1L, .5)
  site <- sample.int(3L, n, replace = TRUE, prob = c(.3, .3, .4))
  C <- cbind(C, sex - .5, (site == 2L) - .3, (site == 3L) - .4)
  SD <- matrix(c(.25, 0, 0, 0, .21, -.12, 0, -.12, .24), 3L)
  SC <- rbind(cbind(SC, matrix(0, 7L, 3L)), cbind(matrix(0, 3L, 7L), SD))
  SGC <- cbind(SGC, matrix(0, 3L, 3L))
  STC <- cbind(STC, matrix(0, 3L, 3L))
  gamma <- c(.15, -.10, .10, .08, -.08, .07, .05, .15, .10, -.10)
  varY <- drop(t(beta) %*% ST %*% beta + t(gamma) %*% SC %*% gamma +
                 2 * t(beta) %*% STC %*% gamma + .9^2)
  Y <- drop(T %*% beta + C %*% gamma + .9 * rnorm(n))
  threshold <- qnorm(.9) * sqrt(varY)

  ravel_example <- data.frame(FID = seq_len(n), Y = Y,
    Y_binary = as.integer(Y > threshold), age = 10 + 1.5 * C[, 1],
    sex = factor(sex, levels = c(0, 1), labels = c("female", "male")),
    site = factor(site, levels = 1:3, labels = c("A", "B", "C")))
  for (j in seq_len(6L)) ravel_example[[paste0("PC", j)]] <- C[, j + 1L]
  for (version in seq_len(2L)) {
    error <- matrix(rnorm(n * 3L), n) %*% chol(SU)
    P <- Graw + sqrt((1 - reliability[version]) / reliability[version]) * error
    for (j in seq_len(3L)) {
      name <- paste0("pgs", version, "_", c("child", "mother", "father")[j])
      ravel_example[[name]] <- P[, j]
    }
  }

  covGY <- SGT %*% beta + SGC %*% gamma
  covCY <- t(STC) %*% beta + SC %*% gamma
  joint <- rbind(cbind(SG, SGC), cbind(t(SGC), SC))
  attr(ravel_example, "simulation") <- list(
    model = "Gaussian trio moments; one generation of assortment on total liability",
    seed = 20260913L, n = n, hidden_fraction = hidden_fraction,
    rho_total = rho_total, n_generations = 1L, reliability_random_mating = reliability,
    total_coefficients_raw = beta, residual_sd = .9,
    phenotype_variance = varY, binary_threshold = threshold,
    covariates = c("age", "sex", "site", paste0("PC", 1:6)),
    categorical_covariates = c("sex", "site"),
    expected_visible = setNames(drop(solve(joint, c(covGY, covCY))[1:3]) / sqrt(varY), names(beta)),
    expected_total_std = beta * sqrt(diag(ST)) / sqrt(varY))
  stopifnot(nrow(ravel_example) == n, ncol(ravel_example) == 18L,
            !anyNA(ravel_example), !anyDuplicated(ravel_example$FID))
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  save(ravel_example, file = destination, compress = "xz", version = 2)
})
