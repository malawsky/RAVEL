# Independent generative reference for the optional mixed-model extension.
# This file must not call RAVEL fitting or correction helpers: its covariance
# algebra, inheritance, scores and outcomes are the validation oracle.

mixed_validation_population <- function(pi = .4, rho_total = .3,
                                        parent_variance = 1.2,
                                        mode = "estimated_state",
                                        n_generations = 3L) {
  stopifnot(length(pi) == 1L, is.finite(pi), pi >= 0, pi < 1,
            abs(rho_total) < 1, mode %in% c("estimated_state", "equilibrium", "n_gen"))
  if (mode == "equilibrium") parent_variance <- 1 / (1 - rho_total)
  if (mode == "n_gen") {
    stopifnot(n_generations >= 1L, n_generations == as.integer(n_generations))
    parent_variance <- 1
    if (n_generations > 1L) for (g in seq_len(n_generations - 1L))
      parent_variance <- .5 * (1 + rho_total) * parent_variance + .5
  }
  stopifnot(is.finite(parent_variance), parent_variance > .5)
  w <- c(sqrt(1 - pi), sqrt(pi))
  S <- diag(2L) + (parent_variance - 1) * tcrossprod(w)
  K <- rho_total * parent_variance * tcrossprod(w)
  parental <- rbind(cbind(S, K), cbind(K, S))
  child_parent <- .5 * (S + K)
  child <- child_parent + diag(.5, 2L)
  component <- rbind(cbind(child, child_parent, child_parent),
                     cbind(child_parent, S, K),
                     cbind(child_parent, K, S))
  roles <- c("child", "mother", "father")
  total_map <- kronecker(diag(3L), matrix(w, nrow = 1L))
  visible_map <- matrix(0, 3L, 6L)
  visible_map[cbind(1:3, c(1, 3, 5))] <- 1 / sqrt(diag(component)[c(1, 3, 5)])
  TT <- total_map %*% component %*% t(total_map)
  GG <- visible_map %*% component %*% t(visible_map)
  GT <- visible_map %*% component %*% t(total_map)
  dimnames(TT) <- dimnames(GG) <- dimnames(GT) <- list(roles, roles)
  list(pi = pi, rho_total = rho_total, parent_variance = parent_variance,
       mode = mode, n_generations = n_generations, w = w,
       S = S, K = K, parental = parental, component = component,
       Sigma_TT = TT, Sigma_GG = GG, Sigma_GT = GT,
       total_to_visible = solve(GG, GT),
       total_sd = sqrt(diag(TT)), visible_raw_sd = sqrt(diag(component)[c(1, 3, 5)]),
       rho_visible = GG["mother", "father"], eta_visible = GG["child", "mother"])
}

.mixed_validation_mvn <- function(n, covariance) {
  k <- nrow(covariance)
  matrix(stats::rnorm(n * k), n, k) %*% chol(covariance)
}

# Direct linear maps of independent primitive innovations. This intentionally
# uses a different construction from the fitter's compact T/kinship identity.
mixed_validation_observed_covariance <- function(
    population = mixed_validation_population(), n_children = 2L,
    child_index = rep(seq_len(n_children), each = 3L),
    age = rep(c(4, 10, 16), n_children), sex = rep(rep(c(0, 1), length.out = n_children), each = 3L),
    beta = c(child = .3, mother = .12, father = -.04),
    age_beta = c(child = .025, mother = -.008, father = 0),
    sex_beta = c(child = .12, mother = -.06, father = .02),
    reliability = c(.65, .40), score_scales = matrix(1, 2L, 3L),
    family_variance = .25, individual_variance = .2,
    slope_variance = .003, intercept_slope_covariance = .012,
    residual_variance = .75, age_center = 10) {
  stopifnot(length(child_index) == length(age), length(age) == length(sex),
            all(child_index %in% seq_len(n_children)))
  np <- n_children + 2L
  primitive <- matrix(0, 2L * np, 2L * np)
  primitive[1:4, 1:4] <- population$parental
  diag(primitive)[-(1:4)] <- .5
  component_map <- diag(2L * np)
  for (j in seq_len(n_children)) {
    r <- 2L * (j + 1L) + 1:2
    component_map[r, 1:2] <- component_map[r, 3:4] <- diag(.5, 2L)
  }
  Vmap <- component_map[seq.int(1L, 2L * np, by = 2L), , drop = FALSE]
  Tmap <- kronecker(diag(np), matrix(population$w, nrow = 1L)) %*% component_map
  roles <- c("mother", "father", rep("child", n_children))
  role_index <- match(roles, c("child", "mother", "father"))
  scales <- lapply(1:2, function(k) score_scales[k, role_index])
  Pmap <- rbind(scales[[1]] * Vmap, scales[[2]] * Vmap)
  # Error innovations: two unrelated parents and one independent child
  # segregation error per sibling, separately for each GWAS.
  emap <- diag(np)
  emap[-(1:2), 1:2] <- .5
  diag(emap)[-(1:2)] <- sqrt(.5)
  noise_cov <- matrix(0, 2L * np, 2L * np)
  for (k in 1:2) {
    rows <- (k - 1L) * np + seq_len(np)
    A <- scales[[k]] * emap * sqrt((1 - reliability[k]) / reliability[k])
    noise_cov[rows, rows] <- tcrossprod(A)
  }
  age0 <- age - age_center
  beta_rows <- outer(rep(1, length(age)), beta[c("child", "mother", "father")]) +
    outer(age0, age_beta[c("child", "mother", "father")]) +
    outer(sex, sex_beta[c("child", "mother", "father")])
  Ymap <- beta_rows[, "mother"] * matrix(Tmap[1, ], length(age), ncol(Tmap), byrow = TRUE) +
    beta_rows[, "father"] * matrix(Tmap[2, ], length(age), ncol(Tmap), byrow = TRUE) +
    beta_rows[, "child"] * Tmap[child_index + 2L, , drop = FALSE]
  same <- outer(child_index, child_index, "==")
  residual <- diag(residual_variance, length(age)) + family_variance +
    same * (individual_variance + slope_variance * outer(age0, age0) +
      intercept_slope_covariance * outer(age0, age0, "+"))
  PP <- Pmap %*% primitive %*% t(Pmap) + noise_cov
  PY <- Pmap %*% primitive %*% t(Ymap)
  YY <- Ymap %*% primitive %*% t(Ymap) + residual
  rbind(cbind(PP, PY), cbind(t(PY), YY))
}

mixed_validation_data <- function(
    N_families = 500L, sib_sizes = 2L, visits = 1L, seed = 31415L,
    pi = .4, rho_total = .3, parent_variance = 1.2,
    mode = "estimated_state", n_generations = 3L,
    beta = c(child = .3, mother = .12, father = -.04),
    age_beta = c(child = 0, mother = 0, father = 0),
    sex_beta = c(child = 0, mother = 0, father = 0),
    reliability = c(.65, .40), score_scales = NULL,
    family_variance = .25, individual_variance = 0,
    slope_variance = 0, intercept_slope_covariance = 0,
    residual_variance = .75, families_per_FID = 1L, extended_fraction = 0,
    age_center = 10, age_range = c(4, 16), age_pattern = c("random", "scheduled"),
    age_effect = .03, sex_effect = .12, intercept = NULL,
    outcome = c("continuous", "binary_liability"), prevalence_reference = .1,
    missing_outcome_fraction = 0, genetic_sex_association = 0,
    genetic_age_association = 0) {
  outcome <- match.arg(outcome)
  age_pattern <- match.arg(age_pattern)
  stopifnot(N_families >= 1L, N_families == as.integer(N_families),
            length(reliability) == 2L, all(reliability > 0 & reliability < 1),
            family_variance >= 0, individual_variance >= 0, slope_variance >= 0,
            residual_variance > 0, families_per_FID >= 1L,
            extended_fraction >= 0, extended_fraction <= 1,
            missing_outcome_fraction >= 0, missing_outcome_fraction < 1,
            length(age_range) == 2L, diff(age_range) > 0)
  if (outcome == "binary_liability" && residual_variance != 1)
    stop("Binary validation fixes the conditional probit residual variance to 1.")
  previous_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (is.null(previous_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", previous_seed, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  population <- mixed_validation_population(pi, rho_total, parent_variance, mode, n_generations)
  roles <- c("child", "mother", "father")
  normalize_beta <- function(x) {
    stopifnot(length(x) == 3L, setequal(names(x), roles), all(is.finite(x)))
    x[roles]
  }
  beta <- normalize_beta(beta)
  age_beta <- normalize_beta(age_beta)
  sex_beta <- normalize_beta(sex_beta)
  if (length(sib_sizes) == 1L) sib_sizes <- rep(sib_sizes, N_families)
  stopifnot(length(sib_sizes) == N_families, all(sib_sizes >= 1L),
            all(sib_sizes == as.integer(sib_sizes)))
  family <- rep(seq_len(N_families), sib_sizes)
  N_children <- length(family)
  if (length(visits) == 1L) visits <- rep(visits, N_children)
  stopifnot(length(visits) == N_children, all(visits >= 1L),
            all(visits == as.integer(visits)))
  subject <- rep(seq_len(N_children), visits)
  n <- length(subject)
  fid <- (seq_len(N_families) - 1L) %/% as.integer(families_per_FID) + 1L
  parents <- .mixed_validation_mvn(N_families, population$parental)
  child <- .5 * (parents[family, 1:2, drop = FALSE] + parents[family, 3:4, drop = FALSE]) +
    matrix(stats::rnorm(2L * N_children, sd = sqrt(.5)), N_children, 2L)
  latent_T <- cbind(child = drop(child %*% population$w),
                    mother = drop(parents[family, 1:2, drop = FALSE] %*% population$w),
                    father = drop(parents[family, 3:4, drop = FALSE] %*% population$w))
  latent_G <- sweep(cbind(child = child[, 1], mother = parents[family, 1],
                         father = parents[family, 3]), 2L, population$visible_raw_sd, "/")
  scales <- if (is.null(score_scales)) matrix(1, 2L, 3L,
    dimnames = list(c("PGS1", "PGS2"), roles)) else as.matrix(score_scales)
  stopifnot(identical(dim(scales), c(2L, 3L)), all(is.finite(scales)), all(scales > 0))
  pgs <- list()
  for (k in 1:2) {
    eparents <- matrix(stats::rnorm(2L * N_families), N_families, 2L)
    echild <- .5 * rowSums(eparents[family, , drop = FALSE]) +
      sqrt(.5) * stats::rnorm(N_children)
    raw <- cbind(child = child[, 1], mother = parents[family, 1], father = parents[family, 3])
    errors <- cbind(child = echild, mother = eparents[family, 1], father = eparents[family, 2])
    score <- sweep(raw + sqrt((1 - reliability[k]) / reliability[k]) * errors,
                   2L, scales[k, ], "*")
    for (role in roles) pgs[[paste0("PGS", k, "_", role)]] <- score[, role][subject]
  }
  sex <- stats::rbinom(N_children, 1L,
                       stats::plogis(genetic_sex_association * latent_G[, "child"]))
  age <- if (age_pattern == "random") stats::runif(n, age_range[1], age_range[2]) else
    unlist(lapply(visits, function(k) if (k == 1L) age_center else
      seq(age_range[1], age_range[2], length.out = k)), use.names = FALSE)
  age <- age + genetic_age_association * latent_G[subject, "child"]
  centered_age <- age - age_center
  familial <- sqrt(family_variance) * (
    sqrt(1 - extended_fraction) * stats::rnorm(N_families) +
      sqrt(extended_fraction) * stats::rnorm(max(fid))[fid])
  random_covariance <- matrix(c(individual_variance, intercept_slope_covariance,
                               intercept_slope_covariance, slope_variance), 2L)
  eig <- eigen(random_covariance, symmetric = TRUE)
  if (min(eig$values) < -1e-12) stop("Individual intercept/slope covariance must be positive semidefinite.")
  random_individual <- matrix(stats::rnorm(2L * N_children), N_children, 2L) %*%
    diag(sqrt(pmax(eig$values, 0)), 2L) %*% t(eig$vectors)
  if (is.null(intercept)) intercept <- if (outcome == "continuous") 0 else {
    variance_reference <- drop(t(beta) %*% population$Sigma_TT %*% beta) +
      family_variance + individual_variance + 1
    stats::qnorm(prevalence_reference) * sqrt(variance_reference)
  }
  genetic <- drop(latent_T[subject, , drop = FALSE] %*% beta) +
    centered_age * drop(latent_T[subject, , drop = FALSE] %*% age_beta) +
    sex[subject] * drop(latent_T[subject, , drop = FALSE] %*% sex_beta)
  linear_predictor <- intercept + genetic + age_effect * centered_age + sex_effect * sex[subject] +
    familial[family[subject]] + random_individual[subject, 1L] +
    centered_age * random_individual[subject, 2L]
  latent_Y <- linear_predictor + stats::rnorm(n, sd = sqrt(residual_variance))
  observed_Y <- if (outcome == "continuous") latent_Y else as.numeric(latent_Y > 0)
  if (missing_outcome_fraction > 0)
    observed_Y[stats::runif(n) < missing_outcome_fraction] <- NA_real_
  result <- data.frame(M_ID = sprintf("M%06d", family[subject]),
    F_ID = sprintf("P%06d", family[subject]),
    FID = sprintf("E%06d", fid[family[subject]]),
    IID = sprintf("I%07d", subject), as.data.frame(pgs),
    Y = observed_Y, age = age, sex = sex[subject], stringsAsFactors = FALSE)
  coefficient_labels <- c(child = "beta_D", mother = "beta_M", father = "beta_F")
  coefficients <- function(b, label = "") {
    names_out <- paste0(unname(coefficient_labels), label)
    list(corrected_raw = stats::setNames(unname(b), names_out),
         corrected_std = stats::setNames(unname(b * population$total_sd), names_out),
         visible = stats::setNames(drop(population$total_to_visible %*% b), names_out))
  }
  main <- coefficients(beta)
  age_truth <- coefficients(age_beta, ":age")
  sex_truth <- coefficients(sex_beta, ":sex")
  combined <- lapply(names(main), function(sc) c(main[[sc]], age_truth[[sc]], sex_truth[[sc]]))
  names(combined) <- names(main)
  attr(result, "truth") <- list(population = population, coefficients = combined,
    beta_raw = beta, age_beta_raw = age_beta, sex_beta_raw = sex_beta,
    family_variance = family_variance, individual_variance = individual_variance,
    slope_variance = slope_variance, intercept_slope_covariance = intercept_slope_covariance,
    residual_variance = residual_variance, extended_fraction = extended_fraction,
    fixed_effects = c(intercept = intercept, age = age_effect, sex = sex_effect),
    age_center = age_center, sex_center = 0, outcome = outcome,
    prevalence_reference = if (outcome == "binary_liability") prevalence_reference else NA_real_,
    expected_prevalence = if (outcome == "binary_liability") mean(stats::pnorm(linear_predictor)) else NA_real_,
    N_families = N_families, N_children = N_children, N_rows = n, N_FIDs = max(fid),
    seed = seed, reliability_founder = reliability,
    omitted_covariate_independence = c(sex = genetic_sex_association, age = genetic_age_association))
  attr(result, "latent") <- data.frame(IID = sprintf("I%07d", seq_len(N_children)),
    M_ID = sprintf("M%06d", family), FID = sprintf("E%06d", fid[family]),
    stats::setNames(as.data.frame(latent_T), paste0("T_", roles)),
    stats::setNames(as.data.frame(latent_G), paste0("G_", roles)),
    family_effect = familial[family], individual_intercept = random_individual[, 1L],
    individual_slope = random_individual[, 2L], stringsAsFactors = FALSE)
  result
}
