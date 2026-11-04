.mixed_covariate_preparation <- function(dat, covariates = c("PC1", "PC2", "sex", "age"),
    interactions = NULL, person_mean = "age") {
  p <- .ravel_mixed_prepare(dat, "continuous", "Y",
    c(child1 = "PGS1_child", child2 = "PGS2_child", mother1 = "PGS1_mother",
      mother2 = "PGS2_mother", father1 = "PGS1_father", father2 = "PGS2_father"),
    "M_ID", "IID", "FID", covariates, interactions, NULL,
    FALSE, FALSE, NULL, "estimated_state", 1L, .4, FALSE,
    covariate_model = "joint")
  p$spec$person_mean_covariates <- intersect(person_mean, p$spec$covariates)
  .ravel_mixed_covariate_prepare(p)
}

test_that("child covariate preparation preserves person and family information", {
  dat <- mixed_covariate_data(N_families = 80L, visits = 3L, seed = 811)
  # An entirely missing outcome does not remove the child's covariate/genetic block.
  id <- dat$IID[1L]; dat$Y[dat$IID == id] <- NA_real_
  p <- .mixed_covariate_preparation(dat)
  expect_equal(p$n_subjects, length(unique(dat$IID)))
  meta <- p$spec$covariate
  expect_equal(meta$q, 4L)
  expect_lt(max(abs(colSums(meta$moment_scores))), 1e-10)
  expect_equal(meta$features$aggregation[meta$features$column == "age"],
    "person_mean_of_supplied_visits")
  for (z in p$groups) {
    rows <- dat[dat$M_ID == z$families, ]
    ids <- sort(unique(rows$IID))
    raw <- t(vapply(ids, function(i) vapply(meta$columns,
      function(v) mean(rows[rows$IID == i, v]), numeric(1)), numeric(meta$q)))
    f <- meta$features[match(meta$columns, meta$features$column), ]
    expected <- sweep(sweep(raw, 2L, f$center, "-"), 2L, f$scale, "/")
    expect_equal(unname(z$C[1L, ]), unname(c(t(expected))), tolerance = 1e-12)
  }
  bad <- dat; bad$PC1[1L] <- bad$PC1[1L] + .1
  expect_error(.mixed_covariate_preparation(bad), "PC1.*changes across visits")
  expect_error(.mixed_covariate_preparation(dat, person_mean = character()),
    "age.*changes across visits")
  bad <- dat; bad$PC1[bad$IID == id] <- NA_real_
  expect_error(.mixed_covariate_preparation(bad), "PC1.*finite person-level")
  partial <- dat; partial$PC1[1L] <- NA_real_
  pp <- .mixed_covariate_preparation(partial)
  expect_equal(pp$spec$covariate$eta, meta$eta, tolerance = 1e-12)
  expect_equal(pp$covariate_features$missing_repeat_values[1L], 1L)
})

test_that("covariate raw-moment influence and unpacking derivatives agree with direct changes", {
  for (siblings in c(1L, 2L)) {
    d <- mixed_covariate_data(N_families = 80L, sib_sizes = siblings, seed = 812)
    p <- .mixed_covariate_preparation(d)
    cmeta <- p$spec$covariate
    physical <- function(eta) {
      x <- .ravel_mixed_covariate_unpack(cmeta, eta)
      c(x$mean, x$Qself, x$Qsib)
    }
    numerical <- numDeriv::jacobian(physical, cmeta$eta)
    expect_equal(unname(.ravel_mixed_covariate_unpack_jacobian(cmeta)),
      unname(numerical), tolerance = 1e-7)
    expect_equal(cmeta$has_siblings, siblings > 1L)
    expect_equal(length(cmeta$parameter_names), if (siblings == 1L) 3L else 4L)
    family_C <- lapply(p$groups, function(z) matrix(z$C[1L, ], siblings, cmeta$q, byrow = TRUE))
    names(family_C) <- vapply(p$groups, function(z) z$families, character(1))
    family_C <- family_C[cmeta$family_ids]
    weighted_moments <- function(w) {
      count <- vapply(family_C, nrow, integer(1)); n <- sum(w * count)
      first <- Reduce(`+`, Map(function(x, a) a * colSums(x), family_C, w)) / n
      second <- Reduce(`+`, Map(function(x, a) a * crossprod(x), family_C, w)) / n
      other <- if (siblings > 1L) Reduce(`+`, Map(function(x, a)
        a * (tcrossprod(colSums(x)) - crossprod(x)), family_C, w)) /
          sum(w * count * (count - 1L)) else NULL
      c(first, second[cmeta$triangle], if (siblings > 1L) other[cmeta$triangle])
    }
    w <- rep(1, length(family_C)); plus <- minus <- w
    plus[3L] <- plus[3L] + 1e-5; minus[3L] <- minus[3L] - 1e-5
    derivative <- (weighted_moments(plus) - weighted_moments(minus)) / 2e-5
    expect_equal(unname(cmeta$moment_scores[3L, ] / length(w)), derivative, tolerance = 1e-8)
  }
})

test_that("exchangeable conditioning matches a dense full-family Schur complement", {
  d <- mixed_covariate_data(N_families = 100L, sib_sizes = 3L, visits = 2L, seed = 813)
  p <- .mixed_covariate_preparation(d, interactions = c("age", "sex"))
  par <- p$start
  dnames <- unlist(p$spec$covariate$parameter_names, use.names = FALSE)
  par[dnames] <- .025 * sin(seq_along(dnames))
  par[p$spec$beta_names] <- seq(.02, .18, length.out = length(p$spec$beta_names))
  z <- p$groups[[1L]]
  b <- .ravel_mixed_covariate_blocks(par, p$spec, z$nchild)
  moments <- .ravel_mixed_covariate_unpack(p$spec$covariate)
  Q <- kronecker(diag(z$nchild), moments$Qself - moments$Qsib) +
    kronecker(matrix(1, z$nchild, z$nchild), moments$Qsib)
  expect_equal(b$Qinv, solve(Q), tolerance = 1e-11)
  got <- .ravel_mixed_covariate_joint(par, z, p$spec, blocks = b)
  marginal <- .ravel_mixed_covariance(par, z, p$spec)
  CovXC <- rbind(b$L, got$B %*% b$R) %*% b$D
  expected <- marginal - CovXC %*% solve(Q, t(CovXC))
  expect_equal(got$S, expected, tolerance = 1e-11)
  centered <- z$C - rep(moments$mean, z$nchild)
  expect_equal(got$offset, centered %*% solve(Q, t(CovXC)), tolerance = 1e-11)
  expect_equal(.ravel_mixed_covariate_joint(par, z, p$spec, blocks = b,
    include_offset = FALSE)$H, got$H)
  # The parent's score mean is represented once even though children have visits.
  expect_equal(nrow(got$H), 2L * (z$nchild + 2L) + length(z$child))
  zero <- par; zero[dnames] <- 0
  expect_equal(.ravel_mixed_covariate_joint(zero, z, p$spec)$S,
    .ravel_mixed_covariance(zero, z, p$spec), tolerance = 1e-12)
  expect_equal(unname(.ravel_mixed_covariate_joint(zero, z, p$spec)$offset),
    matrix(0, z$n, ncol(z$observed)), tolerance = 1e-12)
})

test_that("redundant person features and invalid exchangeability are explicit", {
  dat <- mixed_covariate_data(N_families = 80L, visits = 3L, seed = 814)
  p <- .mixed_covariate_preparation(dat,
    covariates = c("PC1", "PC2", "sex", "baseline_age", "age"))
  expect_equal(p$spec$covariate$q, 4L)
  expect_equal(sum(p$covariate_features$status == "linearly_redundant_person_feature"), 1L)
  dat$age <- rep(c(9, 10, 11), length.out = nrow(dat))
  p <- .mixed_covariate_preparation(dat, covariates = "age")
  expect_null(p$spec$covariate)
  expect_equal(p$covariate_features$status, "constant_between_persons")
  moments <- list(mean = c(0, 0), Qself = diag(2), Qsib = -.6 * diag(2))
  expect_error(.ravel_mixed_covariate_inverse(moments, 3L), "not positive definite")
})

.mixed_covariate_audit_fixture <- function(siblings = 2L, binary = FALSE) {
  dat <- mixed_covariate_data(N_families = 80L, sib_sizes = siblings, seed = 815)
  p <- .mixed_covariate_preparation(dat, covariates = c("PC1", "PC2"))
  par <- p$start; d <- unlist(p$spec$covariate$parameter_names, use.names = FALSE)
  par[d] <- .04 * cos(seq_along(d))
  if (binary) {
    p$spec$outcome <- "binary_liability"; p$spec$outcome_scale <- 1
    p$spec$family_random <- TRUE
    par <- c(par[names(par) != "residual_sd"], family_variance = 0)
  }
  eta <- p$spec$covariate$eta; all <- c(par, eta)
  v <- seq(.01, .03, length.out = length(all))
  V <- diag(.002, length(all)) + tcrossprod(v)
  dimnames(V) <- list(names(all), names(all))
  # These retained objects deliberately omit the large first-stage score array.
  p$spec$covariate$moment_scores <- NULL
  list(settings = p$spec, parameters = par, prepared = p,
    score_scaling = p$score_scaling, vcov = V[names(par), names(par), drop = FALSE],
    uncertainty = list(nuisance_vcov = V[names(eta), names(eta), drop = FALSE],
      param_nuisance_cov = V[names(par), names(eta), drop = FALSE]), joint_vcov = V)
}

test_that("covariate audit propagates cross-stage covariance and works without raw moment scores", {
  f <- .mixed_covariate_audit_fixture()
  a <- .ravel_mixed_audit_covariates(f, 1L)
  expect_equal(a$covariate_audit_status$status, "ok")
  expect_true(all(is.finite(a$conditional_implied_covariances$se)))
  expect_true(all(is.finite(a$genetic_covariate_covariances$se)))
  par <- f$parameters; eta <- f$settings$covariate$eta; x <- c(par, eta)
  label <- "implied|conditional_input_PGS_outcome_covariance|PGS1_mother|PGS1_mother"
  fun <- function(v) .ravel_mixed_audit_covariate_values(v[names(par)], v[names(eta)],
    f, f$prepared$groups[[1L]])[[label]]
  gradient <- numDeriv::grad(fun, x)
  expected <- sqrt(drop(gradient %*% f$joint_vcov %*% gradient))
  row <- a$conditional_implied_covariances$lhs == "PGS1_mother" &
    a$conditional_implied_covariances$rhs == "PGS1_mother"
  expect_equal(a$conditional_implied_covariances$se[row], expected, tolerance = 1e-7)
  no_cross <- f; no_cross$uncertainty$param_nuisance_cov[,] <- 0
  b <- .ravel_mixed_audit_covariates(no_cross, 1L)
  expect_gt(abs(a$conditional_implied_covariances$se[row] -
    b$conditional_implied_covariances$se[row]), 1e-7)
  f$prepared <- NULL
  without_data <- .ravel_mixed_audit_covariates(f, 1L)
  expect_equal(without_data$covariate_audit_status$status, "moments_only")
  expect_equal(without_data$covariate_moments$estimate, a$covariate_moments$estimate)
  expect_equal(without_data$covariate_moments$se, a$covariate_moments$se)
  expect_equal(nrow(without_data$conditional_implied_covariances), 0L)
})

test_that("covariate audit does not invent sibling moments and retains binary boundary estimates", {
  f <- .mixed_covariate_audit_fixture(siblings = 1L)
  a <- .ravel_mixed_audit_covariates(f, 1L)
  expect_false(any(grepl("sibling", a$covariate_moments$block)))
  expect_false(any(grepl("child2", a$genetic_covariate_covariances$rhs)))
  binary <- .mixed_covariate_audit_fixture(binary = TRUE)
  # Regular Wald inference is unavailable at this variance boundary.
  binary$vcov <- NULL; binary$uncertainty <- NULL
  a <- .ravel_mixed_audit_covariates(binary, 1L)
  expect_equal(a$covariate_audit_status$status, "estimates_only")
  expect_true(all(is.finite(a$conditional_implied_covariances$estimate)))
  expect_true(all(is.na(a$conditional_implied_covariances$se)))
  expect_true(all(a$conditional_implied_covariances$block ==
    "conditional_input_PGS_probit_liability_covariance"))
})
