gxe_cache <- new.env(parent = emptyenv())

# Gaussian latent trios measured at 1..max_occ ages (long format, one row per
# child-occasion), with child and maternal genetic x age effects.
gxe_data <- function(N = 4000L, max_occ = 4L, seed = 11L, rho = .2, rel = c(.5, .3),
                     gD = .02, gM = .01) {
  key <- paste(N, max_occ, seed, gD, gM, sep = "_")
  if (exists(key, gxe_cache, inherits = FALSE)) return(gxe_cache[[key]])
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv))
  set.seed(seed)
  z <- MASS::mvrnorm(N, c(0, 0), matrix(c(1, rho, rho, 1), 2))
  Gm <- z[, 1]; Gf <- z[, 2]; Gc <- (Gm + Gf) / 2 + sqrt(.5 * (1 - rho)) * rnorm(N)
  P <- list()
  for (k in 1:2) {
    em <- rnorm(N); ef <- rnorm(N); ec <- (em + ef) / 2 + sqrt(.5) * rnorm(N)
    P[[paste0("PGS", k, "_mother")]] <- sqrt(rel[k]) * Gm + sqrt(1 - rel[k]) * em
    P[[paste0("PGS", k, "_father")]] <- sqrt(rel[k]) * Gf + sqrt(1 - rel[k]) * ef
    P[[paste0("PGS", k, "_child")]] <- sqrt(rel[k]) * Gc + sqrt(1 - rel[k]) * ec
  }
  nocc <- sample(seq_len(max_occ), N, replace = TRUE)
  fam <- rep(seq_len(N), nocc)
  a <- stats::runif(length(fam), 4, 16) - 10
  Y <- (.30 + gD * a) * Gc[fam] + (.10 + gM * a) * Gm[fam] + .05 * Gf[fam] + .05 * a +
    rnorm(N, sd = sqrt(.3))[fam] + rnorm(N, sd = .03)[fam] * a + rnorm(length(fam), sd = sqrt(.6))
  out <- data.frame(FID = fam, as.data.frame(P)[fam, ], Y = Y, age = a + 10)
  rownames(out) <- NULL
  gxe_cache[[key]] <- out
  out
}

gxe_fit <- function(...) suppressWarnings(ravel_fit(...))

test_that("fits without interactions are unchanged", {
  d <- gxe_data()[!duplicated(gxe_data()$FID), ]
  fit <- gxe_fit(d)
  expect_null(fit$interactions)
  expect_null(fit$cluster)
  expect_false(grepl("GX", fit$model))
  expect_false(any(grepl(":", as.data.frame(fit)$parameter)))
})

test_that("product indicators share the genetic loadings", {
  d <- gxe_data()
  spec <- .ravel_interaction_spec("age", c("child", "mother"), NULL, d, "continuous")
  syn <- .ravel_interaction_syntax(spec, "age", c("mother", "father", "child"))
  expect_match(syn, "GXc_1 =~ lc1*PGS1_child_x1 + lc2*PGS2_child_x1", fixed = TRUE)
  expect_match(syn, "GXm_1 =~ lm1*PGS1_mother_x1 + lm2*PGS2_mother_x1", fixed = TRUE)
  expect_match(syn, "GXc_1 ~~ NA*GXc_1", fixed = TRUE)
  expect_match(syn, "Y ~ gD_1*GXc_1", fixed = TRUE)
  expect_match(syn, "PGS1_father ~~ PGS1_child_x1", fixed = TRUE)
  prepared <- .ravel_interaction_data(d, spec)
  expect_equal(mean(prepared$age), 0, tolerance = 1e-10)
  expect_equal(prepared$PGS1_child_x1,
    (d$PGS1_child - mean(d$PGS1_child)) * (d$age - mean(d$age)) / sd(d$age), tolerance = 1e-12)
})

test_that("genetic x age effects are recovered in long format with cluster-robust SEs", {
  d <- gxe_data(N = 6000L, seed = 21L)
  fit <- gxe_fit(d, interactions = "age", interaction_roles = c("child", "mother", "father"),
    interaction_center = c(age = 10), cluster = "FID")
  expect_true(fit$converged)
  expect_identical(lavaan::lavInspect(fit$fit, "options")$se, "robust.cluster")
  expect_equal(fit$n_clusters, 6000)
  tab <- fit$interactions$table
  expect_identical(tab$parameter, c("beta_D:age", "beta_M:age", "beta_F:age"))
  truth <- c(.02, .01, 0)
  expect_true(all(abs(tab$estimate - truth) < 4 * tab$se))
  expect_true(all(tab$uncertainty_source == "lavaan cluster-robust"))
  # Main effects are at the chosen centre (age 10).
  expect_lt(abs(fit$beta_vis[["bD"]] - .30), 4 * fit$lavaan_inference$visible$se[["bD"]])
  df <- as.data.frame(fit, include_vanilla = FALSE)
  expect_true(all(c("beta_D:age", "beta_M:age", "beta_F:age") %in% df$parameter))
  expect_match(df$uncertainty_source[df$parameter == "beta_D"], "cluster-robust")
  expect_output(print_simple(fit), "Main effects are at age = 10")
  expect_output(print_ravel_results(fit), "Genetic x covariate interactions")
})

test_that("changing the centre only re-expresses the main effects", {
  d <- gxe_data()
  at10 <- gxe_fit(d, interactions = "age", interaction_center = c(age = 10), cluster = "FID")
  at_mean <- gxe_fit(d, interactions = "age", cluster = "FID")
  g <- at_mean$interactions$table$estimate[1]
  expect_equal(at10$interactions$table$estimate[1], g, tolerance = 1e-4)
  expect_equal(at10$beta_vis[["bD"]],
    at_mean$beta_vis[["bD"]] + g * (10 - at_mean$interactions$centers[["age"]]), tolerance = 1e-4)
})

test_that("several moderators and a quadratic growth term are supported", {
  d <- gxe_data()
  d$age2 <- (d$age - 10)^2
  fit <- gxe_fit(d, interactions = c("age", "age2"), interaction_center = c(age = 10, age2 = 0),
    cluster = "FID")
  expect_true(fit$converged)
  expect_identical(fit$interactions$table$parameter, c("beta_D:age", "beta_D:age2"))
  # The data have no quadratic genetic effect.
  tab <- fit$interactions$table
  expect_lt(abs(tab$estimate[2]), 4 * tab$se[2])
})

test_that("interactions combine with constrained modes and parental siblings", {
  d <- gxe_data()
  fit <- gxe_fit(d, interactions = "age", interaction_center = c(age = 10), cluster = "FID",
    fit_constraints = "equilibrium")
  expect_true(fit$converged)
  expect_true(.ravel_constraint_fit_valid(fit$fit))
  expect_false(is.null(fit$interactions$table))

  s <- ravel_simulate_trios(N = 3000, m = 100, rho_MF = .2, parental_sibs = "maternal", seed = 4)
  s$age <- stats::runif(nrow(s), 4, 16)
  fs <- gxe_fit(s, interactions = "age", parental_sibs = "maternal")
  expect_true(fs$converged)
  df <- as.data.frame(fs, include_vanilla = FALSE)
  expect_true(all(c("delta_M", "beta_D:age") %in% df$parameter))
})

test_that("cluster-robust SEs widen family-level effects in long data", {
  d <- gxe_data()
  naive <- gxe_fit(d, fixed_effects = "age")
  clustered <- gxe_fit(d, fixed_effects = "age", cluster = "FID")
  expect_equal(clustered$beta_vis, naive$beta_vis, tolerance = 1e-6)
  expect_gt(clustered$lavaan_inference$visible$se[["bD"]], naive$lavaan_inference$visible$se[["bD"]])
  expect_null(clustered$interactions)
})

test_that("bootstrap resamples whole clusters", {
  d <- gxe_data(N = 1500L, seed = 5L)
  fit <- gxe_fit(d, interactions = "age", cluster = "FID", n_boot = 10, seed = 2)
  expect_true(fit$converged)
  expect_false(is.null(fit$bootstrap))
})

test_that("invalid interaction and cluster inputs are rejected", {
  d <- gxe_data(N = 1500L, seed = 5L)
  expect_error(gxe_fit(d, interactions = "height"), "Missing interaction covariate")
  expect_error(gxe_fit(d, interactions = "age", interaction_roles = "sibling"), "interaction_roles")
  expect_error(gxe_fit(d, interaction_center = c(age = 10)), "requires interactions")
  expect_error(gxe_fit(d, interactions = "age", interaction_center = c(height = 1)), "named by interaction")
  db <- d; db$Y <- as.integer(db$Y > 0)
  expect_error(gxe_fit(db, interactions = "age", outcome = "binary_liability"), "continuous outcome")
  expect_error(gxe_fit(d, cluster = "family"), "cluster must name")
  dn <- d; dn$FID[1] <- NA
  expect_error(gxe_fit(dn, cluster = "FID"), "missing values")
  dc <- d; dc$age <- as.character(dc$age)
  expect_error(gxe_fit(dc, interactions = "age"), "numeric")
})

# ---------------------------------------------------------------------------
# Total-liability correction of interactions
# ---------------------------------------------------------------------------

# Hidden variation (pi = 0.4) and two constant rounds of assortment; the
# phenotype has a genetic x age effect on the standardized total liabilities.
gxe_hidden <- function(N = 8000L, seed = 3L, sibs = "none") {
  key <- paste("hidden", N, seed, sibs, sep = "_")
  if (exists(key, gxe_cache, inherits = FALSE)) return(gxe_cache[[key]])
  z <- function(x) as.numeric(scale(x))
  d <- ravel_simulate_trios(N = N, m = 300, m_hidden = 300, hidden_fraction = .4, rho_MF = .4,
    grandparent_rho_MF = .4, parental_sibs = sibs, seed = seed)
  Tz <- cbind(z(d$G_child), z(d$G_mother), z(d$G_father))
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv))
  set.seed(seed)
  fam <- rep(seq_len(N), sample(1:4, N, replace = TRUE))
  a <- stats::runif(length(fam), 4, 16) - 10
  Y <- drop(Tz[fam, ] %*% c(.30, .10, .05)) + drop((Tz[fam, ] * a) %*% c(.02, .01, 0)) + .05 * a +
    rnorm(N, sd = sqrt(.3))[fam] + rnorm(length(fam), sd = sqrt(.6))
  out <- data.frame(d[fam, grep("^PGS|^FID", names(d), value = TRUE)], Y = Y, age = a + 10)
  rownames(out) <- NULL
  gxe_cache[[key]] <- out
  out
}

test_that("corrected interactions are the main-effect map applied to the interaction vector", {
  d <- gxe_hidden()
  fit <- gxe_fit(d, pi = .4, am_n_generations = 2L, interactions = "age",
    interaction_center = c(age = 10), cluster = "FID")
  expect_identical(fit$interactions$roles, c("child", "mother", "father"))
  g <- fit$interactions$theta
  lc <- fit$latent_cov_vis
  mapped <- .compute_scale_estimates(g, lc[["eta_m"]], lc[["eta_f"]], lc[["rhoMF"]], .4, "n_gen", 2L)
  ct <- fit$interactions$corrected
  s <- fit$interactions$scales[["age"]]
  expect_equal(ct$estimate[ct$scale == "corrected_std"], unname(mapped$beta_total_assumed_std) / s, tolerance = 1e-10)
  expect_equal(ct$estimate[ct$scale == "corrected_raw"], unname(mapped$beta_total_assumed) / s, tolerance = 1e-10)
  # Recovery of the total-liability truth; the visible coefficients are attenuated.
  std <- ct[ct$scale == "corrected_std", ]
  expect_true(all(abs(std$estimate - c(.02, .01, 0)) < 4 * std$se))
  expect_gt(std$estimate[1], fit$interactions$table$estimate[1])
  expect_true(all(abs(fit$interactions$independence$correlation) < .05))
  df <- as.data.frame(fit, include_vanilla = FALSE)
  expect_true(all(c("beta_D:age", "beta_M:age") %in% df$parameter[df$scale == "corrected_std"]))
  expect_false(any(df$scale == "corrected_raw"))
  expect_true("beta_D:age" %in% as.data.frame(fit, include_raw = TRUE)$parameter[
    as.data.frame(fit, include_raw = TRUE)$scale == "corrected_raw"])
  expect_output(print_simple(fit), "standardized total liability")
})

test_that("pi > 0 requires interaction terms for every family member", {
  d <- gxe_hidden()
  expect_error(gxe_fit(d, pi = .4, interactions = "age", interaction_roles = "child"),
    "every family member")
  expect_null(gxe_fit(gxe_data(N = 1500L, seed = 5L), interactions = "age")$interactions$corrected)
})

test_that("external pi uncertainty widens corrected interaction intervals", {
  d <- gxe_hidden()
  fixed <- gxe_fit(d, pi = .4, am_n_generations = 2L, interactions = "age", cluster = "FID")
  fit <- gxe_fit(d, pi = .4, pi_se = .05, am_n_generations = 2L, interactions = "age", cluster = "FID")
  a <- fixed$interactions$corrected; b <- fit$interactions$corrected
  expect_equal(b$estimate, a$estimate, tolerance = 1e-10)
  expect_true(all(b$se[b$scale == "corrected_std"] > a$se[a$scale == "corrected_std"]))
  expect_match(b$uncertainty_source[1], "pi")
})

test_that("parental-sibling fits correct all five interaction terms and contrasts", {
  d <- gxe_hidden(N = 6000L, seed = 4L, sibs = "both")
  fit <- gxe_fit(d, pi = .4, am_n_generations = 2L, interactions = "age", interaction_center = c(age = 10),
    cluster = "FID", parental_sibs = "both")
  expect_identical(fit$interactions$roles, c("child", "mother", "father", "mat_sib", "pat_sib"))
  ct <- fit$interactions$corrected
  std <- ct[ct$scale == "corrected_std", ]
  expect_identical(std$parameter, c("beta_D:age", "beta_M:age", "beta_F:age", "beta_MS:age",
    "beta_FS:age", "delta_M:age", "delta_F:age"))
  truth <- c(.02, .01, 0, 0, 0, .01, 0)
  expect_true(all(abs(std$estimate - truth) < 4 * std$se))
  vt <- fit$interactions$table
  expect_equal(vt$estimate[vt$parameter == "delta_M:age"],
    vt$estimate[vt$parameter == "beta_M:age"] - vt$estimate[vt$parameter == "beta_MS:age"], tolerance = 1e-10)
})

test_that("a moderator correlated with the genetic factors triggers a warning", {
  d <- gxe_hidden()
  d$age <- d$age + 3 * as.numeric(scale(d$PGS1_child))
  expect_warning(ravel_fit(d, pi = .4, am_n_generations = 2L, interactions = "age", cluster = "FID"),
    "assumes independence")
})

test_that("bootstrap resampling unit is reconciled with cluster", {
  d <- gxe_data(N = 600L, seed = 5L)
  # An absent or finer bootstrap unit must fail before any bootstrap fit.
  expect_error(gxe_fit(d, fixed_effects = "age", cluster = "FID", cluster_var = "typo",
                       n_boot = 3, seed = 2), "cluster_var must name a column")
  dr <- d; dr$row_id <- seq_len(nrow(dr))
  expect_error(gxe_fit(dr, fixed_effects = "age", cluster = "FID", cluster_var = "row_id",
                       n_boot = 3, seed = 2), "equal or nest")

  # Record the data handed to every lavaan fit.
  rec <- new.env(); rec$data <- list()
  trace("sem", where = asNamespace("lavaan"), print = FALSE,
        tracer = bquote(assign("data", c(get("data", .(rec)), list(data)), envir = .(rec))))
  on.exit(untrace("sem", where = asNamespace("lavaan")), add = TRUE)
  whole_families <- function(fit_data) {
    boot <- Filter(function(x) !identical(x, fit_data), rec$data)
    n0 <- table(fit_data$FID)
    length(boot) >= 3L && all(vapply(boot, function(x) {
      n <- table(x$FID); all(n %% n0[names(n)] == 0L)
    }, logical(1)))
  }
  # Default: cluster_var inherits cluster and whole families are resampled.
  gxe_fit(d, fixed_effects = "age", cluster = "FID", n_boot = 3, seed = 2)
  expect_true(whole_families(d))
  # A coarser unit that nests families is allowed and keeps families whole.
  ds <- d; ds$site <- (ds$FID - 1L) %/% 20L
  rec$data <- list()
  gxe_fit(ds, fixed_effects = "age", cluster = "FID", cluster_var = "site", n_boot = 3, seed = 2)
  expect_true(whole_families(ds))
})

gxe_cache <- new.env(parent = emptyenv())

# Gaussian latent trios measured at 1..max_occ ages (long format, one row per
# child-occasion), with child and maternal genetic x age effects.
gxe_data <- function(N = 4000L, max_occ = 4L, seed = 11L, rho = .2, rel = c(.5, .3),
                     gD = .02, gM = .01) {
  key <- paste(N, max_occ, seed, gD, gM, sep = "_")
  if (exists(key, gxe_cache, inherits = FALSE)) return(gxe_cache[[key]])
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv))
  set.seed(seed)
  z <- MASS::mvrnorm(N, c(0, 0), matrix(c(1, rho, rho, 1), 2))
  Gm <- z[, 1]; Gf <- z[, 2]; Gc <- (Gm + Gf) / 2 + sqrt(.5 * (1 - rho)) * rnorm(N)
  P <- list()
  for (k in 1:2) {
    em <- rnorm(N); ef <- rnorm(N); ec <- (em + ef) / 2 + sqrt(.5) * rnorm(N)
    P[[paste0("PGS", k, "_mother")]] <- sqrt(rel[k]) * Gm + sqrt(1 - rel[k]) * em
    P[[paste0("PGS", k, "_father")]] <- sqrt(rel[k]) * Gf + sqrt(1 - rel[k]) * ef
    P[[paste0("PGS", k, "_child")]] <- sqrt(rel[k]) * Gc + sqrt(1 - rel[k]) * ec
  }
  nocc <- sample(seq_len(max_occ), N, replace = TRUE)
  fam <- rep(seq_len(N), nocc)
  a <- stats::runif(length(fam), 4, 16) - 10
  Y <- (.30 + gD * a) * Gc[fam] + (.10 + gM * a) * Gm[fam] + .05 * Gf[fam] + .05 * a +
    rnorm(N, sd = sqrt(.3))[fam] + rnorm(N, sd = .03)[fam] * a + rnorm(length(fam), sd = sqrt(.6))
  out <- data.frame(FID = fam, as.data.frame(P)[fam, ], Y = Y, age = a + 10)
  rownames(out) <- NULL
  gxe_cache[[key]] <- out
  out
}

gxe_fit <- function(...) suppressWarnings(ravel_fit(...))

test_that("post-fit pi correction refuses incomplete interaction roles and flags uncorrected interactions", {
  d <- gxe_data(N = 1500L, seed = 5L)
  # Child-only interaction fitted at pi = 0, main effects reported away from the mean age.
  child_only <- gxe_fit(d, interactions = "age", interaction_center = c(age = 6), cluster = "FID")
  expect_identical(child_only$interactions$roles, "child")
  expect_error(ravel_correct_pi(child_only, pi = .4, pi_se = .05), "every family member")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(child_only, path)
  expect_error(ravel_correct_pi(path, pi = .4, pi_ci = c(.3, .5)), "every family member")
  # pi = 0 is the identity map and remains available.
  expect_s3_class(ravel_correct_pi(child_only, pi = 0, pi_se = 0), "data.frame")

  # Complete roles: main effects equal the direct fit at pi; interactions are flagged, not silently dropped.
  full <- gxe_fit(d, interactions = "age", interaction_roles = c("child", "mother", "father"),
    interaction_center = c(age = 6), cluster = "FID")
  expect_warning(cp <- ravel_correct_pi(full, pi = .4, pi_se = .05), "not corrected")
  expect_false(attr(cp, "interactions")$corrected)
  direct <- gxe_fit(d, pi = .4, interactions = "age", interaction_center = c(age = 6), cluster = "FID")
  main <- cp[cp$scale == "corrected_std" & cp$parameter %in% c("beta_D", "beta_M", "beta_P"), ]
  expect_equal(main$estimate, unname(direct$beta_total_assumed_std[c("bD", "bm", "bf")]), tolerance = 1e-8)

  # ravel_fit's own external-pi path and its display re-correction do not warn.
  msgs <- character()
  withCallingHandlers({
    with_se <- ravel_fit(d, pi = .4, pi_se = .05, interactions = "age", cluster = "FID")
    as.data.frame(with_se, include_vanilla = FALSE)
  }, warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("not corrected", msgs)))
  expect_identical(attr(with_se$pi_correction, "interactions")$corrected_in, "fit$interactions$corrected")
})

test_that("correlated external pi does not zero-fill the interaction/pi covariance", {
  d <- gxe_hidden()
  args <- list(d, pi = .4, pi_se = .05, am_n_generations = 2L, interactions = "age",
               interaction_center = c(age = 4), cluster = "FID")
  base <- do.call(gxe_fit, args)
  V <- lavaan::vcov(base$fit)
  core <- .report_core_labels
  g <- names(base$interactions$theta)
  # pi correlated with the core estimates along Cov(core, gD): a valid core/pi block (R^2 = .9)
  dir <- V[core, "gD_1"]
  cross <- setNames(dir * sqrt(.9 * .05^2 / drop(crossprod(dir, solve(V[core, core], dir)))), core)
  lab <- c(core, g, "pi")
  Vx <- matrix(0, length(lab), length(lab), dimnames = list(lab, lab))
  Vx[c(core, g), c(core, g)] <- V[c(core, g), c(core, g)]
  Vx["pi", "pi"] <- .05^2
  Vx[core, "pi"] <- Vx["pi", core] <- cross
  expect_gte(min(eigen(Vx[c(core, "pi"), c(core, "pi")], TRUE, TRUE)$values), 0)
  # ...whose zero-filled expansion is indefinite (the reviewer's counterexample)
  expect_lt(min(eigen(cov2cor(Vx), TRUE, TRUE)$values), -.01)
  expect_warning(fit <- ravel_fit(d, pi = .4, pi_se = .05, am_n_generations = 2L, interactions = "age",
    interaction_center = c(age = 4), cluster = "FID", core_pi_cov = cross),
    "interaction/pi covariance is not supplied")
  expect_true(all(fit$pi_correction$uncertainty_status == "ok"))  # main effects remain supported
  ct <- fit$interactions$corrected
  expect_equal(ct$estimate, base$interactions$corrected$estimate, tolerance = 1e-10)
  expect_true(all(is.na(ct$se) & is.na(ct$ci_low) & is.na(ct$p_value)))
  expect_match(ct$uncertainty_source[1], "unavailable")
  # An all-zero core_pi_cov is the independent case and keeps inference
  zero <- do.call(gxe_fit, c(args, list(core_pi_cov = setNames(numeric(6), core))))
  expect_equal(zero$interactions$corrected$se, base$interactions$corrected$se, tolerance = 1e-10)
})

test_that("an invalid or unavailable SEM covariance gives NA, not clamped, corrected interaction SEs", {
  fit <- gxe_fit(gxe_hidden(), pi = .4, am_n_generations = 2L, interactions = "age", cluster = "FID")
  est <- list(beta_vis = fit$beta_vis, latent_cov_vis = fit$latent_cov_vis, parental_sibs = NULL)
  call <- function(object, pi_se = NULL) .ravel_interaction_corrected(object, est, fit$interactions, .4,
    fit$am_correction, fit$am_n_generations, character(0), .95, pi_se = pi_se)
  expect_true(all(is.finite(call(fit$fit)$se)))
  bad <- fit$fit
  V <- lavaan::vcov(bad)
  i <- match(c("eta_m", "gD_1"), rownames(V))
  bad@vcov$vcov[i[1], i[2]] <- bad@vcov$vcov[i[2], i[1]] <- 3 * sqrt(V[i[1], i[1]] * V[i[2], i[2]])
  expect_lt(min(eigen(lavaan::vcov(bad), TRUE, TRUE)$values), 0)
  for (pse in list(NULL, .05)) {
    ct <- call(bad, pse)
    expect_true(all(is.na(ct$se)))
    expect_match(ct$uncertainty_source[1], "unavailable")
  }
  bad@Options$se <- "none"   # lavaan::vcov() now errors
  for (pse in list(NULL, .05)) {
    ct <- call(bad, pse)
    expect_true(all(is.na(ct$se)))
    expect_equal(ct$estimate, fit$interactions$corrected$estimate, tolerance = 1e-10)
  }
})

# Regression tests: interaction reporting honours the uncertainty selection.
# (Append to tests/testthat/test-interactions.R; uses gxe_hidden()/gxe_fit()/gxe_cache from that file.)

gxe_boot_fit <- function() {
  if (exists("boot_fit", gxe_cache, inherits = FALSE)) return(gxe_cache$boot_fit)
  gxe_cache$boot_fit <- gxe_fit(gxe_hidden(N = 2000L, seed = 6L), pi = .4, am_n_generations = 2L,
    interactions = "age", interaction_center = c(age = 10), cluster = "FID", n_boot = 12, seed = 3)
}
ix_lines <- function(text) grep(":age ", text, value = TRUE)

test_that("uncertainty = 'none' suppresses interaction intervals in both printers", {
  fit <- gxe_boot_fit()
  for (printer in list(print_simple, print_ravel_results)) {
    rows <- ix_lines(capture.output(printer(fit, uncertainty = "none")))
    expect_gt(length(rows), 0L)
    expect_false(any(grepl("[", rows, fixed = TRUE)))
  }
  # print_simple also shows the corrected interaction rows, without intervals.
  text <- capture.output(print_simple(fit, uncertainty = "none"))
  expect_true(any(grepl("standardized total liability", text)))
  expect_length(ix_lines(text), 6L)
  # Intervals are still shown when uncertainty is requested.
  expect_true(all(grepl("[", ix_lines(capture.output(print_simple(fit, uncertainty = "lavaan"))), fixed = TRUE)))
})

test_that("bootstrap selection uses joint bootstrap draws for interactions", {
  fit <- gxe_boot_fit()
  skip_if(is.null(fit$bootstrap))
  draws <- fit$bootstrap$interaction_draws
  expect_identical(colnames(draws), names(fit$interactions$theta))
  ok <- stats::complete.cases(draws)
  expect_identical(ok, stats::complete.cases(fit$bootstrap$core_draws))
  boot <- as.data.frame(fit, uncertainty = "bootstrap", include_raw = TRUE, include_vanilla = FALSE)
  lav <- as.data.frame(fit, uncertainty = "lavaan", include_raw = TRUE, include_vanilla = FALSE)
  gb <- grepl(":", boot$parameter); gl <- grepl(":", lav$parameter)
  expect_equal(boot$estimate[gb], lav$estimate[gl])
  expect_true(all(startsWith(boot$uncertainty_source[gb], "bootstrap")))
  expect_true(all(boot$p_value_method[gb] == "normal_approx_from_bootstrap_se"))
  expect_true(all(grepl("lavaan", lav$uncertainty_source[gl])))
  expect_true(all(lav$p_value_method[gl] == "normal_approx_from_lavaan_se"))

  s <- fit$interactions$scales[["age"]]
  gD <- draws[ok, "gD_1"] / s
  v <- boot[gb & boot$scale == "visible" & boot$parameter == "beta_D:age", ]
  expect_equal(v$se, stats::sd(gD), tolerance = 1e-10)
  expect_equal(c(v$ci_low, v$ci_high), unname(stats::quantile(gD, c(.025, .975))), tolerance = 1e-10)

  # Corrected rows: the main-effect map applied to each joint (core, g) replicate.
  core <- fit$bootstrap$core_draws[ok, , drop = FALSE]
  g <- draws[ok, c("gD_1", "gM_1", "gF_1"), drop = FALSE]
  mapped <- vapply(seq_len(nrow(g)), function(b) .compute_scale_estimates(g[b, ], core[b, "eta_m"],
    core[b, "eta_f"], core[b, "rhoMF"], .4, "n_gen", 2L)$beta_total_assumed_std[["bD"]] / s, numeric(1))
  cs <- boot[gb & boot$scale == "corrected_std" & boot$parameter == "beta_D:age", ]
  expect_equal(cs$se, stats::sd(mapped), tolerance = 1e-8)
  expect_equal(c(cs$ci_low, cs$ci_high), unname(stats::quantile(mapped, c(.025, .975))), tolerance = 1e-8)
  text <- capture.output(print_simple(fit, uncertainty = "bootstrap"))
  expect_true(any(grepl("visible scale; per SD of G per unit of covariate; bootstrap", text, fixed = TRUE)))
  expect_true(any(grepl("standardized total liability; bootstrap", text, fixed = TRUE)))
})

test_that("a bootstrap request without interaction draws is labelled as lavaan", {
  fit <- gxe_boot_fit()
  skip_if(is.null(fit$bootstrap))
  fit$interactions$bootstrap <- NULL
  fit$bootstrap$interaction_draws <- NULL
  boot <- as.data.frame(fit, uncertainty = "bootstrap", include_vanilla = FALSE)
  gx <- grepl(":", boot$parameter)
  expect_true(all(grepl("lavaan.*no bootstrap draws for interactions", boot$uncertainty_source[gx])))
  expect_true(all(boot$p_value_method[gx] == "normal_approx_from_lavaan_se"))
})

test_that("corrected interaction intervals are withdrawn when the external-pi correction failed", {
  fit <- gxe_boot_fit()
  fit$pi_correction_error <- "synthetic"  # set by ravel_fit for incompatible sibling/pi states
  df <- as.data.frame(fit, uncertainty = "lavaan", include_vanilla = FALSE)
  gx <- grepl(":", df$parameter)
  expect_true(all(is.na(df$se[gx & df$scale == "corrected_std"])))
  expect_true(all(is.na(df$ci_low[gx & df$scale == "corrected_std"])))
  expect_false(anyNA(df$se[gx & df$scale == "visible"]))
})

test_that("fixed parental paths enter corrected interaction SEs as constants", {
  d <- gxe_data(N = 3000L, seed = 5L)
  fit <- gxe_fit(d, pi = .4, fix_parent_paths = TRUE, interactions = "age", cluster = "FID")
  V <- lavaan::vcov(fit$fit)
  expect_false(any(c("bm", "bf") %in% rownames(V)))
  std <- fit$interactions$corrected
  std <- std[std$scale == "corrected_std", ]
  expect_true(all(is.finite(std$se)))
  expect_identical(std$uncertainty_source[1], "lavaan cluster-robust")
  # Independent delta method over the free parameters only: the map depends on
  # the interaction vector and the latent covariances, not on bm or bf.
  g <- names(fit$interactions$theta)
  free <- c("rhoMF", "eta_m", "eta_f", g)
  s <- fit$interactions$scales[["age"]]
  f <- function(v) unname(.compute_scale_estimates(setNames(v[g], c("bD", "bm", "bf")), v[["eta_m"]],
    v[["eta_f"]], v[["rhoMF"]], .4, fit$am_correction, fit$am_n_generations)$beta_total_assumed_std) / s
  v0 <- c(fit$latent_cov_vis[c("rhoMF", "eta_m", "eta_f")], fit$interactions$theta[g])
  J <- sapply(seq_along(v0), function(i) {
    h <- 1e-6 * max(1, abs(v0[[i]])); up <- dn <- v0
    up[i] <- up[i] + h; dn[i] <- dn[i] - h
    (f(up) - f(dn)) / (2 * h)
  })
  expect_equal(std$estimate, f(v0), tolerance = 1e-10)
  expect_equal(std$se, sqrt(diag(J %*% V[free, free] %*% t(J))), tolerance = 1e-5)
  # External pi uncertainty and the finite-generation refit sensitivity still apply.
  for (extra in list(list(), list(fit_constraints = "n_gen", am_n_generations = 2L))) {
    base <- do.call(gxe_fit, c(list(d, pi = .4, fix_parent_paths = TRUE, interactions = "age", cluster = "FID"), extra))
    withpi <- do.call(gxe_fit, c(list(d, pi = .4, pi_se = .05, fix_parent_paths = TRUE, interactions = "age",
      cluster = "FID"), extra))
    b <- base$interactions$corrected; w <- withpi$interactions$corrected
    expect_true(all(is.finite(w$se)))
    expect_true(all(w$se >= b$se - 1e-12))
    expect_match(w$uncertainty_source[1], "pi \\(delta-Wald\\)")
  }
})
