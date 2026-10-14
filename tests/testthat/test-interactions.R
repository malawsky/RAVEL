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
