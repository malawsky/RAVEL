sib_cache <- new.env(parent = emptyenv())

sib_data <- function(parental_sibs = "both", N = 30000L, seed = 2024L, rho = .2,
                     lineage_mat = .10, lineage_pat = 0) {
  key <- paste(parental_sibs, N, seed, rho, lineage_mat, lineage_pat, sep = "_")
  if (!exists(key, sib_cache, inherits = FALSE)) {
    sib_cache[[key]] <- ravel_simulate_trios(
      N = N, m = 200, h2_direct = .09, h2_mat = .0225, h2_pat = .0025,
      rho_MF = rho, parental_sibs = parental_sibs,
      lineage_effect_mat = lineage_mat, lineage_effect_pat = lineage_pat, seed = seed)
  }
  sib_cache[[key]]
}

sib_fit <- function(...) suppressWarnings(ravel_fit(...))

test_that("default trio model and simulator output are unchanged by the sibling option", {
  expect_identical(.build_model(), .build_model(parental_sibs = "none"))
  expect_false(grepl("Gms|Gfs", .build_model()))

  a <- ravel_simulate_trios(N = 500, m = 40, seed = 9)
  b <- ravel_simulate_trios(N = 500, m = 40, seed = 9, parental_sibs = "none")
  expect_identical(a, b)
  expect_false(any(grepl("_sib$|_lineage$", names(a))))

  fit <- sib_fit(a)
  expect_null(fit$parental_sibs)
  expect_false(any(c("beta_MS", "delta_M") %in% as.data.frame(fit)$parameter))
})

test_that("simulator adds the requested siblings and lineage liabilities", {
  d <- ravel_simulate_trios(N = 2000, m = 60, seed = 3, parental_sibs = "maternal")
  expect_true(all(c("PGS1_mat_sib", "PGS2_mat_sib", "G_mat_sib", "T_mat_lineage") %in% names(d)))
  expect_false(any(grepl("pat_sib", names(d))))
  expect_identical(attr(d, "sim_params")$parental_sibs, "maternal")

  d <- sib_data("both")
  # Full siblings share half their genome; mother-to-child is also one half.
  expect_equal(cor(d$G_mother, d$G_mat_sib), .5, tolerance = .05)
  expect_equal(cor(d$G_father, d$G_pat_sib), .5, tolerance = .05)
  # Child and aunt/uncle share a quarter (plus a little through assortment).
  expect_equal(cor(d$G_child, d$G_mat_sib), .25, tolerance = .06)
  # The two parental siblings correlate only through the parents' assortment:
  # rho * (1 - 1/(2V))^2 = rho / 4 for founder parents.
  expect_equal(cor(d$G_mat_sib, d$G_pat_sib), .2 / 4, tolerance = .03)

  dh <- ravel_simulate_trios(N = 1500, m = 40, m_hidden = 40, hidden_fraction = .3,
    seed = 4, parental_sibs = "paternal")
  expect_true(all(c("PGS1_pat_sib", "G_pat_sib", "G_pat_sib_vis") %in% names(dh)))
  dp <- ravel_simulate_trios(N = 1500, m = 40, seed = 5, pop_structure = TRUE,
    parental_sibs = "both")
  expect_true(all(c("PGS1_mat_sib", "PGS1_pat_sib", "pop") %in% names(dp)))
})

test_that("each sibling option fits with the expected degrees of freedom", {
  d <- sib_data("both")
  for (opt in c("maternal", "paternal", "both")) {
    fit <- sib_fit(d, parental_sibs = opt)
    expect_true(fit$converged)
    expect_equal(as.numeric(lavaan::fitMeasures(fit$fit, "df")), if (opt == "both") 10 else 6)
    expect_identical(fit$parental_sibs$mode, opt)
    labels <- if (opt == "maternal") c("bms", "delta_m") else if (opt == "paternal")
      c("bfs", "delta_f") else c("bms", "delta_m", "bfs", "delta_f")
    expect_named(fit$parental_sibs$coefficients, labels)
    co <- fit$parental_sibs$coefficients
    if ("delta_m" %in% labels)
      expect_equal(co[["delta_m"]], fit$beta_vis[["bm"]] - co[["bms"]], tolerance = 1e-8)
    if ("delta_f" %in% labels)
      expect_equal(co[["delta_f"]], fit$beta_vis[["bf"]] - co[["bfs"]], tolerance = 1e-8)
  }
})

test_that("parent-specific contrasts remove maternal lineage confounding", {
  d <- sib_data("both")
  fit <- sib_fit(d, parental_sibs = "both")
  trio <- sib_fit(d)
  tab <- fit$parental_sibs$lavaan_table
  get <- function(p, col = "estimate") tab[tab$parameter == p, col]

  # Parents are founders at h = 0, so true visible coefficients are the raw ones.
  expect_lt(abs(get("delta_M") - .15), 4 * get("delta_M", "se"))
  expect_lt(abs(get("delta_F") - .05), 4 * get("delta_F", "se"))
  expect_lt(abs(get("beta_FS")), 4 * get("beta_FS", "se"))
  expect_gt(get("beta_MS"), 2 * get("beta_MS", "se"))
  # Without siblings the lineage effect loads onto the maternal coefficient.
  expect_gt(trio$beta_vis[["bm"]] - .15, 2 * trio$lavaan_inference$visible$se[["bm"]])

  rc <- fit$parental_sibs$resid_corr
  expect_equal(unname(rc[c("e1_m_ms", "e2_m_ms", "e1_f_fs", "e2_f_fs")]), rep(.5, 4), tolerance = .1)
  expect_equal(unname(rc[c("e1_c_ms", "e2_c_ms", "e1_c_fs", "e2_c_fs")]), rep(.25, 4), tolerance = .1)
  expect_equal(unname(rc[c("e1_f_ms", "e1_m_fs", "e1_ms_fs")]), rep(0, 3), tolerance = .1)
  expect_equal(unname(fit$parental_sibs$latent_cov[c("kap_m", "kap_f")]), c(.5, .5), tolerance = .08)
})

test_that("families without a sibling are retained by FIML", {
  d <- sib_data("both")
  set.seed(77)
  drop_m <- sample(nrow(d), round(.4 * nrow(d)))
  drop_f <- sample(nrow(d), round(.6 * nrow(d)))
  d[drop_m, c("PGS1_mat_sib", "PGS2_mat_sib")] <- NA
  d[drop_f, c("PGS1_pat_sib", "PGS2_pat_sib")] <- NA

  full <- sib_fit(sib_data("both"), parental_sibs = "both")
  fit <- sib_fit(d, parental_sibs = "both")
  expect_true(fit$converged)
  expect_equal(fit$N, nrow(d))
  expect_equal(unname(fit$parental_sibs$n_with_sib), c(.6, .4) * nrow(d))
  se_full <- full$parental_sibs$lavaan_table$se
  se_miss <- fit$parental_sibs$lavaan_table$se
  expect_true(all(se_miss > se_full))
  expect_lt(abs(fit$parental_sibs$coefficients[["delta_m"]] - .15),
    4 * fit$parental_sibs$lavaan_table$se[fit$parental_sibs$lavaan_table$label == "delta_m"])
})

test_that("custom sibling column names are mapped", {
  d <- sib_data("both")[1:6000, ]
  names(d)[names(d) == "PGS1_mat_sib"] <- "aunt_score_a"
  names(d)[names(d) == "PGS2_mat_sib"] <- "aunt_score_b"
  fit <- sib_fit(d, parental_sibs = "maternal",
    col_pgs1_mat_sib = "aunt_score_a", col_pgs2_mat_sib = "aunt_score_b")
  expect_true(fit$converged)
  expect_true("beta_MS" %in% as.data.frame(fit)$parameter)
})

test_that("unsupported combinations and bad sibling data are rejected", {
  d <- sib_data("both")[1:4000, ]
  expect_error(sib_fit(d[setdiff(names(d), "PGS2_pat_sib")], parental_sibs = "paternal"),
    "PGS2_pat_sib")
  expect_error(sib_fit(d, parental_sibs = "both", fix_parent_paths = TRUE), "fix_parent_paths")
  expect_error(sib_fit(d, parental_sibs = "both", pi = .3, pi_se = .05,
    core_pi_cov = c(bD = 0, bm = 0, bf = 0, rhoMF = 0, eta_m = 0, eta_f = 0)), "all core parameters")
  sparse <- d
  sparse$PGS1_mat_sib[-(1:30)] <- NA
  expect_error(sib_fit(sparse, parental_sibs = "maternal"), "at least 50 families")
  chr <- d
  chr$PGS1_mat_sib <- as.character(chr$PGS1_mat_sib)
  expect_error(sib_fit(chr, parental_sibs = "maternal"), "numeric")

})

test_that("sibling results are reported and bootstrapped", {
  d <- sib_data("both")[1:5000, ]
  fit <- sib_fit(d, parental_sibs = "both", n_boot = 12, seed = 5)
  expect_false(is.null(fit$bootstrap$parental_sibs))
  expect_named(fit$bootstrap$parental_sibs$se, c("bms", "delta_m", "bfs", "delta_f"))

  boot_df <- as.data.frame(fit, uncertainty = "bootstrap", include_vanilla = FALSE)
  lav_df <- as.data.frame(fit, uncertainty = "lavaan", include_vanilla = FALSE)
  rows <- c("beta_MS", "delta_M", "beta_FS", "delta_F")
  expect_true(all(rows %in% boot_df$parameter))
  expect_true(all(grepl("bootstrap", boot_df$uncertainty_source[boot_df$parameter %in% rows])))
  expect_equal(boot_df$estimate[boot_df$parameter %in% rows],
    lav_df$estimate[lav_df$parameter %in% rows])
  expect_identical(attr(lav_df, "settings")$parental_sibs, "both")

  expect_output(print_simple(fit), "Parental siblings")
  expect_output(print_simple(fit), "delta_M")
  out <- capture.output(print_ravel_results(fit, uncertainty = "lavaan"))
  expect_true(any(grepl("conditional on the parental siblings", out)))
  expect_true(any(grepl("e1_m_ms", out)))
  expect_true(any(grepl("mat_sib", out)))
  expect_true(all(c("mat_sib", "pat_sib") %in% fit$reliability$role))
  expect_true(all(c("mat_sib", "pat_sib") %in% names(fit$determinacy)))
})

# ---------------------------------------------------------------------------
# Total-liability correction with parental siblings
# ---------------------------------------------------------------------------

# Exact five-person raw covariance of Z = (G*, H*) under assortment on T*
# (rank-one spouse block), founder scale, parents after g_par rounds.
sib_theory <- function(pi, rhoT, g_par) {
  w <- c(sqrt(1 - pi), sqrt(pi)); I <- diag(2)
  V <- 1; for (g in seq_len(g_par)) V <- (1 + rhoT) / 2 * V + .5
  S <- I + (V - 1) * tcrossprod(w); K <- rhoT * V * tcrossprod(w)
  A <- S - I / 2; k <- 1 - 1 / (2 * V); X <- k * K
  B <- list(c_c = (S + K + I) / 2, c_m = (S + K) / 2, c_f = (S + K) / 2, c_ms = (A + X) / 2,
    c_fs = (A + X) / 2, m_m = S, m_f = K, m_ms = A, m_fs = X, f_f = S, f_ms = X, f_fs = A,
    ms_ms = S, ms_fs = k^2 * K, fs_fs = S)
  r <- c("c", "m", "f", "ms", "fs"); Sig <- matrix(0, 10, 10)
  for (i in 1:5) for (j in i:5) {
    blk <- B[[paste0(r[i], "_", r[j])]]
    Sig[2 * i - 1:0, 2 * j - 1:0] <- blk; Sig[2 * j - 1:0, 2 * i - 1:0] <- t(blk)
  }
  nm <- paste0(rep(c("G", "H"), 5), rep(r, each = 2)); dimnames(Sig) <- list(nm, nm)
  Sig
}

test_that("the sibling visible-hidden cross-block matches the analytic five-person covariance", {
  for (cfg in list(c(pi = .4, rhoT = .4, g = 1), c(pi = .4, rhoT = .4, g = 2), c(pi = .25, rhoT = .3, g = 6))) {
    Sg <- sib_theory(cfg[["pi"]], cfg[["rhoT"]], cfg[["g"]])
    gi <- paste0("G", c("c", "m", "f", "ms", "fs")); hi <- paste0("H", c("c", "m", "f", "ms", "fs"))
    R <- cov2cor(Sg[gi, gi]); dimnames(R) <- list(c("c", "m", "f", "ms", "fs"), c("c", "m", "f", "ms", "fs"))
    GH <- diag(1 / sqrt(diag(Sg)[gi])) %*% Sg[gi, hi] %*% diag(1 / sqrt(diag(Sg)[hi]))
    sib <- list(roles = c("ms", "fs"), b = c(bms = .02, bfs = 0), cov = c(kap_m = R["m", "ms"],
      chi_ms = R["f", "ms"], eta_ms = R["c", "ms"], kap_f = R["f", "fs"], chi_fs = R["m", "fs"],
      eta_fs = R["c", "fs"], zeta = R["ms", "fs"]))
    est <- .compute_scale_estimates(c(bD = .3, bm = .1, bf = .05), R["c", "m"], R["c", "f"],
      R["m", "f"], cfg[["pi"]], "n_gen", as.integer(cfg[["g"]] + 1), sib = sib)
    expect_equal(unname(est$sib_total$Sigma_xh), unname(GH), tolerance = 1e-10)
    expect_equal(est$hidden_correction_details$rho_T_corr_assumed, cfg[["rhoT"]], tolerance = 1e-8)
    # The note's equilibrium shortcut (q for every cross-lineage entry) is not used.
    expect_lt(est$sib_total$Sigma_xh["Gm", "Hfs"], est$sib_total$Sigma_xh["Gm", "Hf"])
  }
  # At pi = 0 the corrected standardized coefficients equal the visible ones.
  est0 <- .compute_scale_estimates(c(bD = .3, bm = .1, bf = .05), .6, .6, .2, 0, "n_gen", 1L,
    sib = list(roles = "ms", b = c(bms = .03), cov = c(kap_m = .5, chi_ms = .1, eta_ms = .29)))
  expect_equal(unname(est0$sib_total$std), c(.03, .07), tolerance = 1e-10)
  expect_equal(unname(est0$beta_total_assumed_std), c(.3, .1, .05), tolerance = 1e-10)
})

sib_hidden_data <- function(N = 30000L, seed = 31L) {
  key <- paste("hidden", N, seed, sep = "_")
  if (!exists(key, sib_cache, inherits = FALSE)) {
    sib_cache[[key]] <- ravel_simulate_trios(N = N, m = 300, m_hidden = 300, hidden_fraction = .4,
      h2_direct = .09, h2_mat = .0225, h2_pat = .0025, rho_MF = .3, grandparent_rho_MF = .3,
      parental_sibs = "both", lineage_effect_mat = .10, seed = seed)
  }
  sib_cache[[key]]
}

test_that("grandparental assortment gives parents one prior round of assortment", {
  d <- sib_hidden_data()
  expect_equal(var(d$G_mother), 1 + .3 / 2, tolerance = .04)
  expect_identical(attr(d, "sim_params")$grandparent_rho_MF, .3)
})

test_that("corrected parent-specific contrasts recover the total-liability truth", {
  d <- sib_hidden_data()
  truth_M <- .15 * sd(d$G_mother)
  fit <- sib_fit(d, pi = .4, parental_sibs = "both", am_n_generations = 2L)
  trio <- sib_fit(d, pi = .4, am_n_generations = 2L)
  df <- as.data.frame(fit, include_vanilla = FALSE, include_raw = TRUE)
  row <- df[df$scale == "corrected_std" & df$parameter == "delta_M", ]
  expect_equal(nrow(row), 1L)
  expect_lt(abs(row$estimate - truth_M), 4 * row$se)
  expect_true(all(c("beta_MS", "delta_M", "beta_FS", "delta_F") %in%
    df$parameter[df$scale == "corrected_raw"]))
  # The trio-only maternal coefficient absorbs the maternal lineage effect.
  expect_gt(trio$beta_total_assumed_std[["bm"]] - truth_M,
    2 * trio$lavaan_inference$total_std$se[["bm"]])
  # Raw and standardized contrasts differ by the parents' total-liability SD.
  co <- fit$parental_sibs$corrected
  expect_equal(co$std[["delta_m"]], co$raw[["delta_m"]] * sqrt(co$var_total[["mother"]]), tolerance = 1e-10)
  # Variance decomposition includes the sibling terms.
  expect_true(all(c("mat_sib", "pat_sib", "cross_mother_mat_sib") %in% names(fit$components_vis)))
  expect_true(all(c("mat_sib", "cross_father_pat_sib") %in% names(fit$components_total_assumed)))
  expect_output(print_simple(fit), "standardized total liability")
})

test_that("constrained sibling models satisfy the relatedness restrictions", {
  d <- sib_hidden_data()
  fit <- sib_fit(d, pi = .4, parental_sibs = "both", fit_constraints = "estimated_state")
  expect_true(fit$converged)
  expect_true(.ravel_constraint_fit_valid(fit$fit))
  rc <- fit$parental_sibs$resid_corr
  expect_equal(unname(rc[c("e1_m_ms", "e2_f_fs")]), c(.5, .5), tolerance = 1e-6)
  expect_equal(unname(rc[c("e1_c_ms", "e2_c_fs")]), c(.25, .25), tolerance = 1e-6)
  expect_equal(unname(rc[c("e1_f_ms", "e2_m_fs", "e1_ms_fs")]), c(0, 0, 0))
  pe <- lavaan::parameterEstimates(fit$fit)
  val <- function(l) pe$est[!is.na(pe$label) & pe$label == l][1]
  B <- fit$hidden_correction_details$estimated_state$B
  lc <- fit$parental_sibs$latent_cov
  expect_equal(lc[["kap_m"]], (3 + val("rhoMF") - 2 * B) / 2, tolerance = 1e-5)
  expect_equal(lc[["kap_f"]], lc[["kap_m"]], tolerance = 1e-6)
  expect_equal(lc[["chi_ms"]], lc[["chi_fs"]], tolerance = 1e-6)
  expect_equal(lc[["eta_ms"]], (lc[["kap_m"]] + lc[["chi_ms"]]) / (2 * sqrt(B)), tolerance = 1e-5)
  rel <- fit$reliability
  expect_equal(rel$reliability[rel$role == "mat_sib"], rel$reliability[rel$role == "mother"], tolerance = 1e-5)
  truth_M <- .15 * sd(d$G_mother)
  co <- fit$parental_sibs$corrected$std
  se <- fit$lavaan_inference$sib_total_std$se
  expect_lt(abs(co[["delta_m"]] - truth_M), 4 * se[["delta_m"]])
})

test_that("pi uncertainty propagates to the sibling contrasts", {
  d <- sib_hidden_data()[1:12000, ]
  fixed <- sib_fit(d, pi = .4, parental_sibs = "maternal", am_n_generations = 2L)
  fit <- sib_fit(d, pi = .4, pi_se = .05, parental_sibs = "maternal", am_n_generations = 2L)
  tab <- fit$pi_correction
  expect_true(all(c("beta_MS", "delta_M") %in% tab$parameter[tab$scale == "corrected_std"]))
  se_fixed <- fixed$lavaan_inference$sib_total_std$se[["delta_m"]]
  se_pi <- tab$se[tab$scale == "corrected_std" & tab$parameter == "delta_M"]
  expect_gt(se_pi, se_fixed)
  df <- as.data.frame(fit, include_vanilla = FALSE)
  row <- df[df$scale == "corrected_std" & df$parameter == "delta_M", ]
  expect_equal(row$se, se_pi, tolerance = 1e-10)
  expect_identical(row$inference, "lavaan_pi_delta_wald")
  # The same propagation is available after fitting, with matching estimates.
  again <- ravel_correct_pi(fixed, pi = .4, pi_se = .05)
  expect_equal(again$estimate, tab$estimate, tolerance = 1e-10)
  expect_true(all(c("bms", "kap_m", "chi_ms", "eta_ms") %in% rownames(attr(again, "joint_input_covariance"))))
})

test_that("bootstrap draws carry the sibling core for corrected inference", {
  d <- sib_hidden_data()[1:5000, ]
  fit <- sib_fit(d, pi = .4, parental_sibs = "maternal", am_n_generations = 2L, n_boot = 12, seed = 3)
  expect_true(all(c("bms", "kap_m", "chi_ms", "eta_ms") %in% colnames(fit$bootstrap$core_draws)))
  expect_false(is.null(fit$bootstrap$parental_sibs_total_std))
  boot_df <- as.data.frame(fit, uncertainty = "bootstrap", include_vanilla = FALSE)
  rows <- boot_df[boot_df$scale == "corrected_std" & boot_df$parameter %in% c("beta_MS", "delta_M"), ]
  expect_equal(nrow(rows), 2L)
  expect_true(all(is.finite(rows$se)))
  corrected <- ravel_correct_pi(fit, pi = .4, pi_se = .05, uncertainty = "bootstrap")
  expect_true(all(c("beta_MS", "delta_M") %in% corrected$parameter))
})

test_that("pi-dependent finite-generation sibling fits refit for pi sensitivity", {
  skip_on_cran()
  d <- sib_hidden_data()[1:8000, ]
  fit <- sib_fit(d, pi = .4, pi_se = .05, parental_sibs = "maternal",
    fit_constraints = "n_gen", am_n_generations = 2L)
  expect_true(fit$converged)
  sens <- attr(fit$pi_correction, "pi_uncertainty")$core_pi_sensitivity
  expect_true(all(c("bD", "bms", "kap_m", "chi_ms", "eta_ms") %in% names(sens)))
  expect_true(all(is.finite(sens)))
  expect_identical(fit$constraint_refit$parental_sibs, "maternal")
})
