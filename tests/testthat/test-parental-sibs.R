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
  expect_error(sib_fit(d, parental_sibs = "both", fit_constraints = "equilibrium"),
    "fit_constraints = 'none'")
  expect_error(sib_fit(d, parental_sibs = "both", pi = .3), "pi = 0")
  expect_error(sib_fit(d, parental_sibs = "both", pi = 0, pi_se = .05), "pi = 0")
  expect_error(sib_fit(d, parental_sibs = "both", fix_parent_paths = TRUE), "fix_parent_paths")
  sparse <- d
  sparse$PGS1_mat_sib[-(1:30)] <- NA
  expect_error(sib_fit(sparse, parental_sibs = "maternal"), "at least 50 families")
  chr <- d
  chr$PGS1_mat_sib <- as.character(chr$PGS1_mat_sib)
  expect_error(sib_fit(chr, parental_sibs = "maternal"), "numeric")

  fit <- sib_fit(d, parental_sibs = "maternal")
  expect_error(ravel_correct_pi(fit, pi = .3, pi_se = .05), "parental siblings")
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
