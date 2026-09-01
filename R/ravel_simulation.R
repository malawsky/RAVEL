# =============================================================================
# ravel_simulation.R
# =============================================================================
# RAVEL Simulation Framework for Trio Genetic Data
#
# Author: Daniel Malawsky (dm22@sanger.ac.uk)
# =============================================================================

#' Simulate Trio Data for RAVEL Analysis
#'
#' @description
#' Generates simulated family trio data following a realistic genetic
#' architecture including assortative mating, Mendelian inheritance, and
#' measurement error in polygenic scores.
#'
#' @details
#' The simulation follows this data generating process:
#'
#' 1. **SNP Effects**: m causal variants with effects gamma ~ N(0,1)
#'
#' 2. **Parental Genotypes**: Gm_geno, Gf_geno ~ Binomial(2, AF) for each SNP
#'
#' 3. **Assortative Mating**: Parents are paired via a copula method to achieve
#'    the target standardized correlation \eqn{\rho_T} between
#'    \eqn{T_m^*} and \eqn{T_p^*}.
#'
#' 4. **True Genetic Liabilities**: SNP effects are scaled once to the
#'    founder/base-population genotype variance. Raw total liabilities \eqn{T^*}
#'    are then
#'    formed before any role-specific standardization, so assortative-mating
#'    variance inflation is retained on the total-liability scale.
#'
#' 5. **PGS Weights**: gamma_PGS = gamma + epsilon, where epsilon ~ N(0, sigma^2_w)
#'    sigma^2_w calibrated to achieve target reliability.
#'    Same weights applied to all family members (shared GWAS).
#'
#' 6. **Mendelian Inheritance**: Child inherits one allele per locus from each parent
#'    P(allele=1 | parent genotype=g) = g/2
#'
#' 7. **Phenotype**: \eqn{Y = (\beta^{raw})^T T^* + \epsilon_Y}, where the
#'    legacy returned variables \code{G_mother}, \code{G_father}, and
#'    \code{G_child} contain the raw total liabilities \eqn{T_m^*},
#'    \eqn{T_p^*}, and \eqn{T_c^*}.
#'    - epsilon_Y ~ N(0, sigma^2_e) with sigma^2_e chosen so Var(Y) ~ 1
#'
#' @param N Integer. Number of trios to simulate.
#' @param m Integer. Number of causal SNPs.
#' @param allele_freq Numeric or vector. Allele frequency(ies) for SNPs.
#' @param h2_direct Nonnegative coefficient-setting parameter. The simulator sets
#'   \eqn{\beta_D^{raw}=\sqrt{h2_direct}}. This is not generally an isolated
#'   variance contribution when trio liabilities are correlated or their raw
#'   variances differ from one.
#' @param h2_mat Nonnegative coefficient-setting parameter with
#'   \eqn{\beta_M^{raw}=\sqrt{h2_mat}}; not generally an isolated maternal
#'   variance contribution.
#' @param h2_pat Nonnegative coefficient-setting parameter with
#'   \eqn{\beta_P^{raw}=\sqrt{h2_pat}}; not generally an isolated paternal
#'   variance contribution.
#' @param rho_MF Numeric. Target standardized spouse correlation \eqn{\rho_T}
#'   between parental raw total liabilities. The achieved finite-sample Pearson
#'   correlation is recorded in the simulation attributes and can differ
#'   slightly from the requested target.
#' @param PGS1_reliability Numeric. Target reliability of PGS1 (R^2 with true G).
#' @param PGS2_reliability Numeric. Target reliability of PGS2 (R^2 with true G).
#' @param hidden_fraction Numeric [0,1). The simulator's \eqn{\pi}: the founder
#'   or random-mating fraction assigned to hidden additive variants. When
#'   \code{hidden_fraction > 0}, the simulation generates two sets
#'   of causal variants: "visible" variants (included in PGS with noise) and
#'   "hidden" variants (contribute to true G but completely absent from PGS).
#'   \eqn{T^* = \sqrt{1-\pi}G^* + \sqrt{\pi}H^*}, and the phenotype depends
#'   on \eqn{T^*} while the PGS only measures the visible component.
#' @param m_hidden Integer. Number of hidden SNPs (default: same as m).
#' @param pop_structure Logical. If TRUE, simulate population structure with
#'   Fst-based allele frequency differentiation between populations. Assortative
#'   mating (if rho_MF != 0) is implemented *within each population* to avoid
#'   inducing artificial cross-pop couples.
#' @param n_pops Integer. Number of populations (default: 2).
#' @param Fst Numeric. Fixation index controlling allele frequency differentiation.
#'   Typical values: 0.001 (minimal), 0.01 (moderate), 0.05 (strong).
#' @param pop_effect_on_Y Numeric scalar or vector. Creates population-level
#'   differences in the phenotype (confounding).
#'
#'   * If scalar: the same mean shift is applied to populations 2..n_pops relative
#'     to population 1 (reference = 0).
#'   * If length (n_pops - 1): interpreted as the mean shifts for pops 2..n_pops.
#'   * If length n_pops: treated as absolute effects and re-referenced so pop 1 = 0.
#' @param seed Integer. Random seed for reproducibility.
#' @param return_true_G Logical. Include true G values in output?
#' @param return_visible_genotypes Logical. If TRUE, attach the visible-variant
#'   child, mother, and father genotype matrices as an attribute named
#'   \code{"visible_genotypes"}. This is primarily useful for downstream GRM-
#'   based analyses such as genetic relationship matrix comparisons.
#'
#' @return A data.frame with columns:
#'   - FID: Family ID
#'   - PGS1_mother, PGS2_mother: Mother's polygenic scores
#'   - PGS1_father, PGS2_father: Father's polygenic scores
#'   - PGS1_child, PGS2_child: Child's polygenic scores
#'   - Y: Phenotype
#'   - pop: Population indicator (if pop_structure = TRUE)
#'   - G_mother, G_father, G_child: Legacy names for raw total liabilities
#'     \eqn{T_m^*}, \eqn{T_p^*}, and \eqn{T_c^*} (if return_true_G = TRUE)
#'   - G_mother_vis, G_father_vis, G_child_vis: Role-standardized visible
#'     liabilities used for oracle comparisons when present
#'
#'   Attributes include simulation parameters for reference.
#'
#' @examples
#' \dontrun{
#' dat <- ravel_simulate_trios(
#'   N = 30000, m = 100,
#'   h2_direct = 0.20, h2_mat = 0.05, h2_pat = 0.00,
#'   rho_MF = 0.30,
#'   PGS1_reliability = 0.50, PGS2_reliability = 0.30,
#'   seed = 12345
#' )
#' print_simulation_summary(dat)
#' }
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm rbinom var cor
ravel_simulate_trios <- function(
    N = 30000,
    m = 100,
    allele_freq = 0.25,
    h2_direct = 0.20,
    h2_mat = 0.05,
    h2_pat = 0.00,
    rho_MF = 0.30,
    PGS1_reliability = 0.50,
    PGS2_reliability = 0.30,
    hidden_fraction = 0,
    m_hidden = NULL,
    pop_structure = FALSE,
    n_pops = 2,
    Fst = 0.01,
    pop_effect_on_Y = 0,
    seed = NULL,
    return_true_G = TRUE,
    return_visible_genotypes = FALSE
) {

  # Input validation
  stopifnot(N > 0, m > 0)
  stopifnot(h2_direct >= 0, h2_mat >= 0, h2_pat >= 0)
  stopifnot(h2_direct + h2_mat + h2_pat <= 1)
  stopifnot(rho_MF >= -1 && rho_MF <= 1)
  stopifnot(PGS1_reliability > 0 && PGS1_reliability < 1)
  stopifnot(PGS2_reliability > 0 && PGS2_reliability < 1)
  stopifnot(hidden_fraction >= 0 && hidden_fraction < 1)

  # Default: same number of hidden SNPs as visible
  if (is.null(m_hidden)) m_hidden <- m
  if (hidden_fraction > 0) stopifnot(m_hidden > 0)

  if (!is.null(seed)) set.seed(seed)

  # Expand allele frequency if scalar
  if (length(allele_freq) == 1) allele_freq <- rep(allele_freq, m)
  stopifnot(length(allele_freq) == m)
  stopifnot(all(allele_freq > 0 & allele_freq < 1))

  # Allele frequencies for hidden variants (use same default)
  allele_freq_hidden <- rep(mean(allele_freq), m_hidden)

  # True effect sizes
  beta_D <- sqrt(h2_direct)
  beta_M <- sqrt(h2_mat)
  beta_F <- sqrt(h2_pat)

  # Convenience
  h <- hidden_fraction
  sqrt_1mh <- sqrt(1 - h)
  sqrt_h <- sqrt(h)

  # Keep genetic effects on one fixed base-population scale. This avoids
  # restandardizing transmitted child liabilities after assortative mating,
  # which would otherwise erase AM-induced total-liability variance inflation.
  scale_effects_to_unit_var <- function(gamma, af) {
    genotype_var <- 2 * af * (1 - af)
    gamma / sqrt(sum(genotype_var * gamma^2))
  }

  centered_genetic_score <- function(geno, gamma, af) {
    as.numeric(geno %*% gamma) - sum(2 * af * gamma)
  }

  draw_weight_noise <- function(af, target_score_var) {
    z <- rnorm(length(af))
    genotype_var <- 2 * af * (1 - af)
    z / sqrt(sum(genotype_var * z^2)) * sqrt(target_score_var)
  }

  # ---------------------------------------------------------------------------
  # POPULATION STRUCTURE SETUP
  # ---------------------------------------------------------------------------
  if (pop_structure) {
    stopifnot(n_pops >= 2, Fst >= 0 && Fst < 0.5)

    # Assign trios to populations (roughly equal sizes, random order)
    pop_assignment <- rep(seq_len(n_pops), length.out = N)
    pop_assignment <- sample(pop_assignment, size = N, replace = FALSE)
    pop_sizes <- table(pop_assignment)

    # Population effects on Y (relative to pop 1 reference = 0)
    # - scalar: applies same effect to pops 2..n_pops
    # - length (n_pops - 1): effects for pops 2..n_pops
    # - length n_pops: treated as absolute effects and re-referenced to pop 1
    if (length(pop_effect_on_Y) == 1) {
      pop_effect_vec <- rep(0, n_pops)
      if (n_pops >= 2) pop_effect_vec[2:n_pops] <- pop_effect_on_Y
    } else if (length(pop_effect_on_Y) == (n_pops - 1)) {
      pop_effect_vec <- c(0, pop_effect_on_Y)
    } else if (length(pop_effect_on_Y) == n_pops) {
      pop_effect_vec <- pop_effect_on_Y - pop_effect_on_Y[1]
      pop_effect_vec[1] <- 0
    } else {
      stop("pop_effect_on_Y must be a scalar, length (n_pops - 1), or length n_pops.")
    }

    # Generate population-specific allele frequencies using Fst-based drift
    # Visible variants
    af_by_pop_vis <- matrix(NA, nrow = n_pops, ncol = m)
    af_by_pop_vis[1, ] <- allele_freq  # Pop 1 = ancestral

    for (pop in 2:n_pops) {
      drift_sd <- sqrt(Fst * allele_freq * (1 - allele_freq))
      af_by_pop_vis[pop, ] <- allele_freq + rnorm(m, sd = drift_sd)
      # Bound to valid range
      af_by_pop_vis[pop, ] <- pmin(pmax(af_by_pop_vis[pop, ], 0.01), 0.99)
    }

    # Hidden variants (if applicable)
    if (h > 0) {
      af_by_pop_hid <- matrix(NA, nrow = n_pops, ncol = m_hidden)
      af_by_pop_hid[1, ] <- allele_freq_hidden

      for (pop in 2:n_pops) {
        drift_sd <- sqrt(Fst * allele_freq_hidden * (1 - allele_freq_hidden))
        af_by_pop_hid[pop, ] <- allele_freq_hidden + rnorm(m_hidden, sd = drift_sd)
        af_by_pop_hid[pop, ] <- pmin(pmax(af_by_pop_hid[pop, ], 0.01), 0.99)
      }
    }
  } else {
    pop_assignment <- rep(1, N)
    n_pops <- 1
    pop_effect_vec <- 0
  }

  # ---------------------------------------------------------------------------
  # STEP 1: Generate SNP effect sizes
  # ---------------------------------------------------------------------------
  # Visible variants (contribute to both G* and PGS)
  gamma_vis <- scale_effects_to_unit_var(rnorm(m), allele_freq)

  # Hidden variants (contribute to G* but NOT to PGS)
  if (h > 0) {
    gamma_hid <- scale_effects_to_unit_var(rnorm(m_hidden), allele_freq_hidden)
  }

  # ---------------------------------------------------------------------------
  # STEP 2: Generate parental genotypes (before assortment)
  # ---------------------------------------------------------------------------
  if (pop_structure) {
    # Generate genotypes with population-specific allele frequencies
    Gm_geno_vis_orig <- matrix(NA, nrow = N, ncol = m)
    Gf_geno_vis_orig <- matrix(NA, nrow = N, ncol = m)

    for (pop in 1:n_pops) {
      idx <- which(pop_assignment == pop)
      n_pop <- length(idx)

      for (j in 1:m) {
        Gm_geno_vis_orig[idx, j] <- rbinom(n_pop, 2, af_by_pop_vis[pop, j])
        Gf_geno_vis_orig[idx, j] <- rbinom(n_pop, 2, af_by_pop_vis[pop, j])
      }
    }

    # Hidden variant genotypes (if h > 0)
    if (h > 0) {
      Gm_geno_hid_orig <- matrix(NA, nrow = N, ncol = m_hidden)
      Gf_geno_hid_orig <- matrix(NA, nrow = N, ncol = m_hidden)

      for (pop in 1:n_pops) {
        idx <- which(pop_assignment == pop)
        n_pop <- length(idx)

        for (j in 1:m_hidden) {
          Gm_geno_hid_orig[idx, j] <- rbinom(n_pop, 2, af_by_pop_hid[pop, j])
          Gf_geno_hid_orig[idx, j] <- rbinom(n_pop, 2, af_by_pop_hid[pop, j])
        }
      }
    }
  } else {
    # No population structure - original code
    Gm_geno_vis_orig <- matrix(rbinom(N * m, 2, rep(allele_freq, each = N)), nrow = N)
    Gf_geno_vis_orig <- matrix(rbinom(N * m, 2, rep(allele_freq, each = N)), nrow = N)

    # Hidden variant genotypes (if h > 0)
    if (h > 0) {
      Gm_geno_hid_orig <- matrix(rbinom(N * m_hidden, 2, rep(allele_freq_hidden, each = N)), nrow = N)
      Gf_geno_hid_orig <- matrix(rbinom(N * m_hidden, 2, rep(allele_freq_hidden, each = N)), nrow = N)
    }
  }

  # Compute TOTAL genetic value for assortment ranking
  Gm_vis_temp <- centered_genetic_score(Gm_geno_vis_orig, gamma_vis, allele_freq)
  Gf_vis_temp <- centered_genetic_score(Gf_geno_vis_orig, gamma_vis, allele_freq)

  if (h > 0) {
    Gm_hid_temp <- centered_genetic_score(Gm_geno_hid_orig, gamma_hid, allele_freq_hidden)
    Gf_hid_temp <- centered_genetic_score(Gf_geno_hid_orig, gamma_hid, allele_freq_hidden)

    # Total genetic value (for assortment)
    Gm_total_temp <- sqrt_1mh * Gm_vis_temp + sqrt_h * Gm_hid_temp
    Gf_total_temp <- sqrt_1mh * Gf_vis_temp + sqrt_h * Gf_hid_temp
  } else {
    Gm_total_temp <- Gm_vis_temp
    Gf_total_temp <- Gf_vis_temp
  }

  # ---------------------------------------------------------------------------
  # STEP 3: Assortative mating via Gaussian copula
  # ---------------------------------------------------------------------------
  if (abs(rho_MF) > 1e-10) {

    Sigma_assort <- matrix(c(1, rho_MF, rho_MF, 1), 2, 2)

    # Initialize (random mating), then reorder mothers/fathers to induce assortment
    Gm_geno_vis <- Gm_geno_vis_orig
    Gf_geno_vis <- Gf_geno_vis_orig
    if (h > 0) {
      Gm_geno_hid <- Gm_geno_hid_orig
      Gf_geno_hid <- Gf_geno_hid_orig
    }

    if (pop_structure) {
      # IMPORTANT: mate *within* population to preserve allele-frequency differentiation.
      # This avoids creating artificial cross-pop couples (which would blur stratification)
      # and keeps the trio's pop label meaningful for both parents.
      for (pop in 1:n_pops) {
        idx <- which(pop_assignment == pop)
        n_pop <- length(idx)
        if (n_pop <= 1) next

        z <- MASS::mvrnorm(n_pop, mu = c(0, 0), Sigma = Sigma_assort)

        # Sort within population by total genetic value
        mom_order <- idx[order(Gm_total_temp[idx])]
        dad_order <- idx[order(Gf_total_temp[idx])]

        # Copula ranks
        r1 <- rank(z[, 1], ties.method = "first")
        r2 <- rank(z[, 2], ties.method = "first")

        # Reorder visible genotypes within pop
        Gm_geno_vis[idx, ] <- Gm_geno_vis_orig[mom_order, , drop = FALSE][r1, , drop = FALSE]
        Gf_geno_vis[idx, ] <- Gf_geno_vis_orig[dad_order, , drop = FALSE][r2, , drop = FALSE]

        # Reorder hidden genotypes within pop (if applicable)
        if (h > 0) {
          Gm_geno_hid[idx, ] <- Gm_geno_hid_orig[mom_order, , drop = FALSE][r1, , drop = FALSE]
          Gf_geno_hid[idx, ] <- Gf_geno_hid_orig[dad_order, , drop = FALSE][r2, , drop = FALSE]
        }
      }
    } else {
      # Single population: global mating is equivalent to within-pop mating
      z <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma_assort)

      mom_order <- order(Gm_total_temp)
      dad_order <- order(Gf_total_temp)

      r1 <- rank(z[, 1], ties.method = "first")
      r2 <- rank(z[, 2], ties.method = "first")

      Gm_geno_vis <- Gm_geno_vis_orig[mom_order, , drop = FALSE][r1, , drop = FALSE]
      Gf_geno_vis <- Gf_geno_vis_orig[dad_order, , drop = FALSE][r2, , drop = FALSE]

      if (h > 0) {
        Gm_geno_hid <- Gm_geno_hid_orig[mom_order, , drop = FALSE][r1, , drop = FALSE]
        Gf_geno_hid <- Gf_geno_hid_orig[dad_order, , drop = FALSE][r2, , drop = FALSE]
      }
    }

  } else {
    # Random mating
    Gm_geno_vis <- Gm_geno_vis_orig
    Gf_geno_vis <- Gf_geno_vis_orig
    if (h > 0) {
      Gm_geno_hid <- Gm_geno_hid_orig
      Gf_geno_hid <- Gf_geno_hid_orig
    }
  }

  # ---------------------------------------------------------------------------
  # STEP 4: Compute parental genetic liabilities
  # ---------------------------------------------------------------------------
  # Visible component on the fixed raw scale, plus role-standardized versions
  # for the visible latent scale targeted by the SEM measurement model.
  Gm_vis_raw <- centered_genetic_score(Gm_geno_vis, gamma_vis, allele_freq)
  Gf_vis_raw <- centered_genetic_score(Gf_geno_vis, gamma_vis, allele_freq)
  Gm_vis <- as.numeric(scale(Gm_vis_raw))
  Gf_vis <- as.numeric(scale(Gf_vis_raw))

  # Hidden component and total liability on the fixed raw scale.
  if (h > 0) {
    Gm_hid_raw <- centered_genetic_score(Gm_geno_hid, gamma_hid, allele_freq_hidden)
    Gf_hid_raw <- centered_genetic_score(Gf_geno_hid, gamma_hid, allele_freq_hidden)

    # Total genetic liability: T = sqrt(1-h)*G_vis + sqrt(h)*G_hid.
    # This is intentionally not role-standardized after AM.
    Gm <- sqrt_1mh * Gm_vis_raw + sqrt_h * Gm_hid_raw
    Gf <- sqrt_1mh * Gf_vis_raw + sqrt_h * Gf_hid_raw
  } else {
    Gm <- Gm_vis_raw
    Gf <- Gf_vis_raw
  }

  achieved_rho <- cor(Gm, Gf)
  achieved_rho_vis <- cor(Gm_vis, Gf_vis)

  # ---------------------------------------------------------------------------
  # STEP 5: Generate PGS weights with GWAS estimation noise
  # ---------------------------------------------------------------------------
  sigma_w1 <- sqrt((1 - PGS1_reliability) / PGS1_reliability)
  sigma_w2 <- sqrt((1 - PGS2_reliability) / PGS2_reliability)

  gamma1 <- gamma_vis + draw_weight_noise(allele_freq, target_score_var = sigma_w1^2)
  gamma2 <- gamma_vis + draw_weight_noise(allele_freq, target_score_var = sigma_w2^2)

  # Compute and standardize parental PGS (only from visible variants)
  Gm1 <- as.numeric(scale(centered_genetic_score(Gm_geno_vis, gamma1, allele_freq)))
  Gm2 <- as.numeric(scale(centered_genetic_score(Gm_geno_vis, gamma2, allele_freq)))
  Gf1 <- as.numeric(scale(centered_genetic_score(Gf_geno_vis, gamma1, allele_freq)))
  Gf2 <- as.numeric(scale(centered_genetic_score(Gf_geno_vis, gamma2, allele_freq)))

  # ---------------------------------------------------------------------------
  # STEP 6: Mendelian inheritance to offspring
  # ---------------------------------------------------------------------------
  # Visible variants
  mat_alleles_vis <- matrix(rbinom(N * m, 1, as.vector(Gm_geno_vis) / 2), nrow = N)
  pat_alleles_vis <- matrix(rbinom(N * m, 1, as.vector(Gf_geno_vis) / 2), nrow = N)
  Gc_geno_vis <- mat_alleles_vis + pat_alleles_vis

  # Child's visible genetic liability on the fixed raw scale, plus the
  # role-standardized visible latent scale for oracle comparisons.
  Gc_vis_raw <- centered_genetic_score(Gc_geno_vis, gamma_vis, allele_freq)
  Gc_vis <- as.numeric(scale(Gc_vis_raw))

  # Hidden variants (if applicable)
  if (h > 0) {
    mat_alleles_hid <- matrix(rbinom(N * m_hidden, 1, as.vector(Gm_geno_hid) / 2), nrow = N)
    pat_alleles_hid <- matrix(rbinom(N * m_hidden, 1, as.vector(Gf_geno_hid) / 2), nrow = N)
    Gc_geno_hid <- mat_alleles_hid + pat_alleles_hid

    Gc_hid_raw <- centered_genetic_score(Gc_geno_hid, gamma_hid, allele_freq_hidden)

    # Total child genetic liability
    Gc <- sqrt_1mh * Gc_vis_raw + sqrt_h * Gc_hid_raw
  } else {
    Gc <- Gc_vis_raw
  }

  # Child's PGS (only from visible variants, using same GWAS weights)
  Gc1 <- as.numeric(scale(centered_genetic_score(Gc_geno_vis, gamma1, allele_freq)))
  Gc2 <- as.numeric(scale(centered_genetic_score(Gc_geno_vis, gamma2, allele_freq)))

  # ---------------------------------------------------------------------------
  # STEP 7: Generate phenotype
  # ---------------------------------------------------------------------------
  genetic_component <- beta_D * Gc + beta_M * Gm + beta_F * Gf

  # Add population effect if specified (creates confounding)
  if (pop_structure) {
    # pop_effect_vec is defined in the population-structure setup block
    if (any(abs(pop_effect_vec) > 1e-10)) {
      pop_component <- pop_effect_vec[pop_assignment]
      genetic_component <- genetic_component + pop_component
    }
  }

  genetic_var <- var(genetic_component)

  # Residual variance to achieve Var(Y) ~ 1
  var_e <- max(0.01, 1 - genetic_var)
  Y <- genetic_component + rnorm(N, sd = sqrt(var_e))

  # ---------------------------------------------------------------------------
  # STEP 8: Assemble output
  # ---------------------------------------------------------------------------
  dat <- data.frame(
    FID = seq_len(N),
    PGS1_mother = Gm1,
    PGS2_mother = Gm2,
    PGS1_father = Gf1,
    PGS2_father = Gf2,
    PGS1_child = Gc1,
    PGS2_child = Gc2,
    Y = Y,
    stringsAsFactors = FALSE
  )

  # Add population indicator if pop_structure is enabled
  if (pop_structure) {
    dat$pop <- pop_assignment
  }

  if (return_true_G) {
    # Total genetic liability (what affects phenotype)
    dat$G_mother <- Gm
    dat$G_father <- Gf
    dat$G_child <- Gc

    # Visible genetic liability (what PGS measures, if h > 0)
    if (h > 0) {
      dat$G_mother_vis <- Gm_vis
      dat$G_father_vis <- Gf_vis
      dat$G_child_vis <- Gc_vis
    }
  }

  if (isTRUE(return_visible_genotypes)) {
    attr(dat, "visible_genotypes") <- list(
      child = Gc_geno_vis,
      mother = Gm_geno_vis,
      father = Gf_geno_vis
    )
  }

  # Store simulation parameters as attributes
  attr(dat, "sim_params") <- list(
    N = N,
    m = m,
    m_hidden = if (h > 0) m_hidden else 0,
    h2_direct = h2_direct,
    h2_mat = h2_mat,
    h2_pat = h2_pat,
    beta_D = beta_D,
    beta_M = beta_M,
    beta_F = beta_F,
    rho_MF_target = rho_MF,
    rho_MF_achieved_total = achieved_rho,
    rho_MF_achieved_visible = if (h > 0) achieved_rho_vis else achieved_rho,
    PGS1_reliability_target = PGS1_reliability,
    PGS2_reliability_target = PGS2_reliability,
    PGS1_reliability_achieved_vs_visible = cor(Gm1, Gm_vis)^2,
    PGS2_reliability_achieved_vs_visible = cor(Gm2, Gm_vis)^2,
    hidden_fraction = h,
    return_visible_genotypes = isTRUE(return_visible_genotypes),
    pop_structure = pop_structure,
    n_pops = if (pop_structure) n_pops else 1,
    Fst = if (pop_structure) Fst else NA,
    pop_effect_on_Y = if (pop_structure) pop_effect_vec else 0,
    var_e = var_e,
    seed = seed
  )

  dat
}


#' Print Summary of Simulated Trio Data
#'
#' @param dat Data frame from \code{\link{ravel_simulate_trios}}
#'
#' @export
#' @importFrom stats cov
print_simulation_summary <- function(dat) {
  params <- attr(dat, "sim_params")

  cat("================================================================================\n")
  cat("RAVEL SIMULATION SUMMARY\n")
  cat("================================================================================\n\n")

  cat("Sample size:", params$N, "trios\n")
  cat("Visible SNPs:", params$m, "\n")
  if (params$hidden_fraction > 0) {
    cat("Hidden SNPs:", params$m_hidden, "\n")
    cat("Hidden fraction (h):", params$hidden_fraction, "\n")
  }
  if (params$pop_structure) {
    cat("Populations:", params$n_pops, "\n")
    cat("Fst:", params$Fst, "\n")
    if (any(abs(params$pop_effect_on_Y) > 1e-10)) {
      cat("Pop effects on Y (pop1=0):", paste(round(params$pop_effect_on_Y, 4), collapse = ", "), "(confounding)
")
    }
  }
  cat("\n")

  cat("True effect sizes (on TOTAL G):\n")
  cat("  beta_D (direct):", round(params$beta_D, 4),
      "(h2 =", round(params$h2_direct, 3), ")\n")
  cat("  beta_M (maternal):", round(params$beta_M, 4),
      "(h2 =", round(params$h2_mat, 3), ")\n")
  cat("  beta_F (paternal):", round(params$beta_F, 4),
      "(h2 =", round(params$h2_pat, 3), ")\n\n")

  cat("Assortative mating:\n")
  cat("  Target rho:", params$rho_MF_target, "\n")
  cat("  Achieved rho (total G):", round(params$rho_MF_achieved_total, 4), "\n")
  if (params$hidden_fraction > 0) {
    cat("  Achieved rho (visible G):", round(params$rho_MF_achieved_visible, 4), "\n")
  }
  cat("\n")

  cat("PGS reliability (vs visible G):\n")
  cat("  PGS1: target =", params$PGS1_reliability_target,
      ", achieved =", round(params$PGS1_reliability_achieved_vs_visible, 4), "\n")
  cat("  PGS2: target =", params$PGS2_reliability_target,
      ", achieved =", round(params$PGS2_reliability_achieved_vs_visible, 4), "\n\n")

  if ("G_mother" %in% names(dat)) {
    cat("Covariance structure (empirical, total G):\n")
    cat("  Cov(Gc, Gm):", round(cov(dat$G_child, dat$G_mother), 4),
        "(expected ~", round(0.5 * (1 + params$rho_MF_achieved_total), 4), ")\n")
    cat("  Cov(Gc, Gf):", round(cov(dat$G_child, dat$G_father), 4), "\n")

    if (params$hidden_fraction > 0 && "G_mother_vis" %in% names(dat)) {
      cat("\nCovariance structure (empirical, visible G only):\n")
      cat("  Cov(Gc_vis, Gm_vis):", round(cov(dat$G_child_vis, dat$G_mother_vis), 4), "\n")
      cat("  Cov(Gc_vis, Gf_vis):", round(cov(dat$G_child_vis, dat$G_father_vis), 4), "\n")
    }
  }

  # Population structure diagnostics
  if (params$pop_structure && "pop" %in% names(dat)) {
    cat("\nPopulation structure diagnostics:\n")
    pop_table <- table(dat$pop)
    cat("  Pop sizes:", paste(names(pop_table), "=", pop_table, collapse = ", "), "\n")

    # Mean PGS by population (shows differentiation)
    pgs_by_pop <- tapply(dat$PGS1_mother, dat$pop, mean)
    cat("  Mean PGS1 by pop:", paste(round(pgs_by_pop, 3), collapse = ", "), "\n")

    # Correlation between pop and PGS (proxy for stratification)
    if (params$n_pops == 2) {
      r_pop_pgs <- cor(dat$pop, dat$PGS1_mother)
      cat("  Cor(pop, PGS1):", round(r_pop_pgs, 4), "\n")
    }
  }

  cat("\n================================================================================\n")
}
