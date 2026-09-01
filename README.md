# RAVEL
## Recovering Assortment-Adjusted Variance Explained by Genetic Liabilities

RAVEL is a structural equation framework for estimating direct-offspring and
parental genetic-liability-correlated coefficients from polygenic scores in
parent-offspring trios.

**Author**: Daniel Malawsky (dm22@sanger.ac.uk), Wellcome Sanger Institute

---

## Overview

RAVEL is a trio structural equation framework for estimating direct offspring and parental genetic-
liability-correlated variance explained by leveraging two independently constructed polygenic scores
per individual as multiple noisy indicators of a common latent genetic liability.
Two independently constructed PGSs serve as noisy indicators of a
role-standardized visible genetic liability for each mother, father, and child.
The fitted SEM returns visible-scale coefficients
**b** = (b_D, b_M, b_P). Under explicit assumptions about additive transmission,
hidden genetic variation, and the common assortment pattern across visible and
hidden components, a post-estimation matrix transformation returns coefficients
on both the raw total-liability scale beta^raw and the role-standardized
total-liability scale beta.

- **b_D**: child visible-liability coefficient conditional on parental visible liabilities
- **b_M / b_P**: maternal and paternal visible-liability coefficients
- **rho_G**: visible-scale spouse correlation
- **eta_m / eta_p**: child-mother and child-father visible covariances
- **t_m / t_p**: standardized child-on-parent projection coefficients recovered from rho_G, eta_m, and eta_p

The measurement model corrects for PGS unreliability via ORIV-as-SEM, while
cross-person residual covariances account for shared GWAS estimation noise across
family members.

### Notation used in the theory note and R API

The supplementary theory note uses p/P for the paternal role, while the R API
retains f/F in field names for backward compatibility. Thus `bf`, `eta_f`, and
`tf` correspond to b_P, eta_p, and t_p, respectively. The legacy `rhoMF` field
is the mother-father entry.

The API argument `h` corresponds to the note's pi: the fraction of founder- or
random-mating additive variance assigned to the hidden component.

---

## Installation

```r
remotes::install_github("malawsky/RAVEL")

# or from a local checkout
devtools::install()
```

---

## Usage

### Required data format

A data frame with two PGS per trio member, an outcome, and optionally covariates. By default RAVEL expects columns named `PGS1_mother`, `PGS2_mother`, `PGS1_father`, `PGS2_father`, `PGS1_child`, `PGS2_child`, `Y` — but you can map your own column names:

```r
# Default column names
fit <- ravel_fit(dat, h = 0.5, n_boot = 1000)

# Custom column names
fit <- ravel_fit(dat,
  col_pgs1_mother = "EA1_PGS_mom",  col_pgs2_mother = "EA2_PGS_mom",
  col_pgs1_father = "EA1_PGS_dad",  col_pgs2_father = "EA2_PGS_dad",
  col_pgs1_child  = "EA1_PGS_kid",  col_pgs2_child  = "EA2_PGS_kid",
  col_outcome     = "phenotype",
  h = 0.5, n_boot = 1000
)

print_simple(fit)
```

### Choosing uncertainty source

RAVEL now returns two uncertainty summaries:

- `lavaan_inference`: SEM-based standard errors with normal-approximation CIs, including delta-method propagation to the hidden-corrected raw and standardized total-liability scales
- `bootstrap`: empirical bootstrap summaries when `n_boot > 0`

The print methods can report either source:

```r
fit <- ravel_fit(dat, h = 0.5, n_boot = 500)

print_simple(fit, uncertainty = "bootstrap")
print_simple(fit, uncertainty = "lavaan")
print_ravel_results(fit, uncertainty = "bootstrap")
```

### Fixed effects / covariates

Covariates can be included via `fixed_effects`, which residualizes Y on the specified columns within the SEM. Alternatively, you can residualize Y yourself before fitting (recommended):

```r
# Option 1: let RAVEL handle it
fit <- ravel_fit(dat, fixed_effects = c("sex", "PC1", "PC2"))

# Option 2: residualize beforehand (recommended)
dat$Y <- residuals(lm(Y ~ sex + PC1 + PC2 + PC3 + PC4, data = dat))
fit <- ravel_fit(dat)
```

### Fixing cross-person residual correlations

If you want to impose values on the observed cross-person residual correlations rather than estimate them freely, use `fixed_crossperson_corr` with `crossperson = "free"`:

```r
fit <- ravel_fit(
  dat,
  crossperson = "free",
  fixed_crossperson_corr = c(
    e1_mf = 0, e1_mc = 0.5, e1_fc = 0.5,
    e2_mf = 0, e2_mc = 0.5, e2_fc = 0.5
  )
)
```

---

## `ravel_fit` options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dat` | — | Data frame with 6 PGS columns + `Y` |
| `h` | `0` | Hidden additive-variance fraction pi on the founder or random-mating reference scale. At `h = 0` there is no hidden component. |
| `outcome` | `"continuous"` | `"continuous"` or `"binary_liability"` (probit/WLSMV) |
| `measurement` | `"free_loadings"` | Only `"free_loadings"` is currently supported |
| `crossperson` | `"free"` | `"free"` (6 residual covariances) or `"equal_pc"` (constrain parent-child) |
| `fixed_crossperson` | `NULL` | Optional named vector fixing selected cross-person residual covariances |
| `fixed_crossperson_corr` | `NULL` | Optional named vector fixing selected cross-person residual correlations (`crossperson = "free"` only) |
| `am_correction` | `"n_gen"` | Raw assortative-mating state used for the total-liability mapping: finite n-generation recursion or equilibrium. |
| `am_n_generations` | `1L` | Number of mating rounds in the constant-strength finite-generation model. For `n = 1` the mating parents are on the founder scale and the child is the first generation produced under assortment. |
| `fixed_effects` | `NULL` | Character vector of covariate column names to residualize Y |
| `n_boot` | `0` | Number of bootstrap replicates for 95% CIs |
| `cluster_var` | `"FID"` | Column name for cluster bootstrap; if missing from `dat`, bootstrap falls back to row resampling |
| `ci_level` | `0.95` | Confidence level for bootstrap intervals |
| `seed` | `NULL` | Optional random seed |
| `verbose` | `FALSE` | Print model syntax and fitting details |
| `col_pgs1_mother` | `"PGS1_mother"` | Column name for mother's PGS1 |
| `col_pgs2_mother` | `"PGS2_mother"` | Column name for mother's PGS2 |
| `col_pgs1_father` | `"PGS1_father"` | Column name for father's PGS1 |
| `col_pgs2_father` | `"PGS2_father"` | Column name for father's PGS2 |
| `col_pgs1_child` | `"PGS1_child"` | Column name for child's PGS1 |
| `col_pgs2_child` | `"PGS2_child"` | Column name for child's PGS2 |
| `col_outcome` | `"Y"` | Column name for outcome variable |

### Output fields

| Field | Description |
|-------|-------------|
| `beta_vis` | Named vector `bD`, `bm`, `bf`, corresponding to (b_D,b_M,b_P) on the role-standardized visible-liability scale. |
| `beta_total_assumed` | Coefficients beta^raw on the raw founder-scaled total liabilities T*. |
| `beta_total_assumed_std` | Coefficients beta = D_T beta^raw on role-standardized total liabilities T. |
| `latent_cov_vis` | Entries `rhoMF`, `eta_m`, `eta_f`, corresponding to rho_G, eta_m, and eta_p. Because G_c, G_m, and G_p have unit variance, these covariances are numerically correlations. |
| `latent_cov_total_assumed` | Raw covariances among T_c*, T_m*, and T_p*. The spouse entry is kappa_T, despite retaining the legacy field name `rhoMF`. |
| `latent_corr_total_assumed` | Correlations among standardized total liabilities T_c, T_m, and T_p. The spouse entry is rho_T. |
| `crossperson_resid` | Cross-person residual covariances on the observed PGS scale |
| `crossperson_resid_corr` | Cross-person residual correlations on the observed PGS scale |
| `var_total_assumed` | Role-specific V_T,r = Var(T_r*) for child, mother, and father. Standardized T has variance one. |
| `reliability` | PGS reliability per role and PGS |
| `determinacy` | Factor score determinacy per role |
| `vanilla_regression` | Plain standardized trio regressions on observed `PGS1` and `PGS2` (`lm` for continuous outcomes, probit `glm` for binary outcomes) |
| `R2_model` | Variance explained by the visible latent linear predictor. |
| `R2_total_assumed` | Variance explained by the assumed total-liability linear predictor, conditional on the supplied hidden-variance and AM model. |
| `lavaan_inference` | lavaan-based SEs and delta-method CIs for visible and corrected estimates |
| `bootstrap` | CIs and SEs if `n_boot > 0` |
| `converged` | Logical convergence flag |


### Assortative mating correction for the total-liability scale

RAVEL assumes that the aggregate cross-mate covariance projected into the
modeled visible-hidden component space has a common direction aligned with the
components' covariances with the total liability. Other assortment axes may
exist, provided that their projections onto the modeled visible and hidden
components are proportional to their projection onto T*.

The `am_correction` parameter controls the raw covariance state used in the
visible-to-total mapping, including raw AM variance inflation when `h = 0`:

- **`"n_gen"`** (default): Propagates the raw visible-hidden covariance state for a finite number of mating rounds under a constant standardized total-liability spouse correlation. With `am_n_generations = 1`, the mating parents are founders, so the parental within-person visible-hidden covariance is zero. The child generated by that mating round can nevertheless have nonzero visible-hidden covariance.

- **`"equilibrium"`**: Uses the fixed-point raw covariance state after indefinitely repeated assortment at constant strength.

```r
# One mating round; founder parents have d_parent = 0
fit <- ravel_fit(dat, h = 0.5, am_correction = "n_gen", am_n_generations = 1)

# 15 generations of AM
fit <- ravel_fit(dat, h = 0.5, am_correction = "n_gen", am_n_generations = 15)

# Equilibrium (infinite generations)
fit <- ravel_fit(dat, h = 0.5, am_correction = "equilibrium")
```

### Analytical expected values

`ravel_expected_values()` computes the population quantities analytically implied by a set of total-liability parameters under a specified AM history. This is useful for verifying what RAVEL should recover in expectation, checking simulation calibration, and understanding the scale transformations between visible and total liabilities.

```r
ev <- ravel_expected_values(
  beta_total = c(bD = 0.45, bm = 0.20, bf = 0.00),
  beta_total_scale = "raw",
  h = 0.4,
  rho_T_corr = 0.3,
  tm = 0.5, tf = 0.5,
  am_correction = "equilibrium"
)
ev$beta_vis_expected          # Expected visible-scale b
ev$beta_total_recovered       # Recovered raw beta^raw
ev$beta_total_recovered_std   # Recovered standardized beta
ev$Sigma_xh                   # Code name for the note's Sigma_GH; X is visible G
ev$var_total_expected         # Var(T*) by role, not Var(T)
```

With `beta_total_scale = "raw"`, `beta_total` is beta^raw and multiplies T*.
With `beta_total_scale = "standardized"`, it is beta and multiplies the
role-standardized T. The argument `rho_T_corr` is rho_T, the standardized
total-liability spouse correlation. In this analytical function, `tm` and `tf`
are raw maternal and paternal transmission coefficients; `tf` corresponds to
the note's paternal t_p notation. Both remain `0.5` under the maintained
additive model.

`ravel_expected_pgs_ols()` extends this to the observed-score scale. By
default, `reliability` is interpreted on the random-mating visible-scale,
matching the simulation convention where PGS noise is calibrated before
assortative mating inflates role-specific visible signal variance. Use
`reliability_scale = "observed"` if the values are already final reliabilities
for standardized PGSs (`loading^2 = reliability`). By default it assumes
spousal residual correlation `0` and parent-child residual correlations `0.5`.

```r
pgs_ev <- ravel_expected_pgs_ols(
  beta_total = c(bD = 0.20, bm = 0.05, bf = 0.00),
  beta_total_scale = "raw",
  h = 0.8,
  rho_T_corr = 0.7,
  reliability = c(child = 0.20, mother = 0.20, father = 0.20),
  reliability_scale = "random_mating",
  am_correction = "equilibrium"
)

pgs_ev$beta_pgs_expected   # expected standardized trio OLS on observed PGS
pgs_ev$Sigma_pp            # expected observed-score covariance block
```

### Print functions

- `print_ravel_results(fit)` — full detailed output
- `print_simple(fit)` — compact table with visible + corrected estimates and CIs
- `ravel_to_dataframe(fit, uncertainty = "lavaan")` — tidy table with SEM rows plus the vanilla `PGS1`/`PGS2` comparator regressions

---

## `ravel_simulate_trios` options

The native simulator generates the phenotype from raw founder-scaled total
liabilities:

Y = (beta^raw)^T T* + epsilon_Y.

The corresponding role-standardized target is beta = D_T beta^raw.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N` | — | Number of trios |
| `m` | — | Number of visible causal SNPs |
| `h2_direct` | — | Nonnegative coefficient-setting parameter. The simulator sets beta_D^raw = sqrt(h2_direct). This is not generally an isolated variance contribution when trio liabilities are correlated or their raw variances differ from one. |
| `h2_mat` | `0` | Nonnegative coefficient-setting parameter with beta_M^raw = sqrt(h2_mat); not generally an isolated maternal variance contribution. |
| `h2_pat` | `0` | Nonnegative coefficient-setting parameter with beta_P^raw = sqrt(h2_pat); not generally an isolated paternal variance contribution. |
| `rho_MF` | `0` | Target standardized spouse correlation rho_T between the parental raw total liabilities. The achieved finite-sample Pearson correlation is recorded in the simulation attributes and can differ slightly from the requested target. |
| `PGS1_reliability` | `0.50` | Target reliability of PGS1 |
| `PGS2_reliability` | `0.30` | Target reliability of PGS2 |
| `hidden_fraction` | `0` | The simulator's pi: the founder/random-mating fraction assigned to hidden additive variants. |
| `pop_structure` | `FALSE` | Simulate population stratification |
| `Fst` | `0.01` | Fixation index for allele frequency differentiation |
| `pop_effect_on_Y` | `0` | Direct population effect on phenotype (confounding) |
| `seed` | `NULL` | Random seed |
| `return_true_G` | `FALSE` | Include true genetic values in output |
| `return_visible_genotypes` | `FALSE` | Attach visible genotype matrices as a `visible_genotypes` attribute for downstream GRM analyses |

`G_mother`, `G_father`, and `G_child` are legacy variable names for the raw
total liabilities T_m*, T_p*, and T_c* that enter phenotype generation. When
present, `G_mother_vis`, `G_father_vis`, and `G_child_vis` are role-standardized
visible liabilities used for oracle comparisons. These legacy returned-column
names are retained for backward compatibility.

---

## Requirements

- R >= 4.0
- [lavaan](https://lavaan.ugent.be/) >= 0.6-12
- MASS

---

## License

MIT

---

## Citation

```
Malawsky D. (2025). RAVEL: Recovering Assortment-Adjusted Variance Explained by Genetic Liabilities. R package.
https://github.com/malawsky/RAVEL
```
