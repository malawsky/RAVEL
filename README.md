# RAVEL
## Recovering Assortment-Adjusted Variance Explained by Genetic Liabilities

RAVEL is a structural equation framework for estimating direct-offspring and
parental genetic-liability-correlated coefficients from polygenic scores in
parent-offspring trios.

**Author**: Daniel Malawsky (d.malawsky@gmail.com), Wellcome Sanger Institute

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

### Optional inheritance constraints

The default fit remains unchanged. To impose stronger measurement and
inheritance assumptions during fitting, use:

```r
fit_eq <- ravel_fit(dat, h = 0.4, fit_constraints = "equilibrium")
fit_g3 <- ravel_fit(dat, h = 0.4, fit_constraints = "n_gen",
                   am_n_generations = 3L)
```

All constrained modes retain a separate unknown reliability for each PGS. They fix
within-PGS residual spouse correlation to zero and child-parent correlations
to 0.5, and assume independent errors across PGSs. The scores measure the same
visible genetic component and their errors are orthogonal to that component
and the outcome. Scores should be oriented toward the same genetic target.
Unknown marginal score variances remain nuisance parameters, so separately
standardizing the scores does not supply known population variances.

At equilibrium, each score's reliability is equal across family roles and
both latent child-parent correlations equal `(1 + rho_visible) / 2`. The two
scores' reliabilities are not equated to each other. With two scores and no
covariates there are 13 covariance parameters, compared with 19 when only the
residual correlations are fixed.

For finite generations, child reliability changes with the modeled genetic
variance. The model assumes random-mating founders, Mendelian transmission,
and a constant total-liability spouse correlation across mating rounds. The
generation count is supplied, not estimated; generation one means founder
parents and their first offspring generation. This is a specified finite
assortment history, not arbitrary time-varying assortment. Fits after generation
one depend on the supplied founder hidden fraction `h`. Changing that fraction
or the assumed history requires refitting.

With the same GWAS weights and constant raw error variance, reliability evolves
with the visible signal variance. Let `R_j,n` be the reliability of score j in
generation n, and `B_n = V_G,n+1 / V_G,n` the visible-variance growth factor. Then

```
R_j,n+1 = B_n * R_j,n / (1 + (B_n - 1) * R_j,n)
R_j,n+1 - R_j,n = (B_n - 1) * R_j,n * (1 - R_j,n) /
                  (1 + (B_n - 1) * R_j,n)
```

The two scores share the growth factor, but retain their own reliabilities.
With generation zero denoting random-mating founders, constant total-liability
spouse correlation r, and founder hidden fraction pi, the model gives
`V_G,n = pi + (1 - pi) * (1 - r * ((1 + r)/2)^n) / (1 - r)`.
Positive assortment increases reliability toward equilibrium; negative
assortment decreases it. Changes tend to zero at equilibrium. Separately
standardizing each observed score to variance one leaves reliability unchanged.
The API's `am_n_generations = n + 1` compares parents in generation n with
children in generation n + 1.

The same option is available in `ravel_dataframe_report()`. Constrained modes
select the fixed residual model when `residual_models` is omitted and reject
an explicitly requested free residual model. Report summaries remain on the
visible scale, with the fitting assumptions recorded in the correction packet.

```r
out_eq <- ravel_dataframe_report(df, fit_constraints = "equilibrium")
out_g3 <- ravel_dataframe_report(df, fit_constraints = "n_gen",
                                h = 0.4, am_n_generations = 3L)
```

These restrictions can improve precision when their assumptions hold. The
package's correction uses a supplied hidden fraction and estimates the total
spouse correlation. The earlier pooled-score comparison instead conditioned
on a known total spouse correlation; its numerical SE gains are not a promise
for this different uncertainty budget.

### Estimating the variance state without specifying generations

Use `estimated_state` when the generation count is unknown:

```r
fit_state <- ravel_fit(dat, h = 0.4, fit_constraints = "estimated_state")
fit_state$hidden_correction_details$estimated_state
corrected_state <- ravel_correct_pi(fit_state, pi = 0.4, pi_se = 0.02)
attr(corrected_state, "estimated_state")$summary
```

This model estimates one more covariance parameter than the equilibrium or
specified-generation modes: 14 parameters with two scores and no covariates.
The two PGSs retain separate unknown reliabilities. Parent-child correlations
identify the signal variance growth, even if each observed PGS was separately
standardized. Let `rho` be the fitted visible spouse correlation and `eta` the
common fitted child-parent correlation. On the scale where founder visible
variance is one,

```
B = ((1 + rho) / (2 * eta))^2
V_G,parent = 1 / (2 * B - 1 - rho)
V_G,child = B * V_G,parent
V_T,parent = (V_G,parent - pi) / (1 - pi)
rho_total = rho * V_G,parent / (V_G,parent - pi)
```

The founder anchor comes from Mendelian segregation variance equal to half the
founder variance. Correction also assumes random-mating founders with independent
visible and hidden components, and assortment along the same total genetic
target over time. Founder pi remains external. Pre-existing founder LD, a changing
genetic target, or a different segregation model are not identified by these data.

The estimated-state model is a smooth relaxation of a constant-strength history:
it estimates the accumulated variance and current assortment rather than selecting
an integer number of generations. It also accommodates changing assortment strength
when the target and segregation assumptions remain valid. Diagnostics indicate
whether the fitted state lies on a constant-strength founder-to-equilibrium
trajectory and, when defined, its effective duration. A fractional effective
duration is not an integer-generation estimate. Near equilibrium or random mating,
duration can be unidentifiable even when correction is estimable. The fit does not
force noisy estimates onto the founder or equilibrium boundary.

Fitting this mode is independent of pi; an incompatible supplied pi is rejected
by post-estimation correction. The full fitted covariance propagates uncertainty
in growth and founder scaling, and `ravel_correct_pi()` adds external pi uncertainty
without refitting. Its `estimated_state` attribute includes a table of estimates,
SEs and intervals for growth, parental visible and total variance, and the total
spouse correlation, together with their joint covariance. It does not report a
Wald SE for effective duration. Precision can decrease relative to a correctly specified
generation count. Low score signal and nearly singular reconstructed states can
make inference unstable. The model does not promise identifiable separate
reliabilities when both assortment and all score-outcome associations vanish.
Estimated-state fits explicitly check the rank of the observed-moment derivative
along the fitted constraint surface and the full parameter covariance before
allowing inference. Nonidentified fits and active inequality boundaries are
flagged; their correction packets cannot be used for ordinary Wald inference.

`ravel_dataframe_report(..., fit_constraints = "estimated_state")` supports the
same mode and retains it in the aggregate packet. Correct these packets with
`ravel_correct_packet()` or `ravel_correct_pi()`, or sweep `hidden_fraction` using
`ravel_correct_grid()` with `n_generations` omitted. The grid then uses the
estimated state for every cell and records the generation count as `NA`.
Supplying a generation grid for an estimated-state fit is rejected because
that would change the fitting assumptions.

### Corrected estimates with uncertainty in external pi

Use the separate `ravel_correct_pi()` extension when `h` is estimated externally,
supplying its standard error directly on the raw pi scale:

```r
# pi is the same founder hidden fraction called h by ravel_fit().
corrected <- ravel_correct_pi(fit_eq, pi = 0.4, pi_se = 0.05)
corrected[, c("scale", "parameter", "estimate", "se", "ci_low", "ci_high")]
attr(corrected, "covariance")

# The finite-generation case also propagates how changing pi changes the fit.
corrected_g3 <- ravel_correct_pi(fit_g3, pi = 0.4, pi_se = 0.05)
```

The external variance is `pi_se^2`. No pi Monte Carlo draws or CI-to-SE
conversion are needed. Supply exactly one of `pi_se` or `pi_ci`; `pi_se = 0`
adds no external uncertainty. `ci_level` controls the returned effect intervals.
The method uses a local delta approximation and does not clip confidence
intervals to parameter bounds.

For compatibility, `pi_ci = c(lower, upper)` is also accepted. By default it
is interpreted as a 95% normal-Wald interval and must be symmetric about the
estimate. For a logit-Wald interval, use
`pi_ci_scale = "logit"`; the interval must be symmetric on that scale. An
arbitrary percentile or asymmetric interval does not uniquely specify a
sampling variance and is rejected. `pi_ci_level` describes the supplied
interval. These CI-specific settings apply only with `pi_ci`. Zero external
variance reduces model-based inference to the fixed-pi delta method.

The external pi estimate is assumed independent of the target-trio estimator
unless `core_pi_cov` supplies their covariance. For pi-dependent fits, this
must describe covariance with SEM estimation noise at fixed pi; supplying the
already marginalized plug-in core covariance would double-count pi uncertainty.
The joint propagation retains covariance among all six visible core parameters.
For a live finite-generation fit it also uses perturb-and-refit sensitivity,
so the result includes uncertainty in both fitting and correction. The pi
point estimate must match the value used in that finite fit.

Visible-transfer packets are accepted too. Pi-dependent aggregate packets
need a separately supplied, validated `core_pi_sensitivity`, because they do
not contain individual data for refitting; use the live fit for automatic
sensitivity. Bootstrap inputs use the covariance of complete joint core draws
and return explicitly labeled delta-Wald intervals, not percentile intervals.
Existing fixed-h fitting, correction and reporting APIs retain their meaning.

### Choosing a fitting mode

| `fit_constraints` | Assortment state | External inputs for fitting | External inputs for total correction |
|---|---|---|---|
| `"none"` (default) | Existing unconstrained measurement and latent correlations | None | `h` and a finite history or equilibrium |
| `"equilibrium"` | Equilibrium transmission and equal within-score reliability across roles | None | Founder `h` |
| `"n_gen"` | Constant strength from random-mating founders | Generation count; `h` after generation one | Matching generation count; fitted `h` after generation one |
| `"estimated_state"` | Current growth and accumulated variance inferred from the trio | None | Founder `h` |

The three constrained modes share the residual-error assumptions described above
and estimate separate PGS reliabilities. `h` and `pi` name the same **founder
hidden genetic variance fraction**, not GWAS prediction error and not the
fraction of the outcome variance that is genetic. Larger `pi` means more
uncaptured genetic variation. These are opt-in models: imposing incorrect
restrictions can introduce bias, and more constraints do not guarantee a better
result in a misspecified model.

### Reporting and exporting the extensions

All fitting modes use the existing fit, table, print, report and bootstrap
interfaces. Fit summaries name the constraint and correction mode; estimated-state
summaries also show growth and founder-scaled variances. Table conversion retains
its existing column names (including `beta_F` for the paternal effect), with
assumptions in `attr(table, "settings")`, `attr(table, "constraint_spec")`, and
state diagnostics in `attr(table, "estimated_state")`.

```r
print_simple(fit_state)
fit_table <- as.data.frame(fit_state, include_raw = TRUE)
attr(fit_table, "settings")

# Correction is a separate step: a fixed-h fit table does not include pi SE.
corrected <- ravel_correct_pi(fit_state, pi = .4, pi_se = .02)
print_simple(corrected)                 # print_ravel_results(corrected) also works
print(ravel_plot_report(corrected))     # corrected_std by default
print(ravel_plot_report(corrected, scale = "corrected_raw"))

# A report uses the same constrained fitting mode and selected uncertainty.
state_report <- ravel_dataframe_report(df, fit_constraints = "estimated_state",
                                      h = .4, out_dir = "ravel_state_report")
state_corrected <- ravel_correct_pi(state_report$packet, pi = .4, pi_se = .02)
state_grid <- ravel_correct_grid(state_report, hidden_fraction = c(.3, .4, .5),
                                 pi_se = .02)  # no generation count
print(ravel_plot_report(state_grid))

saveRDS(state_corrected, "corrected_with_pi_uncertainty.rds")
write.csv(state_corrected, "corrected_with_pi_uncertainty.csv", row.names = FALSE)
write.csv(attr(state_corrected, "estimated_state")$summary,
          "estimated_variance_state.csv", row.names = FALSE)
```

`ravel_dataframe_report()` summaries and its automatic plots remain visible-scale
comparisons with the two ordinary PGS regressions. Corrected tables use `beta_P`
for the paternal effect; plotting accepts both paternal labels. Corrected plots
select one predictor scale at a time, retain scenario labels, and distinguish
π-adjusted delta-Wald intervals from fixed-π bootstrap percentile intervals.
Binary outcomes retain their SEM liability units; plain probit regression rows
are displayed separately. Live continuous fits retain the supplied outcome units;
the report workflow uses its prepared, standardized phenotype.

`include_raw = TRUE` returns founder-scale raw effects even at `h = 0`.
With assortment these generally differ from effects per current genetic SD,
although visible and standardized total effects coincide when `h = 0`.
Bootstrap aliases now follow their named scale in this case as well.

Omitting `n_generations` in a correction grid inherits a constrained report's
fitted mode and generation count. For an unconstrained report it retains the
legacy 15-generation correction default. Explicit incompatible cells are flagged;
a pi-dependent finite-generation fit cannot be reused for a sweep over pi.
Supplying `pi_se` asks for independent external pi uncertainty in each grid cell.
Finite-generation packets after generation one additionally require validated
`core_pi_sensitivity`; use `ravel_correct_pi()` on the live fit for automatic
perturb-and-refit sensitivity. Use that function directly if external pi and
trio estimation errors are correlated and `core_pi_cov` is needed.

CSV preserves row-level estimates, uncertainty and scenario labels. RDS also
preserves joint covariance matrices, derivative matrices, state diagnostics and
assumptions stored as attributes. Use RDS for further inference. The report's
`model_status` and `diagnostics` retain identification and fit failure reasons;
failed cells and missing intervals must be inspected before interpretation.

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

Covariates can be included jointly via `fixed_effects`. The SEM allows genetic
factors to correlate with these covariates. Residualizing only Y before fitting
is not equivalent when covariates and genetic predictors correlate; binary Y
must remain coded 0/1 rather than residualized.

With a binary outcome and covariates, a joint probit regression on all six
scores and the encoded covariates supplies marginal variance-one liability
moments. RAVEL fits the same SEM to those moments with positive working
weights and a full joint sandwich covariance for robust SEs. This avoids
singular first-stage information from dummy covariates. It assumes a linear,
homoscedastic conditional-probit liability model; it does not require the
covariates themselves to be Gaussian. Continuous fits and binary fits without
covariates retain their existing estimators.
For this binary-covariate path, baseline-relative fit indices such as CFI/TLI
are marked unavailable rather than interpreted using an unvalidated
independence baseline. See the report diagnostics for fit-statistic limitations.

```r
fit <- ravel_fit(dat, fixed_effects = c("sex", "PC1", "PC2"))
```

### One-call dataframe report

Supply raw scores, outcome and covariates, with one independent trio per row.
Default score names are `pgs1_child`, `pgs1_mother`, `pgs1_father` and their
`pgs2_` counterparts. Categorical covariates are treatment-dummy coded and
incomplete trios are excluded. Both fixed (parent-child 0.5, spouse 0) and
free within-version score-residual correlations are fitted by default.

```r
out <- ravel_dataframe_report(
  df, outcome = "Y", covariates = c("age", "sex", "PC1", "PC2"),
  out_dir = "results", report_name = "analysis_1"
)
out$summary
print(out$plot)
```

RAVEL inputs use each score version's child mean/SD across all roles; parental
variance ratios are preserved. Vanilla trio regressions use each role's own
score SD. Continuous Y is standardized and covariates are fitted jointly.
For binary outcomes use `outcome_type = "binary_liability"`; vanilla probit
and SEM liability units are labelled and plotted separately.

The default reports model-based SEs, CIs and p-values for coefficients and the
visible spouse correlation. Add `n_boot = 500, seed = 123` to select whole-trio
bootstrap SEs and percentile CIs, retaining the original point estimates.
Bootstrap p-values use a normal approximation with the empirical SE. Failed
draws are recorded and insufficient successful draws produce NA uncertainty,
not a silent model-based fallback. Original model-based inference remains in
`out$summary_model_based`. Descriptive moments, reliability and fit indices do
not all carry SEs. See `?ravel_dataframe_report` for the complete specification.

The optional new output folder contains CSV tables, an aggregate-only RDS,
and a PDF/PNG diagnostic plot; existing report folders are not overwritten.
Send `ravel_report.rds`, not just the summary CSV, for later corrections.
`ravel_correct_packet(out$packets$fixed, h = 0.4, ...)` propagates the full
joint parameter covariance or stored bootstrap draws to total-liability
coefficients. Joint-covariate corrections require explicit confirmation that
the mapping applies conditional on those covariates.

Inspect `out$model_status` and `out$diagnostics` before interpretation. Binary
covariate fits remain subject to rare numerical/SE failures; successful checks
do not establish calibration for every covariate mechanism.

To evaluate every combination of assumed hidden fraction and AM history:

```r
total_grid <- ravel_correct_grid(
  out,
  hidden_fraction = seq(.4, .5, .02),
  n_generations = c(1, 3, 15, Inf),
  include_raw = FALSE,
  covariate_mapping_confirmed = TRUE
)
```

Only set the confirmation flag after reviewing the conditional mapping
assumption. Scalars work too; `Inf` denotes equilibrium. Both recorded residual
models are included unless `residual_models` selects one. The full `out`, not
`out$summary`, is needed for corrected SEs. The returned table labels every
scenario and preserves the selected uncertainty method without refitting.
With `include_raw = TRUE` (default), raw total-liability coefficients are also
returned. Incompatible cells remain flagged NA rows. Intervals are conditional
on each assumed scenario, not simultaneous bands over the grid.

### Try the bundled example

`ravel_example` contains 4,000 simulated independent trios with raw scores,
continuous and binary outcomes, and numeric and categorical covariates. It uses
the Gaussian trio-moment model from the validation simulations: hidden fraction
0.4, one generation of assortment with total-liability mate correlation 0.4,
random-mating PGS reliabilities 0.5 and 0.3, and direct plus maternal effects.
No real participant data or simulated SNP genotypes are included.

```r
data("ravel_example", package = "RAVEL")
df <- ravel_example
head(df)

out <- ravel_dataframe_report(
  df, outcome = "Y",
  covariates = c("age", "sex", "site", paste0("PC", 1:6)),
  categorical_covariates = c("sex", "site")
)
out$model_status
out$summary
```

Use `outcome = "Y_binary", outcome_type = "binary_liability"` for the binary
example. The complete installed example fits both outcomes and demonstrates a
total-liability correction:

```r
source(system.file("examples", "dataframe_workflow.R", package = "RAVEL"))
```

See `?ravel_example` for the data-generating assumptions and
`attr(df, "simulation")` for the parameters and population expectations. The
repository script `data-raw/generate_ravel_example.R` regenerates the dataset.

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
| `am_correction` | `"n_gen"` | Total-liability mapping from finite generations, equilibrium, or `"estimated_state"` with matching fitting constraints. |
| `am_n_generations` | `1L` | Number of mating rounds in the constant-strength finite-generation model. For `n = 1` the mating parents are on the founder scale and the child is the first generation produced under assortment. |
| `fixed_effects` | `NULL` | Character vector of jointly adjusted covariate column names |
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
