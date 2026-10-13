# RAVEL API and workflow reference

The [README](../README.md) gives the shortest route to fitting and reporting.
This reference retains advanced options, data preparation details and output
field definitions. The [methods note](methods.md) records model assumptions,
coefficient scales, variance recurrences and uncertainty propagation.

## Column mapping

A data frame with two PGS per trio member, an outcome, and optionally covariates. By default RAVEL expects columns named `PGS1_mother`, `PGS2_mother`, `PGS1_father`, `PGS2_father`, `PGS1_child`, `PGS2_child`, `Y` — but you can map your own column names:

```r
# Default column names
fit <- ravel_fit(dat, pi = 0.5, n_boot = 1000)

# Custom column names
fit <- ravel_fit(dat,
  col_pgs1_mother = "EA1_PGS_mom",  col_pgs2_mother = "EA2_PGS_mom",
  col_pgs1_father = "EA1_PGS_dad",  col_pgs2_father = "EA2_PGS_dad",
  col_pgs1_child  = "EA1_PGS_kid",  col_pgs2_child  = "EA2_PGS_kid",
  col_outcome     = "phenotype",
  pi = 0.5, n_boot = 1000
)

print_simple(fit)
```

## Fitting constraints and correction inputs

| `fit_constraints` | External inputs for fitting | External inputs for total correction |
|---|---|---|
| `"none"` (default) | None | pi and finite history or equilibrium |
| `"equilibrium"` | None | Founder pi |
| `"n_gen"` | Generation count; pi after generation one | Matching generation count; fitted pi after generation one |
| `"estimated_state"` | None | Founder pi |

`pi` and `h` are aliases in `ravel_fit()`. `ravel_dataframe_report()` and older
correction/expectation functions retain `h`; `ravel_correct_pi()` uses `pi`.
At pi zero there is no hidden component. None of these arguments specifies
PGS prediction error or the fraction of outcome variance that is genetic.
The models are opt-in: incorrect restrictions can introduce bias, and more
constraints do not guarantee a better result under misspecification.

`ravel_dataframe_report()` accepts the same `fit_constraints` modes.
Constrained modes select the fixed residual model when `residual_models` is
omitted and reject an explicitly requested free residual model. Summaries
remain on the visible scale, with assumptions recorded in the correction packet.

```r
out_eq <- ravel_dataframe_report(df, fit_constraints = "equilibrium", h = .4)
out_g3 <- ravel_dataframe_report(df, fit_constraints = "n_gen", h = .4,
                                am_n_generations = 3L)
```

## Reporting and exporting

All fitting modes use the existing fit, table, print, report and bootstrap
interfaces. Fit summaries name the constraint and correction mode; estimated-state
summaries also show growth and founder-scaled variances. Table conversion retains
its existing column names (including `beta_F` for the paternal effect), with
assumptions in `attr(table, "settings")`, `attr(table, "constraint_spec")`, and
state diagnostics in `attr(table, "estimated_state")`.

```r
fit_state <- ravel_fit(dat, pi = .4, pi_se = .02,
                       fit_constraints = "estimated_state")
print_simple(fit_state)
fit_table <- as.data.frame(fit_state, include_raw = TRUE)
attr(fit_table, "settings")

# One-call fitting can include pi uncertainty; separate correction remains available.
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

For an estimated-state fit, inspect
`fit_state$hidden_correction_details$estimated_state`; the pi-adjusted
state estimates, SEs and covariance are in
`attr(fit_state$pi_correction, "estimated_state")`.
An estimated-state grid uses the fitted variance state at every pi and records
its generation count as `NA`. Supplying a generation grid for that mode is
rejected because it would change the fitting assumptions.

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

## Choosing an uncertainty source

The fit retains two uncertainty summaries conditional on the supplied pi:

- `lavaan_inference`: SEM-based standard errors with normal-approximation CIs, including delta-method propagation to the hidden-corrected raw and standardized total-liability scales
- `bootstrap`: empirical bootstrap summaries when `n_boot > 0`

The print methods can report either source. When the fit includes `pi_se`,
corrected rows additionally propagate external pi uncertainty and use
delta-Wald intervals; visible bootstrap rows retain their percentile intervals:

```r
fit <- ravel_fit(dat, pi = 0.5, pi_se = 0.05, n_boot = 500)

print_simple(fit, uncertainty = "bootstrap")
print_simple(fit, uncertainty = "lavaan")
print_ravel_results(fit, uncertainty = "bootstrap")
```

## Covariates and binary outcomes

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

## Dataframe reports

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

## Bundled example

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

## Fixing cross-person residual correlations

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

## `ravel_fit` options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dat` | — | Data frame with 6 PGS columns + `Y` |
| `pi` | `NULL` | Preferred name for the founder hidden additive-variance fraction; must be in `[0, 1)`. |
| `h` | `0` | Legacy name for `pi`. If both are supplied explicitly they must agree. With neither supplied, the hidden fraction is zero. |
| `pi_se` | `NULL` | External SE on the raw pi scale. Omit to treat pi as fixed; supply to propagate its uncertainty in the same call. |
| `core_pi_cov` | `NULL` | Optional covariance between the six SEM core-parameter estimation errors at fixed pi and external pi. The default assumes independence. |
| `covariate_mapping_confirmed` | `FALSE` | Confirm that the visible-to-total mapping applies conditional on jointly fitted covariates before adding pi uncertainty. |
| `fit_constraints` | `"none"` | `"none"`, `"equilibrium"`, `"n_gen"`, or `"estimated_state"`; see the README mode table. |
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
| `parental_sibs` | `"none"` | Add a full sibling of the mother (`"maternal"`), father (`"paternal"`) or both. Requires `fit_constraints = "none"` and `pi = 0`. |
| `col_pgs1_mat_sib`, `col_pgs2_mat_sib` | `"PGS1_mat_sib"`, `"PGS2_mat_sib"` | Column names for the maternal sibling's scores; missing values allowed |
| `col_pgs1_pat_sib`, `col_pgs2_pat_sib` | `"PGS1_pat_sib"`, `"PGS2_pat_sib"` | Column names for the paternal sibling's scores; missing values allowed |

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
| `lavaan_inference` | SEM-based inference conditional on pi. Printing and table conversion use `pi_correction` for corrected rows when external uncertainty was supplied. |
| `bootstrap` | Fixed-pi CIs and SEs if `n_boot > 0` |
| `pi_correction` | Corrected estimates with propagated external pi uncertainty when `pi_se` was supplied, including covariance and derivative attributes. |
| `parental_sibs` | With parental siblings: `coefficients` (`bms`, `delta_m`, `bfs`, `delta_f`), `lavaan_table` (rows `beta_MS`, `delta_M`, `beta_FS`, `delta_F`), `latent_cov`, `resid_cov`, `resid_corr`, `reliability`, `determinacy` and `n_with_sib`. `NULL` otherwise. |
| `converged` | Logical convergence flag |

## Analytical expected values

See the [methods note](methods.md#analytical-expected-values) for runnable
`ravel_expected_values()` and `ravel_expected_pgs_ols()` examples, including
raw versus standardized total coefficients and the reliability reference scale.

## Print functions

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
| `parental_sibs` | `"none"` | Also simulate a full sibling of the mother, father or both. Parents and siblings are then children of a random-mating grandparent generation. |
| `lineage_effect_mat`, `lineage_effect_pat` | `0` | Effect on Y per SD of the maternal (paternal) grandparents' midparent total liability; lineage confounding shared by a parent and their sibling. |

`G_mother`, `G_father`, and `G_child` are legacy variable names for the raw
total liabilities T_m*, T_p*, and T_c* that enter phenotype generation. When
present, `G_mother_vis`, `G_father_vis`, and `G_child_vis` are role-standardized
visible liabilities used for oracle comparisons. These legacy returned-column
names are retained for backward compatibility.
