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

`parent_contrasts = TRUE` in `print_simple()`, `print_ravel_results()`,
`as.data.frame()` and `ravel_to_dataframe()` adds three rows per scale:

- `beta_mid` (bm + bf), the coefficient on the parental midpoint (G_m + G_p)/2;
- `delta_M_mid` (bm - bf), the maternal deviation;
- `delta_P_mid` (bf - bm), the paternal deviation.

They use the fit's joint covariance or bootstrap draws (see
[Methods](methods.md#parental-midpoint-and-deviation-contrasts)).

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

## Score-specific phenotype association test

Run the diagnostic on a live fitted result and save it separately from the
scientific effect tables:

```r
fit <- ravel_fit(dat, fit_constraints = "equilibrium", pi = .4)
signal_test <- ravel_test_score_signal(fit)
print(signal_test)
saveRDS(signal_test, "score_signal_test.rds")
```

The null retains a free outcome path from every included genetic liability.
The alternative replaces all those paths with outcome regressions on the six,
eight or ten observed scores. The joint score measurement model and all its
constraints stay the same, and nuisance parameters are refitted in both models.
This tests agreement of the two scores' normalized phenotype associations;
it does not test whether parental or sibling effects are zero.

The initial scope is a continuous phenotype, complete observations, independent
families and no covariates. Scores must be observed rather than imputed, with
both versions oriented to the same genetic target. Repeated or overlapping
families, sampling weights, multigroup models, binary outcomes, incomplete
data and covariate extensions are not validated for this diagnostic. Ordinary
RAVEL fitting remains available under its existing scope. Imputation or
dependence not recorded in the data cannot be discovered from the scores alone;
the analyst must establish that the diagnostic's sampling assumptions hold.

The default `fit_constraints = "none"` model with free residual covariances is
unsupported because the score-only model is unidentified. Merely changing the
reference degrees of freedom would not repair its alternative fit. Identified
constrained modes and the fixed-residual benchmark are checked for full local
rank and numerical regularity. With validated ranks, the reference difference
is three for trios, four for one parental sibling and five for both. Near-zero
assortment, weak identification and active boundaries can prevent inference.
Supplying a fit with fixed parental outcome paths is also unsupported: the
diagnostic does not silently replace the user's null model.
The `crossperson = "equal_pc"` residual-covariance option and manually altered
SEM syntax are also outside the initial diagnostic scope.

The result records status and reason, test statistic, reference degrees of
freedom, p value, the actual adjustment, both log likelihoods and model degrees
of freedom, sample size, included roles, fitting constraints, software version,
and fit/rank diagnostics. Unsupported or invalid comparisons return unavailable
inference. For MLR fits the comparison uses a robust nested-model difference
test, not subtraction of separately scaled fit statistics.

Useful fields are `status`, `reason`, `statistic`, `df`, `p_value`,
`raw_statistic`, `test`, `adjustment`, `log_likelihood`, `model_df`, `N`,
`roles`, `mode`, `measurement_constraints`, `external_pi`, `history`, `software`, `diagnostics` and
`fits`. The last two retain the refitted models, constructive nesting check,
warnings and identification details for review. They are not ordinary
coefficient tables and can be larger than the compact printed result.

Both diagnostic fits use sample-SD units for the scores and phenotype, with
fixed covariance values converted consistently. Observed marginal variances
and unknown score reliabilities remain estimated under the selected fitting
constraints. Public `log_likelihood` values are restored
to the original input units; `diagnostics$canonical_log_likelihood` and
`diagnostics$input_standard_deviations` retain the transformation. The models
in `fits` are on the diagnostic units. The robust test uses copies with the
final equality-constraint tangent supplied explicitly; the original fit is
left unchanged.

Rank checks use naturally scaled parameter and moment coordinates, a relative
singular-value tolerance of `1e-7`, stability checks at `1e-8` and `1e-6`, and
a minimum smallest-to-largest singular-value ratio of `1e-5`. The scaled
projected gradient must be at most `1e-4` and equality residuals at most
`1e-6`. Covariance correlation matrices must have minimum eigenvalue above
`1e-8`; standardized loadings must exceed `1e-5` in magnitude and inequality
margins must exceed `1e-7`. Each observed residual variance must also exceed
`1e-8` of that variable's total fitted variance. These numerical safeguards do not guarantee good
finite-sample behavior near a population singularity. The returned diagnostic
records the values and thresholds, rather than only a pass/fail flag.

Unresolved fitting or test warnings prevent inference. Two specific numerical
warnings from optimizer trial evaluations can be resolved by one restart from
the final admissible optimum: both fits must pass every regularity check, the
restart must emit no warnings, and their log likelihoods must agree within
`1e-5`. Original warnings, retry warnings, validation results and likelihood
agreement remain in `diagnostics$refit_retries`. Other warnings are not waived.

Use the original live fit; a CSV effect table or aggregate correction packet
does not preserve the likelihood and model needed for this test. An explicit
`data` argument must represent the same observations and preprocessing used in
the fitted model; it is not an invitation to switch the analysis sample.

The diagnostic leaves the fit and its scientific coefficients unchanged.
When finite-generation fitting constraints depend on pi, inference is
conditional on that same fixed pi in both models. `pi_se` propagation into
corrected coefficients does not propagate pi uncertainty into this test.
Changing only the post-fit total correction has no effect on the comparison.

Rejection indicates disagreement with the shared-factor outcome restriction
under the retained measurement assumptions. It does not identify the source:
different direct/indirect mixtures and other measurement or covariance
misspecification can also cause rejection. Non-rejection does not establish
genetic correlation one or unbiased effects. A score-specific outcome signal
proportional to each score's loading can bias scientific effects while remaining
invisible to this test; see the [derivation and blind spot](methods.md#score-specific-phenotype-association-test).

## Parental-sibling workflow

The option is off by default. Supply at most one full sibling per selected
parent, scored with the same two weight sets. The report workflow assumes
independent extended families; relatives shared across rows require a sampling
model that accounts for that dependence. `ravel_fit()` expects
`PGS1_mat_sib`, `PGS2_mat_sib`, `PGS1_pat_sib`, `PGS2_pat_sib` as appropriate;
its `col_pgs1_mat_sib` / `col_pgs2_mat_sib` and paternal equivalents map custom
column names. The dataframe report uses lowercase defaults:

```r
out <- ravel_dataframe_report(df, outcome = "Y", parental_sibs = "both",
                              fit_constraints = "estimated_state", h = .4)
out$summary
print(out$plot)  # adds sibling coefficients and parent-specific contrasts
corrected <- ravel_correct_pi(out$packet, pi = .4, pi_se = .05)
subset(corrected, parameter %in% c("delta_M", "delta_F"))
print(ravel_plot_report(corrected))
```

For custom report names, include the requested sibling entries in the named
`pgs_columns` vector, alongside the six usual trio entries. For example, add
`PGS1_mat_sib = "aunt_score_1"` and `PGS2_mat_sib = "aunt_score_2"` for a
maternal sibling. Only the selected side(s) are required. When covariates are present, the
correction assumes the mapping applies conditional on them (see the README
covariate section).

Continuous reports require complete core trio/outcome/covariate data and retain
missing sibling scores for FIML, assuming ignorable missingness conditional on
the modeled variables. Binary reports exclude families with missing
requested sibling scores before fitting every comparison model, warn, and
record the exclusion count. `ravel_fit()` also uses complete families for binary
outcomes. At least 50 families with both sibling scores are required per selected
side; this is an input check, not a guarantee of identification or precision.

Tables, print methods, correction packets and grids include `beta_MS`, `delta_M`
and/or `beta_FS`, `delta_F`. The ordinary parent coefficients are conditional on
the included sibling liabilities. The deltas subtract parent and sibling
coefficients jointly, with their covariance, on each reported liability scale;
they are not standardized by the SD of the parent-minus-sibling difference.
The paternal contrast is named `delta_F` in both live-fit and report tables.
The report's `residual_models = "fixed"` fixes all included within-score
residual correlations: 1/2 for parent-child and full-sibling pairs, 1/4 for
child-aunt/uncle pairs, and zero for unrelated pairs. `"free"` estimates these
covariances. Ordinary PGS comparator rows remain trio regressions and therefore estimate a
different conditional association when siblings enter RAVEL. See the
[methods](methods.md#parental-siblings) before giving a causal interpretation.

Packets preserve ten core parameters for one sibling and fifteen for both,
rather than the usual six. All fitting/correction modes and uncertainty sources
use that full covariance. Save the RDS packet for later inference; CSV retains
results but not their joint covariance. An incompatible hidden covariance makes
the total correction unavailable and leaves visible estimates available. If
`pi_se` was requested, `fit$pi_correction_error` and the exported table's matching
attribute retain the reason; corrected intervals remain unavailable.

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

With a binary outcome and covariates, a joint probit regression on all requested
PGS indicators (including parental siblings) and the encoded covariates supplies marginal variance-one liability
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
  include_raw = FALSE
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

With parental siblings, the same arguments accept `e1_m_ms`, `e1_f_ms`,
`e1_c_ms` for a maternal sibling; `e1_f_fs`, `e1_m_fs`, `e1_c_fs` for a paternal
sibling; and `e1_ms_fs` between the siblings. Replace `e1` with `e2` for PGS2.
Only labels for included relatives are allowed. For example,
`fixed_crossperson_corr = c(e1_m_ms = .5, e2_m_ms = .5)` fixes the maternal
parent-sibling correlations. Omitted pairs remain free unless a constrained
fitting mode fixes them. A residual cannot be fixed on both the covariance
and correlation scales, and settings that conflict with the chosen constrained
mode are rejected.

## `ravel_fit` options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dat` | — | Data frame with 6 trio PGS columns + `Y`, and optional sibling scores |
| `pi` | `NULL` | Preferred name for the founder hidden additive-variance fraction; must be in `[0, 1)`. |
| `h` | `0` | Legacy name for `pi`. If both are supplied explicitly they must agree. With neither supplied, the hidden fraction is zero. |
| `pi_se` | `NULL` | External SE on the raw pi scale. Omit to treat pi as fixed; supply to propagate its uncertainty in the same call. |
| `core_pi_cov` | `NULL` | Optional covariance between the SEM core-parameter estimation errors at fixed pi and external pi. The default assumes independence. |
| `covariate_mapping_confirmed` | `FALSE` | Deprecated and ignored; retained so existing calls run. Corrections with covariates assume the mapping applies conditional on them (README). |
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
| `parental_sibs` | `"none"` | Add a full sibling of the mother (`"maternal"`), father (`"paternal"`) or both. Works with every constraint mode and with `pi`/`pi_se`. |
| `col_pgs1_mat_sib`, `col_pgs2_mat_sib` | `"PGS1_mat_sib"`, `"PGS2_mat_sib"` | Column names for the maternal sibling's scores; continuous fits use FIML for missing values |
| `interactions` | `NULL` | Numeric covariate columns interacting with the genetic liabilities (e.g. `"age"`); each is also added as a main effect. Continuous outcomes only. |
| `interaction_roles` | `"child"` | Liabilities that interact: any of `"child"`, `"mother"`, `"father"`, `"mat_sib"`, `"pat_sib"` (siblings in the model). With `pi > 0` it defaults to every family member, which the total-liability correction requires. |
| `interaction_center` | `NULL` | Named values at which main effects are reported, e.g. `c(age = 10)`; default is each covariate's mean. |
| `cluster` | `NULL` | Column of independent units (e.g. `"FID"` for long-format repeated measures); gives cluster-robust SEs and cluster bootstrap. |
| `col_pgs1_pat_sib`, `col_pgs2_pat_sib` | `"PGS1_pat_sib"`, `"PGS2_pat_sib"` | Column names for the paternal sibling's scores; continuous fits use FIML for missing values |

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
| `parental_sibs` | With parental siblings: `coefficients` (`bms`, `delta_m`, `bfs`, `delta_f`), `lavaan_table` (rows `beta_MS`, `delta_M`, `beta_FS`, `delta_F`), `corrected` (`raw` and `std` total-scale coefficients and contrasts, plus the four- or five-role covariance blocks), `latent_cov`, `resid_cov`, `resid_corr`, `reliability`, `determinacy` and `n_with_sib`. `NULL` otherwise. |
| `interactions` | With `interactions`: moderators, roles, centres, internal scales, `table` of visible-scale coefficients (`beta_D:age`, ..., `delta_M:age` with siblings), `corrected` total-liability rows when `pi > 0`, `independence` (latent correlation of each G with each moderator), and `theta` (model-scale coefficients). `NULL` otherwise. |
| `cluster`, `n_clusters` | Cluster column and number of clusters when `cluster` is supplied. |
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
| `grandparent_rho_MF` | `0` | Spouse correlation among grandparents. Equal to `rho_MF` gives two constant rounds of assortment (`am_n_generations = 2`). |

`G_mother`, `G_father`, and `G_child` are legacy variable names for the raw
total liabilities T_m*, T_p*, and T_c* that enter phenotype generation. When
present, `G_mother_vis`, `G_father_vis`, and `G_child_vis` are role-standardized
visible liabilities used for oracle comparisons. These legacy returned-column
names are retained for backward compatibility.
