# RAVEL 0.2.0.9000 (development)

## Genetic x covariate interactions and growth models

- `ravel_fit(..., interactions = "age")` estimates how the visible genetic
  coefficients change with an observed covariate. Products of the centred
  scores and covariate measure the latent interaction with the same loadings
  as the genetic factor, so the interaction is corrected for PGS measurement
  error. `interaction_roles` selects child, mother and/or father liabilities,
  and `interaction_center` sets where main effects are reported.
- `cluster = "FID"` gives cluster-robust standard errors and family-level
  bootstrap resampling. Long-format data (one row per child-occasion) with an
  age interaction give a population-averaged growth model; a squared-age
  column adds curvature.
- With `pi > 0`, interaction coefficients are corrected to the total-liability
  scale by the main-effect map, for trios and parental-sibling fits
  (`delta_M:age` and so on). This requires interaction terms for every family
  member (the default when `pi > 0`) and a moderator independent of genetic
  liability; a latent-correlation check warns otherwise. External `pi_se` and
  finite-generation refit sensitivity are propagated.
- Interactions require a continuous outcome. Fits without these options are
  unchanged.
- Constraint checks now evaluate lavaan's internal equality labels, which
  appear when parameter labels are shared across rows.

## Score-specific phenotype association diagnostic

- `ravel_test_score_signal()` compares the shared-factor outcome model with
  unrestricted observed-score outcome associations, preserving the score
  measurement and covariance assumptions. It returns a separate result and
  does not replace scientific coefficients or alter total-liability correction.
- Continuous, complete-data comparisons use a robust nested-model difference
  test for MLR fits, with explicit convergence, admissibility, nesting and
  constrained-Jacobian identification checks. Unsupported inputs and invalid
  comparisons return a reason and unavailable inference.
- The fully free residual model is explicitly unsupported because its
  score-only measurement model has an unidentified scale. Non-rejection does
  not establish equal genetic targets or unbiased scientific effects.

## Parental siblings

- `ravel_fit(..., parental_sibs = "maternal" | "paternal" | "both")` adds a full
  sibling of the mother and/or father, scored with the same two PGS weight sets.
  The outcome is also regressed on each sibling's visible liability.
  `delta_M = bm - bMS` and `delta_F = bf - bFS` are the parent-specific
  coefficients, with lavaan delta-method or bootstrap uncertainty.
- Families without a sibling stay in the fit with missing sibling scores (FIML
  for continuous outcomes). Sibling scores use new `col_pgs*_mat_sib` and
  `col_pgs*_pat_sib` column arguments.
- With `pi > 0`, the total-liability correction maps the child, parent and
  sibling coefficients jointly. Its visible-hidden cross-block follows from the
  parents' assortment state: same-lineage entries equal the within-person
  covariance, and cross-lineage entries are scaled by k = 1 - 1/(2 V_T) (k^2
  between the two siblings). Corrected sibling coefficients and contrasts are
  reported on the raw and standardized total scales.
- The constrained modes (`"equilibrium"`, `"n_gen"`, `"estimated_state"`) fix
  sibling residual correlations at 1/2 (own parent), 1/4 (child) and 0, give
  siblings their parent's reliability, and tie the parent-sibling and
  child-sibling latent correlations to the variance state.
- `pi_se`, `ravel_correct_pi()`, bootstrap and pi-sensitivity refits cover the
  sibling parameters. Trio fits are numerically unchanged.
- `ravel_simulate_trios()` gains `parental_sibs`, maternal/paternal lineage
  effects and `grandparent_rho_MF`, drawing parents and their siblings from a
  grandparent generation that can itself assort.
- Dataframe reports now accept the same sibling option and custom column
  mappings. Tables, plots, saved packets and assumption grids include the
  sibling coefficients and parent-specific contrasts, with the full enlarged
  core covariance and bootstrap draws retained for total corrections.
- A report's fixed residual model now fixes sibling correlations to the pedigree
  values as well. Explicit sibling residual covariance or correlation settings
  are also supported by `ravel_fit()` without imposing the other constraints.
- Binary sibling fits use their retained analysis sample for sample-size
  reporting; dataframe comparisons use a common complete sample. Continuous
  reports retain partially missing sibling data under FIML. `ravel_fit()`
  rejects infinite sibling scores explicitly.
- Total corrections reject incompatible implied visible/hidden covariance
  matrices while retaining visible estimates. Documentation distinguishes the
  conditional parent-minus-sibling contrast from a causal effect and records
  the additional assumptions required for the latter.

## Supply pi and its standard error during fitting

- `ravel_fit()` now accepts `pi` as the preferred name for the founder hidden
  fraction and `pi_se` for its external standard error. Existing `h` calls
  remain supported; conflicting `h` and `pi` values are rejected.
- A fit supplied with `pi_se` includes external uncertainty in its corrected
  effect tables and printed results. The stored `pi_correction` table retains
  the joint covariance and propagation details. Conditional SEM and bootstrap
  covariance remain separate so later corrections never count pi uncertainty
  twice.
- `ravel_correct_pi()` remains available for sensitivity analyses. Changing pi
  requires refitting when finite-generation constraints depend on its value.
- The README now focuses on usage; detailed methods and advanced examples are
  retained in separate documentation.

## Estimated assortment state

- New opt-in `fit_constraints = "estimated_state"` identifies genetic variance
  growth and founder scaling from spouse and parent-child covariances without
  supplying generations. Separate unknown PGS reliabilities are retained.
- This adds one covariance parameter relative to specified-history constraints.
  It is a smooth state model, with diagnostics for constant-strength history
  compatibility, rather than integer-generation selection.
- Fitting is independent of pi. Physical compatibility is checked during
  correction, and full covariance propagation includes estimated state and
  external pi uncertainty. Existing modes and defaults are preserved.
- Estimated-state fits verify local identification along the constraint surface
  and the full parameter covariance. Nonidentified or active-boundary fits are
  excluded from regular Wald inference; aggregate packets retain this check.

## Optional constrained fitting

- `ravel_fit()` and `ravel_dataframe_report()` accept `fit_constraints =
  "equilibrium"` or `"n_gen"`. The default `"none"` retains the existing model.
- Each PGS retains its own unknown reliability. Restrictions link reliability
  within each score across family roles, preserve unknown marginal scales, and
  fix the Mendelian residual correlations. Finite-generation restrictions allow
  child reliability to change with the modeled genetic variance.
- Finite histories assume constant total-liability spousal correlation from
  random-mating founders, conditional on the supplied hidden fraction and
  generation count. Correction packets record these assumptions and reject
  incompatible post-estimation changes that would require a new fit.

## External hidden-fraction uncertainty

- New `ravel_correct_pi()` accepts an external pi estimate and its raw-scale
  standard error through `pi_se`. Normal-Wald or logit-Wald confidence intervals
  remain available through `pi_ci`; supply exactly one uncertainty input.
  Joint delta propagation includes the full SEM
  parameter covariance and pi variance for raw and standardized corrected
  estimates. The returned table retains the full transformed covariance.
- Independence of external pi is the default; an explicit cross-covariance can
  represent sample overlap. For finite-generation fits that depend on pi,
  automatic live-fit sensitivity refits propagate the additional uncertainty.
- Bootstrap inputs use their joint core-parameter covariance with delta-method
  propagation and are labeled accordingly; no percentile-interval claim is made.

## Reporting integration

- Correction grids inherit constrained histories, support estimated-state pi
  sweeps, and optionally propagate independent external pi SEs.
- Print helpers accept corrected tables. Plots accept corrected tables and grids,
  select one predictor scale, and preserve scenario and inference labels.
- Fit tables retain fitting assumptions and variance-state attributes. Report
  exports preserve identification failure reasons and omit fictitious generation
  counts for equilibrium and estimated-state modes.
- Explicit raw table requests and bootstrap raw-scale aliases now work at h = 0.
- PNG export tests device availability at runtime, including the native macOS
  fallback when a configured Cairo library cannot load.
