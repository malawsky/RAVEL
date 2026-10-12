# RAVEL 0.2.0.9000 (development)

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
