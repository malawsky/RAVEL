# Family regressions with RAVEL

This directory contains the original parental-deviation analysis script and a
companion for adult outcomes. The original script retains its analysis defaults;
its installation check now requires the family-model functions in the new bundle.
The companion only defines a function when sourced, so it is safe to load before
configuring an analysis.

## Install the bundled package

Keep `RAVEL_0.2.0.9001.tar.gz` alongside these scripts. With this directory as your
working directory, install into a local library:

```r
dir.create("library", showWarnings = FALSE)
.libPaths(c(normalizePath("library"), .libPaths()))
install.packages("RAVEL_0.2.0.9001.tar.gz", repos = NULL, type = "source",
                 lib = "library")
library(RAVEL)
source("ravel_family_regressions.R")
```

Existing dependencies, including lavaan, must be installed. Restart R first if
an older RAVEL namespace is already loaded.

For a runnable demonstration, source `run_synthetic_example.R`. It uses the
included artificial `synthetic_family_example.rds`, runs all eight maternal
regressions and writes tables and fitted objects locally. These observations are
simulated and are not a sample of your data.

## Maternal height with trio-assisted identification

The following uses default score-column names (`PGS1_mother`, `PGS2_mother`,
etc.). The pi values are illustrative; use the calibration for the scored trait.

```r
result <- ravel_family_regressions(
  dat = my_data,
  outcomes = "mom_height",
  genetic_roles = c("child", "mother", "father", "mat_sib", "pat_sib"),
  focal = "mother",
  pi = .4, pi_se = .02,
  pi_grid = c(.3, .4, .5),
  models = list(estimated_state = list(fit_constraints = "estimated_state")),
  output_prefix = "maternal_height"
)
result$table
result$failures
result$sample_counts
```

This fits one joint genetic model per pi assumption and derives all eight adult
regressions containing maternal liability. The child remains available for
identification without being included in the reported outcome regressions.
Both scores of every member in `genetic_roles` need to be supplied; partially
missing observations are handled by the fitter. Binary models use complete
observations and report that restriction.

Use `pgs_columns = c(mother1 = "your_column", mother2 = "your_other_column", ...)`
to map custom names. Only selected genetic members need entries. Supply
`outcome_models = list(sib = c("mother", "mat_sib"), full = c("mother",
"mat_sib", "father", "pat_sib"))` to request particular regressions.

Outcomes are kept in their input units by default. Set `standardize_outcome = TRUE`
for phenotype-SD units. `total_std` standardizes the genetic liability to its
current-generation SD; it does not itself standardize the phenotype.

## Child-free examples

Change `genetic_roles` to `c("mother", "mat_sib")` for a sibling pair, to
`c("mother", "father", "mat_sib")` for partners and the mother's sibling, or
to the four adults for partners and both siblings. No child columns are required.

Partners alone cannot identify a free accumulated variance state. Use an
appropriate history model, or supply `parent_visible_variance` with the
estimated-state model. A supplied state is treated as fixed, so its uncertainty
is not included automatically.

## Covariates, clustering and interpretation

Use `covariates` for numerical adjustment variables and `cluster` for the family
ID used for robust inference. Further supported RAVEL arguments, including
`n_boot`, `interactions` and `interaction_roles`, can be passed through.
Corrected reduced regressions with covariates require
`covariate_mapping_confirmed = TRUE`, asserting the specific hidden-residual
orthogonality assumption described in [family-models.md](family-models.md). Without that
confirmation, visible projections are retained and corrected rows explain why
they are unavailable.

For a cross-sectional adult outcome, prepare one record per focal adult. An
optional `subject_id` checks that rule. Repeating the same mother once per child
changes the weighting even if standard errors are clustered; the function does
not silently choose a child or discard duplicates. With missing scores,
`complete_cases = TRUE` provides a complete-observation sensitivity analysis.
Bootstrap resampling requires at least two observed independent family IDs.
A constant dummy `FID` is rejected. If rows truly are independent, explicitly
pass `cluster_var = NULL, cluster = NULL` to request row resampling.

Inspect row-level correction and inference status as well as `$failures`.
Regression subsets share one joint fit, so they are correlated and are not
independent tests. Sibling contrasts are conditional association contrasts; a
causal own-genetic interpretation requires additional assumptions.

## Reuse a fit directly

```r
fit <- ravel_family_fit(my_data,
  genetic_roles = c("child", "mother", "father", "mat_sib", "pat_sib"),
  outcome_roles = c("mother", "mat_sib", "father", "pat_sib"),
  col_outcome = "mom_height", pi = .4, fit_constraints = "estimated_state")

fit$regression$table
sibling_projection <- ravel_family_regression(fit, c("mother", "mat_sib"))
sibling_projection$table
saveRDS(fit, "maternal_height_fit.rds")
```

Reusing the fit holds genetic-model information constant. Removing a member from
`genetic_roles` is a separate analysis and requires fitting again. Pi-dependent
history models also require refitting when pi changes; the companion always
refits its pi grid to respect this requirement.
