# Family regressions with RAVEL

This directory contains the original parental-deviation script, a matching
parental-phenotype workflow, and a function-based companion for adult outcomes.
The original script retains its fitting defaults. All three use the family-model
functions in the bundled package.

## Parental phenotypes: the editable-settings workflow

Open **`ravel_parental_phenotypes.R`**, edit its SETTINGS block, and run:

```r
source("ravel_parental_phenotypes.R")
```

Keep it alongside `ravel_parental_deviation.R`, which supplies the existing data
preparation and fitting code, `ravel_sem_audit.R`, and the package archive below. With no data file
set, it runs the included synthetic maternal-height example. To use your data,
set `data_file`, `outcome`, the six trio `pgs_columns`, optional `sib_columns`,
and your calibrated `pi`/`pi_se`. The workflow standardizes continuous outcomes
and produces the same main estimates, pi-grid analyses, CSV and Excel-friendly
CSV files as the original script, plus an RDS file retaining fits and diagnostics.

The child's two PGSs always remain in the genetic model. The reported outcome
regressions always contain the focal parent's liability: mother alone, mother
plus father, mother plus her sibling, and the other combinations with any
included parental siblings. With both siblings there are **eight regressions**.
The default is `focal = "mother"`; use `focal = "father"` for paternal outcomes.
You can supply a named list in `outcome_models` to select specific combinations,
but each must include the focal parent.

By default, `genetic_configurations = "available"` uses all supplied relatives
to identify one genetic model per outcome, assortment assumption and pi value.
Every outcome subset uses that same source fit. Set
`genetic_configurations = "nested"` to also repeat the original trio-only,
maternal-sibling, paternal-sibling and both-sibling genetic-block analyses,
where scores are available. Each block runs its own adult predictor subsets,
always retaining the focal parent.

The source fit retains nuisance outcome associations with every genetic member,
including the child, to estimate the joint moments. Selected regressions are
then derived from those moments: excluding the child from the reported
regression does not assert that its source coefficient is zero. Plain PGS1/PGS2
comparisons are separately fitted to each selected adult subset too. Midpoint,
parent-difference and sibling-difference SEs use joint coefficient covariance.

Use one row per focal adult for cross-sectional outcomes. Set `subject_id` to
check for duplicate adults; the script does not choose one child automatically.
For covariate-adjusted corrected subsets, review the assumption below before
setting `covariate_mapping_confirmed = TRUE`. Without it, visible results remain
available and corrected rows state why they are unavailable.

To load either editable-settings script's functions without running its example:

```r
options(ravel.workflow_functions_only = TRUE)
source("ravel_parental_phenotypes.R")
options(ravel.workflow_functions_only = NULL)
```

## Detailed SEM tables

The parental workflow keeps its primary results CSV and adds seven separate
audit CSVs beside it. Each row identifies the outcome, genetic block, model,
main/grid assumption and pi. These describe each shared source fit once, rather
than repeating the same SEM statistics for every adult predictor subset.

| Filename suffix | Contents |
|---|---|
| `_sem_fit_statistics.csv` | Available lavaan fit measures, including chi-square, df, CFI/TLI, RMSEA, SRMR, likelihood/AIC/BIC where defined, convergence, admissibility and N. |
| `_sem_parameters.csv` | Loadings, regressions, intercepts/thresholds, variances, covariances and defined parameters, with SEs, intervals, p-values and fixed/free indicators. |
| `_sem_covariances.csv` | Latent and residual/exogenous variance/covariance parameters with their SEs and intervals. |
| `_sem_standardized_parameters.csv` | Fully standardized estimates and their own SEs/intervals, with extraction warnings retained. |
| `_sem_implied_covariances.csv` | Model-implied observed and latent covariance matrices in long form, with SEs where a validated derivative or direct parameter is available. |
| `_sem_parameter_vcov.csv` | Lower triangle of the sampling covariance of the estimated SEM parameters, for reconstructing joint uncertainty. |
| `_sem_extraction_status.csv` | Extraction status and the full warning/error text for each table-producing step. |

SEM values use the prepared data's units: continuous outcomes are standardized,
each PGS version uses the original workflow's pooled trio scaling, and genetic
latent factors have unit variance by role. The SEM SEs use the fitted estimator
(including robust/clustered uncertainty when applicable) **conditional on the
supplied pi**. They do not add external `pi_se`; the corrected regression tables
continue to propagate that uncertainty. Fixed latent variances have SE zero
because they set the factor scale, not because their population variance was
estimated without error. An unavailable SE or fit index is left missing, with a
status explaining why. Warnings and inadmissible fits remain visible.

All tables, source fits and one-factor objects are also in the saved RDS:

```r
res <- readRDS("ravel_parental_phenotypes_mom_height_fits.rds")
attr(res, "sem_audit")$sem_covariances
attr(res, "diagnostics")[["both|Estimated state"]]
attr(res, "fits")  # main and pi-grid fit records
```

Set `include_sem_audit = FALSE` to skip the extra tables during a quick trial.

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

`ravel_family_regressions.R` defines a function when sourced and never starts an
analysis automatically. It also supports genetic models with no children; the
editable parental-phenotype workflow above deliberately retains the trio.

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

## When the one-factor check is unavailable

The check compares the two PGS versions' phenotype associations under the full
shared-factor model. It is run once per main source fit, not once per projected
adult regression or pi-grid value. It does **not** test that the child's outcome
coefficient is zero. Degrees of freedom are 3 for a trio, 4 with one parental
sibling and 5 with both siblings; finite-generation inference is conditional on
the pi supplied to that fit.

Both editable workflows now record `diagnostic_status`, `diagnostic_df` and
`diagnostic_N` in the CSV and preserve the full diagnostic object in the RDS.
Removing covariates alone does not guarantee that this comparison is available:

| Status | Meaning |
|---|---|
| `ok` | The validated likelihood comparison succeeded. |
| `unsupported_free_measurement` | The fully free residual score model has an unidentified score-only alternative; no valid test is implemented for that model. |
| `unsupported_missingness` | An included score or outcome is missing; the diagnostic does not silently drop observations from the source fit. |
| `unsupported_design` | Repeated or clustered families fall outside the independent-row diagnostic. |
| `unsupported_covariates` / `unsupported_outcome` | Covariate or binary-outcome comparisons are outside its validated scope. |
| `invalid_score_model`, `invalid_H0`, `invalid_H1` | A required fit failed identification or admissibility checks. Inspect the saved reason and rank diagnostics. |
| `fit_warning` / `fit_error` | Inspect the saved warnings or error from the diagnostic refits. |

A constrained sensitivity analysis can be checked separately, but changing the
measurement restrictions changes the model being tested. A missing p-value is
not evidence for or against one factor. A non-significant valid test also does
not establish unbiased genetic effects. The original Desktop bundle and the
current checkout contain the same score-diagnostic implementation; the workflow
now exposes and saves its results instead of retaining only a short message.

[Parental-workflow validation and reproduction](parental-workflow-validation.md)
records the independent regression checks, one-factor reproductions, covariance
SE checks and full default run.
