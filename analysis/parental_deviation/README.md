# Family regressions with RAVEL

This directory contains the original parental-deviation script, a matching
parental-phenotype workflow, and a function-based companion for adult outcomes.
The child script has a short editable settings block; all three use the
functions in the bundled package.

## Parental phenotypes: the editable-settings workflow

Open **`ravel_parental_phenotypes.R`**, edit its SETTINGS block, and run:

```r
source("ravel_parental_phenotypes.R")
```

Keep it alongside `ravel_parental_deviation.R`, which supplies the existing data
preparation, fitting and reporting code, and the package archive below. With no data file
set, it runs the included synthetic maternal-height example. To use your data,
set `data_file`, `outcome`, the six trio `pgs_columns`, optional `sib_columns`,
and your calibrated `pi`/`pi_se`. The workflow standardizes continuous outcomes
and produces the same main estimates and pi-grid analyses, with Excel-ready
CSV tables and an RDS file retaining fits and diagnostics.

The child's two PGSs always remain in the genetic model. The reported outcome
regressions cover every nonempty combination of the available adult liabilities,
including father-only and father-plus-paternal-sibling models for a maternal
phenotype, and vice versa. There are **15 regressions** with both parental
siblings, **seven** with one sibling, or **three** with the parents only.
Use the same settings for maternal and paternal outcomes; no focal parent is
required. Supply a named list in `outcome_models` to select specific combinations.
Older function calls that supply `focal` still run, with a warning that it no
longer restricts the predictors.

By default, `genetic_configurations = "available"` uses all supplied relatives
to identify one genetic model per outcome, assortment assumption and pi value.
Every outcome subset uses that same source fit. Set
`genetic_configurations = "nested"` to also repeat the original trio-only,
maternal-sibling, paternal-sibling and both-sibling genetic-block analyses,
where scores are available. Each block runs its own adult predictor subsets,
using every nonempty adult combination.

The source fit retains nuisance outcome associations with every genetic member,
including the child, to estimate the joint moments. Selected regressions are
then derived from those moments: excluding the child from the reported
regression does not assert that its source coefficient is zero. Plain PGS1/PGS2
comparisons are separately fitted to each selected adult subset too. Midpoint,
parent-difference and sibling-difference SEs use joint coefficient covariance.

Use one row per phenotyped adult for cross-sectional outcomes. Set `subject_id` to
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

## Reporting tables: the same layout for both workflows

Each phenotype gets the existing coefficient table plus five supporting tables.
Every CSV uses the Excel-friendly format (a `sep=` first line); there is no
second plain-CSV copy. Main/grid fits are identified by outcome, genetic block,
assortment model and pi. For parental regressions, fit statistics describe the
shared **source SEM**, once per fit, rather than each derived adult subset.

| Filename suffix | Contents |
|---|---|
| `_excel.csv` | Existing effects, contrasts, SEs, confidence intervals, p-values and pi sensitivity. |
| `_fit_statistics_excel.csv` | Sample counts, likelihood, parameter count and available fit indices, with reasons for unavailable values. Lavaan models include ordinary/scaled/robust chi-square, CFI/TLI, RMSEA with its interval, SRMR and AIC/BIC where supplied. |
| `_parameters_excel.csv` | Measurement and regression parameters, residual/random-effect variances, and relevant covariances. Standardized SEM rows have their own SEs and intervals; mixed-model loadings and variance components use scientific rather than log-optimizer scales. |
| `_genetic_structure_excel.csv` | Genetic state and all pairwise implied liability correlations, including parent–child, spouse and parental-sibling pairs, with SEs/intervals. Visible and pi-corrected total liabilities are distinguished. |
| `_diagnostics_excel.csv` | Auxiliary one-factor statistic, df, p-value, sample size and warnings. |
| `_status_excel.csv` | Convergence, admissibility/inference status, estimation/SE method, warnings and extraction failures. |
| `_fits.rds` | Full fitted objects, detailed diagnostics and the reporting tables. |

**Mixed-model AIC:** random effects and clustered SEs do not themselves prevent
AIC reporting. The joint-covariate model estimates covariate moments separately,
so ordinary `AIC` is unavailable. Its numerical
`AIC_conditional_moments_fixed` is shown as **descriptive**, treating those
moments as fixed; random effects are integrated out. Compare only fits using
the same observations and scaling. The custom mixed likelihood does not fit
the reference models needed for SEM chi-square/CFI/TLI/RMSEA/SRMR, so those
entries state why they are unavailable. Zero-df SEM fit indices are flagged as
uninformative.

The genetic table retains every pair among included liabilities. With child,
mother, father and both parental siblings this gives ten visible correlations
and ten total-liability correlations per source fit, including mother–her
sibling, father–his sibling and cross-family pairs. Mixed child-sibling models
label the reference children `child` and `child_sib`. Correlation SEs use joint
parameter uncertainty conditional on pi; unavailable corrections or SEs have
explicit statuses. These are marginal genetic correlations, not correlations
of PGS measurements or standardized residual-covariance parameters.

Continuous workflow outcomes are standardized. Parameter tables identify
source units and standardized estimates separately. Source-SEM and genetic
structure SEs condition on the supplied pi; the primary corrected effects
retain their external `pi_se` adjustment. Fixed parameters, boundaries and
unavailable uncertainty remain labelled.

The large sampling-covariance and implied-observation matrices are not exported
by default. Use the saved fit when those details are needed:

```r
res <- readRDS("ravel_parental_phenotypes_mom_height_fits.rds")
attr(res, "reporting_tables")$fit_statistics
attr(res, "fits")  # parental source main/grid fit records

child <- readRDS("ravel_parental_deviation_Y_fits.rds")
attr(child, "report_fit_records")  # child source main/grid fit records
```

The optional parental setting `include_sem_audit = TRUE` additionally retains
the full SEM audit in the RDS and requires `ravel_sem_audit.R` beside the scripts.
It is off by default. Start with a fresh output prefix/folder when switching
from older scripts so earlier CSVs are not mixed with the new exports.

## Install the bundled package

Keep `RAVEL_0.2.0.9003.tar.gz` alongside these scripts. With this directory as your
working directory, install into a local library:

```r
dir.create("library", showWarnings = FALSE)
.libPaths(c(normalizePath("library"), .libPaths()))
install.packages("RAVEL_0.2.0.9003.tar.gz", repos = NULL, type = "source",
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
included artificial `synthetic_family_example.rds`, runs all 15 adult
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
  pi = .4, pi_se = .02,
  pi_grid = c(.3, .4, .5),
  models = list(estimated_state = list(fit_constraints = "estimated_state")),
  output_prefix = "maternal_height"
)
result$table
result$failures
result$sample_counts
```

This fits one joint genetic model per pi assumption and derives all 15 adult
regressions. The child remains available for
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

For a cross-sectional adult outcome, prepare one record per phenotyped adult. An
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

## Optional random effects for child outcomes

Edit the **SETTINGS** block in `ravel_parental_deviation.R`. Set your data file,
design, phenotype names, PGS columns, covariates and pi. Then choose:

```r
se_cluster <- "FID"               # Sandwich SE groups; also used for diagnostic sampling.
iid <- "IID"                      # Person ID for diagnostics and repeated measurements.
family_random_effect <- "M_ID"     # Nuclear-family intercept; NULL = off.
individual_random_effect <- FALSE # TRUE adds an IID intercept in longitudinal data.
interactions <- NULL              # c("age", "sex") fits/reports genetic interactions; NULL = none.
nuclear_family_id <- "M_ID"        # Shared-parent ID when the family effect is off.
```

For repeated continuous outcomes, set `design <- "longitudinal"`, specify `age`
and `age_center`, and turn on `individual_random_effect` if wanted. Age is
included as a fixed covariate. `interactions` independently controls genetic
slopes; `random_slope <- "age"` optionally adds an individual age slope and
requires the individual intercept. Otherwise leave `random_slope <- NULL`.

`family_random_effect` identifies nuclear families with the same two parents;
`se_cluster` can group several such families into an extended FID. Mixed fits
still need `nuclear_family_id` when the family intercept is off. Child PCs,
sex and age can stay in `covariates`; no confirmation flag is needed.

Cross-sectional runs with effects/interactions off retain the original
estimator. The editable longitudinal trio workflow uses the mixed fitter,
including when `interactions = NULL`. Existing direct function calls and their
advanced settings remain available. Parental-sibling longitudinal models retain
their existing age-only workflow: set `interactions <- age` and leave random
effects off. Binary outcomes support cross-sectional sibling models only.

The usual coefficient-table format is retained. Both workflows use the
Excel-ready reporting files described above. The RDS preserves every main/grid
fit and diagnostic:

```r
res <- readRDS("ravel_parental_deviation_Y_fits.rds")
attr(res, "reporting_tables")$fit_statistics
attr(res, "mixed_fit_records", exact = TRUE)
```

Keep the script beside `RAVEL_0.2.0.9003.tar.gz`, restart R and run it as before.
See [model assumptions and advanced API](../../docs/mixed-family-models.md) and
[validation and limitations](../../docs/mixed-family-validation.md) for details.

## The automatic one-factor check

The check compares the two PGS versions' phenotype associations in a separate
unadjusted shared-factor model. It is attempted once per main source fit, not once per projected
adult regression or pi-grid value. It does **not** test that the child's outcome
coefficient is zero. Covariates produce the explicit warning **"Covariates were
not used to compute the one-factor test."** Interactions are omitted too. A
valid result does not validate the adjusted model or its growth slopes.

The child script uses `iid` for the person and `se_cluster` for independent
family sampling. In the parental script and direct function calls, set
`subject_id` to the phenotyped person's ID. There, `diagnostic_cluster = "auto"` uses
the original `FID` when present, otherwise `family_id`. One individual per FID
and one target-age occasion are selected reproducibly before missing-score or
outcome exclusions. Missing selected observations are not replaced by another
person or visit. The selected rows and original/auxiliary sample sizes are saved.
Use `diagnostic_cluster = NULL` only for independent individuals.
In the mixed workflow, a custom diagnostic grouping must keep each primary
extended family together; splitting it returns an `invalid_design` diagnostic
while preserving the primary model fit.

Continuous outcomes use the established likelihood comparison. Binary outcomes
attempt a separate exploratory conditional-probit Wald calculation. Its
**nominal p-value is reported with a calibration warning**: null simulations
rejected 40/499 tests (8.0%) at a nominal 5% threshold. A successful calculation
has `status = "ok"` and `calibration_status = "uncalibrated"`; the warning appears
in the output notes and saved diagnostic. Primary binary mixed-model inference
is separate. Degrees of freedom are 3 for a
trio, 4 with one parental sibling and 5 with both; inspect the returned test
name and status. External pi is held fixed for this auxiliary test.

Both editable workflows now record `diagnostic_status`, `diagnostic_df` and
`diagnostic_N` in the CSV and preserve the full diagnostic object in the RDS.
An auxiliary attempt can still be unavailable for substantive or numerical reasons:

| Status | Meaning |
|---|---|
| `ok` | The named auxiliary calculation succeeded. Read its warnings; binary p-values remain nominal and uncalibrated. |
| `unsupported_free_measurement` | The fully free residual score model has an unidentified score-only alternative; no valid test is implemented for that model. |
| `needs_subject_id`, `needs_occasion`, `ambiguous_occasion` | Supply enough ID/visit information to select independent observations unambiguously. |
| `invalid_design`, `invalid_arguments` | Check supplied IDs, their nesting and the requested columns. |
| `insufficient_data` | Fewer than 50 complete independently selected observations remain. |
| `invalid_score_model`, `invalid_H0`, `invalid_H1` | A required fit failed identification or admissibility checks. Inspect the saved reason and rank diagnostics. |
| `fit_warning` / `fit_error` | Inspect the saved warnings or error from the diagnostic refits. |
| `invalid_probit_fit`, `invalid_influence`, `invalid_null_scale` | The binary diagnostic cannot provide regular inference; no p-value is substituted. |

A constrained sensitivity analysis can be checked separately, but changing the
measurement restrictions changes the model being tested. A missing p-value is
not evidence for or against one factor. A non-significant valid test also does
not establish unbiased genetic effects. Transfer the newly built package archive
with the updated scripts; older archives do not contain the auxiliary wrapper.

For primary clustering, large extended-family IDs are retained: a 70-person FID
is not a reason to discard clustering. Supplied family IDs must be observed,
nonempty and include at least two families. Binary fits on the original estimator
retain their stated clustering limitation; the optional mixed path has its own
full FID sandwich implementation.

[Parental-workflow validation and reproduction](parental-workflow-validation.md)
records the independent regression checks, one-factor reproductions, covariance
SE checks and full default run.
