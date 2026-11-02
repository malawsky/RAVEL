# Optional family random effects

`ravel_mixed_fit()` is a separate development interface. The established
`ravel_fit()` models retain their existing behavior, and the analysis scripts
retain that estimator by default.
Use the new interface when the residual family or individual variation is
itself part of the model you want to fit.

## Start with your identifiers

Use one row per child for a cross-sectional outcome, or one row per child visit
for a continuous longitudinal outcome. Each row carries the two child PGSs and
the two PGSs for each parent, using the same six column names as `ravel_fit()`.
Scores must be complete and constant across visits. Parental scores must agree
within a maternal family. Shared scores are represented once per maternal-family
likelihood block. Dependence between blocks is handled through FID inference.

| Column | Purpose |
| --- | --- |
| `M_ID` | Maternal family, defining the family random intercept |
| `IID` | Child, defining repeated measurements and individual random effects |
| `FID` | Extended family, defining clusters for the full estimator's sandwich SEs |
| `F_ID` (optional) | Check that children grouped under the same mother share a father |

Several maternal families may belong to the same FID. Every IID must belong to
one maternal family and every maternal family to one FID. The current genetic
block assumes full siblings; a maternal ID alone does not establish this.
Supply `paternal_id = "F_ID"` to reject families with multiple or missing father
IDs, then apply your own explicit full-sibling restriction before fitting.
The function never silently removes half siblings or changes their family ID.
When a father is shared across different maternal families, his scores must
agree and those families must share an FID.

## Continuous sibling outcomes

```r
fit <- ravel_mixed_fit(dat,
  family_id = "M_ID", subject_id = "IID", cluster = "FID",
  paternal_id = "F_ID", pi = .4,
  fit_constraints = "estimated_state", family_random = TRUE)

as.data.frame(fit)
fit$variance_components
fit$inference_status
fit$one_factor
```

The child, maternal and paternal total genetic liabilities enter together.
The family intercept captures residual variation shared by siblings after these
terms. It is assumed independent of the genetic components and need not be
purely environmental. FID clustering allows residual dependence between the
maternal families within an extended family; it does not change which families
share a fitted random intercept. CR0 sandwich SEs are asymptotic and are not a
small-number-of-clusters correction.

## Repeated continuous outcomes and interactions

```r
fit_long <- ravel_mixed_fit(long_dat,
  family_id = "M_ID", subject_id = "IID", cluster = "FID",
  paternal_id = "F_ID", pi = .4, pi_se = .03,
  family_random = TRUE, individual_random = TRUE,
  random_slope = "age",
  interactions = c("age", "sex"),
  interaction_center = c(age = 10, sex = 0),
  covariate_mapping_confirmed = TRUE)
```

This estimates a family intercept, an individual intercept, an individual age
slope and its covariance with the individual intercept. Genetic-by-age and
genetic-by-sex terms are fitted for all three genetic roles, alongside the
moderators' main effects. A genetic-by-age coefficient is the change in the
genetic effect per year; the random age slope is residual variation in growth
between children. They are different model components.

Main genetic coefficients refer to age 10 and sex 0 in this example. Sex must
be a numeric coding such as 0/1. Continuous coefficients are returned in your
original phenotype units, and interactions per original moderator unit.
Internal score/outcome standardization is only for numerical stability.

`covariate_mapping_confirmed = TRUE` explicitly accepts independence of the
genetic components and score errors from the supplied covariates/moderators.
Do not use it as an automatic switch for ancestry-correlated or genetically
associated covariates. A small sample correlation does not establish the
assumption. Age-dependent selection or missingness may also violate it.

Family and individual intercepts need observed sibling families to separate
their variances. Individual effects require repeated outcomes; random age
slopes additionally require within-child age variation. The fitter checks these
conditions after missing outcome/covariate exclusions. Children with missing
outcomes can still contribute their complete scores. Missing-data mechanisms
and sparse visit schedules require substantive assessment; the simulations do
not establish validity under selective attendance.

## Binary sibling outcomes

```r
fit_binary <- ravel_mixed_fit(dat,
  outcome = "binary_liability", family_random = TRUE,
  paternal_id = "F_ID", pi = .4,
  interactions = "sex", covariate_mapping_confirmed = TRUE)
```

Binary outcomes must be numeric 0/1 and cross-sectional. Install `mvtnorm` for
this path. It fits a joint PGS likelihood and multivariate probit family
probability, with occasion residual variance fixed to one. Coefficients use
this **conditional probit scale**, which differs from the original fitter's
unit marginal liability scale. Binary longitudinal outcomes and individual
random slopes are unsupported. Deterministic integration is limited to at most
20 observed outcomes in a maternal family, with large families potentially
expensive. Very rare outcomes can leave family variance weakly identified.

An estimated variance at zero is a boundary solution. It is reported explicitly;
ordinary Wald SEs are not certified at that boundary. Inspect both
`fit$converged` and `fit$inference_status`, as optimizer completion alone does
not guarantee usable uncertainty. `fit$messages` records failed information,
gradient or pi-sensitivity checks.

## Inspect and save results

```r
audit <- ravel_mixed_audit(fit)
audit$fit_statistics
audit$state
audit$genetic_covariances
audit$implied_covariances
audit$implied_outcome_design
audit$parameters

write.csv(as.data.frame(fit), "effects.csv", row.names = FALSE)
write.csv(audit$fit_statistics, "fit-statistics.csv", row.names = FALSE)
write.csv(audit$genetic_covariances, "genetic-covariances.csv", row.names = FALSE)
write.csv(fit$variance_components, "residual-variances.csv", row.names = FALSE)
saveRDS(fit, "mixed-fit.rds")
```

Coefficient tables retain `visible`, `corrected_std` and `corrected_raw` labels.
The visible coefficients are marginal trio projections of the total model;
they are not regressions on estimated factor scores or on all siblings' scores.
The total scales denote one present-day role SD and one founder SD respectively.
Joint coefficient and variance-component covariance matrices are available for
contrasts. The audit also returns full internal parameter covariance and the
joint covariance of state and reference-family genetic moments.

`implied_covariances` adds the fitted score/outcome covariance matrix and its
SEs for the first retained covariance-design group. Inspect
`implied_outcome_design` for its children, moderator values and ages, or select
another group with `design_group`. This is a particular observed design pattern,
not an average over ages or sexes. Values use the input score and phenotype
units; for binary outcomes, Y denotes the latent probit response. These SEs
condition on pi. Fitting with `keep_data = FALSE` makes this extraction
unavailable without affecting the other audit tables.

Genetic covariance tables describe an exchangeable reference family with two
children by default; `n_children` changes this reference, not the fitted sample.
Their SEs, state SEs and internal parameter SEs condition on the supplied pi.
Coefficient and residual variance tables include independent external pi
uncertainty when `pi_se` is supplied and sensitivity refits succeed. This joint
likelihood depends on pi in every supported assortment mode, so pi sensitivity
requires full refits rather than only a post-fit coefficient transformation.

Log-likelihood and AIC refer to internally standardized data. Compare only fits
using the same observations and scaling. When maternal families within FID are
dependent, these describe a working-independence likelihood; robust SEs do not
make the AIC dependence-adjusted. Chi-square SEM fit, CFI and RMSEA are not
available for this custom family likelihood.

## Automatic one-factor diagnostic

The mixed fitter automatically attempts a separate unadjusted diagnostic.
Covariates produce the warning: **"Covariates were not used to compute the
one-factor test."** Interactions and random effects are omitted as well. This
checks the unadjusted cross-sectional score-association restriction, not the
adjusted growth model or its random effects.

For related families, it selects one child per FID using a reproducible seed;
for repeated outcomes, it selects one prespecified occasion or an observation
nearest the target age. Selection precedes complete-case exclusions and does
not use phenotype values. The result records selected rows and original and
auxiliary counts. Set `diagnostic_control = list(age = "age", target_age = 10,
selection_seed = 1)` to choose this explicitly.

Continuous outcomes use the established nested likelihood comparison. Binary
outcomes use a separate stacked-influence Wald test combining the score
measurement model and an unrestricted probit association model. Inspect the
test name, `status` and `reason`; a failed or unidentified test has no valid
p-value. Non-rejection is not evidence that all model assumptions hold.

The exploratory binary diagnostic has **not passed null calibration**: 40 of
499 valid simulated null tests rejected at a nominal 5% threshold. Its nominal
p-value is reported with a warning, `status = "ok"` for a successful numerical
calculation and `calibration_status = "uncalibrated"`. The warning is retained in
the saved diagnostic and table notes. Treat this p-value cautiously; primary
binary mixed-model estimation and uncertainty are assessed separately.

For an existing legacy fit, the same auxiliary workflow is available without
changing that fit:

```r
diagnostic <- ravel_score_diagnostic(dat, primary_fit = old_fit,
  family_id = "M_ID", subject_id = "IID", cluster = "FID",
  covariates = "sex")
```

The fully free, unidentified score-only model remains unsupported. Omitting
covariates cannot repair that identification problem, and the wrapper does not
silently substitute stronger genetic restrictions.

## Current scope

The mixed path supports `estimated_state`, `equilibrium` and `n_gen`, complete
two-GWAS scores, full-sibling child outcomes, and all three trio liability paths.
Its estimated-state parameterization requires total parental founder-scale
variance greater than 0.5. This is stronger than the legacy bound on visible
parental variance when pi is positive; agreement with the legacy fitter in a
regular state is not equivalence over the entire legacy parameter domain.
It does not yet implement the free genetic model, parental-sibling latent
blocks, arbitrary subsets of outcome liabilities, missing PGSs or binary
longitudinal outcomes. Those limitations do not remove functionality from
`ravel_fit()` or the existing adult-phenotype workflows.

Numerical and simulation validation is recorded separately in
[the validation report](mixed-family-validation.md). The editable
`ravel_parental_deviation.R` script exposes this path through optional
`random_family`, `random_individual` and `random_slope` settings. Leaving them
off and `genetic_interactions = NULL` retains the established estimator;
explicit interaction settings also select the mixed path. `maternal_id` groups shared parents,
`subject_id` identifies children and `family_id` supplies the outer sandwich
clusters. The advanced `mixed_settings` list also remains available. See the
[workflow quick start](../analysis/parental_deviation/README.md#optional-random-effects-for-child-outcomes).
