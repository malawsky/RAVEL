# RAVEL

**Recovering Assortment-Adjusted Variance Explained by Genetic Liabilities**

RAVEL estimates offspring, maternal and paternal genetic-liability coefficients
from parent–offspring trios, using two independently constructed polygenic scores
(PGSs) per person. It accounts for PGS measurement error and can correct estimates
for assortative mating and genetic variation not captured by the scores.

**Author:** Daniel Malawsky (<d.malawsky@gmail.com>), Wellcome Sanger Institute

## Installation

The development branch contains the constrained models and external-pi uncertainty
features described below:

```r
remotes::install_github("malawsky/RAVEL", ref = "dev")
library(RAVEL)
```

Requires R >= 4.0, lavaan >= 0.6-12 and MASS. To install a local checkout, run
`devtools::install()` from its directory.

## Quick start

With your trio data in `dat` and the column names listed below:

```r
fit <- ravel_fit(
  dat,
  pi = 0.40,
  pi_se = 0.05,
  fit_constraints = "estimated_state"
)

print_simple(fit)
results <- as.data.frame(fit, include_raw = TRUE)
```

This fits the constrained model without specifying a number of assortment
generations. It returns visible and total corrected estimates; the corrected
SEs and confidence intervals include the supplied uncertainty in pi.
No second correction call is needed.

`pi` is the fraction of **founder/random-mating genetic variance hidden from the
PGSs**, between 0 and 1 (excluding 1). It is not PGS measurement error or outcome
heritability. `pi_se` is its external standard error on that same scale. Omit
`pi_se` to treat pi as fixed. The older name `h` remains supported; if both `h`
and `pi` are supplied, they must agree.

To try the example without your own data, create `dat` first:

```r
dat <- ravel_simulate_trios(
  N = 10000, m = 200, hidden_fraction = 0.40, seed = 12345
)
```

This simulator produces one generation of assortment. The quick-start model
estimates the variance state rather than being told that history. The example
`pi_se` illustrates an external uncertainty input; the simulator does not estimate it.

## Input data

Use one row per trio, with two PGS versions targeting the same trait and
independent GWAS estimation errors. Apply each version's weights to all three
family members and orient both versions toward the same genetic target.

| Person | PGS columns expected by `ravel_fit()` |
|---|---|
| Mother | `PGS1_mother`, `PGS2_mother` |
| Father | `PGS1_father`, `PGS2_father` |
| Child | `PGS1_child`, `PGS2_child` |

The outcome column defaults to `Y`. Use the `col_pgs1_mother`,
`col_pgs2_mother`, and corresponding father/child arguments to map other names;
use `col_outcome` for the outcome. Raw or separately standardized PGSs can be
used. Their reliability is estimated, including a separate reliability for each
PGS version. Standardizing a score does not make its population variance known.

## Choose the assumptions explicitly

| `fit_constraints` | What is imposed during fitting? | Generation count |
|---|---|---|
| `"none"` | Original unconstrained SEM | Used only for the later total correction |
| `"equilibrium"` | Equilibrium transmission and stable within-score reliability across roles | Not supplied |
| `"n_gen"` | A specified constant-strength history from random-mating founders | Supply `am_n_generations` |
| `"estimated_state"` | Inheritance relations; variance growth estimated from the trio | Not supplied |

All three constrained modes fix within-PGS residual spouse correlation to 0
and child–parent residual correlations to 0.5. They assume independent errors
between the two PGS versions and estimate their reliabilities separately.
They do not assume the two PGSs are equally reliable. More restrictions can
improve precision when these assumptions fit the data.

```r
# Equilibrium
fit_eq <- ravel_fit(dat, pi = .4, pi_se = .05,
                    fit_constraints = "equilibrium")

# Three mating rounds from random-mating founders
fit_g3 <- ravel_fit(dat, pi = .4, pi_se = .05,
                    fit_constraints = "n_gen", am_n_generations = 3L)

# Original SEM, with an equilibrium total correction
fit_free <- ravel_fit(dat, pi = .4, pi_se = .05,
                      fit_constraints = "none", am_correction = "equilibrium")
```

**Defaults have not changed:** bare `ravel_fit(dat)` uses `pi = 0`,
`fit_constraints = "none"`, `am_correction = "n_gen"`,
`am_n_generations = 1L`, and no bootstrap. One generation means founder parents
and the children produced by their first assortative mating round. The default
SEM itself does not impose this generation history. Selecting a constrained mode
automatically selects its matching correction unless an incompatible one is
explicitly supplied, which is rejected.

For the assumptions, identification limits and equations, see
[Methods and assumptions](docs/methods.md).

## Read, save and plot results

```r
print_simple(fit)                         # compact table
print_ravel_results(fit)                  # detailed output
results <- as.data.frame(fit, include_raw = TRUE)
write.csv(results, "ravel_results.csv", row.names = FALSE)
saveRDS(fit, "ravel_fit.rds")

# Full corrected covariance and uncertainty details
fit$pi_correction
attr(fit$pi_correction, "covariance")
print(ravel_plot_report(fit$pi_correction))
```

Visible coefficients describe the genetic component measured by the scores.
Corrected coefficients include the assumed hidden component. Standardized total
coefficients are per current genetic SD; raw total coefficients use the founder
reference scale. Maternal and paternal coefficients are genetic-liability-correlated
associations, whose causal interpretation requires the model assumptions.

The ordinary print and table functions include external-pi uncertainty when
`pi_se` was supplied. Save RDS for further inference: CSV does not preserve joint
covariance matrices and other attributes. Inspect convergence and diagnostics;
failed or nonidentified fits are not suitable for ordinary confidence intervals.

## How uncertainty is handled

The default uses SEM-based SEs. Supplying `pi_se` adds its sampling variance by
the delta method, assuming the external pi estimate is independent of the trio
estimation errors. No Monte Carlo sampling of pi is used. Corrected intervals
are normal-approximation intervals and are not clipped to parameter bounds.
Finite-generation constrained fits after generation one also propagate how pi
affects the fit, using automatic perturb-and-refit sensitivity.

```r
fit_boot <- ravel_fit(dat, pi = .4, pi_se = .05,
                      fit_constraints = "estimated_state",
                      n_boot = 500, seed = 123)
print_simple(fit_boot, uncertainty = "bootstrap")
```

With external-pi uncertainty, corrected bootstrap-based results use the covariance
of joint bootstrap estimates plus external pi variance in a delta-Wald interval.
Visible bootstrap intervals remain percentile intervals. The underlying fixed-pi
SEM and bootstrap results remain available in the fit. See the
[reference](docs/reference.md) for uncertainty selectors and correlated external
pi estimates (`core_pi_cov`).

## Change pi after fitting

`ravel_correct_pi()` remains useful for sensitivity analyses. It recalculates
from the visible estimates and leaves the original fit unchanged:

```r
alternative <- ravel_correct_pi(fit, pi = .45, pi_se = .05)
print_simple(alternative)
```

This works without refitting for unconstrained, equilibrium and estimated-state
fits. A finite-generation constrained fit after generation one depends on pi:
refit with the new value in `ravel_fit()` instead. The separate correction
function rejects a mismatched pi for that case. The
[reference](docs/reference.md) also covers assumption grids and aggregate packets.

## Covariates and binary outcomes

Pass numeric covariates to `fixed_effects` so they are fitted jointly. Do not
residualize only the outcome when it correlates with genetic predictors. For
categorical covariates, the dataframe workflow below handles dummy coding.

```r
fit_adjusted <- ravel_fit(
  dat, pi = .4, pi_se = .05, fit_constraints = "estimated_state",
  fixed_effects = c("age", "PC1", "PC2"),
  covariate_mapping_confirmed = TRUE
)
```

Set `covariate_mapping_confirmed = TRUE` only after checking that the
visible-to-total mapping applies conditional on your covariates. Covariate
adjustment alone does not establish this assumption.

For binary outcomes, keep `Y` coded 0/1 and use `outcome = "binary_liability"`.
Results are on the SEM probit-liability scale. Baseline-relative fit indices
are unavailable for binary fits with covariates; inspect the diagnostic output.

## Dataframe reports and the bundled example

`ravel_dataframe_report()` prepares raw data, compares RAVEL with ordinary PGS
regressions, and can save tables and plots. Its default score names use lowercase
`pgs1_` and `pgs2_`, unlike `ravel_fit()`. It jointly fits numeric and categorical
covariates and standardizes continuous outcomes.

```r
data("ravel_example", package = "RAVEL")
out <- ravel_dataframe_report(
  ravel_example, outcome = "Y",
  covariates = c("age", "sex", "site", paste0("PC", 1:6)),
  categorical_covariates = c("sex", "site"),
  fit_constraints = "estimated_state", h = .4
)
out$model_status
out$summary
print(out$plot)

# Report summaries are visible-scale comparisons; correct the packet separately.
corrected <- ravel_correct_pi(
  out$packet, pi = .4, pi_se = .05, covariate_mapping_confirmed = TRUE
)
print_simple(corrected)
```

The confirmation is appropriate for the bundled simulation's covariate mechanism.
This dataset contains 4,000 simulated trios, founder pi 0.4, one generation of
assortment, and founder PGS reliabilities 0.5 and 0.3. See `?ravel_example` for the
full generating assumptions. Add `out_dir = "results"` to save a report in a new
folder. Keep the aggregate RDS packet for later corrections and inspect
`out$model_status` and `out$diagnostics` before interpretation.

## Further documentation

- [Methods and assumptions](docs/methods.md): theory, equations and limitations.
- [API and workflow reference](docs/reference.md): options, output fields, column
  mapping, reporting, uncertainty and simulation details.
- `?ravel_fit`, `?ravel_correct_pi`, `?ravel_dataframe_report`,
  `?ravel_correct_grid`: installed function documentation.

## License and citation

MIT. Please cite:

> Malawsky D. (2025). RAVEL: Recovering Assortment-Adjusted Variance Explained by
> Genetic Liabilities. R package. https://github.com/malawsky/RAVEL
