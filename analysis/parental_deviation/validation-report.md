# RAVEL family regression validation

The independent validation checks reduced adult regressions while retaining family members as genetic-model auxiliaries. The Monte Carlo results below contain 600 fits across 6 cells, with 0 fits having an estimation, correction or inference failure. These results assess the specified Gaussian pedigree model; they do not establish calibration under arbitrary assortative mating or measurement misspecification.

## Exact population checks

The fixture constructs parents, conditional family-of-origin means, sibling segregation and child segregation from independent Gaussian innovations. Two independent GWAS error processes have the specified pedigree residual correlations. The generator and oracle projections do not call RAVEL correction functions.

Exact covariance samples match the N-denominator covariance used by the likelihood. They check algebraic recovery, not sampling coverage. The test file covers all eight maternal adult subsets, trio-only and one-sibling adapters, child-free sibling/partner configurations, zero pi, negative and changing assortment, raw and standardized coefficients, variance explained, joint contrast covariance and external pi uncertainty.

## Monte Carlo design

Each cell uses 100 random samples of 2,000 families and fits the estimated-state model. Baseline pi is 0.4, current total-liability spouse correlation is 0.3, and parental total-liability variance is 1.15. PGS error variances are 0.4 and 1.6. The changed-assortment cell uses spouse correlation -0.2 and parental variance 1.5. The uncertain-pi cell draws an independent calibration estimate with SD 0.03 and supplies that SE to the model.

The phenotype has nonzero child, parental, sibling and family-of-origin contributions. Consequently, a reduced-regression target differs from a subset of the full-model coefficient vector. Oracle targets use the true total-liability covariance and outcome cross-covariance. The same seeds are shared between configurations for paired comparison; cells are not pooled as independent evidence. Preliminary ten-replicate smoke runs reuse the first ten seeds and are excluded from the counts and tables here.

| Genetic model and scenario | Replicates | Replicates with any failure | Coverage range across reported coefficients and contrasts |
| --- | --- | --- | --- |
| Partners plus both siblings | 100 | 0 | 0.91 to 0.97 |
| Partners plus both siblings; changed assortment | 100 | 0 | 0.93 to 0.96 |
| Partners plus maternal sibling | 100 | 0 | 0.92 to 0.97 |
| Mother and her sibling | 100 | 0 | 0.93 to 0.97 |
| Mother and her sibling; uncertain pi | 100 | 0 | 0.93 to 0.95 |
| Trio plus both parental siblings | 100 | 0 | 0.91 to 0.97 |

## Maternal sibling contrast

This table compares the same outcome regression, Y on maternal and maternal-sibling total liability, across genetic-model configurations. Bias and standard errors use founder-scaled raw liability units. Coverage refers to nominal 95% intervals among valid estimates; all failure counts remain visible above.

| Genetic model and scenario | Bias | Bias MCSE | Empirical SD | Mean SE | Coverage | Binomial 95% CI |
| --- | --- | --- | --- | --- | --- | --- |
| Partners plus both siblings | -0.0038 | 0.0062 | 0.0623 | 0.0608 | 0.97 | [0.915, 0.994] |
| Partners plus both siblings; changed assortment | -0.0037 | 0.0065 | 0.0654 | 0.0624 | 0.95 | [0.887, 0.984] |
| Partners plus maternal sibling | -0.0032 | 0.0062 | 0.0623 | 0.0609 | 0.97 | [0.915, 0.994] |
| Mother and her sibling | -0.0024 | 0.0063 | 0.0625 | 0.0612 | 0.97 | [0.915, 0.994] |
| Mother and her sibling; uncertain pi | -0.0011 | 0.0063 | 0.0633 | 0.0619 | 0.95 | [0.887, 0.984] |
| Trio plus both parental siblings | -0.0041 | 0.0062 | 0.0617 | 0.0606 | 0.97 | [0.915, 0.994] |

With 100 replicates, coverage near 0.95 has Monte Carlo SE about 0.022. The displayed binomial intervals show this uncertainty, including boundary cases. These runs can reveal gross implementation or calibration problems; they are too small to distinguish modest coverage differences or establish small precision advantages.

## Adding the child to the same four adult genetic model

Both configurations below use the same population, sample size and random samples, and target Y on all four adult liabilities. The only additional genetic-model member in the trio-assisted configuration is the child. These compare regression-estimator precision; variance-state SEs were not recorded, so no direct state-precision improvement is claimed.

| Genetic model | Parameter | Bias (MCSE) | Empirical SD | Mean SE | Coverage [binomial 95% CI] |
| --- | --- | --- | --- | --- | --- |
| Four adults only | Mother minus sibling | -0.0049 (0.0067) | 0.0670 | 0.0630 | 0.94 [0.874, 0.978] |
| Four adults only | Maternal coefficient | -0.0058 (0.0042) | 0.0420 | 0.0398 | 0.95 [0.887, 0.984] |
| Child plus four adults | Mother minus sibling | -0.0053 (0.0066) | 0.0664 | 0.0628 | 0.94 [0.874, 0.978] |
| Child plus four adults | Maternal coefficient | -0.0064 (0.0042) | 0.0416 | 0.0393 | 0.94 [0.874, 0.978] |

The baseline oracle own-minus-sibling contrast is 0.253652 in the sibling-pair regression and 0.280000 when both partners and both siblings enter the outcome regression. That difference is expected: these are different conditional associations. It is not evidence of an incorrect pi correction, and neither contrast is automatically a causal direct effect.

## Reproduction and source records

Run from the repository root. Full replicate records include random seeds, fitted pi, estimates, SEs, failures and warnings. Each cell's metadata text records R source hashes and session information. Files are under `dev/family-simulation-results/main/`; summary CSVs contain every coefficient and contrast, bias MCSE, coverage MCSE and binomial intervals.

The recorded hashes identify the code at each simulation's launch. Later review fixes added bootstrap-unit validation, support for fixed-zero legacy paths and reuse of pi-sensitivity refits. The estimated-state covariance formulas used in these simulation cells were unchanged; the final regression tests additionally cover those review fixes.

```sh
Rscript dev/run-extension-tests.R --all
Rscript dev/family-simulations.R --configuration=trio --reps=100 --n=2000 --output=dev/family-simulation-results/main
Rscript dev/family-simulations.R --configuration=siblings --reps=100 --n=2000 --output=dev/family-simulation-results/main
Rscript dev/family-simulations.R --configuration=adult_one --reps=100 --n=2000 --output=dev/family-simulation-results/main
Rscript dev/family-simulations.R --configuration=adult_both --reps=100 --n=2000 --output=dev/family-simulation-results/main
Rscript dev/family-simulations.R --configuration=siblings --pi-se=.03 --reps=100 --n=2000 --output=dev/family-simulation-results/main
Rscript dev/family-simulations.R --configuration=adult_both --scenario=changed --reps=100 --n=2000 --output=dev/family-simulation-results/main
Rscript dev/family-validation-report.R
```

The default seed is 20261029. The simulation source is `dev/family-simulations.R`; independent population construction and oracle calculations are in `tests/testthat/helper-family-validation.R`. Final exact and regression checks should be interpreted alongside the package check results, rather than replacing them.
