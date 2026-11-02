# Mixed-family model validation

Development validation for the optional `ravel_mixed_fit()` path. The legacy
RAVEL estimator and its fitting scripts are not replaced. This report separates
numerical verification, Monte Carlo behavior under the model, and unsupported
applications. The checks establish behavior in the stated designs; they are not
a broad calibration claim across arbitrary pedigrees or missingness mechanisms.

## What is being estimated

The likelihood represents shared parents and each child's scores once per
maternal family. Let `K` have diagonal 1, unrelated-spouse covariance 0 and
parent-child/full-sibling covariance 1/2. In founder genetic units, total
genetic covariance `T` has parental variance `Vp`, spouse covariance `rho*Vp`,
parent-child and sibling covariance `Vp*(1+rho)/2`, and child variance
`Vp*(1+rho)/2 + 1/2`.

The founder-scale visible genetic covariance is `(1-pi)*T + pi*K`, and its
cross-covariance with total liability is `sqrt(1-pi)*T`. Each independent GWAS
adds its own neutral-kinship measurement error, with role-specific loadings.
Outcome paths use role-standardized total liability. The full covariance of
scores `P` and outcomes `Y` is

```
Cov(P, Y) = Cov(P, T_std) B'
Var(Y)   = B Var(T_std) B' + R
```

Each row of `B` selects the observed child's liability and the shared parental
liabilities, with age/sex effects evaluated for that observation. For two
observations `i,j`, the optional residual terms are

```
R[i,j] = residual_variance * I(i=j) + family_variance
       + I(same child) * (individual_variance
         + intercept_slope_covariance * (age_i + age_j)
         + slope_variance * age_i * age_j)
```

Ages here are centered on the reported reference age. Disabled terms are
omitted. Binary outcomes fix occasion residual variance at one and integrate
the conditional multivariate-normal liability over the observed binary pattern.
The binary probability and the Gaussian PGS density are fitted jointly.

This retains the conditional hidden-genetic covariance, including its
dependence on genetic interaction slopes. Plugging estimated factor scores into
an ordinary mixed model would not fit this likelihood.

Continuous means are profiled by generalized least squares. Inference restores
the full mean/covariance information matrix and aggregates the full family
scores at FID before forming the sandwich. Independent external pi uncertainty
is propagated by full sensitivity refits in all three assortment modes.

## Numerical and independent checks

- An independently written simulator draws visible and hidden genetic
  components from primitive covariance equations; it does not call the fitter
  or its correction functions.
- Its full repeated-family covariance, including age and sex interactions,
  correlated individual intercept/age slope and score scales, agrees with the
  fitted covariance equations to `1e-10`.
- With one child per family and no random effects, the new total/visible/raw
  coefficients agree with `ravel_fit()` to `9.3e-6` in a separate fitting check.
  This checks a regular state, not global parameter-domain equivalence: the
  mixed estimated-state model bounds total parental variance above 0.5, whereas
  the legacy constrained syntax bounds visible parental variance.
- Direct unprofiled likelihood calculations reproduce the continuous full
  information matrix (relative error `9.24e-7`), family scores (`3.39e-10`),
  and FID sandwich (`2.68e-5`) on an independent fixture.
- Compressing families that differ only in fixed-effect mean designs preserves
  likelihood, gradient and fitted means. It reduces redundant covariance
  calculations without treating different covariate values as identical.
- Binary orthant derivatives are checked against direct probability
  differences in dimensions 1--4. Full family likelihood gradients are also
  checked by independent differences, including the family-variance derivative
  at zero.
  An independent binary information matrix agrees to relative error `1.75e-5`;
  the manually reconstructed extended-FID sandwich agrees to `1.05e-14`.
- Independent wider-step pi refits agree with the default sensitivity in all
  three assortment modes; the largest discrepancy times `pi_se` is `0.000423`
  of the corresponding coefficient SE in the saved check.
  An interior binary fixture also passes a wider-step pi sensitivity check:
  the largest derivative discrepancy times `pi_se=.02` is `0.0001304` effect
  units, small relative to its coefficient SEs.
- Outcome rescaling and independent PGS rescaling from `1e-5` to `1e3` preserve
  the genetic state and appropriately transform coefficients, SEs and residual
  variances. Internal likelihood values remain invariant to this rescaling.

The deterministic tests also check identifier nesting, score constancy across
visits, design identification after missing outcomes, and diagnostic sampling
before exclusions. They distinguish optimizer success from regular inference.

## Boundary checks

Binary family variance is optimized directly on `[0, infinity)`. A reproducible
100-family fixture reaches exactly zero with a valid one-sided optimum and
explicit `boundary` status; ordinary Wald SEs are withheld.

Continuous fitting retains log-SD coordinates but additionally checks analytic
natural-variance derivatives. A true-zero 200-family fixture otherwise stopped
at internal variance `1.087e-6` with optimizer success, despite a natural-variance
score of `21.68`. The added certification flags this solution and withholds
ordinary SEs. True-zero and near-zero fixtures and correlated random-slope
derivative checks are retained as regression tests. This is honest detection
of nonregular inference, not a boundary confidence-interval method.

## Initial continuous pilot

Two scenarios used 50 independent replicates each: the baseline family model
and three maternal families per FID with additional cross-family residual
dependence. All 100 fits converged with finite inference. Present-day total
coefficient coverage was 94% for all three baseline paths and 92%, 94%, 94%
for the child, father and mother paths under extended-family dependence.
Coverage MCSE is approximately 3.4--3.8 percentage points at this pilot size.
All six coefficient biases were within 1.7 Monte Carlo SEs of zero.

These pilots support larger checks; they do not establish calibration for rare
binary outcomes, sparse longitudinal designs or a small number of FIDs.
The attempted-fit and coefficient-level records are retained in
`dev/mixed-validation/pilot-continuous`.

## Larger simulations and focused pilots

| Scenario | Replicates | Regular inference | 95% total-liability coefficient coverage |
| --- | ---: | ---: | --- |
| Continuous baseline, 1,000 maternal families | 500 | 500 | Child 95.8%, mother 95.2%, father 96.4% |
| Extended-FID dependence, 1,000 maternal families | 500 | 500 | Child 96.2%, mother 95.2%, father 93.6% |
| Age + sex interactions and individual age slope, 500 families | 50 | 49 | 88--98% across nine coefficients, counting every attempt |
| Binary siblings, 500 families | 20 | 20 | Child 95%, mother 95%, father 90% |

All scenarios had two children per maternal family. Continuous cross-sectional
truth used founder coefficients `(0.30, 0.12, -0.04)`, pi 0.4, total spouse
correlation 0.3, parental variance 1.2, family variance 0.25 and occasion residual
variance 0.75. Extended FIDs contained three maternal families, with 60% of each
family effect's variance shared across the FID. This preserves its marginal
variance and induces cross-family residual covariance 0.15. Genotypes of
different maternal families were independent: this tests that specified residual
dependence, not an arbitrary KING-derived genetic pedigree.

All 1,000 continuous fits provided finite inference. Coverage MCSE ranges from
0.83 to 1.09 percentage points. All six genetic coefficient biases were within
one Monte Carlo SE of zero. Family variance showed modest downward bias:
`-0.00612` (MCSE `0.00159`) and `-0.00500` (MCSE `0.00178`) against truth 0.25,
with coverage 95.6% and 93.8%. This finite-sample variance bias is not hidden by
describing the estimator as universally unbiased.

The longitudinal pilot used four scheduled ages (4, 8, 12, 16), individual
intercept variance 0.2, age-slope variance 0.003 and covariance 0.012. Founder
age slopes were `(0.025, -0.008, 0)` and sex slopes `(0.12, -0.06, 0.02)`.
All 50 fits converged; one reached a random-effect correlation boundary and
withheld ordinary inference. Conditional coverage among 49 regular fits was
89.8--100%, with non-extreme coverage MCSE 2.0--4.3 percentage points. The largest
conditional bias was about 2.1 MCSE. The zero father-age coefficient rejected
0/49 times (exact 95% interval 0--7.25%), which is not evidence of zero type-I
error. These are pilot results, not precise interaction calibration.

The binary pilot used 15% reference prevalence and occasion residual variance
one. All 20 fits had regular inference; genetic biases were within 1.12 MCSE.
Coverage MCSE was 4.9--6.7 percentage points, and the exact interval around 90%
coverage was 68.3--98.8%. A separate 300-family binary sex-interaction workflow
check verifies routing, contrasts, scale labels and audits, not calibration.

The large runs retain their starting source hashes. Later natural-variance
gradient certification changed the acceptance gate, not the optimizer or
likelihood. Twenty selected seeds were replayed on the final code: first/last,
boundary, largest-gradient and smallest-variance cases. All statuses were
preserved (19 regular, one boundary) and coefficient differences were below
`5e-16`. This selected recertification is not a replay of every historical fit.
Exact counts, bias including all finite converged point estimates, coverage,
Monte Carlo errors and unavailable-inference denominators are saved in
`dev/mixed-validation/consolidated`.

## Auxiliary diagnostic

The unadjusted diagnostic is fitted on one reproducibly sampled individual per
FID and one selected occasion. Covariates and interactions are omitted with
warnings; non-rejection does not validate the adjusted mixed model.

An initial 500-replicate binary null experiment exposed a real numerical defect:
109 restriction-covariance failures and one score-model convergence failure.
Only 390 tests were valid, and their 1.79% rejection rate was below nominal.
The cause was an unstable choice of score rows for the Wald restrictions,
with condition numbers around `3e7`. A norm-pivoted QR basis removes these
numerical failures, but does **not** resolve inferential calibration: after the
repair, 499/500 calculations were numerically valid and 40 rejected. Their rejection rate
is 8.02% (Monte Carlo SE 1.22 percentage points; exact 95% interval 5.79--10.76%),
or 8.0% across all attempts. One score-model convergence failure remains.

A known-measurement oracle on the same independently regenerated 500 datasets
rejects 26/500 (5.2%; exact 95% interval 3.42--7.53%). The excess tail appears
after measurement estimation/restriction transformation; a coding error versus
finite-sample effect has not been established. No empirical correction has been
tuned to these results. Public binary diagnostic p-values are reported as
nominal asymptotic results with an explicit calibration warning. Numerical
success has `status = "ok"`; `calibration_status = "uncalibrated"` and the saved
warning distinguish it from evidence of calibrated inference. The fitted
models, statistic and calibration results are retained for investigation.

The QR-calibration run loaded the repaired restriction basis before the final
singular-value rejection guard and audit fields were added. Its Wald calculation
matches the current private calculation; it is not represented as an exact
final-source run. Current-source regression tests exercise previously failing
seeds and confirm that the public p-value is returned together with its warning.
Original failed results,
repaired-basis results, oracle results and their provenance are retained in
`dev/mixed-validation/binary-diagnostic-*`.

## Reproduce

Run from the repository root with `pkgload`, `testthat` and optionally `mvtnorm`
installed:

```r
pkgload::load_all()
testthat::test_local(filter = "mixed|score-diagnostic|binary-score-signal")
```

The standalone drivers are `dev/run-mixed-validation.R`,
`dev/check-mixed-nesting-pi.R` and `dev/recheck-mixed-certification.R`.
Saved output includes every attempt, unavailable inference, conditional and
all-attempt coverage, Monte Carlo uncertainty, session information and source
hashes. The large simulations retain their starting source versions; replay
checks distinguish later stricter convergence gates from changes to estimates.

The pre-existing test suite completed 252 tests and 4,316 passing expectations,
with no failures, errors, warnings or skips. It was started before the mixed
modules were added; the original R implementation files remain unchanged.
Separate new-module and integrated-workflow checks are recorded under
`dev/mixed-validation/checks` and the workflow validation directories. The
updated nominal-p reporting passed all 318 assertions in 32 new-module tests;
the current suite summary is in `checks/reported-p-tests.csv`, alongside the
earlier per-test record. The binary child-workflow check verifies the p-value
and calibration warning in printed output, CSV notes and the saved RDS.
All 36 integrated workflow checks and 28 simple-setting checks passed. The
latter includes actual simulated fits with custom ID columns: a family-only
cross-sectional model (500 outcomes, 250 maternal families) and an individual-only
longitudinal model (1,500 observations, 250 maternal families), with automatic
genetic age slopes retained. Both converged with regular inference; these are
workflow checks, not additional Monte Carlo calibration. Results and source
hashes are in `workflow-check` and `simple-settings-check`.
The optional dispatch retains the original estimator by default; the intentional
script changes are the auxiliary diagnostic and stricter supplied-family-ID
handling.

## Remaining limits

This implementation assumes independent external pi, the specified genetic
measurement model, full siblings, and random residual effects independent of
genetic liabilities. Covariates/moderators must satisfy the stated genetic
independence assumption. These simulations do not establish robustness to
violations, phenotype-dependent selection, binary longitudinal data, missing
PGSs, parental-sibling blocks or free genetic models. FID sandwich inference
does not remove bias from a misspecified genetic mean or dependence between
random effects and genotypes. Natural-scale state/genetic audit SEs condition
on pi; only the coefficient and residual-variance outputs currently propagate
external pi uncertainty.
