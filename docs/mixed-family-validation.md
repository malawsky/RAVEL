# Mixed-family model validation

Development validation for the optional `ravel_mixed_fit()` path. The legacy
RAVEL estimator and its fitting scripts are not replaced. This report separates
numerical verification, Monte Carlo behavior under the model, and unsupported
applications. The checks establish behavior in the stated designs; they are not
a broad calibration claim across arbitrary pedigrees or missingness mechanisms.
The earlier Monte Carlo results below concern the original mixed likelihood
with independent genetic covariates. The new default joint-covariate extension
has separate validation records; those earlier runs do not establish its
coverage merely because both paths share genetic covariance components.

## Sparse sibling covariates (0.2.0.9003)

The first-stage self/sibling covariance estimates now come from a coherent
exchangeable working likelihood. The within-sibling and largest-family mean
blocks are estimated together, with exact shared columns represented on their
observed support. Estimated means, observed first-stage curvature and the
cross-stage influence all enter the existing FID sandwich. Marginal AM/pi
mapping and the original RAVEL estimator are unchanged.

A 10,000-family, 20-covariate simulation with 5% sibling pairs reproduced the
previous indefinite block (minimum eigenvalue -0.2176). The new within-sibling
block had minimum eigenvalue 0.0687. Family-weight perturbations matched the
reported influence to about 6e-12. Independent checks covered negative sibling
covariances, one/two/three-child families, and exactly shared covariates.

Three final-code integration runs had regular inference and finite effect SEs:

| Design | Families | Outcome rows | Checks |
| --- | ---: | ---: | ---: |
| Cross-sectional, family effect, factor IDs | 4,000 | 4,200 | 16 passed |
| Cross-sectional, family effect, exactly shared PC | 4,000 | 4,200 | 17 passed |
| Longitudinal, individual effect, shared PC, genetic sex interactions | 2,000 | 6,300 | 19 passed |

Each had 95% singleton families and 5% two-child families. These are integration
checks, not estimates of full-model coverage. The short editable-settings
workflow separately passed all 29 checks, including two fits and one-factor
results, under version 0.2.0.9003.

A separate 100-replicate first-stage covariance experiment used 2,000 families
per replicate, 5% pairs, two highly shared PCs, sex and age. The old empirical
subtraction failed positive definiteness in 98/100 replicates; the new estimator
completed all 100. Nominal 95% coverage was 94% for PC variance, 95% for sibling
PC covariance, 92% for within-sibling PC variance, and 96% for age variance.
The Bernoulli(.5) sex-variance result is nonregular at its maximum and is not
used as evidence of nominal calibration. This small experiment tests covariance
inference only; it does not establish full SEM coverage or robustness to
different covariate distributions by family size.

The regression suite passed 56 tests / 647 assertions, with no failures,
errors, test warnings or skips. Package check had no errors or warnings and
one existing CPU-runtime NOTE; tests were run separately. All 103 archived
R/help/test files matched the source checkout.

Results, every attempt and source hashes are in
`dev/mixed-validation/sparse-sibling-covariates/`. The authoritative final
folders are `full-fit-final`, `shared-fit-final`, `shared-long` and
`moment-mc-final`. An invalid shared-covariate simulation fixture is retained
with its correction note. An earlier run finished the same 100 covariance seeds
but hit a driver-file parsing error afterward; its clean rerun produced
identical estimate tables. These are 100 distinct seeds, not 200 replicates.
After these runs, the sole implementation edit qualified `utils::head` for the
package checker; its operation is unchanged.

Reproduce from the package root:

```sh
Rscript dev/check-sparse-sibling-covariates.R --output=/tmp/ravel-sparse-ordinary
Rscript dev/check-sparse-sibling-covariates.R --fit-scenario=shared --output=/tmp/ravel-sparse-shared
Rscript dev/check-sparse-sibling-covariates.R --fit-scenario=shared_long --output=/tmp/ravel-sparse-long
Rscript dev/check-sparse-sibling-covariates.R --fit=false --moment-replicates=100 --output=/tmp/ravel-sparse-moments
Rscript dev/check-negative-sibling-covariates.R
```

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

## Genetically associated child covariates

The default `covariate_model = "joint"` retains the marginal AM/pi model while
conditioning on all children's person-level covariates in a maternal family.
Let `V` be founder-scale visible genetics, `T` role-standardized total genetics,
and `C` the complete child-covariate stack. Write

```
A = Var(V)                 D = Cov(V, C)       Q = Var(C)
R = Cov(T, V) A^-1         Omega = Var(T) - R A R'
```

The identifying assumption is that `U = T - R V` is independent of `V` and `C`,
and PGS measurement errors are independent of `C`. Thus `Cov(T,C) = R D`.
Marginal pi alone does not identify this allocation of hidden-genetic
association; a different hidden-covariate assumption can produce a different
total-liability correction while having the same visible association.

The implemented conditional moments are

```
E(V | C) = D Q^-1 (C - E(C))
A_C      = A - D Q^-1 D'
E(T | C) = R E(V | C)
Var(T | C) = R A_C R' + Omega
```

The score/outcome mean offset is `[L E(V|C); B(C) E(T|C)]`, where `L` contains
the score loadings and `B(C)` contains the actual observation's genetic slopes.
The interaction-dependent offset is retained, including products induced by
age/sex interactions. An ordinary outcome covariate term alone would not
represent this model. Scores are still counted once per person, and a parent's
score mean depends on the complete family covariate stack rather than a visit.

For `n` children, the covariate covariance is
`Q_n = I_n (x) (Q_self - Q_sibling) + J_n (x) Q_sibling`, where `(x)` denotes
the Kronecker product. Genetic-covariate cross-covariances distinguish
mother-to-child, father-to-child, child-own and child-other-sibling associations.
Singleton-only data omit the unidentifiable sibling blocks. Means and raw
second moments are estimated from unique children and ordered within-family
child pairs; their per-family influence enters the **same FID sandwich** as
the likelihood scores. This propagates first-stage uncertainty and its
cross-covariance with genetic and outcome parameters. It does not treat the
empirical population covariate covariance as known.

PCs/sex must be stable within person. Explicitly allowed changing features
(by default supplied `age`/random-slope columns) use their mean over observed
supplied visits for genetic conditioning, while visit-level values remain in
the outcome design. A missing copy of an observed stable feature is reusable;
a wholly missing person feature is an error. Constant or exactly redundant
features are recorded and removed from the conditioning vector. No nearest-PD
replacement is silently applied to an invalid empirical covariate covariance.
Complete conditional features are required even for a child with missing Y;
that child still contributes scores and family covariate information.

The conditional genetic model assumes linear means and common conditional
covariances, with exchangeable child self/sibling covariate moments across
families. Sex can be binary; the code does not assert that binary sex itself is
Gaussian. These assumptions do not establish validity under arbitrary ancestry
structure, informative visits or selection, or covariance distributions that
change with family size. `covariate_model = "independent"` remains available
with its explicit independence confirmation. The no-covariate path and legacy
`ravel_fit()` are unchanged.

Under the current constrained Mendelian genetic block, adding siblings does
not change the existing rows of the visible-to-total projection: a parent's
projection uses the parental visible liabilities, and a child's projection
uses its own and its parents' visible liabilities. Thus the common marginal
total-covariate associations remain compatible across family sizes even with
free sibling cross-covariances. Algebra and numerical checks are retained in
`dev/mixed-covariate-validation/projectivity-review`. Conditioning on more
observed child covariates can change conditional means and covariances. The
positive-definiteness checks cover the family sizes in the fitted data; they
do not certify larger hypothetical families or an infinite pedigree.

Independent numerical checks are recorded in
`dev/mixed-covariate-validation/derivative-audit/checks.csv`. On that fixture,
relative information error was `4.16e-6`, maximum family-score discrepancy
`1.60e-8`, and maximum moment-gradient discrepancy `1.36e-7`. Direct weighted
family perturbations checked the moment influence to `8.75e-11`. A separately
assembled stacked covariance agreed to relative error `2.49e-5`.
The preparation/audit regression tests passed 124 assertions, including dense
Schur conditioning with interactions, weighted-moment influence, no-sibling
Jacobians, missing outcomes, redundant features, and the contribution of the
parameter--moment covariance to audit SEs. Audit tests also cover retained
moments without prepared data and binary boundary estimates with unavailable
regular SEs.

The final combined preparation and independent-oracle suite passed all 17
tests and 292 assertions. The existing optional mixed-model and diagnostic
regression suite passed all 32 tests and 318 assertions. Neither run had
failures, errors, test warnings or skips. These are separate from the earlier
tests of the unchanged original SEM estimator.

The new joint-covariate path also passed a separate external-pi sensitivity
check in all three constrained modes, using 400 two-child families with two
genetically associated PCs, sex and age, `pi = .4` and `pi_se = .03`. All nine
main/wider-step fits converged away from boundaries; all three main fits
returned regular inference. Across 33 coefficient and residual-variance
comparisons, stored sensitivities from pi steps of 0.001 agreed with full
refits at steps of 0.003: the largest discrepancy multiplied by `pi_se` was
less than 0.15% of the corresponding conditional SE. All nine coefficient
covariance matrices exactly matched the conditional FID covariance plus the
independent-pi variance addition. The sequential driver, every fit status,
derivatives and source hashes are retained in
`dev/mixed-covariate-validation/pi-sensitivity`. This validates the new joint
path separately from the earlier independence-model pi checks below.

The independent correlated-covariate generator preserves marginal AM/pi and
the stated hidden-projection assumption. Initial continuous fits with 1,500
families completed with regular inference in all three constrained modes;
their attempted-fit records are in
`dev/mixed-covariate-validation/continuous-smoke`. These are recovery and
numerical checks, not Monte Carlo calibration.

The subsequent fixed-seed Monte Carlo completed **270 original attempts**:
100 continuous cross-sectional, 100 with additional extended-FID dependence,
50 longitudinal, and 20 binary. All used `estimated_state`, fixed `pi = .4`,
founder-scale parental total variance 1.2 and total spouse correlation .3.
Each maternal family had two children, two genetically associated PCs and
PC-associated binary sex.
The continuous scenarios used 1,000 maternal families (2,000 outcomes); the
extended-FID scenario grouped three maternal families per FID and allocated
60% of the family residual variance to the shared FID component, preserving
the marginal maternal-family variance .25. The longitudinal and binary pilots
used 500 maternal families. The longitudinal design had four visits at ages
7, 9, 11 and 13, genetic-by-age and genetic-by-sex effects for all three roles,
and family/individual intercepts plus a correlated individual age slope.
Its true individual intercept variance, slope variance and covariance were
.2, .003 and .012. Binary outcomes used a probit residual variance of one.

All 200 continuous cross-sectional fits converged and returned regular
inference. The original longitudinal runs had 39/50 certified converged points
and 38/50 regular intervals: nine intercept--slope correlation boundaries,
one final-source gradient failure, and two historical gradient failures.
Binary had 20/20 converged points but only 17/20 regular intervals; seeds
541006, 541009 and 541015 reached the family-variance boundary. There were no
fit exceptions in the 270 original attempts. Convergence and regular inference
are separate quantities: a certified boundary point can lack regular Wald SEs.

The table reports present-day total-liability (`corrected_std`) effects.
Conditional coverage uses regular finite intervals; all-attempt coverage counts
an unavailable interval as uncovered and therefore also reflects inference
availability. These denominators are not interchangeable. MCSE is the observed
binomial Monte Carlo standard error in percentage points. All scales, point
biases, empirical SDs, average SEs, state recovery, and exact binomial 95%
coverage intervals are retained in `dev/mixed-covariate-validation/summary`.

| Scenario | Effect | Covered / regular | Conditional coverage (MCSE, pp) | Covered / attempted | All-attempt coverage (MCSE, pp) |
|---|---|---:|---:|---:|---:|
| Continuous | `beta_D` | 96/100 | 96.0% (2.0) | 96/100 | 96.0% (2.0) |
| Continuous | `beta_M` | 95/100 | 95.0% (2.2) | 95/100 | 95.0% (2.2) |
| Continuous | `beta_F` | 98/100 | 98.0% (1.4) | 98/100 | 98.0% (1.4) |
| Extended FID | `beta_D` | 94/100 | 94.0% (2.4) | 94/100 | 94.0% (2.4) |
| Extended FID | `beta_M` | 94/100 | 94.0% (2.4) | 94/100 | 94.0% (2.4) |
| Extended FID | `beta_F` | 92/100 | 92.0% (2.7) | 92/100 | 92.0% (2.7) |
| Longitudinal | `beta_D` | 36/38 | 94.7% (3.6) | 36/50 | 72.0% (6.3) |
| Longitudinal | `beta_M` | 37/38 | 97.4% (2.6) | 37/50 | 74.0% (6.2) |
| Longitudinal | `beta_F` | 36/38 | 94.7% (3.6) | 36/50 | 72.0% (6.3) |
| Longitudinal | `beta_D:age` | 38/38 | 100.0% (0.0) | 38/50 | 76.0% (6.0) |
| Longitudinal | `beta_M:age` | 37/38 | 97.4% (2.6) | 37/50 | 74.0% (6.2) |
| Longitudinal | `beta_F:age` | 37/38 | 97.4% (2.6) | 37/50 | 74.0% (6.2) |
| Longitudinal | `beta_D:sex` | 36/38 | 94.7% (3.6) | 36/50 | 72.0% (6.3) |
| Longitudinal | `beta_M:sex` | 37/38 | 97.4% (2.6) | 37/50 | 74.0% (6.2) |
| Longitudinal | `beta_F:sex` | 38/38 | 100.0% (0.0) | 38/50 | 76.0% (6.0) |
| Binary | `beta_D` | 14/17 | 82.4% (9.2) | 14/20 | 70.0% (10.2) |
| Binary | `beta_M` | 16/17 | 94.1% (5.7) | 16/20 | 80.0% (8.9) |
| Binary | `beta_F` | 15/17 | 88.2% (7.8) | 15/20 | 75.0% (9.7) |

For 38/38 coverage, the observed MCSE is zero, but the exact 95% interval is
90.75--100%; this does not establish perfect coverage. The 100-attempt
continuous estimates have coverage MCSEs of 1.4--2.7 percentage points.
Their child/maternal/paternal coefficient biases were
-0.0018/0.0020/0.0060 without extended-FID dependence and
-0.0094/0.0080/0.0079 with it. Absolute spouse-correlation and parental-variance
biases were at most .0027 and .0071, respectively. These are finite pilot
results under the stated model, not broad robustness guarantees.

The small random-slope pilot exposes an important limit: the correlation
boundary and gradient failures were retained, with Wald intervals withheld.
Among its 39 converged points, mean individual slope variance was .00328
(truth .003), individual intercept variance .19635 (.2), and family variance
.23355 (.25). These are conditional point-recovery summaries; variance-interval
calibration was not studied. The separate 150-family irregular-age fixture,
with varying person-level mean ages and the same genetic age/sex interactions,
passed the final natural-variance check and full stacked-SE calculation. Its normalized native
and natural-variance scores were 1.35e-7 and 5.97e-7, respectively; records are
in `irregular-long-replay`. One successful irregular fixture does not establish
longitudinal calibration across arbitrary visit schedules.

**The binary pilot does not establish nominal 95% calibration.** Child-effect
coverage was 14/17 (82.35%) among regular fits, with MCSE 9.25 percentage points
and exact 95% interval 56.57--96.20%; all-attempt coverage was 14/20 (70%).
Maternal and paternal coverage were 16/17 and 15/17. The small denominator,
three variance boundaries and observed child undercoverage remain explicit
limitations; successful optimization alone does not validate the intervals.
No boundary or final-source failure was retried to improve these counts.

Source provenance is retained per attempt. The initial 23 continuous,
seven longitudinal and one binary records were collected after the three-stage
optimizer but before the final natural-variance polish callback; they are
labeled `historical_pre_natural_polish`. Their exact launch-time hashes were
not retained by the older driver, as recorded in the immutable initial
checkpoint notes. Resumed batches write source hashes before fitting. On the
final source, continuous/FID inference was regular for all 177/177 attempts;
longitudinal had 34/43 converged points and 33/43 regular intervals, with nine
correlation boundaries and one gradient failure; binary had 19/19 converged
points and 16/19 regular intervals.

Only the two historical longitudinal failures, seeds 531005 and 531007, were
replayed from their saved points under the final callback; both then supplied
regular full SEs. Their separate `longitudinal-replays` records do not replace
the original 50-attempt denominator. Evaluating saved seed 531038 without
optimization gave a native score per family of 9.31e-8 but a natural slope-
variance score of 1.368e-5, above the unchanged 1e-5 gate. That final-source
failure remains unchanged. The earlier 20-attempt optimizer-development pilot
is retained separately in `pilot-continuous` and is not pooled into these
counts. Interrupted in-flight seeds were resumed and completed; interruption
resolution records preserve the original history. Simulation elapsed times
include machine contention and an unexplained 4,130-second binary outlier;
they are not runtime benchmarks.

The actual child-workflow collector passed **42 integration checks** in
`dev/check-mixed-covariate-workflow.R`, with retained records under
`dev/mixed-covariate-validation/workflow`. A cross-sectional fixture used 1,200
nuclear families within 600 extended-family clusters and fitted both the main
pi value and a full pi-grid refit. A separate longitudinal fixture used 200
nuclear families and 1,200 outcomes, with a custom `VisitAge` column,
genetic-by-age and genetic-by-sex slopes, and both family and individual random
intercepts. Both used custom nuclear-family, individual and outer-cluster IDs,
child PCs, and the default joint model without an independence-confirmation
flag. All three fits provided regular inference. The checks exercise automatic
person-mean age mapping, saved-fit serialization, new audit tables with finite
SEs, and export of the nested covariance matrices for every main/grid fit.

Both final workflow fixtures returned an auxiliary one-factor p-value with the
warning that covariates were omitted; the longitudinal result also retains the
interaction warning. An earlier 400-family cross-sectional setup produced
`invalid_H0` in its 200-person auxiliary sample despite regular main/grid mixed
fits. That failed auxiliary result is retained in
`initial-cross-200-auxiliary-diagnostic.csv`; it was not converted into a
p-value. The larger workflow fixture verifies successful reporting, not a
promise that every auxiliary fit converges or a diagnostic calibration claim.

The existing simple-settings driver was also rerun and passed **28 checks**,
including real family-only cross-sectional and individual-only longitudinal
fits. Its routing and status results reproduced the earlier records in
`dev/mixed-validation/simple-settings-check`; session metadata records the
rerun. These integration checks establish that the intended options reach the
estimator and survive reporting. They are not additional coverage simulations.

A fixed-covariate cross-sectional benchmark used 33,000 outcomes in 30,000
maternal families (27,000 single-child and 3,000 two-child families), with 20
child PCs, sex and age: 103 likelihood covariance parameters, 29 mean
coefficients and 528 estimated covariate moments. The complete fit with the
stacked FID sandwich finished in 290.659 seconds with regular inference; the
fit without SEs took 75.772 seconds. Founder-scale child/maternal/paternal
effects were 0.29484/0.13010/-0.05513 (truth 0.30/0.12/-0.04), with SEs
0.01261/0.01136/0.01070; estimated assortment was 0.3070 (truth 0.30) and
parental variance 1.2162 (truth 1.20). The prepared object used 222 MB and the
completed fit with SEs and `keep_data=FALSE` used 77.7 MB. This is a runtime
and recovery check for two shared covariance patterns, not an irregular-age
longitudinal interaction benchmark or a coverage study. Records and launch-time
source hashes are in `dev/mixed-validation/covariate-performance-with-se`;
subsequent no-covariate validation and optimizer-retry metadata/polish changes
did not affect this run's active path, which converged in its initial optimizer
stage with a maximum score per family of `8.58e-7`.

Joint-model audit tables retain both marginal genetic covariances and
conditional genetic/observed covariances for an actual retained family design.
The conditional tables include uncertainty in covariate moments. Marginal
observed covariance at a fixed interaction design is not an unconditional
average over ages, sex or PCs. Audit uncertainty remains conditional on pi.
Ordinary joint-model AIC is not supplied for the two-stage moment procedure:
`AIC` is `NA`, while `AIC_conditional_moments_fixed` is descriptive. Likelihood
parameter counts and first-stage moment counts are reported separately.

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

For the correlated-covariate extension, run for example:

```sh
Rscript dev/check-mixed-covariates.R --families=1500 --se=true --scenarios=continuous_pc,continuous_equilibrium,continuous_n_gen
Rscript dev/check-mixed-covariate-workflow.R
Rscript dev/check-mixed-covariate-performance.R --families=30000 --pcs=20 --no-fit
```

The first driver also exposes longitudinal interaction, binary, external-pi,
extended-FID and 20-PC scenarios; `--replicates` controls repeated attempts.
`--resume=true` retains every completed attempt, including failures, and runs
only pending seeds with new launch-time source metadata. The balanced
longitudinal Monte Carlo uses its retained `longitudinal/driver.R` variant.
`Rscript dev/summarize-mixed-covariates.R` regenerates the combined coverage,
source-status and state summaries. Saved-point score and variance-recovery
checks have separate dev drivers and do not refit or replace failed attempts.
Performance timings establish computational behavior for the specified design,
not inferential calibration. The historical validation below and its source
hashes remain a record of the earlier independence mode.

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

The subsequent editable-settings simplification passed 29 focused checks via
`Rscript dev/check-workflow-editable-settings.R`. Two simulated fits had regular
inference and successful one-factor diagnostics: a family-intercept model with
600 cross-sectional outcomes and an individual-intercept model with 1,800
longitudinal observations and genetic sex interactions. Both used custom IDs
and genetically associated covariates. Checks also cover explicit interaction
selection, fixed age adjustment, and the legacy parental-sibling age units.
Results and source hashes are in `dev/mixed-validation/editable-settings-check`.
These are workflow checks, not calibration simulations. In the new editable
block, `interactions = NULL` means no genetic interactions; longitudinal trios
use the mixed fitter even when all random effects are off. Direct function
interfaces retain their existing defaults.

## Remaining limits

This implementation assumes independent external pi, the specified genetic
measurement model, full siblings, and random residual effects independent of
genetic liabilities. Covariates/moderators must satisfy the selected mode's
assumptions: hidden-projection and measurement-error independence in the joint
mode, or genetic independence in the earlier mode. These simulations do not establish robustness to
violations, phenotype-dependent selection, binary longitudinal data, missing
PGSs, parental-sibling blocks or free genetic models. FID sandwich inference
does not remove bias from a misspecified genetic mean or dependence between
random effects and genotypes. Natural-scale state/genetic audit SEs condition
on pi; only the coefficient and residual-variance outputs currently propagate
external pi uncertainty.
