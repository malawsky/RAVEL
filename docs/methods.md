# RAVEL methods and assumptions

This note records the measurement, inheritance and uncertainty assumptions used
by the implementation. For installation and ordinary use, start with the
[README](../README.md); detailed examples and output fields are in the
[API and workflow reference](reference.md).

## Genetic components and coefficient scales

RAVEL is a trio structural equation framework for estimating direct offspring and parental genetic-
liability-correlated variance explained by leveraging two independently constructed polygenic scores
per individual as multiple noisy indicators of a common latent genetic liability.
Two independently constructed PGSs serve as noisy indicators of a
role-standardized visible genetic liability for each mother, father, and child.
The fitted SEM returns visible-scale coefficients
**b** = (b_D, b_M, b_P). Under explicit assumptions about additive transmission,
hidden genetic variation, and the common assortment pattern across visible and
hidden components, a post-estimation matrix transformation returns coefficients
on both the raw total-liability scale beta^raw and the role-standardized
total-liability scale beta.

- **b_D**: child visible-liability coefficient conditional on parental visible liabilities
- **b_M / b_P**: maternal and paternal visible-liability coefficients
- **rho_G**: visible-scale spouse correlation
- **eta_m / eta_p**: child-mother and child-father visible covariances
- **t_m / t_p**: standardized child-on-parent projection coefficients recovered from rho_G, eta_m, and eta_p

The measurement model corrects for PGS unreliability via ORIV-as-SEM, while
cross-person residual covariances account for shared GWAS estimation noise across
family members.

### Notation used in the theory note and R API

The supplementary theory note uses p/P for the paternal role, while the R API
retains f/F in field names for backward compatibility. Thus `bf`, `eta_f`, and
`tf` correspond to b_P, eta_p, and t_p, respectively. The legacy `rhoMF` field
is the mother-father entry.

The API argument `pi` (also accepted as the legacy name `h`) is the fraction of founder- or
random-mating additive variance assigned to the hidden component.

## Constrained measurement and inheritance models

All constrained modes retain a separate unknown reliability for each PGS. They fix
within-PGS residual spouse correlation to zero and child-parent correlations
to 0.5, and assume independent errors across PGSs. The scores measure the same
visible genetic component and their errors are orthogonal to that component
and the outcome. Scores should be oriented toward the same genetic target.
Unknown marginal score variances remain nuisance parameters, so separately
standardizing the scores does not supply known population variances.

At equilibrium, each score's reliability is equal across family roles and
both latent child-parent correlations equal `(1 + rho_visible) / 2`. The two
scores' reliabilities are not equated to each other. With two scores and no
covariates there are 13 covariance parameters, compared with 19 when only the
residual correlations are fixed.

For finite generations, child reliability changes with the modeled genetic
variance. The model assumes random-mating founders, Mendelian transmission,
and a constant total-liability spouse correlation across mating rounds. The
generation count is supplied, not estimated; generation one means founder
parents and their first offspring generation. This is a specified finite
assortment history, not arbitrary time-varying assortment. Fits after generation
one depend on the supplied founder hidden fraction `h`. Changing that fraction
or the assumed history requires refitting.

With the same GWAS weights and constant raw error variance, reliability evolves
with the visible signal variance. Let `R_j,n` be the reliability of score j in
generation n, and `B_n = V_G,n+1 / V_G,n` the visible-variance growth factor. Then

```
R_j,n+1 = B_n * R_j,n / (1 + (B_n - 1) * R_j,n)
R_j,n+1 - R_j,n = (B_n - 1) * R_j,n * (1 - R_j,n) /
                  (1 + (B_n - 1) * R_j,n)
```

The two scores share the growth factor, but retain their own reliabilities.
With generation zero denoting random-mating founders, constant total-liability
spouse correlation r, and founder hidden fraction pi, the model gives
`V_G,n = pi + (1 - pi) * (1 - r * ((1 + r)/2)^n) / (1 - r)`.
Positive assortment increases reliability toward equilibrium; negative
assortment decreases it. Changes tend to zero at equilibrium. Separately
standardizing each observed score to variance one leaves reliability unchanged.
The API's `am_n_generations = n + 1` compares parents in generation n with
children in generation n + 1.

These restrictions can improve precision when their assumptions hold. The
package's correction uses a supplied hidden fraction and estimates the total
spouse correlation. The earlier pooled-score comparison instead conditioned
on a known total spouse correlation; its numerical SE gains are not a promise
for this different uncertainty budget.

## Estimating the variance state

The `fit_constraints = "estimated_state"` model estimates one more covariance parameter than the equilibrium or
specified-generation modes: 14 parameters with two scores and no covariates.
The two PGSs retain separate unknown reliabilities. Parent-child correlations
identify the signal variance growth, even if each observed PGS was separately
standardized. Let `rho` be the fitted visible spouse correlation and `eta` the
common fitted child-parent correlation. On the scale where founder visible
variance is one,

```
B = ((1 + rho) / (2 * eta))^2
V_G,parent = 1 / (2 * B - 1 - rho)
V_G,child = B * V_G,parent
V_T,parent = (V_G,parent - pi) / (1 - pi)
rho_total = rho * V_G,parent / (V_G,parent - pi)
```

The founder anchor comes from Mendelian segregation variance equal to half the
founder variance. Correction also assumes random-mating founders with independent
visible and hidden components, and assortment along the same total genetic
target over time. Founder pi remains external. Pre-existing founder LD, a changing
genetic target, or a different segregation model are not identified by these data.

The estimated-state model is a smooth relaxation of a constant-strength history:
it estimates the accumulated variance and current assortment rather than selecting
an integer number of generations. It also accommodates changing assortment strength
when the target and segregation assumptions remain valid. Diagnostics indicate
whether the fitted state lies on a constant-strength founder-to-equilibrium
trajectory and, when defined, its effective duration. A fractional effective
duration is not an integer-generation estimate. Near equilibrium or random mating,
duration can be unidentifiable even when correction is estimable. The fit does not
force noisy estimates onto the founder or equilibrium boundary.

Fitting this mode is independent of pi; an incompatible supplied pi is rejected
by post-estimation correction. The full fitted covariance propagates uncertainty
in growth and founder scaling, and `ravel_correct_pi()` adds external pi uncertainty
without refitting. Its `estimated_state` attribute includes a table of estimates,
SEs and intervals for growth, parental visible and total variance, and the total
spouse correlation, together with their joint covariance. It does not report a
Wald SE for effective duration. Precision can decrease relative to a correctly specified
generation count. Low score signal and nearly singular reconstructed states can
make inference unstable. The model does not promise identifiable separate
reliabilities when both assortment and all score-outcome associations vanish.
Estimated-state fits explicitly check the rank of the observed-moment derivative
along the fitted constraint surface and the full parameter covariance before
allowing inference. Nonidentified fits and active inequality boundaries are
flagged; their correction packets cannot be used for ordinary Wald inference.

## Uncertainty in the external hidden fraction

`ravel_fit(..., pi = estimate, pi_se = standard_error)` fits the model and
applies the same uncertainty propagation as `ravel_correct_pi()` in one call.
`pi` and the legacy name `h` denote the same quantity. Supplying different
explicit values for them is an error. Omitting `pi_se` retains fixed-pi
inference. The separate correction function remains useful for subsequent
sensitivity analyses and aggregate packets.

The external variance is `pi_se^2`. No pi Monte Carlo draws or CI-to-SE
conversion are needed. In `ravel_correct_pi()`, supply exactly one of `pi_se`
or `pi_ci`; `pi_se = 0` adds no external uncertainty. `ci_level` controls the returned effect intervals.
The method uses a local delta approximation and does not clip confidence
intervals to parameter bounds.

For compatibility, `ravel_correct_pi()` also accepts `pi_ci = c(lower, upper)`.
By default it is interpreted as a 95% normal-Wald interval and must be symmetric about the
estimate. For a logit-Wald interval, use
`pi_ci_scale = "logit"`; the interval must be symmetric on that scale. An
arbitrary percentile or asymmetric interval does not uniquely specify a
sampling variance and is rejected. `pi_ci_level` describes the supplied
interval. These CI-specific settings apply only with `pi_ci`. Zero external
variance reduces model-based inference to the fixed-pi delta method.

The external pi estimate is assumed independent of the target-trio estimator
unless `core_pi_cov` supplies their covariance. For pi-dependent fits, this
must describe covariance with SEM estimation noise at fixed pi; supplying the
already marginalized plug-in core covariance would double-count pi uncertainty.
The joint propagation retains covariance among all six visible core parameters.
For a live finite-generation fit it also uses perturb-and-refit sensitivity,
so the result includes uncertainty in both fitting and correction. The pi
point estimate must match the value used in that finite fit.

Visible-transfer packets are accepted too. Pi-dependent aggregate packets
need a separately supplied, validated `core_pi_sensitivity`, because they do
not contain individual data for refitting; use the live fit for automatic
sensitivity. Bootstrap inputs use the covariance of complete joint core draws
and return explicitly labeled delta-Wald intervals, not percentile intervals.
Existing fixed-h fitting, correction and reporting APIs retain their meaning.

## Assortative mating correction for the total-liability scale

RAVEL assumes that the aggregate cross-mate covariance projected into the
modeled visible-hidden component space has a common direction aligned with the
components' covariances with the total liability. Other assortment axes may
exist, provided that their projections onto the modeled visible and hidden
components are proportional to their projection onto T*.

The `am_correction` parameter controls the raw covariance state used in the
visible-to-total mapping, including raw AM variance inflation when `h = 0`:

- **`"n_gen"`** (default): Propagates the raw visible-hidden covariance state for a finite number of mating rounds under a constant standardized total-liability spouse correlation. With `am_n_generations = 1`, the mating parents are founders, so the parental within-person visible-hidden covariance is zero. The child generated by that mating round can nevertheless have nonzero visible-hidden covariance.

- **`"equilibrium"`**: Uses the fixed-point raw covariance state after indefinitely repeated assortment at constant strength.

- **`"estimated_state"`**: Uses the founder-scaled variance state inferred during an estimated-state fit; it requires matching fitting constraints. See the variance-state equations above.

```r
# One mating round; founder parents have d_parent = 0
fit <- ravel_fit(dat, pi = 0.5, am_correction = "n_gen", am_n_generations = 1)

# 15 generations of AM
fit <- ravel_fit(dat, pi = 0.5, am_correction = "n_gen", am_n_generations = 15)

# Equilibrium (infinite generations)
fit <- ravel_fit(dat, pi = 0.5, am_correction = "equilibrium")
```

## Analytical expected values

`ravel_expected_values()` computes the population quantities analytically implied by a set of total-liability parameters under a specified AM history. This is useful for verifying what RAVEL should recover in expectation, checking simulation calibration, and understanding the scale transformations between visible and total liabilities.

```r
ev <- ravel_expected_values(
  beta_total = c(bD = 0.45, bm = 0.20, bf = 0.00),
  beta_total_scale = "raw",
  h = 0.4,
  rho_T_corr = 0.3,
  tm = 0.5, tf = 0.5,
  am_correction = "equilibrium"
)
ev$beta_vis_expected          # Expected visible-scale b
ev$beta_total_recovered       # Recovered raw beta^raw
ev$beta_total_recovered_std   # Recovered standardized beta
ev$Sigma_xh                   # Code name for the note's Sigma_GH; X is visible G
ev$var_total_expected         # Var(T*) by role, not Var(T)
```

With `beta_total_scale = "raw"`, `beta_total` is beta^raw and multiplies T*.
With `beta_total_scale = "standardized"`, it is beta and multiplies the
role-standardized T. The argument `rho_T_corr` is rho_T, the standardized
total-liability spouse correlation. In this analytical function, `tm` and `tf`
are raw maternal and paternal transmission coefficients; `tf` corresponds to
the note's paternal t_p notation. Both remain `0.5` under the maintained
additive model.

`ravel_expected_pgs_ols()` extends this to the observed-score scale. By
default, `reliability` is interpreted on the random-mating visible-scale,
matching the simulation convention where PGS noise is calibrated before
assortative mating inflates role-specific visible signal variance. Use
`reliability_scale = "observed"` if the values are already final reliabilities
for standardized PGSs (`loading^2 = reliability`). By default it assumes
spousal residual correlation `0` and parent-child residual correlations `0.5`.

```r
pgs_ev <- ravel_expected_pgs_ols(
  beta_total = c(bD = 0.20, bm = 0.05, bf = 0.00),
  beta_total_scale = "raw",
  h = 0.8,
  rho_T_corr = 0.7,
  reliability = c(child = 0.20, mother = 0.20, father = 0.20),
  reliability_scale = "random_mating",
  am_correction = "equilibrium"
)

pgs_ev$beta_pgs_expected   # expected standardized trio OLS on observed PGS
pgs_ev$Sigma_pp            # expected observed-score covariance block
```

## Parental siblings

`parental_sibs` adds a full sibling of the mother (a_m) and/or the father (a_p).
Each sibling has the same one-factor measurement model and shared-weight
scores as the trio. In the free model (`fit_constraints = "none"`):

- Each sibling adds three latent correlations: kappa (parent with own sibling),
  eta_a (child with aunt/uncle) and chi (parent with the spouse's sibling).
  With both siblings, zeta (between the siblings) is also estimated: seven
  additional correlations in total.
- Within each score version, residual covariances between each sibling and
  every other family member are free. For outcrossed families with unassorted
  GWAS noise, the expected residual correlations are 1/2 with the sibling's own
  parent, 1/4 with the child and 0 with unrelated members. The constrained
  modes fix them at these values.
- The outcome model is
  `Y = b_D G_c + b_M G_m + b_P G_p + b_MS G_am + b_FS G_ap + e`.
  Writing the parent and sibling as a pair mean B = (G_m + G_am)/2 and
  deviation D = (G_m - G_am)/2 gives `delta_M = b_M - b_MS` as the coefficient
  on D while holding B fixed. The same applies on the paternal side. B is a
  pair mean, not an observed grandparental lineage mean.

The contrast is a conditional association: the child's, the other parent's
and any included sibling's liabilities remain in the outcome model. A causal
interpretation as the parent-specific genetic effect beyond transmission
requires all of the following:

1. The sibling has no causal effect on the child (for example, through
   co-rearing). Otherwise the contrast also subtracts that sibling effect.
2. Lineage confounding is exchangeable between the parent and sibling. A
   confounder linear in the grandparents' midparent liability satisfies this
   under the inheritance and total-liability assortment model below.
3. There is no remaining parent-specific genotype-outcome confounding or
   selection bias. Conditioning on the child and spouse's liabilities also
   conditions out pathways through those liabilities; the contrast does not
   automatically identify a total effect of partner choice or spouse nurture.

For example, let P and S be founder-scaled total liabilities of a parent and
sibling, L the grandparents' midparent liability, V = Var(P) = Var(S), and
A = Var(L) = Cov(P,S) = V - 1/2. Under the linear Gaussian inheritance model,
`E[L | P,S] = a(P + S)`, where `a = A/(V + A)`. If
`Y = gamma P + ell L + ...`, the parent and sibling coefficients are
`gamma + ell a` and `ell a`; their difference is gamma. This also shows why
the sibling coefficient alone is not the effect per unit of the true lineage
mean. Additional conditioned family liabilities preserve this argument under
the specified assortment mechanism, which pairs parents through their total
liabilities without selecting siblings separately.

The visible contrast uses unit-variance visible liabilities. `corrected_std`
uses a unit of parental total-liability SD, assuming equal parent/sibling
marginal variances in their shared generation. Neither is rescaled to one SD
of D. `corrected_raw` uses founder units. With hidden variation, the joint
correction must precede taking the contrast; the visible difference is not
in general a simple attenuation of the total-liability difference.

### Sibling covariances under assortment

A parental sibling belongs to the parents' generation. Let S be the parents'
raw covariance of Z = (G*, H*) on the founder scale, w = (sqrt(1 - pi),
sqrt(pi)), V_T = w'Sw, and K = rho_T S w w' S / V_T the spouse block (A3). With
segregation variance I/2:

```
parent and own sibling:          S - I/2           (variance of the grandparental midparent)
parent and spouse's sibling:     X = rho_T S w w'(S - I/2) / V_T
maternal and paternal siblings:  (S - I/2)' S^{-1} X
child and a sibling:             half the sum of the two parents' blocks
```

This is exact for any number of generations or an estimated state; no
equilibrium approximation is needed. On the visible scale it gives
kappa = 1 - 1/(2 V_G), chi = k rho_G and zeta = k^2 rho_G with
k = 1 - 1/(2 V_T). Because k uses the total liability, chi differs from
rho_G kappa when pi > 0.

For the hidden correction, the standardized visible-hidden entries are d for a
parent and their own sibling (the same as within a person), k q between a
parent and the spouse's sibling, and k^2 q between the two siblings, where d and
q are the trio's within-person and spouse entries. Child rows use the fitted
visible and hidden transmission coefficients, as in the trio. Using the spouse
entry q unchanged for a parent and the spouse's sibling would overstate the
magnitude of that cross-lineage term.

The constrained modes add three pi-free restrictions: each sibling score has
the parent's reliability; kappa = (3 + rho_G - 2B)/2, where B is the visible
variance growth factor; and the child-sibling correlation is
(kappa + chi)/(2 sqrt(B)). chi (equal across lineages) and zeta are estimated
because they depend on pi. The two PGS versions retain distinct, unknown
reliabilities; the equality is between each sibling and their own parent
within a score version.

For uncertainty, the core vector includes all fitted outcome coefficients and
latent correlations (ten parameters with one sibling, fifteen with both).
The same joint transformation produces raw and standardized total coefficients
and parent-minus-sibling contrasts. Its Jacobian propagates the full covariance,
including parent/sibling covariance and external pi uncertainty when supplied.
Fixed-pi bootstrap intervals use each replicate's joint correction and contrast.
Correction packets retain the enlarged vector, covariance and joint bootstrap
draws. The implied full visible/hidden covariance must be positive semidefinite;
incompatible correction assumptions leave total-scale results unavailable while
preserving visible-scale estimates.

## Score-specific phenotype association test

This optional diagnostic tests the outcome restriction of the shared genetic
factor model. It is separate from estimating biological effects and from the
visible-to-total correction. With n included family roles and two scores per
role, write

```
P_kr = lambda_kr G_r + epsilon_kr,   C = Cov(G),
Y = b'G + u.
```

Under the null, `Cov(epsilon_kr, Y) = 0` and `Cov(G, u) = 0`, so
`c_kr = Cov(P_kr, Y) = lambda_kr (C b)_r`. For every role this implies

```
lambda_2r c_1r - lambda_1r c_2r = 0.
```

All included latent outcome paths are free under this null. In particular,
this is not a test that maternal or paternal effects vanish. The alternative
replaces every latent outcome regression with a regression on all 2n observed
scores. Both models retain the same joint score model and refit all free
nuisance parameters, including loadings and outcome residual variance. The
alternative score slopes are diagnostic quantities, not corrected genetic
effects. Latent and score outcome paths must not be combined in that model.

### Observable nesting

For an admissible null covariance, let `S = Cov(P)`, `c = Cov(P, Y)` and
`v = Var(Y)`. The alternative exactly reproduces the covariance using

```
q = solve(S, c),
Var(u_alt) = v - c' solve(S, c),
intercept_alt = E[Y] - q' E[P].
```

The residual variance is positive when the joint null covariance is positive
definite. This constructive check uses the full joint score distribution;
treating the scores as fixed regressors would change the likelihood being
compared. The implementation preserves unit latent variances, models the
scores as random, and removes every sibling outcome path and any definition
that refers to a removed effect label.

### Identification and reference degrees of freedom

When the score-only model is locally identified, C is nonsingular, loadings
are nonzero and both fits are regular interior points, the null adds n outcome
associations and one variance; the alternative adds 2n associations and one
variance. The observable dimension difference is n: three for trios, four
with one parental sibling and five with both. These values are checked against
the actual constrained moment-Jacobian ranks, rather than accepted from the
number of roles alone.

For equality constraints `h(theta) = 0`, let A be their Jacobian and N a basis
for its null space. With moment Jacobian J, the identifiable tangent dimension
is `rank(J N)`. Parameter and moment scaling, singular values and conditioning
matter in addition to rank. The score-only model and both full models must
pass these checks. Exact zero assortment without variance growth can be a
singularity even in a constrained model; near-singular fits are unsuitable
for the regular reference distribution.

Changing score units by positive factors preserves the comparison when the
measurement constraints are transformed consistently. Scores must retain the
common trait orientation required by RAVEL. Reversing an individual score's
sign also requires reversing the associated loading and residual-correlation
signs; it is not the same model under the positive-loading constraints.

The default free residual model has an additional scale invariance. Let
`Lambda_k = diag(lambda_kr)`. Its score covariance blocks are

```
S11 = Lambda_1 C Lambda_1 + Theta_1,
S12 = Lambda_1 C Lambda_2,
S22 = Lambda_2 C Lambda_2 + Theta_2.
```

For positive a near one, replace `Lambda_1` by `a Lambda_1` and `Lambda_2`
by `Lambda_2 / a`, and add `(1 - a^2) Lambda_1 C Lambda_1` to `Theta_1`
and `(1 - a^-2) Lambda_2 C Lambda_2` to `Theta_2`. This preserves the score
covariance and local strict admissibility. Unit latent variances do not remove
the invariance. Nonzero null outcome associations can identify the scale,
whereas unrestricted alternative associations do not. With exactly this
defect, an identified observable formulation would have difference n - 1;
changing the reported degrees of freedom does not identify a naive fit.
The current diagnostic therefore returns an unsupported status for this mode.

### Likelihood comparison and interpretation

For regular Gaussian ML, the unadjusted statistic is
`2 * (logLik(alternative) - logLik(null))`. RAVEL's continuous fits use MLR,
so the reported inference uses the Satorra (2000) robust nested-model
difference test with the delta option for covariance nesting. It does not
subtract separately scaled fit statistics. The software's actual adjustment
and fit diagnostics are retained; unresolved warnings, a fallback to unscaled inference,
materially worse alternative likelihood or invalid nesting prevent an ordinary
success result. See the [lavaan difference-test documentation](https://rdrr.io/cran/lavaan/man/lavTestLRT.html).

Nonlinear equality constraints require the tangent at the fitted parameter
values. The diagnostic supplies an orthonormal basis for the null space of
their final Jacobian to copies of the fitted models before calling the public
robust difference-test function. This makes the delta calculation use the
actual constrained moment spaces, including when a software version does not
populate its linear-constraint projection for nonlinear restrictions. The
stored scientific fits are unchanged. This calculation is checked against an
independent constrained-information implementation.

For reproducibility across input units, both diagnostic fits first express
each observed score and the outcome in units of its sample SD. Intercepts and
population marginal variances remain free; fixed residual covariances are
converted by the corresponding pair of SDs. This is an invertible change of
units, not a reliability assumption. It also fixes the metric used by the
finite-sample delta projection between the two fitted moment spaces. Positive
rescaling of an input score therefore gives the same diagnostic problem.

Both models use the same data, mean structure and external fitting constraints.
When a specified finite history makes the measurement model depend on pi,
the test is conditional on that supplied pi. Adding pi sampling error to
corrected-effect standard errors does not make the diagnostic unconditional.
A change confined to post-fit correction does not change its likelihoods.

Rejection shows incompatibility of the score-phenotype associations with the
shared-factor outcome restriction under the retained measurement assumptions.
Different direct/indirect mixtures, measurement misspecification and other
covariance misspecification can all produce it; the test does not distinguish
these explanations. Non-rejection does not establish genetic correlation one
or unbiased RAVEL effects. In particular, if

```
Cov(epsilon_1r, Y) = lambda_1r a_r,
Cov(epsilon_2r, Y) = lambda_2r a_r,
```

the extra association cancels in the tested contrast. For nonsingular C it can
be absorbed by replacing b with `b + solve(C, a)`, potentially biasing the
scientific effects without rejection. The initial implementation is limited
to complete, independently observed families with a continuous phenotype and
no covariates. Imputation, incomplete observations, covariates and complex
sampling need separate validation. A future null bootstrap must simulate
under the null and refit both models in every replicate.

## Genetic x covariate interactions

`interactions` adds terms G_r x x for an observed covariate x and roles r in
`interaction_roles`:

```
Y = b_D G_c + b_M G_m + b_P G_p + g_D G_c x + g_M G_m x + g_P G_p x + c x + Gamma'W + e
```

The covariate is centred at a chosen value x_0 (default: its mean), so b
refers to x = x_0 and g is the change in b per unit of x. Internally x is also
divided by its SD for numerical conditioning, and g is reported per unit of
the original covariate.

**Measurement.** Because x is observed without error,

```
(P_kr - mean P_kr) x = lambda_kr (G_r x) + epsilon_kr x .
```

The product indicators therefore measure the latent G_r x with the same
loadings lambda_kr as G_r, so g is corrected for score unreliability in the
same way as b. Fitting the observed score products directly instead attenuates
g. In the model:

- The latent interaction has a free variance and free covariances with the
  genetic factors and covariates.
- Product-indicator residuals have free within-score covariances with one
  another and with the base indicators. The latter equal theta E[x], which is
  nonzero when x_0 is not the mean.

**Assumptions.**

- x is measured without error.
- Score measurement errors are independent of x.
- The interaction is linear in G given x.

Constrained modes apply unchanged, because the product loadings are the base
loadings.

**Growth models.** Use long-format data, with one row per child-occasion and
the trio scores repeated, plus `interactions = "age"` and `cluster = "FID"`.
This fits a population-averaged growth model: estimating equations treat rows
as independent, and cluster-robust (sandwich) standard errors account for
repeated rows and within-family correlation. Families with more occasions get
more weight. Main effects and interaction coefficients are consistent under
random intercepts and slopes, but their variances are not estimated. A squared
age term (centred, for example, at age2 = 0 for (age - 10)^2) allows curvature.

**Validation.** Gaussian latent trios were simulated: 6,000 families with 1-4
occasions at ages 4-16, random intercepts and slopes, and true
g = (0.020, 0.010, 0) for (child, mother, father). Over 200 replicates:

| | Bias in g | 95% coverage of g | Main-effect coverage |
|---|---|---|---|
| Long format | below 0.001 | 0.95, 0.91, 0.95 | 0.92-0.96 |
| Cross-sectional (one row per child) | below 0.001 | 0.97-0.98 | 0.95-0.96 |

In long format the child main-effect SE was about 10% small (Monte Carlo SD
0.028, mean SE 0.026). The naive observed-score interaction was attenuated
(0.015 for a true 0.020).

**Total-liability scale.** Suppose the phenotype acts on the role-standardized
total liabilities, Y = beta'T + gamma'(T x~) + c x~ + e, with x~ = x - x_0.
Write A = Sigma_GG, B = Cov(G, T), mu = E[x~] and v = E[x~^2]. If x is
independent of every family member's visible and hidden components, then:

```
Var(G) = A          Cov(G, G x~) = mu A        Var(G x~) = v A
Cov(G, x~) = 0      Cov(G x~, x~) = 0
Cov(G, Y) = B(beta + mu gamma)                Cov(G x~, Y) = B(mu beta + v gamma)
```

b = A^{-1} B beta and g = A^{-1} B gamma solve the normal equations for any mu.
So the interaction vector is corrected by exactly the main-effect map,
including the five-role map with parental siblings, whatever the centring.

Two conditions apply:

1. Every family member in the model needs an interaction term. A^{-1}B is not
   diagonal under assortment with hidden variation, so omitted terms load on
   the included ones and, when mu is not 0, on the main effects through
   mu gamma (A^{-1}B e_r - B_rr e_r). With pi > 0, `interaction_roles` therefore
   defaults to all family members, and an incomplete set is rejected.
2. The moderator is independent of the genetic components. This is plausible
   for age but not, for example, for socioeconomic status. The fit reports
   the latent correlation of each G with the moderator and warns when one is
   clearly nonzero.

Corrected interactions (`beta_D:age`, ..., and `delta_M:age` with siblings) use
delta-method standard errors over the joint SEM covariance of the trio core,
the sibling core and the interaction coefficients. That covariance is
cluster-robust when `cluster` is set. `pi_se` adds external pi variance, and
finite-generation constrained fits include the refit sensitivity of the
interaction coefficients.

**Correction validation.** Gaussian latent pedigrees were simulated with
pi = 0.4, two rounds at rho_T = 0.4, 8,000 families and 1-4 occasions, with the
interaction on total liability. Over 100 replicates:

| | Child | Mother | Father |
|---|---|---|---|
| True gamma | 0.020 | 0.010 | 0 |
| Corrected mean | 0.0201 | 0.0097 | -0.0001 |
| Coverage | 0.95 | 0.96 | 0.95 |
| Visible mean | 0.0148 | 0.0088 | 0.0025 |

Corrected main effects in the same fits had 0.96 coverage.

**Scope.** Interactions need a continuous outcome.

## Covariate adjustment and binary outcomes

Covariates can be included jointly via `fixed_effects`. The SEM allows genetic
factors to correlate with these covariates. Residualizing only Y before fitting
is not equivalent when covariates and genetic predictors correlate; binary Y
must remain coded 0/1 rather than residualized.

With a binary outcome and covariates, a joint probit regression on all requested
PGS indicators (including parental siblings) and the encoded covariates supplies marginal variance-one liability
moments. RAVEL fits the same SEM to those moments with positive working
weights and a full joint sandwich covariance for robust SEs. This avoids
singular first-stage information from dummy covariates. It assumes a linear,
homoscedastic conditional-probit liability model; it does not require the
covariates themselves to be Gaussian. Continuous fits and binary fits without
covariates retain their existing estimators.
For this binary-covariate path, baseline-relative fit indices such as CFI/TLI
are marked unavailable rather than interpreted using an unvalidated
independence baseline. See the report diagnostics for fit-statistic limitations.

## Simulation coefficient scales

The native simulator generates the phenotype from raw founder-scaled total
liabilities:

Y = (beta^raw)^T T* + epsilon_Y.

The corresponding role-standardized target is beta = D_T beta^raw.

`G_mother`, `G_father`, and `G_child` are legacy variable names for the raw
total liabilities T_m*, T_p*, and T_c* that enter phenotype generation. When
present, `G_mother_vis`, `G_father_vis`, and `G_child_vis` are role-standardized
visible liabilities used for oracle comparisons. These legacy returned-column
names are retained for backward compatibility.
