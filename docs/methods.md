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

- Seven new latent correlations are estimated: kappa (parent with own sibling),
  eta_a (child with aunt/uncle), chi (parent with the spouse's sibling) and,
  with both siblings, zeta (between the siblings).
- Within each score version, residual covariances between each sibling and
  every other family member are free. For outcrossed families with unassorted
  GWAS noise, the expected residual correlations are 1/2 with the sibling's own
  parent, 1/4 with the child and 0 with unrelated members. The constrained
  modes fix them at these values.
- The outcome model is
  `Y = b_D G_c + b_M G_m + b_P G_p + b_MS G_am + b_FS G_ap + e`.
  Writing the parent and sibling as a lineage mean B = (G_m + G_am)/2 and
  deviation D = (G_m - G_am)/2 gives `delta_M = b_M - b_MS` as the coefficient
  on D. The same applies on the paternal side.

`delta_M` is the parent-specific coefficient. It is a causal effect of the
mother's genotype, beyond transmission to the child, under three conditions:

1. The sibling has no causal effect on the child (for example, through
   co-rearing). Otherwise the contrast estimates `beta_M - beta_MS`.
2. Lineage confounding acts through the grandparental generation, so it is
   shared equally by the parent and the sibling. A confounder linear in the
   grandparents' midparent liability satisfies this exactly under
   assortment on total liability. Parent-specific influences, such as birth
   order or the parent's own environment, are not removed.
3. The mother's genotype can affect the child through her choice of partner,
   for example through the father's non-genetic nurture. That pathway is part
   of `delta_M`.

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
visible and hidden transmission coefficients, as in the trio. An earlier
equilibrium derivation that used q for every cross-lineage entry overstates
these terms. Under assortment with hidden variation the visible-scale contrasts
are not simply attenuated, so use the corrected rows when pi > 0.

The constrained modes add three pi-free restrictions: each sibling score has
the parent's reliability; kappa = (3 + rho_G - 2B)/2, where B is the visible
variance growth factor; and the child-sibling correlation is
(kappa + chi)/(2 sqrt(B)). chi (equal across lineages) and zeta are estimated
because they depend on pi.

## Covariate adjustment and binary outcomes

Covariates can be included jointly via `fixed_effects`. The SEM allows genetic
factors to correlate with these covariates. Residualizing only Y before fitting
is not equivalent when covariates and genetic predictors correlate; binary Y
must remain coded 0/1 rather than residualized.

With a binary outcome and covariates, a joint probit regression on all six
scores and the encoded covariates supplies marginal variance-one liability
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
