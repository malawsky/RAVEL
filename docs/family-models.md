# Adult outcomes and family regression subsets

This extension separates the people used to identify the genetic model from
the liabilities included in a reported outcome regression. Existing `ravel_fit()`
calls retain their interface and estimands.

For selected-regression results, use `fit$regression`, `as.data.frame(fit)`, or
`ravel_family_regression(fit, outcome_roles)`. A trio-assisted family fit also
retains the legacy full-SEM fields. Explicit legacy helpers such as
`print_simple()` and `ravel_correct_pi()` still describe that full source model;
their tables must not be mistaken for the selected-regression tables. Use the
family projection function or analysis companion for the selected model's pi
correction and sensitivity grid.

## Two sets of roles

`genetic_roles` lists the observed family members whose two independent PGSs
enter the joint SEM. The supported pedigree roles are `child`, `mother`,
`father`, `mat_sib` (mother's full sibling), and `pat_sib` (father's full sibling).
`outcome_roles` is a subset of those roles and defines a particular regression.
These names describe pedigree relationships, not the owner of the outcome.
For maternal height, set `col_outcome = "mom_height"`.

The primary use case retains the trio, optionally with either or both parental
siblings, in `genetic_roles`. The child's scores can then identify the parental
variance state under RAVEL's estimated-state restrictions. No child phenotype
is needed. The child does not have to appear in `outcome_roles`.

An omitted outcome predictor is not an absent genetic-model member, and neither
is the same as a member with partially missing observations. The fitted joint
model retains nuisance outcome associations with all genetic members. This
avoids falsely imposing zero effects and preserves outcome information that can
be necessary to identify the two-indicator measurement model.

## Regression estimand and pi correction

Let A be all genetic-model members and R the requested regression predictors.
Write G for their noise-corrected visible liabilities, T for their total
liabilities, and b_A for the full visible regression coefficients. RAVEL's
existing correction constructs

    M = Sigma_GG^{-1} Cov(G, T),
    beta_A = M^{-1} b_A.

The requested total-liability projection is then

    beta_R = Sigma_TT[R,R]^{-1} Sigma_TT[R,A] beta_A.

The analogous visible projection uses Sigma_GG and b_A. These operations must
be performed in this order: correcting a regression restricted to R with a
submatrix of M is generally a different operation, and slicing beta_A reports
coefficients conditional on all A instead of the requested reduced regression.

This construction inherits the assumptions behind the full hidden-component
correction. In particular, identifying the corrected outcome moments requires
the full outcome residual to be orthogonal to the relevant visible and hidden
components, or the corresponding justified lineage projection. It does not
turn an association into a causal effect.

The primary contrast for a mother and her sibling is beta_m - beta_ms.
Writing B = (T_m + T_ms)/2 and D = (T_m - T_ms)/2, this is the coefficient on D
in the equivalent mean/deviation parameterization. With only the mother's
phenotype, it is not a literal regression of phenotype differences. A causal
own-genetic interpretation additionally needs assumptions about sibling effects
and exchangeable family-of-origin influences.
The contrast need not be invariant when partner or child liabilities are removed:
with assortment and nonzero omitted associations, the full-adult and sibling-pair
regressions can have different contrasts. Each is compared with its own population
projection in the validation suite.

## Covariates and interactions

A regression retaining covariates W must use conditional genetic covariance
matrices, rather than marginal covariance matrices, when reducing its set of
genetic predictors. The visible matrix is

    Sigma_GG.W = Sigma_GG - Sigma_GW Sigma_WW^{-1} Sigma_WG.

For a hidden-component correction, the visible data do not identify Cov(T,W)
without an additional mapping assumption. The implemented mapping requires
Cov(T - M'G, W) = 0, implying Cov(T,W) = M' Cov(G,W). This assumption must be
explicitly confirmed before corrected reduced covariate regressions are used.
Visible projections can still be reported without that confirmation.

Corrected interaction projections retain RAVEL's assumptions concerning
moderator independence from the family's visible and hidden genetic components.
These assumptions should not be inferred merely from small observed
score-moderator correlations. Genetic main effects and interactions must be
transformed together, retaining their joint uncertainty.

## Child-free identification

Under the current founder scaling, exchangeability and segregation assumptions,
an own-sibling visible correlation kappa identifies the parental visible
variance as V_G = 1 / (2(1-kappa)). The total variance is
V_T = (V_G-pi)/(1-pi). When both partners are observed, their visible correlation
rho_G identifies current total-liability assortment as
rho_T = rho_G V_G/(V_G-pi).

A sibling pair alone has no information about its current partners' assortment;
its within-lineage correction does not require inventing that parameter.
Partners without siblings or a child do not generally identify an unconstrained
accumulated state. Such a fit needs an appropriate history assumption or a
supplied founder-scaled `parent_visible_variance`.
The supplied-state option is restricted to `estimated_state` and treats that
state as fixed; uncertainty in an external state estimate is not propagated.
For a sibling-only equilibrium or finite-generation analysis, any reported
assortment value is implied by the chosen history, not an observed current spouse
correlation. The state metadata marks this distinction explicitly.

These covariance identities do not by themselves certify measurement-model
identification. Two free PGS loadings per role with unrestricted same-score
residual covariances can leave a loading-ratio invariance in the score-only
model. Joint outcome moments or appropriate residual-kinship constraints can
resolve it. The fit must pass numerical local-identification and covariance
checks; random-mating and weak-outcome-signal limits warrant particular care.

## Inference and compatibility

Reduced estimates use the full fitted-parameter covariance, including the
estimated state and the moments entering the projection. External pi uncertainty
is propagated jointly. If fitting constraints depend on pi, sensitivity requires
refitting; a pi grid must not be treated as a purely post-fit transformation.

Bootstrap fits now retain full free-parameter draws in addition to the existing
core draws. This supports uncertainty in the genetic-covariate and product
moments without changing the legacy bootstrap summaries.

The analysis companion in `analysis/parental_deviation` enumerates the eight
maternal regressions when all four adults are available. Each outcome/assumption
shares one fitted genetic model and sample across its regression subsets. It
records failures explicitly and provides sample counts by role. Cross-sectional
records should contain one row per focal adult; family clustering alone does
not undo duplicate-adult weighting caused by multiple child records.
