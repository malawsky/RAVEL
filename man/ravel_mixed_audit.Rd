\name{ravel_mixed_audit}
\alias{ravel_mixed_audit}
\title{Extract detailed mixed-family model tables}
\usage{ravel_mixed_audit(fit, n_children = 2L, design_group = 1L)}
\arguments{
\item{fit}{A result from \code{ravel_mixed_fit}.}
\item{n_children}{Number of children in an exchangeable reference full-sibling
family, from 1 to 20. This does not select observations or refit the model.}
\item{design_group}{A retained covariance-design group for implied observed
score/outcome covariances and conditional genetic tables, or NULL to omit them.
Design-dependent tables require keep_data=TRUE in the fit. Covariate moment
tables remain available without retained data. The default selects the first
design, not an average over visits or sexes.}
}
\description{Extracts fit statistics, counts, statuses, internal parameters,
natural-scale genetic state, genetic covariance/correlation tables, and their
joint uncertainty matrices without changing the coefficient tables.}
\value{A named list of tables, covariance matrices, the auxiliary diagnostic,
and notes on interpretation. With joint covariates, additional tables are
\code{covariate_moments}, \code{genetic_covariate_covariances},
\code{conditional_genetic_covariances} and \code{conditional_implied_covariances}.
\code{covariate_audit_vcov} contains their covariance matrices;
\code{covariate_features} and \code{covariate_audit_status} record feature handling
and extraction availability. The primary coefficient table is unchanged.}
\details{The genetic tables describe a reference family. Founder covariance
uses founder genetic units; correlations use present-day role standard deviations.
State, genetic covariance and internal parameter SEs are conditional on the
supplied pi, even when the fit includes external pi uncertainty. Coefficient and
residual variance tables retain any external-pi adjustment from the fit.
Unavailable or boundary uncertainty remains unavailable.

Implied covariances include scores, outcomes and their cross-covariances for one
retained design pattern. The outcome-design table identifies the children,
moderators and random-slope ages. Entries use input score/outcome units; for
binary fits, Y denotes the latent conditional-probit response rather than the
observed binary variable. Their uncertainty is conditional on pi.

For joint covariates, \code{implied_covariances} remains marginal over genetics
at the selected outcome design and is also named
\code{marginal_implied_covariances}. This does not average over the covariate
distribution. \code{conditional_implied_covariances} conditions on all child
covariates in that family and is the covariance used by the likelihood.
\code{conditional_genetic_covariances} uses the same actual family size;
\code{n_children} changes only the separate marginal reference-family tables.
Singleton-only data do not acquire an unidentified sibling-covariate block.

Covariate moments and genetic-covariate covariances use input covariate units.
The latter include visible founder units, total present-day role SD units and
total founder units, labelled in the table. Their uncertainty includes estimated
covariate moments and their covariance with fitted genetic parameters through
the stacked FID sandwich, while remaining conditional on supplied pi. With
\code{keep_data=FALSE}, covariate moments and available SEs remain extractable;
conditional design-dependent tables return an explicit unavailable status.
When regular uncertainty is unavailable, estimates remain with missing SEs.

Likelihood criteria use internally standardized measurements. Compare only models
fitted to the same observations with the same numerical scaling. When multiple
maternal families share a dependent FID, the likelihood is a working-independence
likelihood; robust SEs do not make its criteria dependence-adjusted. For the
joint covariate model, ordinary AIC is NA and
\code{AIC_conditional_moments_fixed} is descriptive only: it treats the separately
estimated covariate moments as fixed. Likelihood parameter counts and first-stage
moment counts are reported separately.}
\seealso{\code{\link{ravel_mixed_fit}}}
