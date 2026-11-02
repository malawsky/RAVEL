\name{ravel_mixed_audit}
\alias{ravel_mixed_audit}
\title{Extract detailed mixed-family model tables}
\usage{ravel_mixed_audit(fit, n_children = 2L, design_group = 1L)}
\arguments{
\item{fit}{A result from \code{ravel_mixed_fit}.}
\item{n_children}{Number of children in an exchangeable reference full-sibling
family, from 1 to 20. This does not select observations or refit the model.}
\item{design_group}{A retained covariance-design group for implied observed
score/outcome covariances, or NULL to omit them. Requires keep_data=TRUE in the
fit. The default selects the first design, not an average over visits or sexes.}
}
\description{Extracts fit statistics, counts, statuses, internal parameters,
natural-scale genetic state, genetic covariance/correlation tables, and their
joint uncertainty matrices without changing the coefficient tables.}
\value{A named list of tables, covariance matrices, the auxiliary diagnostic,
and notes on interpretation.}
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

Likelihood/AIC use internally standardized measurements. Compare only models
fitted to the same observations with the same numerical scaling. When multiple
maternal families share a dependent FID, the likelihood is a working-independence
likelihood; robust SEs do not turn its AIC into a dependence-adjusted criterion.}
\seealso{\code{\link{ravel_mixed_fit}}}
