\name{ravel_score_diagnostic}
\alias{ravel_score_diagnostic}
\title{Automatically attempt an unadjusted shared-factor diagnostic}
\usage{
ravel_score_diagnostic(
  dat, ..., primary_fit = NULL, family_id = "M_ID", subject_id = "IID",
  cluster = "FID", covariates = NULL, interactions = NULL, age = NULL,
  target_age = NULL, occasion = NULL, occasion_value = NULL,
  selection_seed = 1L
)
}
\arguments{
\item{dat}{Data frame for the auxiliary analysis. When \code{primary_fit}
supplies fixed residual covariances, its PGS units must match these data.}
\item{...}{Named arguments to \code{\link{ravel_fit}}, including column mappings
and genetic restrictions. They override settings recovered from \code{primary_fit}.}
\item{primary_fit}{Optional \code{ravel_fit} result supplying genetic settings
and column maps. It is not modified or required to have converged.}
\item{family_id}{Maternal-family ID column, used to check nesting when present.}
\item{subject_id}{Individual ID column. Required when \code{cluster} repeats.}
\item{cluster}{Extended-family ID column. One individual is selected per
cluster. Use \code{NULL} only when supplied individuals are independent.}
\item{covariates}{Primary adjustment variables, omitted from this diagnostic.}
\item{interactions}{Primary genetic moderators, omitted from this diagnostic.}
\item{age}{Optional numeric age column. Repeated records are selected at the
age nearest \code{target_age}, or at the earliest age if no target is supplied.}
\item{target_age}{Optional finite scalar specifying the diagnostic age.}
\item{occasion}{Optional visit column, used for a prespecified visit or the
earliest ordered visit without age, and to break age ties.}
\item{occasion_value}{Optional exact visit value. Another individual is not
substituted when the selected individual lacks that visit.}
\item{selection_seed}{Nonnegative integer seed for uniform individual
selection within sorted clusters. The caller's random-number state is restored.}
}
\value{
A \code{ravel_score_signal_test} with additional \code{auxiliary},
\code{auxiliary_fit} and \code{selection} metadata. Unsupported or failed
analyses retain missing inference and an explicit status. Binary outcomes use
a distinct exploratory conditional-probit calculation. Its nominal p-values
are reported with a calibration warning because it has not passed null
calibration; binary outcomes are never treated as continuous outcomes.
}
\description{
Refits an auxiliary model without covariates or genetic interactions and runs
\code{\link{ravel_test_score_signal}} for continuous outcomes. Related or
repeated records are reduced to one reproducibly selected individual and one
occasion per independent cluster. This marginal diagnostic does not test the
adjusted or longitudinal model.
}
\details{
Covariates generate the warning: ``Covariates were not used to compute the
one-factor test.'' Genetic interactions receive a separate omission warning.
The primary model and its results are unchanged.

Successful binary calculations return \code{status = "ok"} for numerical
completion, \code{calibration_status = "uncalibrated"}, and the nominal
asymptotic \code{p_value}. A prominent warning reports that null simulations
rejected 40 of 499 numerical results (8.0\%) at a nominal 5\% threshold.
The warning is emitted and retained in the result, along with the fitted
models, statistic and degrees of freedom. \code{diagnostics$exploratory}
also records \code{calibrated = FALSE}. Interpret the nominal p-value
cautiously; this warning concerns the auxiliary diagnostic, not the primary
binary model.

Selection never uses outcome values. It occurs before removing incomplete
selected score/outcome records, and missing records are not replaced. The
returned selection metadata records input and selected sample sizes, row
indices, selected individuals and clusters, missingness, the random seed and
the visit-selection rule. An explicit age or visit rule is required when an
individual contributes repeated observations. Ambiguous duplicate visits are
reported rather than silently resolved.

The original diagnostic's identification restrictions remain in force. In
particular, the unconstrained score measurement model is not silently replaced
by an equilibrium or other constrained model. An auxiliary marginal check does
not test random effects or interaction slopes. Removing adjustment can change
the marginal associations being tested. External pi is treated as fixed in
this comparison.
}
\seealso{\code{\link{ravel_test_score_signal}}, \code{\link{ravel_fit}}}
