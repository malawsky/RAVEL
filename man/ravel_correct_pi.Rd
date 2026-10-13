\name{ravel_correct_pi}
\alias{ravel_correct_pi}
\title{Propagate External Hidden-Fraction Uncertainty through RAVEL Correction}
\usage{
ravel_correct_pi(object, pi, pi_ci = NULL,
  pi_ci_scale = c("normal", "logit"), pi_ci_level = .95,
  am_correction = c("n_gen", "equilibrium", "estimated_state"), am_n_generations = 15L,
  ci_level = .95, uncertainty = c("auto", "lavaan", "bootstrap"),
  core_pi_cov = NULL, core_pi_sensitivity = NULL,
  covariate_mapping_confirmed = FALSE, pi_se = NULL)
}
\arguments{
\item{object}{An admissible \code{ravel_results} fit, visible-transfer packet,
or path to an RDS containing either.}
\item{pi}{Estimated hidden genetic fraction on the founder/random-mating
variance scale, in \eqn{[0,1)}. This has the same definition as the existing
\code{h} argument, rather than a present-generation hidden variance fraction.}
\item{pi_se}{Recommended uncertainty input: a finite nonnegative standard
error for \code{pi}, on the raw probability scale. Supply this argument by name.
Zero means no external sampling uncertainty. Supply exactly one of
\code{pi_se} or \code{pi_ci}.}
\item{pi_ci}{Alternative uncertainty input: length-two lower and upper confidence bounds for \code{pi}.
Their distributional interpretation is specified by \code{pi_ci_scale}.
Equal bounds at \code{pi} mean zero external uncertainty.}
\item{pi_ci_scale}{Interpret the supplied interval as a normal-Wald interval
on the raw probability scale or as a symmetric Wald interval on the logit scale.
Arbitrary percentile intervals are not supported. Do not specify this argument
when supplying \code{pi_se}.}
\item{pi_ci_level}{Confidence level of the supplied interval; default 95 percent.
Do not specify this argument when supplying \code{pi_se}.}
\item{am_correction}{Finite-generation, equilibrium, or estimated-state correction. When omitted,
inherits the live fit's correction regime or a constrained packet's history;
otherwise defaults to finite-generation correction. Must match the fitted
history when using a constrained SEM. \code{"estimated_state"} requires a
matching estimated-state fit or transfer packet and infers the accumulated
state without supplying a generation count.}
\item{am_n_generations}{Positive integer generation count for finite-generation
correction; when omitted, inherits the live fit's or constrained packet's
generation count, falling back to 15. Must match a finite-generation constrained fit.}
\item{ci_level}{Confidence level for the resulting corrected-effect intervals.}
\item{uncertainty}{Use the object's selected mode, model-based joint core
covariance, or covariance of stored joint core bootstrap draws. All output
intervals from this function use the joint delta method and are Wald intervals.}
\item{core_pi_cov}{Optional named vector of sampling covariances for all core parameters
between the core estimates \code{bD}, \code{bm}, \code{bf}, \code{rhoMF},
\code{eta_m}, \code{eta_f}, evaluated at fixed true \code{pi}, and the external
\code{pi} estimator on the raw probability scale. Defaults to independence.
For a fit that depends on \code{pi}, this is not the covariance of its final
plug-in core estimates with \code{pi}: that marginal covariance additionally
contains \code{core_pi_sensitivity * Var(pi)}.}
\item{core_pi_sensitivity}{Optional named derivative of the entire fitted
core with respect to the \code{pi} supplied during SEM fitting. Computed by
perturbed refitting for a live finite-generation fit when needed; required for
a pi-dependent transferred packet. Applies only to pi-dependent constrained
fits, and uses the raw probability scale.}
\item{covariate_mapping_confirmed}{Explicit confirmation that the correction
mapping applies conditional on jointly fitted covariates.}
}
\description{Adds sampling uncertainty in an externally estimated hidden
genetic fraction to RAVEL's raw and standardized corrected effects and spouse
covariance/correlation. It can also recompute a correction at a new compatible
pi value using the fitted visible core.}
\details{
For a new fit, supply \code{pi} and \code{pi_se} directly to
\code{\link{ravel_fit}} to obtain corrected uncertainty in one call. This
function remains available for existing fits and sensitivity analyses:
\code{ravel_correct_pi(fit, pi = 0.4, pi_se = 0.02)}. The external sampling
variance is \code{pi_se^2}, with no confidence interval inferred or invented.
The standard error must be on the raw probability scale, rather than the logit
scale. A positive standard error at \code{pi = 0} requires a boundary-aware
analysis and is rejected; \code{pi = 0, pi_se = 0} is allowed.

The correction always starts from the fitted visible core and its conditional
covariance, including when the input fit already contains a
\code{pi_correction} table. Previous corrected estimates and external
uncertainty are not used as inputs or counted twice. The original fit object
is not modified. Supply the desired uncertainty input for the new correction;
the function does not implicitly reuse an earlier pi SE.

Parental-sibling fits and transfer packets retain the entire enlarged joint
core covariance. Maternal siblings add \code{bms}, \code{kap_m},
\code{chi_ms}, and \code{eta_ms}; paternal siblings add \code{bfs},
\code{kap_f}, \code{chi_fs}, and \code{eta_fs}; including both also adds
\code{zeta}. Optional covariance and sensitivity vectors must name every
included core parameter. Uncertainty in each parent-specific contrast includes
its covariance with the corresponding sibling coefficient.

Confidence intervals remain supported for compatibility. A confidence interval
alone does not determine a variance, so this input requires a normal-Wald or
logit-Wald interpretation, with bounds symmetric about
the point estimate on that scale (relative tolerance \eqn{10^{-4}} and absolute
tolerance \eqn{10^{-8}}). On the selected scale its standard error is the
half-width divided by the normal quantile; logit variance is converted to raw
pi variance by the first-order derivative \eqn{pi(1-pi)}. Raw normal-Wald
bounds may lie outside [0,1]; they are a local sampling approximation and are
never clipped. The point estimate must be physical. A nondegenerate interval
at \code{pi = 0} requires a boundary-aware method and is rejected.

Let \eqn{G} be the correction derivative with respect to the fitted core
parameters, \eqn{g} its partial derivative with respect to pi, and \eqn{s}
the derivative of the fitted core with respect to an input pi. The total pi
derivative is \eqn{a = g + Gs}. If \eqn{V} is the conditional core covariance,
\eqn{v} the external pi variance and \eqn{c} the supplied conditional cross
covariance, the returned covariance is
\eqn{GVG' + ava' + Gca' + ac'G'}.
The default \eqn{c=0} assumes an independent external estimate; sample overlap
requires appropriate covariance information. No uncertainty in the number of
generations or in the assumed mating history is included.

With \code{am_correction = "estimated_state"}, the accumulated assortment
state is itself inferred from the fitted parent-child and spouse correlations;
its uncertainty is included in the joint core covariance. No fixed generation
count or equilibrium assumption is imposed, and \code{fixed_am_history} in
the returned metadata is \code{FALSE}. The estimated-state measurement fit
is independent of pi; changing pi changes the correction, subject to physical
state compatibility, but does not require refitting. Its maternal and paternal
latent parent-child correlations must agree within \eqn{10^{-5}}.
The reconstructed state must be compatible with the supplied pi, including
full-liability parental variance above one half, positive definite component
covariances, and a full-liability spouse correlation strictly between -1 and 1.

For finite-generation constrained fits that depend on pi, the fitted input
must equal \code{pi}; refit first if it differs. The fit's sensitivity must be
included, rather than treating the fitted visible core as unchanged. At
equilibrium and after one generation the constrained measurement model does
not depend on pi. Covariance matrices may be singular because of exact
constraints, but must be positive semidefinite.

Bootstrap mode requires joint core-parameter draws with at least 80 percent
and at least ten complete draws. It uses their covariance in the joint delta
method; its intervals are not bootstrap percentile intervals. Legacy results
without joint draws can be used with explicit \code{uncertainty = "lavaan"}.
These are local first-order intervals, which may extend beyond parameter
boundaries and are not truncated. Transformation failures are reported.
}
\value{A data frame with corrected estimates, standard errors,
confidence limits, p-values and inference metadata. The
\code{pi_uncertainty_source} column and \code{pi_uncertainty} attribute identify
direct SE input (\code{source = "se"}) or CI input (\code{source = "ci"}).
Direct SE input retains the supplied raw-scale SE and has no CI bounds, level,
or scale in the metadata. The \code{covariance}
attribute is the full transformed covariance, including
dependence between raw and standardized scales. Attributes \code{jacobian},
\code{partial_jacobian}, \code{joint_input_covariance},
\code{variance_components}, and \code{pi_uncertainty} retain the derivatives,
joint input covariance, variance contributions and external-uncertainty assumptions.
Trio-only input produces eight rows; each parental sibling adds four rows
(sibling coefficient and parent-specific contrast on raw and standardized scales).}
\section{Estimated-state output}{For \code{am_correction = "estimated_state"},
the \code{estimated_state} attribute contains point diagnostics plus a
\code{summary} table, joint \code{covariance}, and \code{jacobian} for offspring
visible-variance growth, parental visible variance, parental total-liability
variance, and total-liability spouse correlation. Their uncertainty includes
the fitted core covariance and external pi uncertainty. Effective generations
under constant current assortment are descriptive only; an incompatible
constant-history trajectory is labeled, and no generation-count SE is reported.
Exact random mating is distinguished from a numerically unresolved duration
near random mating; neither is reported as a finite duration or forced to equilibrium.}
\seealso{\code{\link{ravel_fit}}, \code{\link{ravel_correct_packet}}}

\section{Reporting}{The output includes its outcome units, a \code{ci_level} attribute, and a \code{constraint_spec} attribute. It can be passed directly to \code{print_simple}, \code{print_ravel_results}, or \code{ravel_plot_report}. Save as RDS to preserve attributes for further inference.}
