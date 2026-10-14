\name{ravel_correct_packet}
\alias{ravel_correct_packet}
\title{Transform an Aggregate Visible-Scale Packet to Total Liability}
\usage{
ravel_correct_packet(packet, h, am_correction = c("n_gen", "equilibrium", "estimated_state"),
  am_n_generations = 15L, ci_level = .95, covariate_mapping_confirmed = FALSE,
  uncertainty = c("auto", "lavaan", "bootstrap"))
}
\arguments{
\item{packet}{One element of \code{out$packets}, or path to an RDS of that packet.}
\item{h}{Assumed hidden fraction, in [0,1).}
\item{am_correction}{Finite-generation, equilibrium, or estimated-state correction. If omitted,
inherits the history of a constrained packet; legacy packets retain the
finite-generation default.}
\item{am_n_generations}{Assumed number of generations for finite-generation AM.
If omitted, inherits a constrained finite-generation packet's count; otherwise
defaults to 15.}
\item{ci_level}{Confidence level.}
\item{covariate_mapping_confirmed}{Deprecated and ignored. Retained so existing calls
still run; it never changed any estimate. Corrections with covariates assume
the visible-to-total mapping applies conditional on them (see the README).}
\item{uncertainty}{Use the packet's selected mode, model-based delta inference,
or transform stored joint bootstrap draws.}
}
\description{Returns raw and standardized total-liability coefficients and
spouse covariance/correlation, conditional on the specified hidden fraction
and mating history. Verifies compatibility of the correction implementation.}
\details{Model-based SEs propagate the full core-parameter covariance, including sibling parameters when present. Bootstrap
inference transforms each stored joint core-parameter draw, retaining dependence between
coefficients and all family covariances. Inadmissible transformations count as failed
draws; at least 80 percent of requested draws and at least 10 draws must be usable.
Intervals are percentile intervals and p-values normal approximations using the
bootstrap SE. Uncertainty in the externally assumed hidden fraction or AM history
is not included. No original data are needed.

For constrained packets, correction history and generation count must match
the fitting assumptions. For finite-generation fits after generation one,
\code{h} must also match the fitted value; changing it requires refitting.
Estimated-state packets infer the variance state without a generation count and require matching correction mode. Their fitted core is independent of \code{h}, which is checked for physical compatibility during correction.
Use \code{\link{ravel_correct_pi}} to include external hidden-fraction uncertainty.}
\value{Data frame of coefficients and spouse covariance/correlation on two
scales (eight rows without parental siblings), including requested sibling
coefficients and parent-specific contrasts
with estimates, SEs, confidence limits, p-values, inference labels and usable-draw counts.}

\section{Reporting}{The output includes its outcome units, a \code{ci_level} attribute, and a \code{constraint_spec} attribute. It can be passed directly to \code{print_simple}, \code{print_ravel_results}, or \code{ravel_plot_report}. Save as RDS to preserve attributes for further inference.}
