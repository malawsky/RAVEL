\name{ravel_correct_packet}
\alias{ravel_correct_packet}
\title{Transform an Aggregate Visible-Scale Packet to Total Liability}
\usage{
ravel_correct_packet(packet, h, am_correction = c("n_gen", "equilibrium"),
  am_n_generations = 15L, ci_level = .95, covariate_mapping_confirmed = FALSE,
  uncertainty = c("auto", "lavaan", "bootstrap"))
}
\arguments{
\item{packet}{One element of \code{out$packets}, or path to an RDS of that packet.}
\item{h}{Assumed hidden fraction, in [0,1).}
\item{am_correction}{Finite-generation or equilibrium correction.}
\item{am_n_generations}{Assumed number of generations for finite-generation AM.}
\item{ci_level}{Confidence level.}
\item{covariate_mapping_confirmed}{Explicit confirmation that the standard
mapping applies conditional on jointly fitted covariates. Required when present.}
\item{uncertainty}{Use the packet's selected mode, model-based delta inference,
or transform stored joint bootstrap draws.}
}
\description{Returns raw and standardized total-liability coefficients and
spouse covariance/correlation, conditional on the specified hidden fraction
and mating history. Verifies compatibility of the correction implementation.}
\details{Model-based SEs propagate the full six-parameter covariance. Bootstrap
inference transforms each stored six-parameter draw, retaining dependence between
coefficients and all trio covariances. Inadmissible transformations count as failed
draws; at least 80 percent of requested draws and at least 10 draws must be usable.
Intervals are percentile intervals and p-values normal approximations using the
bootstrap SE. Uncertainty in the externally assumed hidden fraction or AM history
is not included. No original data are needed.}
\value{Eight-row data frame of coefficients and spouse covariance/correlation
with estimates, SEs, confidence limits, p-values, inference labels and usable-draw counts.}
