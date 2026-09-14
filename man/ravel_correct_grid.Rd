\name{ravel_correct_grid}
\alias{ravel_correct_grid}
\title{Total-Liability Corrections Across a Sensitivity Grid}
\usage{
ravel_correct_grid(out, hidden_fraction, n_generations = 15L,
  residual_models = NULL, uncertainty = c("auto", "lavaan", "bootstrap"),
  ci_level = NULL, covariate_mapping_confirmed = FALSE, include_raw = TRUE)
}
\arguments{
\item{out}{Full output of \code{ravel_dataframe_report()}, or its saved RDS path.
The summary table alone is insufficient.}
\item{hidden_fraction}{Numeric scalar or vector of assumed hidden fractions
in [0,1), using the same definition as \code{h} in \code{ravel_correct_packet()}.}
\item{n_generations}{Numeric scalar or vector of positive integer generation
counts. Use \code{Inf} for equilibrium. Zero is not supported by the correction.}
\item{residual_models}{NULL selects all recorded residual models; otherwise
select \code{"fixed"}, \code{"free"}, or both, when present in the report.}
\item{uncertainty}{\code{auto} uses each packet's selected inference;
\code{lavaan} propagates its joint covariance; \code{bootstrap} transforms
its stored joint draws. No new bootstrap fits are run.}
\item{ci_level}{Confidence level; defaults to the report's level, or 0.95 for
older reports without that setting.}
\item{covariate_mapping_confirmed}{Explicit confirmation that the correction
mapping applies conditional on the jointly fitted covariates. Required when
the original packet flags this assumption for review.}
\item{include_raw}{Return both standardized and raw genetic-liability scales
by default. FALSE returns standardized coefficients and spouse correlation only.}
}
\description{Forms the Cartesian product of assumed hidden fractions and
assortment-generation counts, separately for each selected residual model,
and applies \code{ravel_correct_packet()} to every cell. No SEM or PGS model
is refitted. Duplicate input values are removed while preserving their order.}
\details{Uses the original six core estimates and their full joint covariance
or stored joint bootstrap draws, not just the coefficient SEs in a summary
table. Original point estimates and uncertainty calculations are identical to
the corresponding scalar packet calls. Conventional PGS regression rows are
not corrected by this function.

The standardized output contains beta_D, beta_M, beta_P and rho_corrected;
the raw output contains beta_D, beta_M, beta_P and cov_corrected. Raw refers
to genetic-liability units; the phenotype units remain those of the original
report, identified by outcome_scale. Neither h nor generation uncertainty
is integrated into the intervals. These are conditional, pointwise sensitivity
results, not simultaneous confidence bands.

Invalid arguments, incompatible correction-code versions and missing requested
bootstrap draws stop before the grid runs. Missing original fits and incompatible
parameter cells are retained as NA rows with status and error fields. Failed
bootstrap transformations remain subject to the original minimum-success rule.
A warning is issued if any corrections or intervals are unavailable; there
is no automatic fallback between residual models or uncertainty methods.}
\value{A data frame with one row per scenario, residual model, scale and
parameter. Columns include outcome_variable, regression, residual_model,
scale, outcome_scale, parameter, estimate, se, ci_low, ci_high, p_value, N,
inference, scenario_id, hidden_fraction, n_generations, am_correction, ci_level,
n_boot_success, uncertainty_status, fit_status, correction_status, error and warning.
There are eight rows per cell/model, or four when include_raw=FALSE.}
\examples{
\dontrun{
total <- ravel_correct_grid(out,
  hidden_fraction = seq(.4, .5, .02), n_generations = c(1, 3, 15, Inf),
  include_raw = FALSE, covariate_mapping_confirmed = TRUE)
total[total$parameter == "beta_D", ]
}
}
