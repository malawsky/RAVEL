\name{ravel_correct_grid}
\alias{ravel_correct_grid}
\title{Total-Liability Corrections Across a Sensitivity Grid}
\usage{
ravel_correct_grid(out, hidden_fraction, n_generations = 15L,
  residual_models = NULL, uncertainty = c("auto", "lavaan", "bootstrap"),
  ci_level = NULL, covariate_mapping_confirmed = FALSE, include_raw = TRUE,
  pi_se = NULL, core_pi_sensitivity = NULL)
}
\arguments{
\item{out}{Full output of \code{ravel_dataframe_report()}, or its saved RDS path.
The summary table alone is insufficient.}
\item{hidden_fraction}{Numeric scalar or vector of assumed hidden fractions
in [0,1), using the same definition as \code{h} in \code{ravel_correct_packet()}.}
\item{n_generations}{Numeric scalar or vector of positive integer generation
counts, or \code{Inf} for equilibrium. When omitted, inherits each constrained
packet's fitted history; unconstrained packets retain the default of 15.
Do not supply this argument for estimated-state packets: their hidden-fraction
sweep has no imposed generation count. Zero is not supported.}
\item{residual_models}{NULL selects all recorded residual models; otherwise
select \code{"fixed"}, \code{"free"}, or both, when present in the report.}
\item{uncertainty}{\code{auto} uses each packet's selected inference;
\code{lavaan} propagates its joint covariance; \code{bootstrap} transforms
its stored joint draws. With \code{pi_se}, bootstrap inference instead uses
their joint covariance in the delta method, as in \code{ravel_correct_pi()}.
No new bootstrap fits are run.}
\item{ci_level}{Confidence level; defaults to the report's level, or 0.95 for
older reports without that setting.}
\item{covariate_mapping_confirmed}{Explicit confirmation that the correction
mapping applies conditional on the jointly fitted covariates. Required when
the original packet flags this assumption for review.}
\item{include_raw}{Return both standardized and raw genetic-liability scales
by default. FALSE returns standardized coefficients and spouse correlation only.}
\item{pi_se}{Optional finite nonnegative scalar standard error of the external
hidden-fraction estimate, on the raw probability scale, applied to every grid
cell. \code{NULL} retains conditional fixed-fraction inference. Supplying it
uses \code{ravel_correct_pi()} and assumes the external estimate is independent
of the fitted core estimates. Use that scalar function with \code{core_pi_cov}
when the estimates are dependent. Zero includes no external uncertainty.}
\item{core_pi_sensitivity}{Optional derivative of the six fitted core estimates
with respect to the hidden fraction supplied during fitting. With one selected
residual model, supply a numeric vector named \code{bD}, \code{bm}, \code{bf},
\code{rhoMF}, \code{eta_m}, \code{eta_f}; alternatively supply a list of such
vectors named by selected residual model. Required with positive \code{pi_se}
for transferred finite-generation constrained fits beyond generation one.
These validated derivatives must come from perturbed refitting; the grid does
not estimate them. Only used with \code{pi_se}.}
}
\description{Forms the Cartesian product of assumed hidden fractions and
assortment-generation counts, or a hidden-fraction sweep for estimated-state
fits, separately for each selected residual model. Applies
\code{ravel_correct_packet()} to each cell, or \code{ravel_correct_pi()} when
\code{pi_se} is supplied. No SEM or PGS model is refitted. Duplicate input
values are removed while preserving their order.}
\details{Uses the original six core estimates and their full joint covariance
or stored joint bootstrap draws, not just the coefficient SEs in a summary
table. Original point estimates and uncertainty calculations are identical to
the corresponding scalar packet calls. Conventional PGS regression rows are
not corrected by this function.

The standardized output contains beta_D, beta_M, beta_P and rho_corrected;
the raw output contains beta_D, beta_M, beta_P and cov_corrected. Raw refers
to genetic-liability units; the phenotype units remain those of the original
report, identified by outcome_scale. With \code{pi_se = NULL}, hidden-fraction
uncertainty is excluded. With \code{pi_se}, that sampling variance is propagated
jointly with the core covariance; intervals use the first-order delta method
and normal-Wald limits. Generation uncertainty is never integrated. These are
pointwise sensitivity results, not simultaneous confidence bands.

Constrained packets retain their fitted assumptions. Omitted generation counts
inherit the fitted finite count or equilibrium. Estimated-state packets produce
\code{am_correction = "estimated_state"} and \code{n_generations = NA}; no
constant duration is invented. They require a successful identification check
recorded by report generation. An explicit generation argument with such a
packet is rejected. For other constrained packets, an explicit mismatched
history becomes a failed cell. Finite-generation fits beyond generation one
depend on the supplied hidden fraction, so changing it requires a new SEM fit
and produces a failed cell here. Their uncertainty at the fitted fraction can
be propagated using a supplied \code{core_pi_sensitivity}.

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
Supplying \code{pi_se} additionally returns \code{pi_se} and
\code{pi_uncertainty_source}. There are eight rows per cell/model, or four when
include_raw=FALSE. Scalar covariance and state attributes are not attached to
the combined table; use the corresponding scalar correction for those details.}
\examples{
\dontrun{
total <- ravel_correct_grid(out,
  hidden_fraction = seq(.4, .5, .02), n_generations = c(1, 3, 15, Inf),
  include_raw = FALSE, covariate_mapping_confirmed = TRUE)
total[total$parameter == "beta_D", ]

# Inherit the constrained fit's history; estimated-state fits need no count.
state_grid <- ravel_correct_grid(out_state, hidden_fraction = c(.3, .4, .5))
state_grid_with_pi_se <- ravel_correct_grid(out_state,
  hidden_fraction = c(.3, .4, .5), pi_se = .02)
}
}
