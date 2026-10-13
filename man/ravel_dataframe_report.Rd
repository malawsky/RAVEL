\name{ravel_dataframe_report}
\alias{ravel_dataframe_report}
\title{Prepare Complete Trios and Report Visible RAVEL and PGS Regressions}
\usage{
ravel_dataframe_report(df, outcome = "Y", covariates = character(),
  categorical_covariates = character(),
  outcome_type = c("continuous", "binary_liability"),
  continuous_adjustment = c("joint", "residualize_Y"),
  pgs_columns = c(PGS1_child = "pgs1_child", PGS1_mother = "pgs1_mother",
    PGS1_father = "pgs1_father", PGS2_child = "pgs2_child",
    PGS2_mother = "pgs2_mother", PGS2_father = "pgs2_father"),
  ci_level = .95, out_dir = NULL, residual_models = c("fixed", "free"),
  uncertainty = c("auto", "lavaan", "bootstrap"), n_boot = 0L,
  seed = NULL, make_plot = TRUE, report_name = NULL,
  fit_constraints = c("none", "equilibrium", "n_gen", "estimated_state"), h = 0,
  am_n_generations = 1L,
  parental_sibs = c("none", "maternal", "paternal", "both"))
}
\arguments{
\item{df}{Data frame with one independent trio per row. Not modified.}
\item{outcome}{Numeric outcome column; binary outcomes must be coded 0/1.}
\item{covariates}{Names of jointly adjusted covariates.}
\item{categorical_covariates}{Numeric categorical columns requiring treatment coding.
Factors, characters, logicals and numeric 0/1 columns are detected automatically.}
\item{outcome_type}{Continuous outcome or binary probit liability.}
\item{continuous_adjustment}{Joint adjustment by default. Residualizing only Y
is an optional, different estimand when covariates and scores correlate.}
\item{pgs_columns}{Named mapping from the six trio model score names to input
names, plus \code{PGS1_mat_sib}/\code{PGS2_mat_sib} and/or
\code{PGS1_pat_sib}/\code{PGS2_pat_sib} when requested. Omitted mappings use
lowercase input names, for example \code{pgs1_mat_sib}.}
\item{ci_level}{Confidence level.}
\item{out_dir}{Optional new output directory, or existing parent directory when
\code{report_name} is supplied. Existing report directories are never overwritten.}
\item{residual_models}{Fit fixed residual correlations, free correlations, or both.}
\item{uncertainty}{\code{auto} selects bootstrap if \code{n_boot > 0}, otherwise
model-based SEs. \code{lavaan} selects SEM delta-method and ordinary lm/glm inference.}
\item{n_boot}{Number of whole-trio bootstrap draws, zero or at least 10.
Use at least 500 for final inference; small values are for software checks only.}
\item{seed}{Optional bootstrap seed. The previous random-number state is restored.}
\item{make_plot}{Return a diagnostic plot object and, if writing output, PDF/PNG files.}
\item{report_name}{Optional single folder name inside \code{out_dir}.}
\item{fit_constraints}{Optional strong measurement/inheritance restrictions;
see \code{\link{ravel_fit}}. Default \code{none} preserves the existing workflow.
Constrained modes require the fixed residual model; an omitted
\code{residual_models} automatically selects \code{fixed}.}
\item{h}{Founder hidden fraction assumed during constrained fitting. Ignored
when \code{fit_constraints = "none"}. Finite fits after the first generation
depend on this value. Summary rows remain on the visible scale.}
\item{am_n_generations}{Number of mating rounds for constrained \code{n_gen}
fitting. Generation one has founder parents. Ignored for equilibrium and estimated-state fits.}
\item{parental_sibs}{Optional full siblings of the parents: \code{none},
\code{maternal}, \code{paternal} or \code{both}. Each relative contributes the
same two scores. Reports include sibling coefficients and parent-specific
contrasts; ordinary PGS comparisons remain trio regressions.}
}
\description{A complete-case, covariate-adjusted workflow for two score versions
in trios. Fits visible-scale RAVEL and conventional trio PGS regressions and
retains aggregate quantities for later hidden-variant corrections.}
\details{
Rows with any missing or nonfinite trio PGS, outcome or covariate are excluded.
For continuous outcomes, missing sibling scores are retained and fitted by FIML;
no sibling values are imputed. Descriptive covariance matrices use pairwise
observations and are not the FIML sufficient statistics. Each requested sibling
needs both scores observed in at least 50 families. Binary reports exclude
families with missing requested sibling scores from all compared models and
record that exclusion.
Repeated recorded family/individual identifiers are rejected. Categorical
covariates use K-1 treatment dummies; a full-rank design is required.

RAVEL uses the child's mean and SD for all roles within each score version,
preserving parental variance ratios. The main PGS regression estimates use
each role's own score SD. Continuous Y is standardized after complete-case
selection. Binary Y is unchanged. Ordinary probit PGS coefficients and SEM
liability coefficients have different outcome scales and are plotted separately.

Both residual models are fitted by default. The fixed model uses within-version
parent-child residual correlations of 0.5 and spouse residual correlations of 0;
the free model estimates them. With parental siblings, the fixed model also
uses residual correlations of 0.5 with the sibling's own parent, 0.25 with the
child, and zero with the other parent or the other parent's sibling. The free
model estimates these sibling residual covariances. Cross-version residual
covariances remain zero.
Invalid fits and violations of requested fixed constraints are flagged, not used.

Opt-in \code{fit_constraints} adds equilibrium, finite-generation or estimated-state restrictions
while retaining distinct unknown reliabilities for the two scores. The full
specification is recorded in each packet. A later correction must use the same
fitting history and, for pi-dependent finite fits, the same \code{h}; changing
these requires refitting. Whole-trio bootstrap fits retain these settings.
Use a live \code{ravel_fit} result with \code{\link{ravel_correct_pi}} for
automatic propagation through a pi-dependent fit. An aggregate-only packet
requires a separately supplied, validated core-parameter sensitivity to pi.
Estimated-state fits estimate growth without a generation count and do not
depend on \code{h}; their packets can be corrected at a physically compatible
external pi without supplying refit sensitivity. This smooth model relaxes an
exact constant-strength integer-generation history; see \code{\link{ravel_fit}}.

Bootstrap draws resample entire families jointly across models, retaining
continuous-outcome sibling missingness.
Scaling and fitting are repeated, with the original dummy-column design retained.
Point estimates are always from the original sample. Bootstrap SEs are empirical
SDs; intervals are percentile intervals. P-values are normal approximations using
bootstrap SEs, not resampling-based null-test p-values. At least 80 percent of draws
and at least 10 usable draws are required. Failures and effective counts are
reported; there is no silent fallback to model-based inference.

All coefficient rows and the visible spouse correlation have uncertainty when
estimable. Raw SEM parameter SEs are retained. Descriptive moments, reliability,
derived residual correlations and model-fit statistics do not all have SEs.
Fixed correlations are assumptions, not estimated quantities. Model-based
inference is retained separately even when bootstrap is selected.

For binary outcomes with covariates, a conditional probit regression on all
requested scores and the encoded covariates supplies marginal variance-one liability
moments. The SEM uses positive working weights and their full joint sandwich
covariance. This assumes a linear, homoscedastic conditional-probit liability
model, not normally distributed covariates. Separation, rank deficiency and
inadmissible SEM solutions can still prevent estimation.
Baseline-relative fit indices are reported as unavailable for this path: the
default independence baseline is not validated for these moments and can
contradict mutually exclusive dummy coding. Unvalidated baseline values are
retained separately in the transfer packet. Adjusted chi-square and RMSEA are
moment-matched approximations; coefficient-SE checks do not establish their
null calibration. See the fit-measure note in \code{diagnostics}.
Successful fits are screened, but this does not prove calibration for arbitrary
covariate mechanisms. Later hidden corrections require confirmation that the
mapping is appropriate conditional on the chosen covariates.
}
\value{A list with \code{summary}, \code{summary_model_based},
\code{summary_bootstrap}, \code{scaling}, \code{fit_statistics},
\code{pgs_covariance}, \code{pgs_variances}, \code{descriptives},
\code{covariate_coding}, \code{reliability}, \code{residual_parameters},
\code{diagnostics}, \code{model_status}, \code{vanilla}, \code{packets},
\code{packet}, \code{bootstrap}, \code{settings}, \code{plot} and \code{files}.
\code{packet} is the fixed-model compatibility alias, NULL if unavailable.
Packets retain all core estimates and their joint covariance, including sibling
coefficients and latent covariances when requested; bootstrap
packets also retain joint draws including transmission covariance parameters.
\code{vanilla} retains raw, child-reference and own-role coefficient tables,
covariance matrices, model-based inference and optional bootstrap tables.
No individual-level rows or fitted model objects are exported.}
\examples{
\dontrun{
out <- ravel_dataframe_report(df, outcome = "Y", covariates = c("age", "sex"),
  n_boot = 500, seed = 123, out_dir = "results", report_name = "analysis_1")
out$summary
print(out$plot)
}
}
