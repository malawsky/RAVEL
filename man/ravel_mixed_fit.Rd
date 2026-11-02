\name{ravel_mixed_fit}
\alias{ravel_mixed_fit}
\title{Fit an optional family random-effects genetic model}
\usage{
ravel_mixed_fit(dat, outcome = c("continuous", "binary_liability"),
  col_outcome = "Y", col_pgs1_child = "PGS1_child",
  col_pgs2_child = "PGS2_child", col_pgs1_mother = "PGS1_mother",
  col_pgs2_mother = "PGS2_mother", col_pgs1_father = "PGS1_father",
  col_pgs2_father = "PGS2_father", family_id = "M_ID",
  subject_id = "IID", cluster = "FID", paternal_id = NULL,
  fixed_effects = NULL, interactions = NULL, interaction_center = NULL,
  family_random = TRUE, individual_random = FALSE, random_slope = NULL,
  fit_constraints = c("estimated_state", "equilibrium", "n_gen"),
  am_n_generations = 1L, pi = 0, pi_se = NULL,
  covariate_mapping_confirmed = FALSE, diagnostic = TRUE,
  diagnostic_control = list(), control = list(), start = NULL,
  se = TRUE, keep_data = TRUE)
}
\arguments{
\item{dat}{Data frame with one row per individual or measurement occasion.}
\item{outcome}{Continuous outcomes, including repeated observations, or
cross-sectional binary outcomes coded 0 and 1. Binary repeated observations
are not supported. Binary fitting requires the optional \pkg{mvtnorm} package.}
\item{col_outcome}{Outcome column.}
\item{col_pgs1_child,col_pgs2_child,col_pgs1_mother,col_pgs2_mother,col_pgs1_father,col_pgs2_father}{The two independent PGS columns for each trio role.}
\item{family_id}{Maternal-family identifier defining the shared random intercept.}
\item{subject_id}{Unique individual identifier, constant across visits.}
\item{cluster}{Extended-family identifier for the CR0 sandwich covariance.
Multiple maternal families may share a cluster, but each maternal family must
lie within one cluster. These IDs need not equal maternal IDs.}
\item{paternal_id}{Optional paternal identifier to verify the full-sibling
assumption within maternal families. Multiple fathers are explicitly rejected;
no family is silently redefined or removed. A father shared across maternal
families must have consistent scores and belong to one FID.}
\item{fixed_effects}{Numeric outcome adjustment columns.}
\item{interactions}{Numeric moderators, such as age and binary sex, interacting
with all three genetic liabilities. Their main effects are included.}
\item{interaction_center}{Named reference values for covariates/moderators.
By default numeric 0/1 columns use zero and other columns use their sample mean.}
\item{family_random}{Include a maternal-family random intercept.}
\item{individual_random}{Include a distinct individual random intercept for
continuous repeated outcomes. Both family and individual intercepts require
some observed sibling families to distinguish their variances.}
\item{random_slope}{Optional continuous time column for an individual random
slope. Includes an intercept--slope covariance when the individual intercept is
also requested. This is separate from genetic-by-time interaction coefficients.}
\item{fit_constraints}{Genetic measurement/history restrictions. This optional
likelihood currently implements the three constrained modes; the original
\code{ravel_fit} retains all its modes and defaults.}
\item{am_n_generations}{Positive integer for the finite-generation model.}
\item{pi}{Externally supplied hidden genetic fraction, in [0,1).}
\item{pi_se}{Optional standard error for independent external pi. Full
likelihood sensitivity refits propagate it into coefficient and variance-component
uncertainty in every mode. Correlation with the fitted sample is not modelled.}
\item{covariate_mapping_confirmed}{Must be TRUE when covariates or interactions
are supplied, asserting that their conditional genetic distribution satisfies
the model's independence assumption. Small observed correlations do not establish
this assumption.}
\item{diagnostic}{Automatically attempt a separate unadjusted score-signal
diagnostic. Defaults to TRUE. Set FALSE explicitly for estimator-only simulations.}
\item{diagnostic_control}{Named controls passed to \code{ravel_score_diagnostic}:
age, target_age, occasion, occasion_value, or selection_seed. Related/repeated
data use a reproducible independent cross-sectional auxiliary sample.}
\item{control}{Named optimizer controls for \code{optim}, method L-BFGS-B.}
\item{start}{Optional named starting values for internal parameters.}
\item{se}{Calculate full joint FID-clustered uncertainty.}
\item{keep_data}{Retain prepared data and model structures for audit/refitting.}
}
\value{A \code{ravel_mixed_results} object with coefficient and variance tables,
their joint covariance matrices, conditional parameter covariance, fitted state,
fit statistics, convergence/inference statuses, and the auxiliary diagnostic.
The \code{table} has visible, corrected_std and corrected_raw scales.}
\description{Fits a separate, optional full-family likelihood without changing
the established \code{ravel_fit} fitting path. Shared parental measurements and
repeated individual scores are represented once per maternal family.}
\details{Genetic measurement error and hidden genetic components are integrated
jointly with outcomes. Complete PGSs are required. Continuous missing outcomes
can contribute score information without an observed phenotype. The genetic
block assumes full siblings with the same two parents and constrained Mendelian
measurement-error covariance; it does not infer siblinghood solely from M_ID.
The estimated-state total parental variance is constrained to exceed 0.5 in
founder units, a stronger domain restriction than the legacy visible-variance
bound when pi is positive.

Continuous coefficients are returned in the input outcome units. Internal
standardization is numerical only. Binary coefficients are in conditional probit
units with occasion residual variance fixed to one, not the original fitter's
unit marginal liability variance. Corrected_std denotes one present-day genetic
liability SD; corrected_raw denotes one founder genetic liability SD. Visible
coefficients are marginal trio projections, not regressions on estimated factor
scores or conditional regressions on every sibling's scores.

The shared family residual need not be purely environmental. The model assumes
residual random effects independent of the genetic components. Robust covariance
accommodates dependence across maternal families within FID, but does not remove
bias from a misspecified conditional genetic mean. Random-effect boundaries or
singular information are reported explicitly without ordinary Wald inference.
Likelihood/AIC are available; omnibus chi-square, CFI and RMSEA are not supplied
for this family likelihood.

The auxiliary diagnostic omits covariates and genetic interactions with warnings.
It is not a test of the fitted random effects or interaction slopes. See
\code{ravel_score_diagnostic} for sampling, scope and failure statuses.}
\seealso{\code{\link{ravel_fit}}, \code{\link{ravel_score_diagnostic}}}
