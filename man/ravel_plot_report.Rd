\name{ravel_plot_report}
\alias{ravel_plot_report}
\alias{print.ravel_report_plot}
\title{Plot PGS, Visible RAVEL, or Corrected RAVEL Coefficients}
\usage{
ravel_plot_report(out, scale = NULL)
\method{print}{ravel_report_plot}(x, ...)
}
\arguments{
\item{out}{Output of \code{ravel_dataframe_report}, or a coefficient data frame from \code{ravel_to_dataframe}, \code{ravel_correct_packet}, \code{ravel_correct_pi}, or \code{ravel_correct_grid}.}
\item{scale}{Optional single predictor-scale label. Direct coefficient tables default to \code{corrected_std} if available, otherwise \code{visible}. Report objects default to their existing visible and PGS comparison. Use \code{corrected_raw} for founder-scaled effects.}
\item{x}{A \code{ravel_report_plot} object.}
\item{...}{Reserved for S3 compatibility.}
}
\description{Constructs a diagnostic comparison of offspring, maternal and
paternal coefficients with the selected confidence intervals. Print to display.}
\details{PGS regressions use own-role score SDs; visible RAVEL uses standardized
latent factors. Corrected tables use either current-role total-liability SDs or
founder-scale raw total liabilities, selected one scale at a time. Scenario
labels distinguish correction modes and hidden-fraction/history assumptions.
Both legacy \code{beta_F} and report/correction \code{beta_P} labels are accepted.
Binary models retain SEM liability units and are separated from ordinary probit
regressions. Live continuous fits retain their supplied outcome units, whereas
reports use their prepared standardized phenotype. External-pi delta-Wald
intervals are distinguished from fixed-pi percentile bootstrap intervals.
Confidence levels are read from metadata; none is invented when that metadata
is absent. Failed fits and unavailable intervals are labelled. No fitting occurs.
Subset tables to one outcome before plotting.}
\value{An aggregate-only plot specification with class \code{ravel_report_plot}.
Its print method returns it invisibly after drawing with base graphics.}
