\name{ravel_plot_report}
\alias{ravel_plot_report}
\alias{print.ravel_report_plot}
\title{Plot PGS and Visible-Scale RAVEL Coefficients}
\usage{
ravel_plot_report(out)
\method{print}{ravel_report_plot}(x, ...)
}
\arguments{
\item{out}{Output of \code{ravel_dataframe_report}.}
\item{x}{A \code{ravel_report_plot} object.}
\item{...}{Reserved for S3 compatibility.}
}
\description{Constructs a diagnostic comparison of offspring, maternal and
paternal coefficients with the selected confidence intervals. Print to display.}
\details{PGS regressions use own-role score SDs; RAVEL uses visible latent factors
fitted with child-reference PGS input scaling. Binary models are shown in separate
rows because ordinary probit and SEM liability outcome units differ. Failed fits
and unavailable confidence intervals are labelled. No fitting occurs here.}
\value{An aggregate-only plot specification with class \code{ravel_report_plot}.
Its print method returns it invisibly after drawing with base graphics.}
