\name{ravel_mixed_results}
\alias{as.data.frame.ravel_mixed_results}
\alias{print.ravel_mixed_results}
\title{Inspect optional family mixed-model results}
\usage{
\method{as.data.frame}{ravel_mixed_results}(x, ...)
\method{print}{ravel_mixed_results}(x, ...)
}
\arguments{
\item{x}{A result from \code{ravel_mixed_fit}.}
\item{...}{Unused additional arguments.}
}
\description{Return the coefficient table or print the fitted sample counts,
convergence status, corrected coefficients and diagnostic messages.}
\value{The coefficient data frame, or invisibly the input result for printing.}
