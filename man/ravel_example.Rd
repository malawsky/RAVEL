\name{ravel_example}
\alias{ravel_example}
\docType{data}
\title{Simulated Complete Trios for the Dataframe Workflow}
\description{
A reproducible example of raw polygenic scores, phenotypes and covariates for
\code{ravel_dataframe_report()}. All observations are simulated; no real
individual-level data are included.
}
\usage{data("ravel_example", package = "RAVEL")}
\format{
A data frame with 4,000 rows and 18 variables:
\describe{
  \item{FID}{Unique simulated trio identifier. Each row is an independent trio.}
  \item{Y}{Unstandardized continuous offspring phenotype.}
  \item{Y_binary}{Numeric 0/1 phenotype obtained by thresholding \code{Y} at
  a fixed normal-theory upper-decile cutoff. Prevalence is approximately 10 percent.}
  \item{age}{Simulated continuous covariate, with population mean 10 and SD 1.5.}
  \item{sex}{Simulated factor with levels \code{female} and \code{male}.}
  \item{site}{Simulated factor with levels \code{A}, \code{B} and \code{C}.}
  \item{PC1, PC2, PC3, PC4, PC5, PC6}{Correlated continuous covariates.
  These illustrate PC-like input columns, not principal components computed from genotypes.}
  \item{pgs1_child, pgs1_mother, pgs1_father}{Raw score version 1 for each
  trio member, with random-mating reliability 0.5.}
  \item{pgs2_child, pgs2_mother, pgs2_father}{Raw score version 2 for each
  trio member, with random-mating reliability 0.3.}
}
The \code{simulation} attribute records the seed, generating parameters,
covariate names and population expected visible and standardized total coefficients.
}
\details{
The data use the Gaussian trio-moment construction from the covariate validation
simulations, rather than explicitly simulated SNP genotypes. The hidden fraction
is 0.4 and the total-liability mate correlation is 0.4, with one generation of
assortment. Total liabilities have parental variance 1 and offspring variance 1.2.
The phenotype uses raw total-liability coefficients 0.45, 0.12 and 0 for the
offspring, mother and father, respectively, plus covariate effects and independent
normal residual noise with SD 0.9. These raw coefficients are not heritabilities.

Seven continuous covariates are generated as linear functions of the visible
liabilities plus independent normal noise. The categorical covariates are
independent of the liabilities. Together they produce ten covariate design columns.
Under this construction, the hidden-liability mapping applies conditional on
these covariates; this is not an automatic guarantee for covariates in real data.

Score errors are independent across versions, with parent-child residual
correlations 0.5 and spouse residual correlation 0. The scores have not been
separately standardized by role; simulated generational variance differences
are preserved. There are no missing values. Use the raw columns directly in
the dataframe workflow. This single dataset demonstrates software use, not
repeated-sampling calibration.

The repository script \code{data-raw/generate_ravel_example.R} regenerates the
dataset with seed 20260913. An installed end-to-end example is available at
\code{system.file("examples", "dataframe_workflow.R", package = "RAVEL")}.
}
\examples{
data("ravel_example", package = "RAVEL")
df <- ravel_example
head(df)
attr(df, "simulation")$expected_total_std
\dontrun{
source(system.file("examples", "dataframe_workflow.R", package = "RAVEL"))
}
}
\seealso{\code{\link{ravel_dataframe_report}}, \code{\link{ravel_correct_grid}}}
\keyword{datasets}
