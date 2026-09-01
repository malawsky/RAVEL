#' @keywords internal
"_PACKAGE"

#' RAVEL: Recovering Assortment-adjusted Variance Explained by genetic Liabilities
#'
#' A structural equation modeling framework for estimating direct and indirect
#' genetic effects on offspring phenotypes using parent-offspring trio data
#' with polygenic scores. RAVEL accounts for measurement error, assortative
#' mating, hidden genetic variants, cross-parent tagging, and
#' population stratification.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{ravel_fit}}: Fit the RAVEL SEM model to trio data
#'   \item \code{\link{ravel_simulate_trios}}: Simulate trio data for validation
#'   \item \code{\link{print_ravel_results}}: Print formatted model results
#'   \item \code{\link{print_simulation_summary}}: Print simulation summary
#' }
#'
#' @section Key Features:
#' \itemize{
#'   \item Multiple PGS as indicators of latent genetic liability (measurement error correction)
#'   \item Cross-person residual covariances for shared GWAS noise
#'   \item Hidden heritability correction for incomplete PGS coverage
#'   \item Cross-parent tagging correction for assortment confounds
#'   \item Population structure simulation with Fst-based differentiation
#'   \item Fixed effects support for ancestry and other covariates
#'   \item Bootstrap standard errors with cluster sampling
#'   \item Binary outcome support via liability threshold model
#' }
#'
#' @name RAVEL-package
#' @aliases RAVEL
#'
#' @importFrom lavaan sem inspect parameterEstimates
#' @importFrom MASS mvrnorm
#' @importFrom stats var cor sd rnorm rbinom runif quantile coef lm cov reshape setNames
NULL
