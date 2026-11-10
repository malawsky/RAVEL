# Run after installing the bundle, with R_LIBS_USER pointing to that library.
# Use the installed namespace: load_all() would hide packaging/lookup failures.
stopifnot(requireNamespace("RAVEL", quietly = TRUE))
options(ravel.workflow_functions_only = TRUE)
e <- new.env(parent = globalenv())
sys.source("analysis/parental_deviation/ravel_parental_deviation.R", e)
stopifnot(e$ravel_is_current())
# A same-version installation without the helper must fail the startup check.
e$exists <- function(x, ...) {
  if (identical(x, ".ravel_pgs_comparators")) return(FALSE)
  base::exists(x, ...)
}
stopifnot(!e$ravel_is_current())
rm("exists", envir = e)
source("tests/testthat/helper-mixed-validation.R")
d <- mixed_validation_data(N_families = 300L, family_variance = .5,
  sex_effect = 0, age_effect = 0, seed = 715L)
roles <- c("child", "mother", "father")
pgs <- setNames(unlist(lapply(roles, function(r) paste0("PGS", 1:2, "_", r))),
  unlist(lapply(roles, function(r) paste0(r, 1:2))))
x <- getFromNamespace(".ravel_pgs_comparators", "RAVEL")(
  d, outcome = "Y", pgs_columns = pgs, family_id = "M_ID",
  subject_id = "IID", cluster = "FID", family_random = TRUE)
stopifnot(length(x$fits) == 2L, all(is.finite(x$table$estimate)),
  all(c("PGS1", "PGS2") %in% x$table$regression))
cat("Installed bundle: startup checks and both family-random-effect PGS comparisons passed.\n")
