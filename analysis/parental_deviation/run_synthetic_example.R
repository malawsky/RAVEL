# Run from this directory after installing the bundled package (see README.md).
# All observations are simulated. No real family data are included.
if (dir.exists("library")) .libPaths(c(normalizePath("library"), .libPaths()))
library(RAVEL)
source("ravel_family_regressions.R")
dat <- readRDS("synthetic_family_example.rds")

results <- ravel_family_regressions(
  dat = dat, outcomes = "mom_height",
  genetic_roles = c("child", "mother", "father", "mat_sib", "pat_sib"),
  focal = "mother", subject_id = "FID", pi = .4,
  models = list(estimated_state = list(fit_constraints = "estimated_state")),
  output_prefix = "synthetic_maternal_height"
)
print(results$table[results$table$scale == "total_std" &
  results$table$kind %in% c("coefficient", "contrast"),
  c("regression", "parameter", "estimate", "se", "ci_low", "ci_high", "inference_status")],
  row.names = FALSE)
stopifnot(nrow(results$failures) == 0L)
