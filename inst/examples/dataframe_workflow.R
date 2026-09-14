library(RAVEL)
data("ravel_example", package = "RAVEL")
df <- ravel_example
covariates <- c("age", "sex", "site", paste0("PC", 1:6))

out <- ravel_dataframe_report(
  df, outcome = "Y", outcome_type = "continuous",
  covariates = covariates, categorical_covariates = c("sex", "site"),
  continuous_adjustment = "joint", residual_models = c("fixed", "free"),
  uncertainty = "lavaan", make_plot = TRUE
)
print(out$model_status)
print(out$summary)

# This simulated covariate mechanism satisfies the conditional mapping assumption.
total <- ravel_correct_grid(
  out, hidden_fraction = .4, n_generations = 1,
  covariate_mapping_confirmed = TRUE, include_raw = TRUE
)
print(total)

out_binary <- ravel_dataframe_report(
  df, outcome = "Y_binary", outcome_type = "binary_liability",
  covariates = covariates, categorical_covariates = c("sex", "site"),
  residual_models = c("fixed", "free"), uncertainty = "lavaan", make_plot = TRUE
)
print(out_binary$model_status)
print(out_binary$summary)

# Inspect out$fit_statistics and out$packets$fixed$settings$model_syntax.
# Add out_dir = "results", report_name = "example" to save the complete report.
# Send the resulting ravel_report.rds for later corrections, not just summary.csv.
