# Diagnose saved nonregular points without running an optimizer or replacing
# any Monte Carlo outcome. Original attempt/source labels remain authoritative.
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-mixed-validation.R")
source("tests/testthat/helper-mixed-covariates.R")
directory <- file.path("dev", "mixed-covariate-validation", "longitudinal")
attempts <- utils::read.csv(file.path(directory, "attempts.csv"), stringsAsFactors = FALSE)
parameters <- utils::read.csv(file.path(directory, "parameters.csv"), stringsAsFactors = FALSE)
attempts <- attempts[attempts$inference_status != "ok" & !attempts$error, ]
out <- list()
for (i in seq_len(nrow(attempts))) {
  a <- attempts[i, ]
  d <- mixed_covariate_data(N_families = a$N_families, seed = a$seed,
    contract = mixed_covariate_contract(sex_pc_association = .5,
      baseline_pc_association = 0, baseline_between = 0, baseline_within = 0),
    visits = 4L, individual_variance = .2, slope_variance = .003,
    intercept_slope_covariance = .012,
    age_beta = c(child = .025, mother = -.008, father = 0),
    sex_beta = c(child = .12, mother = -.06, father = .02))
  p <- RAVEL:::.ravel_mixed_prepare(d, outcome = "continuous", col_outcome = "Y",
    pgs_columns = c(child1 = "PGS1_child", child2 = "PGS2_child",
      mother1 = "PGS1_mother", mother2 = "PGS2_mother",
      father1 = "PGS1_father", father2 = "PGS2_father"),
    family_id = "M_ID", subject_id = "IID", cluster = "FID",
    fixed_effects = c("PC1", "PC2", "sex", "age"), interactions = c("age", "sex"),
    interaction_center = c(age = 10, sex = 0), family_random = TRUE,
    individual_random = TRUE, random_slope = "age", fit_constraints = "estimated_state",
    am_n_generations = 1L, pi = .4, covariate_mapping_confirmed = FALSE,
    covariate_model = "joint")
  p$spec$person_mean_covariates <- "age"
  p <- RAVEL:::.ravel_mixed_covariate_prepare(p)
  p <- RAVEL:::.ravel_mixed_covariate_compress(p)
  z <- parameters[parameters$scenario == a$scenario & parameters$replicate == a$replicate, ]
  par <- p$start; par[z$parameter] <- z$estimate
  ev <- RAVEL:::.ravel_mixed_covariate_evaluator(p, p$spec)
  value <- ev(par, derivatives = TRUE)
  status <- RAVEL:::.ravel_mixed_certification(par, value, p, p$spec)
  natural <- status$projected_natural_variance_gradient
  j <- which.max(abs(natural))
  out[[i]] <- data.frame(replicate = a$replicate, seed = a$seed,
    original_status = a$inference_status, original_source_run_id = a$source_run_id,
    max_score_per_family = max(abs(status$projected_gradient)) / a$N_families,
    max_natural_score_per_family = max(abs(natural)) / a$N_families,
    largest_natural_score_component = names(natural)[j],
    current_gradient_ok = status$gradient_ok,
    current_boundary = paste(status$boundary_parameters, collapse = ";"))
}
out <- do.call(rbind, out)
utils::write.csv(out, file.path(directory, "nonregular-score-audit.csv"), row.names = FALSE)
writeLines(c("Saved-point evaluation only; no optimizer or new fitted estimate.",
  capture.output(sessionInfo()), capture.output(tools::md5sum(c(
    "tests/testthat/helper-mixed-covariates.R", "tests/testthat/helper-mixed-validation.R",
    list.files("R", pattern = "mixed", full.names = TRUE))))),
  file.path(directory, "nonregular-score-session.txt"))
print(out, row.names = FALSE)
