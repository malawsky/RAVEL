# Run from the repository root. Uses synthetic data only.
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-family-validation.R")
source("analysis/parental_deviation/ravel_family_regressions.R")
population <- family_validation_population()
dat <- family_validation_data(population, N = 1000L, seed = 419L)
settings <- list(state = list(fit_constraints = "estimated_state"))
full <- ravel_family_regressions(dat, "Y", population$roles, pi = .4,
  models = settings, subject_id = "FID")
stopifnot(length(full$outcome_models) == 8L, nrow(full$failures) == 0L,
  all(is.finite(full$table$estimate)), all(is.finite(full$table$se)))

columns <- c("PGS1_mother", "PGS2_mother", "PGS1_mat_sib", "PGS2_mat_sib", "Y", "FID")
child_free <- dat[columns]
child_free$PGS1_mat_sib[seq_len(100L)] <- NA_real_
pair <- ravel_family_regressions(child_free, "Y", c("mother", "mat_sib"),
  pi = .4, pi_grid = c(.3, .5), models = settings, subject_id = "FID")
stopifnot(length(pair$outcome_models) == 2L, length(pair$fits) == 3L,
  nrow(pair$failures) == 0L, all(is.finite(pair$table$estimate)),
  pair$sample_counts$n_with_both_scores[pair$sample_counts$role == "mat_sib"] == 900L)

duplicate <- rbind(dat, dat[1L, ])
error <- tryCatch(ravel_family_regressions(duplicate, "Y", population$roles,
  models = settings, subject_id = "FID"), error = identity)
stopifnot(inherits(error, "error"), grepl("subject_id", conditionMessage(error)))

duplicates <- list(list(fit_constraints = "estimated_state"), list(fit_constraints = "equilibrium"))
names(duplicates) <- c("same_name", "same_name")
error <- tryCatch(ravel_family_regressions(dat, "Y", population$roles,
  models = duplicates), error = identity)
stopifnot(inherits(error, "error"), grepl("distinct", conditionMessage(error)))

cat("Family analysis wrapper: eight trio-assisted models, child-free FIML, pi-grid refits, duplicate-adult and duplicate-model guards passed.\n")
