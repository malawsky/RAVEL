# Rebuild the concise report from completed simulation CSVs.
# Run from the repository root: Rscript dev/family-validation-report.R
root <- "dev/family-simulation-results/main"
paths <- list.files(root, pattern = "-summary\\.csv$", full.names = TRUE)
stopifnot(length(paths) > 0L)
summaries <- lapply(paths, utils::read.csv, stringsAsFactors = FALSE)
stats <- do.call(rbind, summaries)
replicate_paths <- sub("-summary\\.csv$", "-replicates.csv", paths)
draws <- do.call(rbind, lapply(replicate_paths, utils::read.csv,
                              stringsAsFactors = FALSE))
cell <- function(x) paste(x$configuration, x$scenario, x$pi_se, sep = "/")
stats$cell <- cell(stats)
draws$cell <- cell(draws)
label <- function(configuration, scenario, pi_se) {
  base <- switch(configuration, trio = "Trio plus both parental siblings",
    siblings = "Mother and her sibling", adult_one = "Partners plus maternal sibling",
    adult_both = "Partners plus both siblings")
  paste0(base, if (scenario == "changed") "; changed assortment" else "",
         if (pi_se > 0) "; uncertain pi" else "")
}
fmt <- function(x, digits = 4L) formatC(x, digits = digits, format = "f")
markdown_table <- function(header, rows) {
  c(paste0("| ", paste(header, collapse = " | "), " |"),
    paste0("| ", paste(rep("---", length(header)), collapse = " | "), " |"),
    vapply(rows, function(x) paste0("| ", paste(x, collapse = " | "), " |"), character(1)))
}
cells <- split(stats, stats$cell)
counts <- lapply(cells, function(s) {
  d <- draws[draws$cell == s$cell[1], , drop = FALSE]
  bad <- unique(d$replicate[d$status != "ok"])
  c(label(s$configuration[1], s$scenario[1], s$pi_se[1]),
    length(unique(d$replicate)), length(bad),
    paste0(fmt(min(s$coverage), 2), " to ", fmt(max(s$coverage), 2)))
})
contrast <- stats[stats$outcome_model == "sibling_pair" & stats$parameter == "delta_m", ]
contrast_rows <- lapply(seq_len(nrow(contrast)), function(i) {
  s <- contrast[i, ]
  c(label(s$configuration, s$scenario, s$pi_se),
    fmt(s$bias), fmt(s$bias_mcse), fmt(s$empirical_sd), fmt(s$mean_se),
    fmt(s$coverage, 2), paste0("[", fmt(s$coverage_ci_low, 3), ", ",
                             fmt(s$coverage_ci_high, 3), "]"))
})
comparison <- stats[stats$outcome_model == "all_adults" &
  stats$configuration %in% c("trio", "adult_both") & stats$scenario == "baseline" &
  stats$pi_se == 0 & stats$parameter %in% c("mother", "delta_m"), ]
comparison_rows <- lapply(seq_len(nrow(comparison)), function(i) {
  s <- comparison[i, ]
  c(if (s$configuration == "trio") "Child plus four adults" else "Four adults only",
    if (s$parameter == "mother") "Maternal coefficient" else "Mother minus sibling",
    paste0(fmt(s$bias), " (", fmt(s$bias_mcse), ")"),
    fmt(s$empirical_sd), fmt(s$mean_se),
    paste0(fmt(s$coverage, 2), " [", fmt(s$coverage_ci_low, 3), ", ",
           fmt(s$coverage_ci_high, 3), "]"))
})
source("tests/testthat/helper-family-validation.R")
p <- family_validation_population()
pair <- family_validation_oracle(p, c("mother", "mat_sib"))
adult <- family_validation_oracle(p, c("mother", "father", "mat_sib", "pat_sib"))
fit_count <- sum(vapply(split(draws, draws$cell), function(x)
  length(unique(x$replicate)), integer(1)))
failed_count <- sum(vapply(split(draws, draws$cell), function(x)
  length(unique(x$replicate[x$status != "ok"])), integer(1)))
lines <- c(
  "# RAVEL family regression validation",
  "",
  paste0("The independent validation checks reduced adult regressions while retaining family members as genetic-model auxiliaries. The Monte Carlo results below contain ",
         fit_count, " fits across ", length(cells), " cells, with ", failed_count,
         " fits having an estimation, correction or inference failure. These results assess the specified Gaussian pedigree model; they do not establish calibration under arbitrary assortative mating or measurement misspecification."),
  "",
  "## Exact population checks",
  "",
  "The fixture constructs parents, conditional family-of-origin means, sibling segregation and child segregation from independent Gaussian innovations. Two independent GWAS error processes have the specified pedigree residual correlations. The generator and oracle projections do not call RAVEL correction functions.",
  "",
  "Exact covariance samples match the N-denominator covariance used by the likelihood. They check algebraic recovery, not sampling coverage. The test file covers all eight maternal adult subsets, trio-only and one-sibling adapters, child-free sibling/partner configurations, zero pi, negative and changing assortment, raw and standardized coefficients, variance explained, joint contrast covariance and external pi uncertainty.",
  "",
  "## Monte Carlo design",
  "",
  "Each cell uses 100 random samples of 2,000 families and fits the estimated-state model. Baseline pi is 0.4, current total-liability spouse correlation is 0.3, and parental total-liability variance is 1.15. PGS error variances are 0.4 and 1.6. The changed-assortment cell uses spouse correlation -0.2 and parental variance 1.5. The uncertain-pi cell draws an independent calibration estimate with SD 0.03 and supplies that SE to the model.",
  "",
  "The phenotype has nonzero child, parental, sibling and family-of-origin contributions. Consequently, a reduced-regression target differs from a subset of the full-model coefficient vector. Oracle targets use the true total-liability covariance and outcome cross-covariance. The same seeds are shared between configurations for paired comparison; cells are not pooled as independent evidence. Preliminary ten-replicate smoke runs reuse the first ten seeds and are excluded from the counts and tables here.",
  "",
  markdown_table(c("Genetic model and scenario", "Replicates", "Replicates with any failure",
                   "Coverage range across reported coefficients and contrasts"), counts),
  "",
  "## Maternal sibling contrast",
  "",
  "This table compares the same outcome regression, Y on maternal and maternal-sibling total liability, across genetic-model configurations. Bias and standard errors use founder-scaled raw liability units. Coverage refers to nominal 95% intervals among valid estimates; all failure counts remain visible above.",
  "",
  markdown_table(c("Genetic model and scenario", "Bias", "Bias MCSE", "Empirical SD",
                   "Mean SE", "Coverage", "Binomial 95% CI"), contrast_rows),
  "",
  "With 100 replicates, coverage near 0.95 has Monte Carlo SE about 0.022. The displayed binomial intervals show this uncertainty, including boundary cases. These runs can reveal gross implementation or calibration problems; they are too small to distinguish modest coverage differences or establish small precision advantages.",
  "",
  "## Adding the child to the same four adult genetic model",
  "",
  "Both configurations below use the same population, sample size and random samples, and target Y on all four adult liabilities. The only additional genetic-model member in the trio-assisted configuration is the child. These compare regression-estimator precision; variance-state SEs were not recorded, so no direct state-precision improvement is claimed.",
  "",
  markdown_table(c("Genetic model", "Parameter", "Bias (MCSE)", "Empirical SD",
                   "Mean SE", "Coverage [binomial 95% CI]"), comparison_rows),
  "",
  paste0("The baseline oracle own-minus-sibling contrast is ", fmt(pair$contrasts_raw[["delta_m"]], 6),
         " in the sibling-pair regression and ", fmt(adult$contrasts_raw[["delta_m"]], 6),
         " when both partners and both siblings enter the outcome regression. That difference is expected: these are different conditional associations. It is not evidence of an incorrect pi correction, and neither contrast is automatically a causal direct effect."),
  "",
  "## Reproduction and source records",
  "",
  "Run from the repository root. Full replicate records include random seeds, fitted pi, estimates, SEs, failures and warnings. Each cell's metadata text records R source hashes and session information. Files are under `dev/family-simulation-results/main/`; summary CSVs contain every coefficient and contrast, bias MCSE, coverage MCSE and binomial intervals.",
  "",
  "The recorded hashes identify the code at each simulation's launch. Later review fixes added bootstrap-unit validation, support for fixed-zero legacy paths and reuse of pi-sensitivity refits. The estimated-state covariance formulas used in these simulation cells were unchanged; the final regression tests additionally cover those review fixes.",
  "",
  "```sh",
  "Rscript dev/run-extension-tests.R --all",
  "Rscript dev/family-simulations.R --configuration=trio --reps=100 --n=2000 --output=dev/family-simulation-results/main",
  "Rscript dev/family-simulations.R --configuration=siblings --reps=100 --n=2000 --output=dev/family-simulation-results/main",
  "Rscript dev/family-simulations.R --configuration=adult_one --reps=100 --n=2000 --output=dev/family-simulation-results/main",
  "Rscript dev/family-simulations.R --configuration=adult_both --reps=100 --n=2000 --output=dev/family-simulation-results/main",
  "Rscript dev/family-simulations.R --configuration=siblings --pi-se=.03 --reps=100 --n=2000 --output=dev/family-simulation-results/main",
  "Rscript dev/family-simulations.R --configuration=adult_both --scenario=changed --reps=100 --n=2000 --output=dev/family-simulation-results/main",
  "Rscript dev/family-validation-report.R",
  "```",
  "",
  "The default seed is 20261029. The simulation source is `dev/family-simulations.R`; independent population construction and oracle calculations are in `tests/testthat/helper-family-validation.R`. Final exact and regression checks should be interpreted alongside the package check results, rather than replacing them.")
writeLines(lines, "dev/family-validation-report.md")
cat("Wrote dev/family-validation-report.md\n")
