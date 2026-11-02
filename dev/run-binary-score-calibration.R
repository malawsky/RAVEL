# Reproducible null calibration for the auxiliary binary score diagnostic.
# Each replicate regenerates both genotypes and outcomes independently.
# Run from the repository root, for example:
# Rscript dev/run-binary-score-calibration.R --replicates=500 --basis=lapack
# Rscript dev/run-binary-score-calibration.R --replicates=500 --basis=legacy
# The legacy switch reproduces the discarded QR basis for failure evidence;
# it is never installed into the package or used by a public fitting function.

pkgload::load_all(".", quiet = TRUE)
RAVEL:::.ensure_lavaan_runtime()
arguments <- commandArgs(trailingOnly = TRUE)
option <- function(key, default) {
  value <- arguments[startsWith(arguments, paste0("--", key, "="))]
  if (!length(value)) default else sub(paste0("^--", key, "="), "", value[[length(value)]])
}
basis <- option("basis", "lapack")
first <- as.integer(option("first", "1"))
replicates <- as.integer(option("replicates", "500"))
n <- as.integer(option("N", "1500"))
cutoff <- as.numeric(option("cutoff", "0.6"))
output <- option("output", file.path("dev", "mixed-validation", paste0("binary-diagnostic-", basis)))
stopifnot(basis %in% c("lapack", "legacy"), is.finite(first), first >= 1L,
          is.finite(replicates), replicates >= 1L, is.finite(n), n >= 50L, is.finite(cutoff))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
files <- c("dev/run-binary-score-calibration.R", list.files("R", pattern = "\\.R$", full.names = TRUE))
utils::write.csv(data.frame(file = files, md5 = unname(tools::md5sum(files))),
                 file.path(output, "source-hashes.csv"), row.names = FALSE)
writeLines(c(paste("basis:", basis), paste("replicate IDs:", first, "through", first + replicates - 1L),
  paste("N per independent sample:", n), paste("binary threshold:", cutoff),
  "seed = 620000 + replicate ID; genetic and outcome draws are new in every replicate.",
  capture.output(sessionInfo())), file.path(output, "session.txt"))

diagnostic <- RAVEL:::.ravel_test_binary_score_signal
if (basis == "legacy") {
  source <- paste(readLines("R/ravel_binary_score_signal.R"), collapse = "\n")
  old <- "qr(t(map$A), LAPACK = TRUE)"
  if (!grepl(old, source, fixed = TRUE)) stop("The current QR implementation was not recognized.")
  source <- sub(old, "qr(t(map$A))", source, fixed = TRUE)
  local <- new.env(parent = asNamespace("RAVEL"))
  eval(parse(text = source), local)
  diagnostic <- local$.ravel_test_binary_score_signal
}
constraints <- RAVEL:::.ravel_prepare_constraints("equilibrium", .4, "equilibrium", 1L,
  TRUE, "free", RAVEL:::.normalize_fixed_crossperson("free"),
  RAVEL:::.normalize_fixed_crossperson("free"))
template <- list(outcome = "binary_liability",
  constraint_spec = list(mode = "equilibrium", h = .4,
    n_generations = NA_integer_, depends_on_pi = FALSE),
  crossperson_structure = "free",
  fixed_crossperson = constraints$fixed_cov[!is.na(constraints$fixed_cov)],
  fixed_crossperson_corr = constraints$fixed_corr[!is.na(constraints$fixed_corr)],
  parental_sibs = NULL, score_signal_input = list(dependent_families = FALSE))
template$model <- RAVEL:::.ravel_score_signal_model(template)
genetic_chol <- chol(matrix(c(1, .6, .6, .6, 1, .2, .6, .2, 1), 3L))
error_chol <- chol(matrix(c(1, .5, .5, .5, 1, 0, .5, 0, 1), 3L))
attempts <- vector("list", replicates)
for (i in seq_len(replicates)) {
  replicate <- first + i - 1L
  set.seed(620000L + replicate)
  genetic <- matrix(stats::rnorm(n * 3L), n, 3L) %*% genetic_chol
  error1 <- matrix(stats::rnorm(n * 3L), n, 3L) %*% error_chol
  error2 <- matrix(stats::rnorm(n * 3L), n, 3L) %*% error_chol
  liability <- drop(genetic %*% c(.4, .1, .05)) + stats::rnorm(n)
  dat <- data.frame(Y = as.numeric(liability > cutoff))
  roles <- c("child", "mother", "father")
  dat[paste0("PGS1_", roles)] <- genetic + sqrt(.6) * error1
  dat[paste0("PGS2_", roles)] <- genetic + sqrt(1.3) * error2
  elapsed <- system.time(answer <- suppressWarnings(diagnostic(template, dat)))[["elapsed"]]
  attempts[[i]] <- data.frame(replicate = replicate, N = n, cases = sum(dat$Y),
    cutoff = cutoff, status = answer$status, p = answer$p_value,
    statistic = answer$statistic, elapsed = elapsed, reason = answer$reason)
  utils::write.table(attempts[[i]], file.path(output, "attempts.csv"), sep = ",",
    row.names = FALSE, col.names = i == 1L, append = i != 1L)
  if (i == 1L || i %% 10L == 0L) {
    cat(replicate, answer$status, answer$p_value, "elapsed", elapsed, "\n")
    flush.console()
  }
}
attempts <- do.call(rbind, attempts)
valid <- attempts$status == "ok" & is.finite(attempts$p)
rejected <- sum(attempts$p[valid] < .05)
rate <- if (any(valid)) rejected / sum(valid) else NA_real_
interval <- if (any(valid)) stats::binom.test(rejected, sum(valid))$conf.int else rep(NA_real_, 2L)
summary <- data.frame(attempts = nrow(attempts), valid = sum(valid), failed = sum(!valid),
  rejected = rejected, rejection_among_valid = rate,
  rejection_per_attempt = rejected / nrow(attempts),
  rejection_mcse = sqrt(rate * (1 - rate) / sum(valid)),
  rejection_ci_low = interval[1L], rejection_ci_high = interval[2L],
  failure_rate = mean(!valid))
utils::write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(table(status = attempts$status)),
                 file.path(output, "status-counts.csv"), row.names = FALSE)
print(summary)
