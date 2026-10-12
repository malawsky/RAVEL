# Run from the repository root: Rscript dev/run-extension-tests.R
# Add --all to run the complete package regression suite instead.
pkgload::load_all(".", quiet = TRUE)
# Match ravel_fit's existing fallback on restricted hosts where detectCores()
# returns NA, including tests that call lavaan directly before ravel_fit.
RAVEL:::.ensure_lavaan_runtime()
filter <- if ("--all" %in% commandArgs(trailingOnly = TRUE)) NULL else
  "constrained|extension-adversarial|pi-uncertainty|estimated-state|feature-"
results <- testthat::test_local(".", filter = filter, reporter = "summary",
                               stop_on_failure = FALSE)
tab <- as.data.frame(results)
cat(sprintf("\nTests: %d; passed expectations: %d; failures: %d; errors: %d; warnings: %d\n",
  nrow(tab), sum(tab$passed), sum(tab$failed), sum(tab$error), sum(tab$warning)))
if (any(tab$failed > 0 | tab$error)) quit(status = 1L)
