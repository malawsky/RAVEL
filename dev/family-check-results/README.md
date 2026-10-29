# Final package validation — 29 October 2026

Package: RAVEL 0.2.0.9001. Base commit: `09383e8`.
Branch: `codex/family-outcome-subsets`.

## Acceptance checks

- `R CMD check --no-manual`: no errors or warnings; one pre-existing NOTE for
  the `unlockBinding("detectCores", ...)` runtime workaround.
- Packaged tests: **4,287 passing expectations, zero failures and zero warnings**.
  Four existing tests were skipped by their CRAN-mode guards.
- The two files containing those four tests were rerun with `NOT_CRAN=true`:
  **18 tests, 137 passing expectations, zero failures/errors/warnings/skips**.
  These expectations overlap the packaged suite and should not be added to its
  count as independent tests.
- The analysis wrapper passed all eight maternal regressions, child-free FIML,
  pi-grid refits, duplicate-adult validation and duplicate-model-name validation.
- The synthetic maternal-height example ran successfully against the installed
  built package and produced finite inference for every requested regression.
- The package archive's R code, help files and tests were byte-compared against
  the checkout: 76 files matched, with no differences.
- `git diff --check` passed after normalizing trailing padding in captured text
  logs. Log values and simulation records were not altered.

The packaged test harness now initializes the same existing CPU-count fallback
as the development runner before tests invoke lavaan directly. This fixes a
host-specific startup failure without changing the statistical models.

## Independent review and simulations

The independent reviewer verified fixes for dummy bootstrap IDs and unlabelled
fixed-zero legacy outcome paths, including reordered two-moderator interactions.
The analysis wrapper also rejects duplicate assumption-model names.

The independent Monte Carlo validation completed **600 fits across six cells**,
with no fitting, correction or inference failures. Nominal 95% interval coverage
ranged from 0.91 to 0.97. Each cell contains only 100 replicates, so these runs
do not establish fine differences in calibration or precision. See
[the simulation report](../family-validation-report.md) for assumptions,
Monte Carlo uncertainty, source hashes and reproduction commands.

## Logs and reproduction

The adjacent files contain the package check, packaged tests, non-CRAN tests,
wrapper smoke check and installed-example output.

```sh
env LC_ALL=C LANG=C R CMD build . --no-build-vignettes --no-manual
env LC_ALL=C LANG=C R CMD check --no-manual RAVEL_0.2.0.9001.tar.gz
env LC_ALL=C LANG=C NOT_CRAN=true Rscript dev/run-extension-tests.R --all
env LC_ALL=C LANG=C Rscript dev/run-family-wrapper-smoke.R
```

The original Desktop checkout and analysis files were not modified. The new
source checkout, package bundle and analysis companion are available locally;
no branch has been pushed to GitHub.
