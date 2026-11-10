Reporting scripts updated **4 November 2026** with consistent Excel-ready
exports and all implied genetic-liability pair correlations. The package stays
at **0.2.0.9003**; this reporting update requires only the revised scripts.
All 247 reporting checks passed; see
[workflow validation](parental-workflow-validation.md) and BUILD-MANIFEST.txt.

Current bundle: **RAVEL 0.2.0.9003** fixes the joint-covariate covariance block
for sparse sibling families and retains exactly shared covariates. Transfer the
updated child script and package together; no analysis settings need changing.
Three final integration simulations passed with regular inference, and all 100
first-stage covariance simulations succeeded. These are not 100 full SEM fits.
See [the current validation report](../../docs/mixed-family-validation.md#sparse-sibling-covariates-0209003),
[archive comparison](checks/archive-compare-9003.json), and BUILD-MANIFEST.txt.
The package archive's 103 R, help and test files match the checkout.
The regression suite passed 56 tests / 647 assertions; package check had no
errors or warnings and one existing CPU-runtime NOTE (tests run separately).

Earlier bundle update (2026-11-03): RAVEL 0.2.0.9002 adds genetically associated
child covariates to the optional mixed model. Sex, age and child PCs can be
included with family/individual random effects and genetic interactions without
asserting blanket genetic-covariate independence. Marginal AM and pi scales are
retained under the documented hidden-residual and score-error assumptions.
The original `ravel_fit()` implementation is unchanged.

The current package check has **0 errors, 0 warnings and 1 pre-existing
CPU-runtime NOTE**, with tests run separately. The archive's 100 R, help and test
files match the checkout byte for byte. See
[the current validation report](../../docs/mixed-family-validation.md),
[package check](checks/package-check-9002.log),
[archive comparison](checks/archive-compare-9002.json), and BUILD-MANIFEST.txt.
The report below is retained as the record of the earlier family-model release;
its test counts and publication status do not describe this update.

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
[the simulation report](validation-report.md) for assumptions,
Monte Carlo uncertainty, source hashes and reproduction commands.

## Logs and reproduction

The `checks/` files contain the package check, packaged tests, non-CRAN tests,
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


## Version 0.2.0.9004: diagnostics and PGS reporting (2026-11-05)

The diagnostic constructor now installs transformed fitted starts before the
first SEM evaluation. The test suite reproduces a constructor failure from
invalid syntax starts and checks a successful same-model refit. A separate
fault-injection test reproduces the reported `NA/NaN/Inf in foreign function
call (arg 1)` message and verifies a validated restart, preserved error history,
and failure when the restarted fit fails validation. The exact user's dataset
was unavailable; the tests distinguish that injected error from the naturally
reproduced constructor failure.

Observed-PGS comparators passed 39 seeded synthetic checks: ordinary regression
and FID CR0, continuous family and longitudinal family/individual random effects,
correlated age slopes and age/sex interactions against independent mixed-model
fits, and joint score/Hessian sandwich calculations. Additional checks cover
binary probit regression, 5% sibling families, tiny raw PGS units, variance
boundaries, separation and insufficient clusters. No new runtime dependency was
added. These validate numerical calculations, not Monte Carlo coverage.

Reporting passed 95 PGS/auxiliary-index checks and all 247 existing reporting
checks. New PGS covariance/correlation SEs matched independent derivatives;
auxiliary fit measures matched native lavaan output. Fresh raw-PGS and family
mixed workflow runs retained both vanilla regressions and successfully wrote
all seven output files. Older synthetic child fits were reused for the existing
reporting-only checks.

`R CMD check --no-tests --no-manual` completed with no errors or warnings and the
existing `unlockBinding("detectCores", ...)` NOTE. A separate comparison matched
106 R/help/test/namespace files in the built archive to the source tree. The
archive and updated child script are required; update the parental script too
when using parental outcomes.

Commands: `dev/check-pgs-comparators.R`, `dev/check-pgs-reporting.R`,
`dev/check-pgs-raw-workflow.R`, `dev/check-pgs-mixed-workflow.R`, and
`dev/check-paper-reporting.R --reuse-fits`. Synthetic caches for reporting checks
can be generated with `dev/check-paper-reporting.R --prepare-only`.

Diagnostic suites passed 420 assertions: 306 score-signal, 61 continuous
auxiliary diagnostic, 25 binary diagnostic and 28 numerical-restart checks.
The independent review found no correctness blocker in the start transformations
or same-model restart. Logs are retained in `checks/*9004.log`.


## Version 0.2.0.9005: final constraint tangent (2026-11-05)

A valid estimated-state fit with a deliberately incompatible cached lavaan
constraint basis reproduces the diagnostic rejection. Validation now uses the
final reconstructed tangent; cache compatibility remains metadata. The complete
robust comparison reproduces the baseline statistic and p-value within 1e-7.
Wrong-sized, nonfinite and rank-deficient cached bases are covered; a nonfinite
actual constraint Jacobian still fails. No fitted coefficients or original
fitted objects are changed.

The tangent-cache suite passed 21 assertions, and the existing numerical-refit
suite passed 28. R CMD check (--no-tests --no-manual) passed with the same existing
NOTE and no errors or warnings. Archive sources match the working R/help/tests.
