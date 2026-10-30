# Parental-phenotype workflow validation

Validated on 29 October 2026 using bundled or synthetic observations only.
The original Desktop analysis files were not modified. No GitHub push was made.

## Regression workflow

All 13 independent check groups in `dev/run-parental-phenotype-checks.R` passed:

- Eight maternal and eight paternal regressions, always containing the focal
  parent, matched independent population calculations on raw and standardized
  total-liability scales.
- The child's nonzero source outcome association remained available for joint
  identification, while no reported regression included child liability.
- Plain PGS regressions used the selected adult predictors; coefficients and
  contrast SEs matched independently fitted regressions.
- Main pi uncertainty and fixed-pi grid uncertainty remained distinct, including
  when their pi values were identical. Finite-generation grid results matched
  a separate direct fit at the grid's pi.
- Nested genetic blocks produced 2, 4, 4 and 8 focal-parent regressions. Custom
  mappings, duplicate-adult guards, missing-trio filtering and sibling FIML
  sample counts behaved as specified.

Additional smoke checks covered binary outcomes, repeated adult measurements
with age interactions, and missing sibling scores. Unsupported one-factor
configurations retained explicit reasons rather than a numeric p-value.

The complete editable-settings example also passed: three assortment models
at the main pi and five grid values produced 18 source fits, eight maternal
regressions per source, and 1,570 primary table rows. Every reported liability
coefficient and SE was finite. The run wrote the primary CSV, Excel-friendly
CSV, seven SEM audit CSVs and the fitted-object RDS.

## One-factor diagnostic

All nine end-to-end diagnostic cases passed. The bundled example produced valid
tests under estimated-state, finite-generation and equilibrium constraints.
A synthetic maternal phenotype with the child's scores and both parental
siblings produced a valid five-degree-of-freedom comparison. A custom family ID
with a dummy unused FID also worked.

The fully free residual model, missing sibling observations, retained missing
outcomes and repeated families returned their expected explicit unsupported
statuses. The existing package diagnostic suite passed **306 expectations in
12 tests**, with no failures, errors or warnings.

The old Desktop package archive contains the same diagnostic implementation as
this checkout. The user's particular unavailable result could not be reproduced
without its status or data. No statistical diagnostic restrictions were relaxed.
The workflow now checks that the diagnostic function is installed and saves its
status, reason, warnings, refits and identification details.

## Detailed SEM exports

All 12 independent audit check groups passed. Raw parameters, covariance
parameters and their SEs, all fit indices, and the parameter sampling covariance
matched lavaan's fitted object. Model-implied observed covariance SEs matched
independent numerical derivatives of the forward covariance model, with maximum
absolute discrepancies of `6.08e-13` for a continuous outcome and `4.83e-13` for
a binary outcome.

Missing fits, nonconvergence and a converged Heywood example were labelled
explicitly. Standardization warnings were retained without affecting raw
parameter SEs. The binary outcome's fixed variance received no invented SE when
it had no matching estimated moment. The standalone CSV exports round-tripped
successfully.

## Reproduction

From the repository root:

```sh
LC_ALL=C LANG=C Rscript dev/run-parental-phenotype-checks.R
LC_ALL=C LANG=C Rscript dev/one-factor-workflow-check.R
LC_ALL=C LANG=C Rscript dev/run-sem-audit-checks.R
```

To exercise the actual editable-settings example, install the bundled package,
change to the analysis directory and run `Rscript ravel_parental_phenotypes.R`.
Its generated observations are synthetic; replace the settings for your data.

The SEM audit tables preserve extraction warnings and unavailable indices.
Source-SEM uncertainty is conditional on the supplied pi; external pi
uncertainty continues to enter the corrected regression tables separately.
