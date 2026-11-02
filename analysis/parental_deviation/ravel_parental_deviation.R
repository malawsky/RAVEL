## =============================================================================
## RAVEL: maternal vs paternal (parental deviation) analysis
## =============================================================================
##
## What this script does
##   1. Fits RAVEL to trio polygenic-score data, either
##        - cross-sectional: one row per child, or
##        - longitudinal:    one row per child per measurement occasion, with a
##                           family ID and age (an average growth model).
##   2. Corrects the effects for an assumed hidden genetic fraction `pi`, with
##      its standard error `pi_se`.
##   3. Repeats that correction over a grid of assumed `pi` values, each treated
##      as exact (a sensitivity analysis: no pi uncertainty is added).
##   Steps 1-3 are run for each assortative-mating model in `models` (by
##   default: estimated state, 3 generations of assortment, and equilibrium),
##   and for each outcome in `outcome`, which gives one table per outcome.
##   With parental-sibling scores (`sib_columns`), everything is also run with
##   the mother's sibling, the father's sibling, and both siblings added; the
##   trio-only analysis is always included for reference.
##
##   It returns every coefficient, with SE, 95% CI and p-value, for
##     - plain PGS regressions (PGS1 and PGS2 separately, for comparison),
##     - RAVEL visible effects (what the PGSs measure, corrected for PGS noise),
##     - RAVEL effects corrected for pi (main value and each grid value):
##   child (beta_D), maternal (beta_M), paternal (beta_F), parental midpoint
##   (beta_M + beta_F), the parental deviation test beta_M - beta_F
##   (H0: beta_M = beta_F), the spouse correlation, and in the longitudinal
##   design the age slopes (beta_D:age, ...). Midpoint and difference tests use
##   the joint covariance of the two coefficients, not their separate SEs.
##
##   With parental siblings it adds the sibling coefficients (beta_MS, beta_FS)
##   and the parent-minus-sibling contrasts delta_M = beta_M - beta_MS and
##   delta_F = beta_F - beta_FS. Families without a sibling stay in the RAVEL
##   fits, which use all available scores (full-information ML).
##
##   For each model it attempts a separate unadjusted one-factor diagnostic.
##   Covariates and genetic interactions are omitted with a warning. One
##   individual per original FID and one prespecified/nearest-age occasion are
##   selected before missing scores/outcomes are removed; IDs and sample counts
##   are saved. This does not validate adjusted effects or longitudinal slopes.
##   Continuous outcomes use the score-signal likelihood test; binary outcomes
##   use a separate conditional-probit Wald test. Fully free score models can
##   remain unidentified. Every attempt retains its status and reason.
##
## What you need
##   Two polygenic scores (PGS1, PGS2) built from independent GWASs of the same
##   trait, computed for child, mother and father, plus the outcome. For the
##   longitudinal design also a family ID and age at each measurement.
##
## Install (no internet needed)
##   Keep this script and RAVEL_0.2.0.9001.tar.gz in the same folder. The first
##   run installs RAVEL from that file. RAVEL needs the R package 'lavaan'
##   (>= 0.6-12), which most psychology R setups already have; check with
##   packageVersion("lavaan").
##
## Run
##   1. Set R's working directory to this folder, e.g. setwd("C:/path/to/folder")
##      (RStudio: Session > Set Working Directory > To Source File Location).
##   2. Edit the SETTINGS block, then run source("ravel_parental_deviation.R"),
##      or from a terminal: Rscript ravel_parental_deviation.R
##   With data_file = NULL it runs on RAVEL's bundled example data.
## =============================================================================

## ---- Load RAVEL, installing the bundled version if needed ------------------
# Older RAVEL releases (e.g. 0.2.0) lack pi, fit_constraints, interactions and
# parental contrasts, so the bundled version replaces any older install.
# When replacing the archive with a newer development build of the SAME version,
# reinstall it explicitly with install.packages("RAVEL_<version>.tar.gz",
# repos = NULL, type = "source"), then restart R. API checks cannot detect every
# reporting change within the same development version.
ravel_is_current <- function() {
  requireNamespace("RAVEL", quietly = TRUE) &&
    all(c("ravel_family_fit", "ravel_family_regression") %in% getNamespaceExports("RAVEL")) &&
    all(c("ravel_score_diagnostic", "ravel_mixed_fit", "ravel_mixed_audit") %in% getNamespaceExports("RAVEL")) &&
    all(c("pi", "pi_se", "fit_constraints", "interactions", "cluster") %in%
        names(formals(RAVEL::ravel_fit))) &&
    "parent_contrasts" %in% names(formals(RAVEL::print_simple))
}
if (!ravel_is_current()) {
  if (!requireNamespace("lavaan", quietly = TRUE) || utils::packageVersion("lavaan") < "0.6.12")
    stop("RAVEL needs the R package 'lavaan' (version 0.6-12 or later), which is not ",
         "installed here. Please install lavaan first, then run this script again.")
  bundle <- list.files(".", pattern = "^RAVEL_.*[.]tar[.]gz$", full.names = TRUE)
  if (!length(bundle))
    stop("RAVEL is missing or too old, and no RAVEL_*.tar.gz file was found in ", getwd(),
         ". Set the working directory to the folder containing this script.")
  if ("package:RAVEL" %in% search()) try(detach("package:RAVEL", unload = TRUE), silent = TRUE)
  if ("RAVEL" %in% loadedNamespaces()) try(unloadNamespace("RAVEL"), silent = TRUE)
  if ("RAVEL" %in% loadedNamespaces())
    stop("An older RAVEL is loaded in this R session. Restart R ",
         "(RStudio: Session > Restart R), then run this script again.")
  message("Installing RAVEL from ", basename(bundle[length(bundle)]), " ...")
  utils::install.packages(bundle[length(bundle)], repos = NULL, type = "source")
  if (!ravel_is_current())
    stop("RAVEL could not be updated. Restart R and run this script again.")
}
library(RAVEL)

## ---- SETTINGS: edit these ---------------------------------------------------

if (!isTRUE(getOption("ravel.workflow_functions_only", FALSE))) {
data_file <- NULL                 # .csv, .tsv/.txt or .rds, e.g. "my_trios.rds"; NULL = bundled example
design    <- "cross_sectional"    # "cross_sectional" or "longitudinal"

outcome   <- "Y"                  # outcome column(s), e.g. c("depression", "adhd", "edu");
                                  # each outcome is analysed separately, with its own table
binary_outcome <- FALSE           # TRUE if all outcomes are 0/1, or name the binary ones,
                                  # e.g. c("adhd"); binary needs cross-sectional data
drop_missing_outcome <- TRUE      # drop rows without an outcome (see note below)
# Only complete trios are analysed: rows missing any child, mother or father
# score are always dropped. Rows with a missing outcome are dropped too, so RAVEL
# and the plain regressions use the same rows. Cross-sectional only: FALSE keeps those trios in the RAVEL fit
# (lavaan full-information ML), where their scores still inform the score model.
# Longitudinal rows without an outcome or age are always dropped (not occasions).

pgs_columns <- c(                 # your column names for the six scores
  child1  = "PGS1_child",  child2  = "PGS2_child",
  mother1 = "PGS1_mother", mother2 = "PGS2_mother",
  father1 = "PGS1_father", father2 = "PGS2_father"
)

# Parental siblings (optional): a full sibling of the mother and/or the father,
# scored with the same two PGSs. Families without a sibling stay in the fits
# (FIML). When set, every analysis is run four ways: trio only, + maternal
# sibling, + paternal sibling, and + both siblings. NULL = trio only.
sib_columns <- NULL
# sib_columns <- c(mat_sib1 = "PGS1_mat_sib", mat_sib2 = "PGS2_mat_sib",
#                  pat_sib1 = "PGS1_pat_sib", pat_sib2 = "PGS2_pat_sib")

covariates <- NULL                # numeric covariates, e.g. c("sex", "PC1", "PC2")

family_id  <- "FID"               # extended-family column for sandwich SEs; change to your ID column.
                                  # Required for longitudinal data. Cross-sectional: families with more
                                  # than one row (e.g. siblings) get SEs clustered by this column;
                                  # NULL = rows are independent trios (no family handling).
subject_id <- "IID"               # individual ID; needed to select one person per FID
                                  # and one occasion for the auxiliary one-factor test.
diagnostic_cluster <- "auto"       # "auto": original FID if present, otherwise family_id;
                                  # name another outer-family column if needed. NULL only
                                  # when the supplied individuals are independent.
# Longitudinal design only
age        <- "age"               # age at measurement
age_center <- 10                  # main effects are reported at this age

# Optional mixed model for child outcomes. All effects off + NULL interactions keeps the existing estimator.
maternal_id <- "M_ID"             # nuclear genetic family, also required for individual-only models.
random_family <- FALSE           # TRUE adds a shared maternal-family intercept.
random_individual <- FALSE       # TRUE adds an individual intercept using subject_id; longitudinal only.
random_slope <- NULL             # e.g. "age" for an individual age slope; longitudinal only.
genetic_interactions <- NULL     # explicit c("age", "sex") opts into mixed fitting even without random effects.
                                  # NULL keeps automatic longitudinal age; character(0) requests no interactions.
covariate_mapping_confirmed <- FALSE # mixed models with covariates/moderators require the stated genetic independence assumption.
# The family effect and the sandwich SEs use different IDs: maternal_id versus family_id.
# Full siblings are assumed. Advanced controls, including a paternal-ID check,
# remain available through mixed_settings; do not combine active simple options with it.
mixed_settings <- NULL

# Hidden-fraction correction
pi         <- 0.528               # assumed fraction of genetic variance the PGSs miss
pi_se      <- 0.021               # its standard error (NULL = treat pi as exact)
pi_grid    <- seq(0.40, 0.60, by = 0.05)   # sensitivity grid: 0.40, 0.45, ..., 0.60
grid_pi_se <- NULL                # grid values are treated as exact (no pi SE added)

# Assortative-mating models: each is fitted and reported separately
models <- list(
  "Estimated state" = list(fit_constraints = "estimated_state"),
  "3 generations"   = list(fit_constraints = "n_gen", am_n_generations = 3L),
  "Equilibrium"     = list(fit_constraints = "equilibrium")
)

include_raw <- FALSE              # TRUE also reports corrected effects on the founder scale
output_csv  <- "ravel_parental_deviation.csv"   # one file per outcome, named
                                  # ravel_parental_deviation_<outcome>.csv, plus an
                                  # Excel copy ravel_parental_deviation_<outcome>_excel.csv
}

## ---- ANALYSIS FUNCTION (no need to edit) ------------------------------------

.ravel_parental_mixed_settings <- function(dat, design, subject_id, maternal_id,
    random_family, random_individual, random_slope, genetic_interactions,
    covariate_mapping_confirmed, mixed_settings) {
  for (flag in list(random_family, random_individual, covariate_mapping_confirmed))
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag))
      stop("random_family, random_individual and covariate_mapping_confirmed must be TRUE or FALSE.")
  if (!is.null(random_slope) && (!is.character(random_slope) || length(random_slope) != 1L ||
      is.na(random_slope) || !nzchar(random_slope))) stop("random_slope must name one numeric time column, or be NULL.")
  if (!is.null(genetic_interactions) && (!is.character(genetic_interactions) ||
      anyNA(genetic_interactions) || any(!nzchar(genetic_interactions)) || anyDuplicated(genetic_interactions)))
    stop("genetic_interactions must name distinct numeric moderator columns, or be NULL.")
  requested <- random_family || random_individual || !is.null(random_slope) || !is.null(genetic_interactions)
  if (!is.null(mixed_settings)) {
    if (requested || !is.null(genetic_interactions) || covariate_mapping_confirmed)
      stop("Use either the simple random-effect/interaction options or advanced mixed_settings, not both; put all mixed controls in one place.")
    if (!identical(maternal_id, "M_ID") && !identical(maternal_id, mixed_settings$maternal_id))
      stop("maternal_id conflicts with advanced mixed_settings$maternal_id.")
    if (!identical(subject_id, "IID") && !identical(subject_id, mixed_settings$subject_id))
      stop("subject_id conflicts with advanced mixed_settings$subject_id.")
    return(mixed_settings)
  }
  if (!requested) return(NULL)
  design <- match.arg(design, c("cross_sectional", "longitudinal"))
  if (design != "longitudinal" && (random_individual || !is.null(random_slope)))
    stop("Individual random effects require design='longitudinal'.")
  if (!is.character(maternal_id) || length(maternal_id) != 1L || is.na(maternal_id) ||
      !nzchar(maternal_id) || !maternal_id %in% names(dat))
    stop("maternal_id must name the nuclear-family column, including when only individual random effects are requested; it must not be replaced by an extended-family FID.")
  answer <- list(maternal_id = maternal_id, subject_id = subject_id,
    family_random = random_family, individual_random = random_individual,
    random_slope = random_slope, covariate_mapping_confirmed = covariate_mapping_confirmed)
  if (!is.null(genetic_interactions)) answer$interactions <- genetic_interactions
  answer
}

ravel_parental_deviation <- function(dat, design = c("cross_sectional", "longitudinal"),
                                     outcome = "Y", pgs_columns, covariates = NULL,
                                     binary_outcome = FALSE, family_id = "FID",
                                     age = "age", age_center = NULL, pi = .4,
                                     pi_se = .05, pi_grid = NULL, grid_pi_se = pi_se,
                                     fit_constraints = "estimated_state",
                                     am_n_generations = NULL, model = fit_constraints,
                                     include_comparison = TRUE, include_raw = FALSE,
                                     drop_missing_outcome = TRUE, parental_sibs = "none",
                                     sib_columns = NULL, subject_id = "IID",
                                     diagnostic_cluster = "auto", mixed_settings = NULL,
                                     maternal_id = "M_ID", random_family = FALSE,
                                     random_individual = FALSE, random_slope = NULL,
                                     genetic_interactions = NULL, covariate_mapping_confirmed = FALSE) {
  resolved_mixed <- .ravel_parental_mixed_settings(dat, design, subject_id, maternal_id,
    random_family, random_individual, random_slope, genetic_interactions,
    covariate_mapping_confirmed, mixed_settings)
  if (!is.null(resolved_mixed)) return(ravel_parental_deviation_mixed(
    dat = dat, design = design, outcome = outcome, pgs_columns = pgs_columns,
    covariates = covariates, binary_outcome = binary_outcome, family_id = family_id,
    age = age, age_center = age_center, pi = pi, pi_se = pi_se, pi_grid = pi_grid,
    grid_pi_se = grid_pi_se, fit_constraints = fit_constraints,
    am_n_generations = am_n_generations, model = model,
    include_comparison = include_comparison, include_raw = include_raw,
    drop_missing_outcome = drop_missing_outcome, parental_sibs = parental_sibs,
    sib_columns = sib_columns, mixed_settings = resolved_mixed, diagnostic_cluster = diagnostic_cluster))
  design <- match.arg(design)
  long <- design == "longitudinal"
  # Preserve original IDs/missing records separately from legacy preprocessing.
  # The auxiliary sample is selected before exclusions, using the same score
  # and outcome units as the main analysis but none of its covariates.
  diagnostic_data <- dat
  diagnostic_covariates <- covariates
  diagnostic_age <- if (long) age else NULL
  if (identical(diagnostic_cluster, "auto")) diagnostic_cluster <-
    if ("FID" %in% names(dat)) "FID" else
      if (!is.null(family_id) && family_id %in% names(dat)) family_id else NULL
  sib_label <- c(none = "Trio only", maternal = "+ maternal sibling", paternal = "+ paternal sibling",
                 both = "+ both siblings")[[parental_sibs]]
  who <- if (parental_sibs == "none") model else paste0(model, ", ", sib_label)
  sib_roles <- switch(parental_sibs, none = character(0), maternal = "mat_sib",
                      paternal = "pat_sib", both = c("mat_sib", "pat_sib"))
  if (length(sib_roles) && is.null(sib_columns)) stop("sib_columns must be set to analyse parental siblings.")
  sib_used <- if (length(sib_roles)) unname(sib_columns[paste0(rep(sib_roles, each = 2), 1:2)]) else character(0)
  if (anyNA(sib_used)) stop("sib_columns needs mat_sib1, mat_sib2, pat_sib1 and pat_sib2 entries.")
  needed <- c(outcome, unname(pgs_columns), sib_used, covariates, if (long) c(family_id, age))
  absent <- setdiff(needed, names(dat))
  if (length(absent)) stop("Column(s) not found in the data: ", paste(absent, collapse = ", "))
  if (long && binary_outcome) stop("The longitudinal design needs a continuous outcome.")
  if (long && is.null(family_id)) stop("The longitudinal design needs family_id (the family ID column).")
  if (long && anyNA(dat[[family_id]])) stop("The family ID column has missing values.")
  if (!is.numeric(pi) || length(pi) != 1L || pi <= 0 || pi >= 1) stop("pi must be one value in (0, 1).")
  if (length(pi_grid) && any(pi_grid <= 0 | pi_grid >= 1)) stop("pi_grid values must lie in (0, 1).")
  if (anyDuplicated(c(pgs_columns, sib_used))) stop("Each PGS column must be a different column.")
  numeric_cols <- c(outcome, unname(pgs_columns), sib_used, covariates, if (long) age)
  not_numeric <- numeric_cols[!vapply(dat[numeric_cols], is.numeric, logical(1))]
  if (length(not_numeric))
    stop("Column(s) must be numeric: ", paste(not_numeric, collapse = ", "),
         ". Check for text entries (e.g. 'NA' strings or decimal commas).")

  diagnostic_target_age <- if (long) {
    if (is.null(age_center)) mean(dat[[age]], na.rm = TRUE) else age_center
  } else NULL

  # Only complete trios are analysed: drop rows missing any of the six scores.
  # Also drop rows without an outcome (and, longitudinally, without age).
  n_input <- nrow(dat)
  roles <- c("child", "mother", "father")
  role_missing <- matrix(vapply(roles, function(r) !stats::complete.cases(dat[pgs_columns[paste0(r, 1:2)]]),
    logical(nrow(dat))), nrow = nrow(dat), dimnames = list(NULL, roles))
  miss_pgs <- rowSums(role_missing) > 0
  miss_y <- is.na(dat[[outcome]]) & (drop_missing_outcome || long)
  miss_age <- if (long) is.na(dat[[age]]) else rep(FALSE, nrow(dat))
  drop <- miss_pgs | miss_y | miss_age
  dropped <- c(incomplete_trio = sum(miss_pgs), missing_outcome = sum(miss_y & !miss_pgs),
               missing_age = sum(miss_age & !miss_pgs & !miss_y))
  if (any(drop)) {
    dat <- dat[!drop, , drop = FALSE]
    message(who, ": dropped ", sum(drop), " row(s): ", sum(miss_pgs), " incomplete trio(s) (missing scores: child ",
      sum(role_missing[, "child"]), ", mother ", sum(role_missing[, "mother"]), ", father ",
      sum(role_missing[, "father"]), "), ", dropped[["missing_outcome"]], " missing outcome",
      if (long) paste0(", ", dropped[["missing_age"]], " missing age"), ". ", nrow(dat), " rows remain.")
  }
  if (nrow(dat) < 50L) stop("Fewer than 50 complete rows remain after dropping incomplete trios and missing outcomes.")
  # Sibling scores may be missing: those families stay in (FIML).
  sib_n <- vapply(sib_roles, function(r) sum(stats::complete.cases(dat[sib_columns[paste0(r, 1:2)]])), numeric(1))
  for (r in sib_roles) if (sib_n[[r]] == 0)
    stop("No family has ", if (r == "mat_sib") "maternal" else "paternal", " sibling scores.")
  sib_info <- if (length(sib_roles)) paste0(paste0(sib_n, " of ", nrow(dat), " rows have a ",
    ifelse(sib_roles == "mat_sib", "maternal", "paternal"), " sibling", collapse = "; "),
    if (any(sib_n < nrow(dat))) "; the others stay in via FIML") else ""

  # Family structure comes only from `family_id`. Longitudinal fits always
  # cluster by family; cross-sectional fits do so when a real family ID repeats
  # (e.g. siblings). RAVEL also inspects any column literally named FID, so such
  # a column is set aside whenever it is not the family ID in use.
  cluster_by <- if (long) family_id else NULL
  family_note <- ""
  n_multi <- 0L
  use_family <- !long && !is.null(family_id) && family_id %in% names(dat)
  if (use_family) {
    ids <- dat[[family_id]]
    if (anyNA(ids) || any(!nzchar(as.character(ids))))
      stop("Supplied family IDs in '", family_id,
        "' must be observed and nonempty; missing IDs cannot be treated as independent families.")
    sizes <- table(dat[[family_id]])
    if (length(sizes) < 2L)
      stop("Family clustering requires at least two distinct IDs in '", family_id,
        "'. A single family is not treated as independent trios; supply valid family IDs,",
        " or set family_id = NULL only when the observations are independent.")
    {
      n_multi <- sum(sizes > 1L)
      if (n_multi > 0L) {
        if (binary_outcome) {
          family_note <- paste(n_multi, "families have more than one row; SEs treat rows as independent",
                               "(clustering needs a continuous outcome)")
        } else {
          cluster_by <- family_id
          family_note <- paste(n_multi, "families have more than one row (e.g. siblings); SEs are clustered by family")
        }
        message(who, ": ", family_note, ".")
      }
    }
  }
  if (!long && "FID" %in% names(dat) && !(use_family && identical(family_id, "FID")))
    names(dat)[names(dat) == "FID"] <- make.unique(c(names(dat), "FID_not_used"))[ncol(dat) + 1L]

  # RAVEL's model formulas need syntactic names (e.g. no spaces or brackets) for
  # covariates and the age moderator, so analyse syntactic copies of those columns.
  syntactic <- function(v) {
    for (i in seq_along(v)) if (make.names(v[i]) != v[i]) {
      new <- make.unique(c(names(dat), make.names(v[i])))[ncol(dat) + 1L]
      dat[[new]] <<- dat[[v[i]]]
      v[i] <- new
    }
    v
  }
  covariates <- syntactic(covariates)
  if (long) age <- syntactic(age)

  # Put variables on comparable scales so the SEM converges (raw PGS sums are
  # often ~1e-5). Each PGS version is centred and scaled by one constant shared
  # by child, mother and father, and covariates are standardized; RAVEL's
  # genetic estimates do not depend on either. A continuous outcome is
  # standardized, so all effects are per SD of the outcome, as in the plain
  # regressions.
  for (k in 1:2) {
    v <- unname(pgs_columns[paste0(roles, k)])
    x <- unlist(dat[v])
    if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) stop("PGS", k, " has no variation.")
    v_all <- c(v, if (length(sib_roles)) unname(sib_columns[paste0(sib_roles, k)]))   # siblings: same constant
    pgs_center <- mean(x); pgs_scale <- stats::sd(x)
    dat[v_all] <- (dat[v_all] - pgs_center) / pgs_scale
    diagnostic_data[v_all] <- (diagnostic_data[v_all] - pgs_center) / pgs_scale
  }
  for (v in covariates) {
    sv <- stats::sd(dat[[v]], na.rm = TRUE)
    if (is.finite(sv) && sv > 0) dat[[v]] <- (dat[[v]] - mean(dat[[v]], na.rm = TRUE)) / sv
  }
  outcome_sd <- NA_real_
  if (!binary_outcome) {
    outcome_sd <- stats::sd(dat[[outcome]], na.rm = TRUE)
    if (!is.finite(outcome_sd) || outcome_sd == 0) stop("The outcome has no variation.")
    outcome_center <- mean(dat[[outcome]], na.rm = TRUE)
    dat[[outcome]] <- (dat[[outcome]] - outcome_center) / outcome_sd
    diagnostic_data[[outcome]] <- (diagnostic_data[[outcome]] - outcome_center) / outcome_sd
  }

  base_args <- list(dat,
    col_outcome = outcome,
    col_pgs1_child = pgs_columns[["child1"]],   col_pgs2_child = pgs_columns[["child2"]],
    col_pgs1_mother = pgs_columns[["mother1"]], col_pgs2_mother = pgs_columns[["mother2"]],
    col_pgs1_father = pgs_columns[["father1"]], col_pgs2_father = pgs_columns[["father2"]],
    outcome = if (binary_outcome) "binary_liability" else "continuous",
    fit_constraints = fit_constraints, fixed_effects = covariates)
  if (!is.null(am_n_generations)) base_args$am_n_generations <- am_n_generations
  if (length(sib_roles)) {
    base_args$parental_sibs <- parental_sibs
    for (r in sib_roles) for (k in 1:2)
      base_args[[paste0("col_pgs", k, "_", r)]] <- sib_columns[[paste0(r, k)]]
  }
  if (long) {
    if (is.null(age_center)) age_center <- mean(dat[[age]], na.rm = TRUE)
    base_args <- c(base_args, list(interactions = age,
      interaction_center = stats::setNames(age_center, age)))
  }
  if (!is.null(cluster_by)) base_args$cluster <- cluster_by

  labels <- c(beta_D = "Child (beta_D)", beta_M = "Maternal (beta_M)", beta_F = "Paternal (beta_F)",
    beta_mid = "Parental midpoint (beta_M + beta_F)",
    delta_M_mid = "Mother - father difference (beta_M - beta_F)",
    delta_P_mid = "Father - mother difference (beta_F - beta_M)",
    beta_MS = "Maternal sibling (beta_MS)", delta_M = "Mother - maternal sibling (delta_M)",
    beta_FS = "Paternal sibling (beta_FS)", delta_F = "Father - paternal sibling (delta_F)",
    rho_visible = "Spouse correlation", rho_corrected = "Spouse correlation")
  label_for <- function(par) {
    out <- unname(labels[par])
    slope <- grepl(":", par, fixed = TRUE)   # e.g. beta_M:age
    out[slope] <- paste0(unname(labels[sub(":.*$", "", par[slope])]), " x ", sub("^.*:", "", par[slope]))
    ifelse(is.na(out), par, out)
  }

  # An auxiliary unadjusted diagnostic; the primary fit is unchanged.
  diagnostic <- list(status = "not_run", reason = "Auxiliary diagnostic has not run.")
  one_factor_row <- function(fit) {
    args <- base_args
    args[[1L]] <- diagnostic_data
    names(args)[1L] <- "dat"
    args[c("fixed_effects", "interactions", "interaction_center", "cluster")] <- NULL
    args$pi <- pi
    args["family_id"] <- list(family_id)
    args["subject_id"] <- list(subject_id)
    args["cluster"] <- list(diagnostic_cluster)
    args$covariates <- diagnostic_covariates
    args$interactions <- if (long) diagnostic_age else NULL
    args$age <- diagnostic_age
    args$target_age <- diagnostic_target_age
    captured_warnings <- character()
    z <- tryCatch(withCallingHandlers(do.call(ravel_score_diagnostic, args),
      warning = function(w) {
        captured_warnings <<- c(captured_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) list(status = "error", reason = conditionMessage(e), N = NA_integer_))
    z$workflow_warnings <- unique(captured_warnings)
    z$workflow <- list(primary_fit_converged = isTRUE(fit$converged),
      N_input = n_input, N_primary = nrow(dat), primary_excluded_rows = which(drop),
      original_cluster = diagnostic_cluster, subject_id = subject_id,
      score_and_outcome_units = "Main analysis scaling constants, applied before auxiliary selection")
    diagnostic <<- z
    valid <- identical(z$status, "ok")
    why <- if (!valid) paste0("not available (", z$status, "): ", z$reason) else ""
    sample_note <- if (!is.null(z$selection$N_used)) paste0(
      "Unadjusted independent auxiliary sample: N = ", z$selection$N_used,
      " from ", z$selection$N_input, " input rows; ", z$selection$rule) else
      "Unadjusted auxiliary diagnostic; sample selection was not completed"
    note <- paste(unique(c(why[nzchar(why)], sample_note, captured_warnings,
      if (fit_constraints == "n_gen") paste("conditional on pi =", pi, "used in fitting"),
      if (binary_outcome) "Binary conditional-probit Wald test; not the continuous likelihood-ratio test")), collapse = "; ")
    data.frame(model = model, analysis = "One-factor test (visible fit)", pi = NA_real_, pi_se = NA_real_,
      effect = if (valid) paste0("One-factor test: chi-square (", z$df, " df)") else "One-factor test",
      parameter = "one_factor_test", estimate = if (valid) z$statistic else NA_real_,
      se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = if (valid) z$p_value else NA_real_, note = note,
      diagnostic_status = z$status, diagnostic_df = if (is.null(z$df)) NA_real_ else z$df,
      diagnostic_N = if (is.null(z$N)) NA_real_ else z$N, stringsAsFactors = FALSE)
  }

  # One RAVEL fit at a given pi; returns the fit, its rows and any warnings.
  fit_at <- function(p, p_se) {
    msgs <- character(0)
    fit <- tryCatch(withCallingHandlers(do.call(ravel_fit, c(base_args, list(pi = p, pi_se = p_se))),
      warning = function(w) { msgs <<- c(msgs, gsub("\\s+", " ", conditionMessage(w))); invokeRestart("muffleWarning") }),
      error = function(e) e)
    if (inherits(fit, "error"))
      return(list(fit = NULL, rows = NULL, note = paste("Fit failed:", conditionMessage(fit))))
    if (!isTRUE(fit$converged))
      return(list(fit = fit, rows = NULL, note = paste(c(paste0("did not converge at pi = ", p),
        if (length(msgs)) paste("first warning:", msgs[1])), collapse = "; ")))
    tab <- tryCatch(as.data.frame(fit, parent_contrasts = TRUE, include_vanilla = !long, include_raw = include_raw),
      error = function(e) e)
    if (inherits(tab, "error"))
      return(list(fit = fit, rows = NULL, note = paste("Results unavailable:", conditionMessage(tab))))
    note <- c(if (!is.null(fit$pi_correction_error)) paste("correction unavailable:", fit$pi_correction_error),
              if (length(msgs)) paste(length(msgs), "warning(s); first:", msgs[1]))
    list(fit = fit, rows = tab, note = if (length(note)) paste(note, collapse = "; ") else "")
  }
  # All rows of one regression and scale; known coefficients first, then any others.
  pick <- function(tab, scale, regression = "RAVEL") {
    out <- tab[tab$regression == regression & tab$scale == scale, , drop = FALSE]
    rename <- c(beta_child = "beta_D", beta_mother = "beta_M", beta_father = "beta_F")
    hit <- out$parameter %in% names(rename)                      # plain-regression names
    out$parameter[hit] <- rename[out$parameter[hit]]
    out[order(grepl(":", out$parameter, fixed = TRUE),          # main effects, then slopes
              match(sub(":.*$", "", out$parameter), names(labels), nomatch = length(labels) + 1L)), , drop = FALSE]
  }
  fail_row <- function(analysis, p, p_se, note)
    data.frame(model = model, analysis = analysis, pi = p, pi_se = if (is.null(p_se)) NA_real_ else p_se,
      effect = NA_character_, parameter = NA_character_, estimate = NA_real_, se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, note = note,
      diagnostic_status = NA_character_, diagnostic_df = NA_real_, diagnostic_N = NA_real_,
      stringsAsFactors = FALSE)
  label_rows <- function(rows, analysis, p, p_se, note, model_label = model) {
    if (is.null(rows) || !nrow(rows)) return(NULL)
    data.frame(model = model_label, analysis = analysis, pi = p, pi_se = if (is.null(p_se)) NA_real_ else p_se,
      effect = label_for(rows$parameter), parameter = rows$parameter,
      estimate = rows$estimate, se = rows$se, ci_low = rows$ci_low, ci_high = rows$ci_high,
      p_value = rows$p_value, note = note,
      diagnostic_status = NA_character_, diagnostic_df = NA_real_, diagnostic_N = NA_real_,
      stringsAsFactors = FALSE, row.names = NULL)
  }

  message(who, ": fitting at pi = ", pi, " ...")
  main <- fit_at(pi, pi_se)
  fit_records <- list(main = list(kind = "main", pi = pi, pi_se = pi_se,
    fit = main$fit, note = main$note))
  if (is.null(main$rows)) {
    message(who, ": ", main$note)
    results <- list(fail_row("Visible (measured by the PGSs)", NA_real_, NULL, main$note),
                    one_factor_row(main$fit),
                    fail_row("Corrected for pi", pi, pi_se, main$note))
  } else {
    # In the n_gen model the measurement fit itself depends on pi.
    vis_note <- paste(c(if (fit_constraints == "n_gen") paste("visible fit uses pi =", pi),
                        if (nzchar(main$note)) main$note), collapse = "; ")
    results <- list(
      label_rows(pick(main$rows, "visible"), "Visible (measured by the PGSs)", NA_real_, NULL, vis_note),
      one_factor_row(main$fit),
      label_rows(pick(main$rows, "corrected_std"), "Corrected for pi", pi, pi_se, main$note),
      if (include_raw) label_rows(pick(main$rows, "corrected_raw"), "Corrected for pi (founder scale)",
        pi, pi_se, main$note))
  }
  for (p in pi_grid) {
    message(who, ": fitting at grid pi = ", p, " ...")
    g <- fit_at(p, grid_pi_se)
    fit_records[[length(fit_records) + 1L]] <- list(kind = "grid", pi = p,
      pi_se = grid_pi_se, fit = g$fit, note = g$note)
    if (is.null(g$rows)) message(who, ": ", g$note)
    results[[length(results) + 1L]] <- if (is.null(g$rows)) fail_row("Grid", p, grid_pi_se, g$note) else rbind(label_rows(pick(g$rows, "corrected_std"), "Grid", p, grid_pi_se, g$note),
        if (include_raw) label_rows(pick(g$rows, "corrected_raw"), "Grid (founder scale)", p, grid_pi_se, g$note))
  }
  if (!long && include_comparison && !is.null(main$rows)) for (v in c("PGS1", "PGS2")) {
    rows_v <- pick(main$rows, "vanilla", v)
    n_v <- if ("N" %in% names(rows_v)) unique(stats::na.omit(rows_v$N)) else integer(0)
    results[[length(results) + 1L]] <- label_rows(rows_v, paste0("Plain ", v, " regression",
      if (length(n_v) == 1L) paste0(" (N = ", n_v, ")")), NA_real_, NULL, "", "Plain PGS regressions (comparison)")
  }
  table <- do.call(rbind, results)
  table <- cbind(siblings = sib_label, table, stringsAsFactors = FALSE)
  attr(table, "fit") <- main$fit
  attr(table, "fit_records") <- fit_records
  attr(table, "analysis_data") <- dat
  attr(table, "one_factor") <- diagnostic
  attr(table, "settings") <- list(design = design,
    covariates = covariates, age = age, cluster_by = cluster_by,
    subject_id = subject_id, diagnostic_cluster = diagnostic_cluster,
    age_center = if (long) age_center else NULL, N = nrow(dat), N_dropped = n_input - nrow(dat),
    dropped = dropped, outcome = outcome, outcome_sd = outcome_sd, family_note = family_note,
    clustered_cross = !long && !is.null(cluster_by), sib_info = sib_info,
    N_families = if (long) length(unique(dat[[family_id]])) else NULL,
    missing_outcome_kept = !long && !drop_missing_outcome && any(is.na(dat[[outcome]])))
  table
}

# Optional family likelihood; activated through resolved simple or advanced controls.
ravel_parental_deviation_mixed <- function(dat, design = c("cross_sectional", "longitudinal"),
    outcome = "Y", pgs_columns, covariates = NULL, binary_outcome = FALSE,
    family_id = "FID", age = "age", age_center = NULL, pi = .4, pi_se = .05,
    pi_grid = NULL, grid_pi_se = pi_se, fit_constraints = "estimated_state",
    am_n_generations = NULL, model = fit_constraints, include_comparison = TRUE,
    include_raw = FALSE, drop_missing_outcome = TRUE, parental_sibs = "none",
    sib_columns = NULL, mixed_settings, diagnostic_cluster = "auto") {
  design <- match.arg(design); long <- design == "longitudinal"
  if (!all(c("ravel_mixed_fit", "ravel_mixed_audit", "ravel_score_diagnostic") %in% getNamespaceExports("RAVEL")))
    stop("The explicit mixed workflow requires a RAVEL build exporting mixed fitting, audits and score diagnostics.")
  if (!identical(parental_sibs, "none") || !is.null(sib_columns))
    stop("The optional mixed workflow does not implement parental-sibling latent blocks.")
  if (!fit_constraints %in% c("estimated_state", "equilibrium", "n_gen"))
    stop("The optional mixed workflow requires estimated_state, equilibrium or n_gen.")
  if (long && binary_outcome) stop("Binary longitudinal mixed outcomes are not supported.")
  allowed <- c("maternal_id", "subject_id", "paternal_id", "family_random", "individual_random",
    "random_slope", "interactions", "interaction_center", "covariate_mapping_confirmed",
    "diagnostic_control", "control")
  if (!is.list(mixed_settings) || is.null(names(mixed_settings)) ||
      anyDuplicated(names(mixed_settings)) || any(!names(mixed_settings) %in% allowed))
    stop("mixed_settings must be a named list containing only documented mixed-model settings.")
  if (!all(c("maternal_id", "subject_id") %in% names(mixed_settings)))
    stop("mixed_settings must explicitly supply maternal_id and subject_id column names.")
  ms <- utils::modifyList(list(paternal_id = NULL, family_random = TRUE,
    individual_random = FALSE, random_slope = NULL, interactions = if (long) age else NULL,
    interaction_center = NULL, covariate_mapping_confirmed = FALSE,
    diagnostic_control = list(), control = list()), mixed_settings, keep.null = TRUE)
  dc <- ms$diagnostic_control
  if (!is.list(dc) || (length(dc) && (is.null(names(dc)) || anyDuplicated(names(dc)) ||
      any(!names(dc) %in% c("age", "target_age", "occasion", "occasion_value", "selection_seed")))))
    stop("diagnostic_control must contain only named age, target_age, occasion, occasion_value or selection_seed settings.")
  if (!long && (isTRUE(ms$individual_random) || !is.null(ms$random_slope)))
    stop("Individual random effects require design = 'longitudinal'.")
  for (id in c("maternal_id", "subject_id"))
    if (!is.character(ms[[id]]) || length(ms[[id]]) != 1L || is.na(ms[[id]]) || !nzchar(ms[[id]]))
      stop("maternal_id and subject_id must name single nonmissing columns.")
  if (!is.character(family_id) || length(family_id) != 1L || is.na(family_id) || !nzchar(family_id))
    stop("family_id must explicitly name the extended-family FID column for mixed inference.")
  if (identical(diagnostic_cluster, "auto")) diagnostic_cluster <-
    if ("FID" %in% names(dat)) "FID" else family_id
  if (!is.null(diagnostic_cluster) && (!is.character(diagnostic_cluster) ||
      length(diagnostic_cluster) != 1L || is.na(diagnostic_cluster) ||
      !nzchar(diagnostic_cluster) || !diagnostic_cluster %in% names(dat)))
    stop("diagnostic_cluster must name an existing outer-family column, be 'auto', or be NULL for independent diagnostic subjects.")
  score_keys <- paste0(rep(c("child", "mother", "father"), each = 2L), 1:2)
  if (length(pgs_columns) != 6L || !setequal(names(pgs_columns), score_keys) ||
      anyNA(pgs_columns) || anyDuplicated(unname(pgs_columns)))
    stop("pgs_columns must name six distinct child1/2, mother1/2 and father1/2 columns.")
  if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi) || pi <= 0 || pi >= 1 ||
      (length(pi_grid) && (!is.numeric(pi_grid) || any(!is.finite(pi_grid) | pi_grid <= 0 | pi_grid >= 1))))
    stop("pi and pi_grid values must lie in (0, 1).")
  needed <- unique(c(outcome, unname(pgs_columns), covariates, family_id,
    ms$maternal_id, ms$subject_id, ms$paternal_id, ms$interactions, ms$random_slope, if (long) age))
  if (length(setdiff(needed, names(dat))))
    stop("Missing column(s): ", paste(setdiff(needed, names(dat)), collapse = ", "))
  numeric_columns <- unique(c(outcome, unname(pgs_columns), covariates,
    ms$interactions, ms$random_slope, if (long) age))
  if (any(!vapply(dat[numeric_columns], is.numeric, logical(1))))
    stop("Outcome, PGSs, covariates and moderators must be numeric.")
  # Keep every original row and ID for diagnostic selection. Neither missing
  # outcomes/scores nor an outcome-dependent fitted sample can choose a replacement.
  diagnostic_data <- dat
  diagnostic_age <- if ("age" %in% names(dc)) dc$age else if (long) age else NULL
  diagnostic_target <- if ("target_age" %in% names(dc)) dc$target_age else
    if (!is.null(diagnostic_age) && diagnostic_age %in% names(dat)) {
      if (diagnostic_age %in% names(ms$interaction_center)) unname(ms$interaction_center[[diagnostic_age]]) else
        if (identical(diagnostic_age, age) && !is.null(age_center)) age_center else
          mean(dat[[diagnostic_age]], na.rm = TRUE)
    } else NULL

  # Match the existing workflow's explicit score/outcome/age row exclusions.
  # Identification and ID nesting are then checked by the new fitter, with no
  # upper-size heuristic or fallback that silently discards extended FIDs.
  n_input <- nrow(dat)
  miss_pgs <- !stats::complete.cases(dat[unname(pgs_columns)])
  miss_y <- is.na(dat[[outcome]]) & (drop_missing_outcome || long)
  miss_age <- if (long) is.na(dat[[age]]) else rep(FALSE, n_input)
  drop <- miss_pgs | miss_y | miss_age
  dropped <- c(incomplete_trio = sum(miss_pgs), missing_outcome = sum(miss_y & !miss_pgs),
    missing_age = sum(miss_age & !miss_pgs & !miss_y))
  dat <- dat[!drop, , drop = FALSE]
  if (any(drop)) message(model, ": mixed workflow dropped ", sum(drop),
    " row(s): ", paste(names(dropped), dropped, collapse = ", "), ".")
  if (nrow(dat) < 50L) stop("Fewer than 50 rows remain after the workflow's exclusions.")
  outcome_sd <- NA_real_
  if (!binary_outcome) {
    outcome_sd <- stats::sd(dat[[outcome]], na.rm = TRUE)
    if (!is.finite(outcome_sd) || outcome_sd <= 0) stop("The outcome has no variation.")
    outcome_center <- mean(dat[[outcome]], na.rm = TRUE)
    dat[[outcome]] <- (dat[[outcome]] - outcome_center) / outcome_sd
    diagnostic_data[[outcome]] <- (diagnostic_data[[outcome]] - outcome_center) / outcome_sd
  }
  if (long && is.null(age_center)) age_center <- mean(dat[[age]], na.rm = TRUE)
  centers <- ms$interaction_center
  if (long && age %in% c(covariates, ms$interactions, ms$random_slope) && !age %in% names(centers))
    centers <- c(centers, stats::setNames(age_center, age))
  if (long && age %in% names(centers)) age_center <- unname(centers[[age]])
  args <- list(dat = dat, col_outcome = outcome,
    outcome = if (binary_outcome) "binary_liability" else "continuous",
    family_id = ms$maternal_id, subject_id = ms$subject_id, cluster = family_id,
    paternal_id = ms$paternal_id, fixed_effects = covariates,
    interactions = ms$interactions, interaction_center = centers,
    family_random = ms$family_random, individual_random = ms$individual_random,
    random_slope = ms$random_slope, fit_constraints = fit_constraints,
    covariate_mapping_confirmed = ms$covariate_mapping_confirmed,
    diagnostic_control = ms$diagnostic_control, control = ms$control, keep_data = TRUE)
  for (role in c("child", "mother", "father")) for (k in 1:2)
    args[[paste0("col_pgs", k, "_", role)]] <- pgs_columns[[paste0(role, k)]]
  if (!is.null(am_n_generations)) args$am_n_generations <- am_n_generations
  base_note <- paste(c("Optional mixed family likelihood; FID-clustered CR0 SEs",
    if (binary_outcome) "conditional probit residual SD = 1" else "effects per workflow-standardized outcome SD",
    if (include_comparison) "Plain PGS mixed comparator is not provided"), collapse = "; ")
  labels <- c(beta_D = "Child (beta_D)", beta_M = "Maternal (beta_M)", beta_F = "Paternal (beta_F)",
    beta_mid = "Parental midpoint (beta_M + beta_F)",
    delta_M_mid = "Mother - father difference (beta_M - beta_F)",
    delta_P_mid = "Father - mother difference (beta_F - beta_M)",
    rho_visible = "Spouse correlation", rho_corrected = "Spouse correlation",
    cov_corrected = "Spouse covariance (founder scale)")
  label_for <- function(parameters) {
    stem <- sub(":.*$", "", parameters); answer <- unname(labels[stem])
    answer[is.na(answer)] <- parameters[is.na(answer)]
    slope <- grepl(":", parameters, fixed = TRUE)
    answer[slope] <- paste0(answer[slope], " x ", sub("^[^:]*:", "", parameters[slope]))
    answer
  }
  format_rows <- function(rows, analysis, p, p_se, note) {
    if (is.null(rows)) rows <- data.frame(parameter = NA_character_, estimate = NA_real_, se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)
    data.frame(siblings = "Trio only", model = model, analysis = analysis, pi = p,
      pi_se = if (is.null(p_se)) NA_real_ else p_se, effect = label_for(rows$parameter),
      parameter = rows$parameter, estimate = rows$estimate, se = rows$se,
      ci_low = rows$ci_low, ci_high = rows$ci_high, p_value = rows$p_value,
      note = note, diagnostic_status = NA_character_, diagnostic_df = NA_real_, diagnostic_N = NA_real_,
      stringsAsFactors = FALSE, row.names = NULL)
  }
  contrasts <- function(fit, scale) {
    tab <- fit$table[fit$table$scale == scale, c("parameter", "estimate", "se", "ci_low", "ci_high", "p_value")]
    V <- fit$coefficient_vcov[[scale]]
    suffixes <- unique(sub("^[^:]+", "", tab$parameter))
    for (suffix in suffixes) {
      parents <- paste0(c("beta_M", "beta_F"), suffix)
      if (!all(parents %in% tab$parameter)) next
      for (name in c("beta_mid", "delta_M_mid", "delta_P_mid")) {
        weight <- switch(name, beta_mid = c(1, 1), delta_M_mid = c(1, -1), delta_P_mid = c(-1, 1))
        estimate <- sum(weight * tab$estimate[match(parents, tab$parameter)])
        variance <- if (is.matrix(V) && all(parents %in% rownames(V)))
          drop(t(weight) %*% V[parents, parents, drop = FALSE] %*% weight) else NA_real_
        se <- if (is.finite(variance) && variance >= -1e-10) sqrt(max(0, variance)) else NA_real_
        tab <- rbind(tab, data.frame(parameter = paste0(name, suffix), estimate = estimate, se = se,
          ci_low = estimate - stats::qnorm(.975) * se, ci_high = estimate + stats::qnorm(.975) * se,
          p_value = 2 * stats::pnorm(-abs(estimate / se))))
      }
    }
    tab[order(grepl(":", tab$parameter, fixed = TRUE),
      match(sub(":.*$", "", tab$parameter), names(labels))), , drop = FALSE]
  }
  run_diagnostic <- function(fit) {
    d <- diagnostic_data
    invalid_design <- function(reason) list(status = "invalid_design", reason = reason,
      statistic = NA_real_, df = NA_real_, p_value = NA_real_, N = NA_integer_,
      workflow = list(N_input = n_input, N_primary = nrow(dat),
        primary_excluded_rows = which(drop), primary_fit_converged = isTRUE(fit$converged),
        original_cluster = diagnostic_cluster, subject_id = ms$subject_id))
    primary_ids <- as.character(d[[family_id]])
    diagnostic_ids <- as.character(d[[if (is.null(diagnostic_cluster)) ms$subject_id else diagnostic_cluster]])
    if (anyNA(primary_ids) || any(!nzchar(primary_ids)) ||
        anyNA(diagnostic_ids) || any(!nzchar(diagnostic_ids)))
      return(invalid_design("Original primary and diagnostic family/individual IDs must be observed and nonempty before diagnostic selection."))
    mapping <- unique(data.frame(primary = primary_ids, diagnostic = diagnostic_ids))
    if (anyDuplicated(mapping$primary))
      return(invalid_design(paste0("Diagnostic grouping would split known extended families in '", family_id,
        "'. Choose diagnostic_cluster='", family_id,
        "' or a coarser family column; NULL requires at most one distinct subject per known extended family.")))
    scaling <- fit$score_scaling
    if (is.null(scaling)) {
      # Identical role/unique-person scaling convention, available even when
      # a primary fit fails before returning its numerical scaling metadata.
      scaling <- data.frame(column = unname(pgs_columns), center = 0, scale = 1)
      for (j in seq_along(pgs_columns)) {
        ids <- dat[[if (grepl("child", names(pgs_columns)[j])) ms$subject_id else ms$maternal_id]]
        x <- dat[[pgs_columns[j]]][!duplicated(ids)]
        scaling$center[j] <- mean(x); scaling$scale[j] <- stats::sd(x)
      }
    }
    if (any(!is.finite(scaling$center) | !is.finite(scaling$scale) | scaling$scale <= 0))
      return(list(status = "invalid_data", reason = "PGS scaling is undefined in the primary analysis sample.", N = n_input))
    for (j in seq_len(nrow(scaling))) d[[scaling$column[j]]] <-
      (d[[scaling$column[j]]] - scaling$center[j]) / scaling$scale[j]
    diagnostic_args <- list(dat = d, col_outcome = outcome,
      outcome = if (binary_outcome) "binary_liability" else "continuous",
      family_id = ms$maternal_id, subject_id = ms$subject_id, cluster = diagnostic_cluster,
      covariates = unique(c(covariates, ms$interactions, ms$random_slope)), interactions = ms$interactions,
      age = diagnostic_age, target_age = diagnostic_target,
      fit_constraints = fit_constraints, pi = pi, pi_se = NULL)
    if (!is.null(am_n_generations)) diagnostic_args$am_n_generations <- am_n_generations
    for (role in c("child", "mother", "father")) for (k in 1:2)
      diagnostic_args[[paste0("col_pgs", k, "_", role)]] <- pgs_columns[[paste0(role, k)]]
    diagnostic_args <- utils::modifyList(diagnostic_args, dc, keep.null = TRUE)
    warnings <- character()
    result <- tryCatch(withCallingHandlers(do.call(RAVEL::ravel_score_diagnostic, diagnostic_args),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) list(status = "error", reason = conditionMessage(e), N = NA_integer_))
    result$workflow_warnings <- unique(warnings)
    result$workflow <- list(N_input = n_input, N_primary = nrow(dat),
      primary_excluded_rows = which(drop), primary_fit_converged = isTRUE(fit$converged),
      original_cluster = diagnostic_cluster, subject_id = ms$subject_id,
      score_and_outcome_units = "Primary role PGS constants and workflow outcome constants, applied before auxiliary selection")
    result
  }
  fit_at <- function(p, p_se, main = FALSE) {
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(do.call(RAVEL::ravel_mixed_fit,
      c(args, list(pi = p, pi_se = p_se, diagnostic = FALSE))), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = identity)
    failed <- inherits(fit, "error")
    failure <- if (failed) conditionMessage(fit) else NULL
    if (failed) fit <- NULL
    diagnostic <- if (main) run_diagnostic(fit) else NULL
    if (!is.null(fit) && main) fit$one_factor <- diagnostic
    note <- paste(unique(c(base_note, if (failed) paste("Fit failed:", failure) else
      paste0("inference status: ", fit$inference_status), fit$messages, warnings)), collapse = "; ")
    audit <- if (failed) list(status = data.frame(converged = FALSE, optimizer_converged = FALSE,
      inference_status = "fit_error"), extraction_status = data.frame(stage = "source_fit", status = "error",
      note = failure), notes = note, one_factor = diagnostic) else
      tryCatch(RAVEL::ravel_mixed_audit(fit), error = function(e)
        list(extraction_status = data.frame(stage = "audit", status = "error", note = conditionMessage(e)),
          notes = note, one_factor = diagnostic))
    attr(audit, "source_metadata") <- data.frame(model = model, outcome = outcome,
      assumption = if (main) "main" else "grid", pi = p,
      supplied_pi_se = if (is.null(p_se)) 0 else p_se, stringsAsFactors = FALSE)
    list(fit = fit, audit = audit, note = note, diagnostic = diagnostic)
  }
  rows_for <- function(record, scale, analysis, p, p_se) {
    available <- !is.null(record$fit) && isTRUE(record$fit$converged)
    rows <- format_rows(if (available) contrasts(record$fit, scale) else NULL,
      analysis, p, p_se, record$note)
    state <- record$audit$state
    key <- if (scale == "visible") "rho_visible" else "rho"
    rho <- if (is.data.frame(state)) state[state$lhs == key, , drop = FALSE] else NULL
    if (scale == "corrected_raw") {
      genetic <- record$audit$genetic_covariances
      rho <- if (is.data.frame(genetic)) genetic[genetic$block == "total_founder_covariance" &
        ((genetic$lhs == "mother" & genetic$rhs == "father") |
         (genetic$lhs == "father" & genetic$rhs == "mother")), , drop = FALSE] else NULL
    }
    if (available && is.data.frame(rho) && nrow(rho)) {
      rho_rows <- data.frame(parameter = switch(scale, visible = "rho_visible",
        corrected_std = "rho_corrected", corrected_raw = "cov_corrected"),
        estimate = rho$estimate, se = rho$se, ci_low = rho$ci_low, ci_high = rho$ci_high,
        p_value = 2 * stats::pnorm(-abs(rho$estimate / rho$se)))
      rows <- rbind(rows, format_rows(rho_rows, analysis, p, p_se,
        paste(record$note, if (scale == "corrected_raw")
          "Founder-scale spouse-covariance SE and CI condition on supplied pi; external pi uncertainty is not included." else
          "Spouse-correlation SE and CI condition on supplied pi; external pi uncertainty is not included.", sep = "; ")))
    }
    rows
  }
  main <- fit_at(pi, pi_se, TRUE)
  records <- list(main = list(kind = "main", pi = pi, pi_se = pi_se, fit = main$fit, note = main$note))
  audits <- list(main = main$audit)
  diagnostic <- if (!is.null(main$diagnostic)) main$diagnostic else
    list(status = "needs_fit", reason = main$note, N = nrow(dat))
  diagnostic_notes <- unique(c(diagnostic$test, diagnostic$reason, diagnostic$auxiliary$warnings))
  test <- format_rows(NULL, "One-factor test (visible fit)", NA_real_, NULL,
    paste(diagnostic_notes[!is.na(diagnostic_notes) & nzchar(diagnostic_notes)], collapse = "; "))
  test$parameter <- "one_factor_test"; test$effect <- "One-factor test"
  test$diagnostic_status <- diagnostic$status
  test$diagnostic_df <- if (is.null(diagnostic$df)) NA_real_ else diagnostic$df
  test$diagnostic_N <- if (is.null(diagnostic$N)) NA_real_ else diagnostic$N
  if (identical(diagnostic$status, "ok")) {
    test$effect <- paste0("One-factor test: chi-square (", diagnostic$df, " df)")
    test$estimate <- diagnostic$statistic; test$p_value <- diagnostic$p_value
  }
  results <- list(rows_for(main, "visible", "Visible (measured by the PGSs)", pi, pi_se), test,
    rows_for(main, "corrected_std", "Corrected for pi", pi, pi_se))
  if (include_raw) results[[length(results) + 1L]] <- rows_for(main, "corrected_raw",
    "Corrected for pi (founder scale)", pi, pi_se)
  for (j in seq_along(pi_grid)) {
    p <- pi_grid[j]; grid <- fit_at(p, grid_pi_se)
    key <- paste0("grid_", j)
    records[[key]] <- list(kind = "grid", pi = p, pi_se = grid_pi_se, fit = grid$fit, note = grid$note)
    audits[[key]] <- grid$audit
    results[[length(results) + 1L]] <- rows_for(grid, "corrected_std", "Grid", p, grid_pi_se)
    if (include_raw) results[[length(results) + 1L]] <- rows_for(grid, "corrected_raw",
      "Grid (founder scale)", p, grid_pi_se)
  }
  table <- do.call(rbind, results); rownames(table) <- NULL
  attr(table, "fit") <- main$fit; attr(table, "fit_records") <- records
  attr(table, "mixed_audits") <- audits; attr(table, "analysis_data") <- dat
  attr(table, "one_factor") <- diagnostic
  attr(table, "settings") <- list(engine = "mixed", design = design, covariates = covariates,
    age = age, age_center = if (long) age_center else NULL, cluster_by = family_id,
    diagnostic_cluster = diagnostic_cluster,
    maternal_id = ms$maternal_id, subject_id = ms$subject_id, mixed_settings = ms,
    N = nrow(dat), N_dropped = sum(drop), dropped = dropped, outcome = outcome,
    outcome_sd = outcome_sd, family_note = base_note, clustered_cross = FALSE, sib_info = "",
    N_families = length(unique(dat[[ms$maternal_id]])),
    N_extended_families = length(unique(dat[[family_id]])),
    missing_outcome_kept = !long && !drop_missing_outcome && any(is.na(dat[[outcome]])))
  table
}

# Flatten audit data frames and sampling-covariance matrices for every main/grid
# fit. Full nested objects, diagnostic details and notes also remain in the RDS.
write_ravel_mixed_workflow_audits <- function(audits, prefix) {
  tables <- list()
  add <- function(kind, frame, source, scale = NA_character_) {
    if (!is.data.frame(frame) || !nrow(frame)) return(invisible(NULL))
    tables[[kind]][[length(tables[[kind]]) + 1L]] <<-
      cbind(source_fit = source, covariance_scale = scale,
        metadata[rep(1L, nrow(frame)), , drop = FALSE], frame, stringsAsFactors = FALSE)
  }
  matrix_table <- function(x) {
    if (is.null(x)) return(NULL)
    ij <- which(lower.tri(x, diag = TRUE), arr.ind = TRUE)
    data.frame(parameter_i = rownames(x)[ij[, 1]], parameter_j = colnames(x)[ij[, 2]],
      covariance = x[ij], stringsAsFactors = FALSE)
  }
  for (source in names(audits)) {
    audit <- audits[[source]]
    metadata <- attr(audit, "source_metadata")
    if (is.null(metadata)) metadata <- data.frame(row.names = 1L)
    for (name in names(audit)) {
      value <- audit[[name]]
      if (is.data.frame(value)) add(name, value, source)
      else if (is.matrix(value)) add(name, matrix_table(value), source)
      else if (name == "coefficient_vcov") for (scale in names(value))
        add(name, matrix_table(value[[scale]]), source, scale)
    }
    if (length(audit$notes)) add("audit_notes", data.frame(note = audit$notes), source)
  }
  paths <- character()
  for (name in names(tables)) {
    file <- paste0(prefix, "_mixed_", name, ".csv")
    utils::write.csv(do.call(rbind, tables[[name]]), file, row.names = FALSE)
    paths[name] <- file
  }
  invisible(paths)
}

# Write a copy of the table that Excel opens straight into columns. Excel's
# separator follows the computer's regional settings (';' where the decimal
# mark is a comma), so the decimal mark is detected and a 'sep=' first line
# tells Excel the separator.
write_excel_csv <- function(table, path, decimal = excel_decimal_mark()) {
  sep <- if (decimal == ",") ";" else ","
  con <- file(path, "w"); on.exit(close(con))
  writeLines(paste0("sep=", sep), con)
  utils::write.table(table, con, sep = sep, dec = decimal, row.names = FALSE, qmethod = "double", na = "")
  invisible(path)
}
excel_decimal_mark <- function() {
  dec <- tryCatch({
    if (.Platform$OS.type == "windows") {
      utils::readRegistry("Control Panel\\International", hive = "HCU")$sDecimal
    } else if (Sys.info()[["sysname"]] == "Darwin") {
      read_default <- function(key) suppressWarnings(system2("defaults", c("read", "-g", key),
        stdout = TRUE, stderr = FALSE))
      custom <- grep("^\\s*0 = ", read_default("AppleICUNumberSymbols"), value = TRUE)
      if (length(custom)) gsub('^\\s*0 = "?|"?;\\s*$', "", custom[1]) else {
        old <- Sys.getlocale("LC_MONETARY"); on.exit(Sys.setlocale("LC_MONETARY", old), add = TRUE)
        loc <- paste0(read_default("AppleLocale")[1], ".UTF-8")
        if (nzchar(suppressWarnings(Sys.setlocale("LC_MONETARY", loc))))
          Sys.localeconv()[["mon_decimal_point"]] else "."
      }
    } else "."
  }, error = function(e) ".")
  if (identical(dec, ",")) "," else "."
}

# Read a data file, choosing the reader from its extension (.rds, .csv, .tsv/.txt).
read_trio_data <- function(path) {
  if (!file.exists(path)) stop("Data file not found: ", path, " (working directory: ", getwd(), ")")
  dat <- switch(tolower(tools::file_ext(path)),
    rds = readRDS(path),
    csv = utils::read.csv(path),
    tsv = , txt = utils::read.delim(path),
    tryCatch(readRDS(path), error = function(e) utils::read.csv(path)))  # unknown extension
  if (!is.data.frame(dat))
    stop("The file must contain a data frame (one row per child, or per child-occasion).")
  dat <- as.data.frame(dat)                                   # e.g. tibbles
  dat[] <- lapply(dat, function(x)                            # e.g. SPSS-labelled columns
    if (inherits(x, "haven_labelled")) as.vector(unclass(x)) else x)
  dat
}

# Console summary: one line per analysis with every beta as "estimate (SE)",
# the difference test, and the age slopes in the longitudinal design.
# The CSV has every coefficient with SE, 95% CI and p-value.
print_parental_deviation <- function(table, digits = 3) {
  f <- function(x) if (is.na(x)) "NA" else formatC(round(x, digits) + 0, digits = digits, format = "f")
  fp <- function(x) if (is.na(x)) "NA" else if (x < .001) "<0.001" else formatC(x, digits = 3, format = "f")
  s <- attr(table, "settings")
  cat("\nRAVEL results (", s$design, ", N rows = ", s$N,
      if (!is.null(s$N_families)) paste0(" in ", s$N_families, " families"),
      if (isTRUE(s$N_dropped > 0)) paste0("; rows dropped: ", s$dropped[["incomplete_trio"]],
        " incomplete trios, ", s$dropped[["missing_outcome"]], " missing outcome",
        if (s$design == "longitudinal") paste0(", ", s$dropped[["missing_age"]], " missing age")),
      if (!is.null(s$age_center)) paste0(", main effects at age ", s$age_center), ")\n", sep = "")
  if (!is.null(s$family_note) && nzchar(s$family_note)) cat(s$family_note, ".", if (isTRUE(s$clustered_cross))
    " The plain regressions keep ordinary SEs (rows treated as independent).", "\n", sep = "")
  if (isTRUE(s$missing_outcome_kept))
    cat(if (identical(s$engine, "mixed"))
      "Missing outcomes retain their score information in the family likelihood.\n" else
      "Trios with a missing outcome are kept in the RAVEL fits (FIML) but not in the plain regressions.\n")
  if (is.finite(s$outcome_sd)) {
    cat("Outcome '", s$outcome, "' standardized (SD = ", signif(s$outcome_sd, 4),
        "): effects are per SD of the outcome.\n", sep = "")
  } else cat("Outcome '", s$outcome, "' is binary: effects are on ",
    if (identical(s$engine, "mixed")) "the conditional probit scale (occasion residual SD = 1).\n" else
      "RAVEL's probit liability scale.\n", sep = "")
  cat("Each cell is estimate (SE). Corrected rows use the standardized total-liability scale",
      if (any(grepl("founder", table$analysis))) " unless marked founder scale", ".\n", sep = "")
  line <- "  %-34s %5s %5s %15s %15s %15s %15s %15s %7s\n"
  if (is.null(table$siblings)) table$siblings <- "Trio only"
  sib_info <- attr(table, "sib_info")
  keys <- unique(table[, c("siblings", "model", "analysis", "pi", "pi_se")])
  for (i in seq_len(nrow(keys))) {
    new_sibs <- i == 1L || keys$siblings[i] != keys$siblings[i - 1L]
    if (new_sibs && length(unique(table$siblings)) > 1L) {
      info <- if (!is.null(sib_info)) sib_info[[keys$siblings[i]]] else ""
      cat("\n---------- ", keys$siblings[i], " ----------\n", if (length(info) && nzchar(info)) paste0(info, "\n"), sep = "")
    }
    if (new_sibs || keys$model[i] != keys$model[i - 1L]) {
      cat("\n", keys$model[i], "\n", sep = "")
      cat(sprintf(line, "Analysis", "pi", "SE", "Child", "Mother", "Father", "Midpoint", "M - F", "p(M-F)"))
    }
    r <- table[table$siblings == keys$siblings[i] & table$model == keys$model[i] &
               table$analysis == keys$analysis[i] & identical_or_na(table$pi, keys$pi[i]), , drop = FALSE]
    cell <- function(par) {
      k <- which(r$parameter == par)[1]
      if (is.na(k)) "" else paste0(f(r$estimate[k]), " (", f(r$se[k]), ")")
    }
    if (identical(r$parameter[1], "one_factor_test")) {
      cat(if (is.na(r$estimate[1])) sprintf("  One-factor test: %s\n", r$note[1]) else
        sprintf("  One-factor test (visible fit): chi-square = %s, df = %s, p = %s%s\n", f(r$estimate[1]),
          sub("^.*\\((\\d+) df\\).*$", "\\1", r$effect[1]), fp(r$p_value[1]),
          if (nzchar(r$note[1])) paste0(" (", r$note[1], ")") else ""))
      next
    }
    if (all(is.na(r$estimate))) {
      cat(sprintf("  %-34s %5s %5s   no results: see Notes\n", substr(keys$analysis[i], 1, 34),
        if (is.na(keys$pi[i])) "" else f(keys$pi[i]), if (is.na(keys$pi_se[i])) "" else f(keys$pi_se[i])))
      next
    }
    pd <- r$p_value[r$parameter == "delta_M_mid"]
    cat(sprintf(line, substr(keys$analysis[i], 1, 34),
      if (is.na(keys$pi[i])) "" else f(keys$pi[i]), if (is.na(keys$pi_se[i])) "" else f(keys$pi_se[i]),
      cell("beta_D"), cell("beta_M"), cell("beta_F"), cell("beta_mid"), cell("delta_M_mid"),
      if (length(pd)) fp(pd[1]) else ""))
    sib_par <- intersect(c("beta_MS", "delta_M", "beta_FS", "delta_F"), r$parameter)
    if (length(sib_par)) {
      pv <- function(par) { x <- r$p_value[r$parameter == par][1]; if (grepl("^delta", par)) paste0(" p=", fp(x)) else "" }
      short <- c(beta_MS = "Mat sib", delta_M = "Mother - mat sib", beta_FS = "Pat sib", delta_F = "Father - pat sib")
      cat("      siblings: ", paste0(short[sib_par], " ", vapply(sib_par, cell, ""), vapply(sib_par, pv, ""),
        collapse = ", "), "\n", sep = "")
    }
    slopes <- r$parameter[grepl(":", r$parameter, fixed = TRUE)]
    if (length(slopes))
      cat("      slopes: ", paste0(slopes, " ", vapply(slopes, cell, ""), collapse = ", "), "\n", sep = "")
  }
  has_note <- !is.na(table$note) & nzchar(table$note) & !(table$parameter %in% "one_factor_test")  # printed inline
  notes <- unique(paste0(ifelse(table$siblings == "Trio only", "", paste0(table$siblings, ", ")),
                         table$model, ": ", table$note)[has_note])
  if (length(notes)) cat("\nNotes:\n", paste0("  - ", notes, "\n"), sep = "")
  invisible(table)
}
identical_or_na <- function(x, y) if (is.na(y)) is.na(x) else (!is.na(x) & abs(x - y) < 1e-12)

## ---- RUN --------------------------------------------------------------------

if (!isTRUE(getOption("ravel.workflow_functions_only", FALSE))) {
if (is.null(data_file)) {
  data("ravel_example", package = "RAVEL")
  dat <- ravel_example
  names(dat) <- sub("^pgs", "PGS", names(dat))   # example uses lowercase names
  message("No data_file set: running on RAVEL's bundled example (cross-sectional).")
  design <- "cross_sectional"
} else {
  dat <- read_trio_data(data_file)
  message("Read ", nrow(dat), " rows and ", ncol(dat), " columns from ", data_file)
}

# Analyse one outcome with every model and sibling configuration; returns its table.
if ((!is.null(mixed_settings) || random_family || random_individual ||
     !is.null(random_slope) || !is.null(genetic_interactions)) && !is.null(sib_columns))
  stop("The optional mixed workflow does not implement parental-sibling latent blocks; leave sib_columns=NULL.")
sib_configs <- if (is.null(sib_columns)) "none" else c("none", "maternal", "paternal", "both")
run_outcome <- function(y) {
  is_binary <- if (is.logical(binary_outcome)) isTRUE(binary_outcome) else y %in% binary_outcome
  tables <- list(); fits <- list(); diagnostics <- list(); sib_info <- character(0); settings <- NULL
  mixed_fit_records <- list(); mixed_audits <- list()
  for (ps in sib_configs) {
    parts <- tryCatch(lapply(names(models), function(m) {
      do.call(ravel_parental_deviation, c(list(dat, design = design, outcome = y,
        pgs_columns = pgs_columns, covariates = covariates, binary_outcome = is_binary,
        family_id = family_id, subject_id = subject_id, diagnostic_cluster = diagnostic_cluster,
        mixed_settings = mixed_settings, maternal_id = maternal_id,
        random_family = random_family, random_individual = random_individual,
        random_slope = random_slope, genetic_interactions = genetic_interactions,
        covariate_mapping_confirmed = covariate_mapping_confirmed,
        age = age, age_center = age_center, pi = pi, pi_se = pi_se,
        pi_grid = pi_grid, grid_pi_se = grid_pi_se, model = m, include_raw = include_raw,
        drop_missing_outcome = drop_missing_outcome, parental_sibs = ps, sib_columns = sib_columns,
        include_comparison = TRUE), models[[m]]))
    }), error = function(e) e)
    if (inherits(parts, "error")) {
      message("Sibling configuration '", ps, "' could not be analysed: ", conditionMessage(parts))
      next
    }
    is_comp <- function(p) p$model == "Plain PGS regressions (comparison)"
    result <- do.call(rbind, lapply(parts, function(p) p[!is_comp(p), ]))
    comparison <- Filter(nrow, lapply(parts, function(p) p[is_comp(p), ]))   # same in every model
    if (length(comparison)) result <- rbind(result, comparison[[1]])        # one copy, last
    tables[[ps]] <- result
    fits[[ps]] <- stats::setNames(lapply(parts, attr, "fit", exact = TRUE), names(models))
    diagnostics[[ps]] <- stats::setNames(lapply(parts, attr, "one_factor", exact = TRUE), names(models))
    for (j in seq_along(parts)) {
      audits <- attr(parts[[j]], "mixed_audits", exact = TRUE)
      records <- attr(parts[[j]], "fit_records", exact = TRUE)
      if (length(audits)) for (key in names(records)) {
        source_key <- paste(ps, names(models)[j], key, sep = "|")
        mixed_fit_records[[source_key]] <- records[[key]]
        mixed_audits[[source_key]] <- audits[[key]]
      }
    }
    st <- attr(parts[[1]], "settings")
    sib_info[[result$siblings[1]]] <- st$sib_info
    if (is.null(settings)) settings <- st
  }
  if (!length(tables)) stop("no sibling configuration could be analysed")
  result <- do.call(rbind, tables)
  rownames(result) <- NULL
  result <- cbind(outcome = y, result, stringsAsFactors = FALSE)
  attr(result, "settings") <- settings
  attr(result, "sib_info") <- sib_info
  attr(result, "fits") <- fits      # fits[[sibling configuration]][[model]]
  attr(result, "diagnostics") <- diagnostics
  attr(result, "mixed_fit_records") <- mixed_fit_records
  attr(result, "mixed_audits") <- mixed_audits
  result
}

if (anyDuplicated(outcome)) stop("Each outcome should be listed once.")
results <- list()
for (y in outcome) {
  message("\n==== Outcome: ", y, " ====")
  res <- tryCatch(run_outcome(y), error = function(e) e)
  if (inherits(res, "error")) {
    message("Outcome ", y, " could not be analysed: ", conditionMessage(res))
    next
  }
  results[[y]] <- res
  cat("\n================ Outcome: ", y, " ================\n", sep = "")
  print_parental_deviation(res)
  base <- sub("[.]csv$", "", output_csv)
  file_y <- paste0(base, "_", gsub("[^A-Za-z0-9_.-]+", "_", y), ".csv")
  utils::write.csv(res, file_y, row.names = FALSE)
  excel_y <- sub("[.]csv$", "_excel.csv", file_y)
  write_excel_csv(res, excel_y)
  rds_y <- sub("[.]csv$", ".rds", file_y)
  saveRDS(res, rds_y)
  if (length(attr(res, "mixed_audits", exact = TRUE))) {
    audit_paths <- write_ravel_mixed_workflow_audits(
      attr(res, "mixed_audits", exact = TRUE), sub("[.]csv$", "", file_y))
    message("Mixed-model audit tables written: ", length(audit_paths),
      "; all main/grid fits remain in the result RDS.")
  }
  message("Table for ", y, " written to ", normalizePath(file_y),
          "\nExcel copy (opens in columns): ", normalizePath(excel_y),
          "\nFitted objects and one-factor diagnostics: ", normalizePath(rds_y))
}

message("\nDone. Analysed: ", if (length(results)) paste(names(results), collapse = ", ") else "none",
        if (length(setdiff(outcome, names(results))))
          paste0("; could not analyse: ", paste(setdiff(outcome, names(results)), collapse = ", ")))

# `results` holds one table per outcome; each keeps its fits by sibling
# configuration ("none", "maternal", "paternal", "both") and model, e.g.:
#   print_simple(attr(results[["Y"]], "fits")[["none"]][["Estimated state"]], parent_contrasts = TRUE)
}
