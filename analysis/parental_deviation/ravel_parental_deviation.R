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
##   For each model it also runs RAVEL's one-factor test on the main fit: do
##   the PGS1 and PGS2 associations with the outcome agree with one shared
##   genetic factor? (chi-square, 3 df; a small p-value is evidence against
##   one factor.) It concerns the visible fit, so it is not repeated over the
##   grid. It needs a cross-sectional continuous outcome without covariates.
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
ravel_is_current <- function() {
  requireNamespace("RAVEL", quietly = TRUE) &&
    all(c("ravel_family_fit", "ravel_family_regression") %in% getNamespaceExports("RAVEL")) &&
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

family_id  <- "FID"               # family identifier (rows of a family share it). Required
                                  # for longitudinal data. Cross-sectional: families with more
                                  # than one row (e.g. siblings) get SEs clustered by family;
                                  # NULL = rows are independent trios (no family handling).
# Longitudinal design only
age        <- "age"               # age at measurement
age_center <- 10                  # main effects are reported at this age

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

## ---- ANALYSIS FUNCTION (no need to edit) ------------------------------------

ravel_parental_deviation <- function(dat, design = c("cross_sectional", "longitudinal"),
                                     outcome = "Y", pgs_columns, covariates = NULL,
                                     binary_outcome = FALSE, family_id = "FID",
                                     age = "age", age_center = NULL, pi = .4,
                                     pi_se = .05, pi_grid = NULL, grid_pi_se = pi_se,
                                     fit_constraints = "estimated_state",
                                     am_n_generations = NULL, model = fit_constraints,
                                     include_comparison = TRUE, include_raw = FALSE,
                                     drop_missing_outcome = TRUE, parental_sibs = "none",
                                     sib_columns = NULL) {
  design <- match.arg(design)
  long <- design == "longitudinal"
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
    if (anyNA(ids)) {
      k <- sum(is.na(ids))
      ids <- as.character(ids)
      ids[is.na(ids)] <- paste0(".no_family_id_", seq_len(k))
      dat[[family_id]] <- ids
      message(who, ": ", k, " row(s) have no family ID; each is treated as its own family.")
    }
    sizes <- table(dat[[family_id]])
    if (length(sizes) < 2L || max(sizes) > 20L) {   # e.g. a constant PLINK-style FID
      family_note <- paste0("Column '", family_id, "' does not look like a family ID (one value covers ",
        max(sizes), " rows), so it is ignored and rows are treated as independent trios")
      message(who, ": ", family_note, ".")
      use_family <- FALSE
    } else {
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
    dat[v_all] <- (dat[v_all] - mean(x)) / stats::sd(x)
  }
  for (v in covariates) {
    sv <- stats::sd(dat[[v]], na.rm = TRUE)
    if (is.finite(sv) && sv > 0) dat[[v]] <- (dat[[v]] - mean(dat[[v]], na.rm = TRUE)) / sv
  }
  outcome_sd <- NA_real_
  if (!binary_outcome) {
    outcome_sd <- stats::sd(dat[[outcome]], na.rm = TRUE)
    if (!is.finite(outcome_sd) || outcome_sd == 0) stop("The outcome has no variation.")
    dat[[outcome]] <- (dat[[outcome]] - mean(dat[[outcome]], na.rm = TRUE)) / outcome_sd
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

  # RAVEL's one-factor (score-signal) test on a fitted model, as one table row.
  one_factor_row <- function(fit) {
    why <- if (long) "not available for the longitudinal design (repeated rows)" else
      if (binary_outcome) "not available for a binary outcome" else
      if (length(covariates)) "not available with covariates" else
      if (n_multi > 0L) "not available: some families have more than one row (needs one independent trio per row)" else
      if (length(sib_used) && any(!stats::complete.cases(dat[sib_used])))
        "not available: some families lack sibling scores (the test needs every family's scores)" else NULL
    z <- if (is.null(why)) tryCatch(suppressWarnings(ravel_test_score_signal(fit)), error = function(e) e)
    if (is.null(why) && inherits(z, "error")) why <- paste("failed:", conditionMessage(z))
    if (is.null(why) && !identical(z$status, "ok")) why <- paste0("not available (", z$status, "): ", z$reason)
    note <- if (!is.null(why)) why else
      if (fit_constraints == "n_gen") paste("conditional on pi =", pi, "used in fitting") else ""
    data.frame(model = model, analysis = "One-factor test (visible fit)", pi = NA_real_, pi_se = NA_real_,
      effect = if (is.null(why)) paste0("One-factor test: chi-square (", z$df, " df)") else "One-factor test",
      parameter = "one_factor_test", estimate = if (is.null(why)) z$statistic else NA_real_,
      se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = if (is.null(why)) z$p_value else NA_real_, note = note, stringsAsFactors = FALSE)
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
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, note = note, stringsAsFactors = FALSE)
  label_rows <- function(rows, analysis, p, p_se, note, model_label = model) {
    if (is.null(rows) || !nrow(rows)) return(NULL)
    data.frame(model = model_label, analysis = analysis, pi = p, pi_se = if (is.null(p_se)) NA_real_ else p_se,
      effect = label_for(rows$parameter), parameter = rows$parameter,
      estimate = rows$estimate, se = rows$se, ci_low = rows$ci_low, ci_high = rows$ci_high,
      p_value = rows$p_value, note = note, stringsAsFactors = FALSE, row.names = NULL)
  }

  message(who, ": fitting at pi = ", pi, " ...")
  main <- fit_at(pi, pi_se)
  if (is.null(main$rows)) {
    message(who, ": ", main$note)
    results <- list(fail_row("Visible (measured by the PGSs)", NA_real_, NULL, main$note),
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
  attr(table, "settings") <- list(design = design,
    age_center = if (long) age_center else NULL, N = nrow(dat), N_dropped = n_input - nrow(dat),
    dropped = dropped, outcome = outcome, outcome_sd = outcome_sd, family_note = family_note,
    clustered_cross = !long && !is.null(cluster_by), sib_info = sib_info,
    N_families = if (long) length(unique(dat[[family_id]])) else NULL,
    missing_outcome_kept = !long && !drop_missing_outcome && any(is.na(dat[[outcome]])))
  table
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
    cat("Trios with a missing outcome are kept in the RAVEL fits (FIML) but not in the plain regressions.\n")
  if (is.finite(s$outcome_sd)) {
    cat("Outcome '", s$outcome, "' standardized (SD = ", signif(s$outcome_sd, 4),
        "): effects are per SD of the outcome.\n", sep = "")
  } else cat("Outcome '", s$outcome, "' is binary: effects are on RAVEL's probit liability scale.\n", sep = "")
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
sib_configs <- if (is.null(sib_columns)) "none" else c("none", "maternal", "paternal", "both")
run_outcome <- function(y) {
  is_binary <- if (is.logical(binary_outcome)) isTRUE(binary_outcome) else y %in% binary_outcome
  tables <- list(); fits <- list(); sib_info <- character(0); settings <- NULL
  for (ps in sib_configs) {
    parts <- tryCatch(lapply(names(models), function(m) {
      do.call(ravel_parental_deviation, c(list(dat, design = design, outcome = y,
        pgs_columns = pgs_columns, covariates = covariates, binary_outcome = is_binary,
        family_id = family_id, age = age, age_center = age_center, pi = pi, pi_se = pi_se,
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
    fits[[ps]] <- stats::setNames(lapply(parts, attr, "fit"), names(models))
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
  message("Table for ", y, " written to ", normalizePath(file_y),
          "\nExcel copy (opens in columns): ", normalizePath(excel_y))
}

message("\nDone. Analysed: ", if (length(results)) paste(names(results), collapse = ", ") else "none",
        if (length(setdiff(outcome, names(results))))
          paste0("; could not analyse: ", paste(setdiff(outcome, names(results)), collapse = ", ")))

# `results` holds one table per outcome; each keeps its fits by sibling
# configuration ("none", "maternal", "paternal", "both") and model, e.g.:
#   print_simple(attr(results[["Y"]], "fits")[["none"]][["Estimated state"]], parent_contrasts = TRUE)
