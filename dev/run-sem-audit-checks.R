# Independent checks of the separate source-SEM audit tables.
# Run from the repository root: Rscript dev/run-sem-audit-checks.R
pkgload::load_all(".", quiet = TRUE)
RAVEL:::.ensure_lavaan_runtime()
source("analysis/parental_deviation/ravel_sem_audit.R")
dir.create("dev/sem-audit-check", recursive = TRUE, showWarnings = FALSE)
checks <- list()
check <- function(label, expression) {
  force(expression)
  checks[[length(checks) + 1L]] <<- data.frame(check = label, status = "passed")
  cat("PASS:", label, "\n")
}
equal <- function(observed, expected, tolerance = 1e-7) {
  comparison <- all.equal(as.numeric(observed), as.numeric(expected), tolerance = tolerance,
    check.attributes = FALSE)
  if (!isTRUE(comparison)) stop(paste(comparison, collapse = "; "))
}
wrap <- function(f, key = "main", pi_se = .02) {
  answer <- data.frame()
  attr(answer, "settings") <- list(outcome = f$outcome_variable)
  attr(answer, "fits") <- stats::setNames(list(list(fit = f, kind = "main", pi = f$h,
    pi_se = pi_se, outcome = f$outcome_variable, parental_sibs = f$parental_sibs,
    model = f$constraint_spec$fit_constraints, genetic_roles = c("child", "mother", "father"))), key)
  answer
}
data("ravel_example", package = "RAVEL")
dat <- ravel_example
names(dat) <- sub("^pgs", "PGS", names(dat))
continuous <- ravel_fit(dat, fit_constraints = "estimated_state", pi = .4)
audit <- ravel_sem_audit(wrap(continuous))
lav <- continuous$fit
check("all requested source tables are present and retain supplied-pi uncertainty scope", {
  stopifnot(length(audit) == 7L, all(vapply(audit, nrow, integer(1)) > 0L),
    all(audit$sem_parameters$supplied_pi_se == .02),
    all(grepl("external pi uncertainty is not included", audit$sem_parameters$uncertainty_scope)))
})
check("raw parameter estimates and covariance-parameter SEs match lavaan directly", {
  pe <- lavaan::parameterEstimates(lav, ci = TRUE)
  equal(audit$sem_parameters$est, pe$est)
  equal(audit$sem_parameters$se, pe$se)
  equal(audit$sem_covariances$se, pe$se[pe$op == "~~"])
  stopifnot(all(audit$sem_parameters$status == "ok"),
    all(audit$sem_parameters$se[audit$sem_parameters$is_fixed] == 0))
})
check("all fit measures and the complete lower-triangle parameter sampling covariance are retained", {
  fm <- suppressWarnings(lavaan::fitMeasures(lav))
  found <- audit$sem_fit_statistics
  equal(found$value[match(names(fm), found$statistic)], as.numeric(fm))
  stopifnot(all(found$status[match(names(fm)[!is.finite(fm)], found$statistic)] == "unavailable"))
  V <- lavaan::vcov(lav)
  equal(audit$sem_parameter_vcov$covariance, V[lower.tri(V, diag = TRUE)])
})

# Derivatives below are independent numerical perturbations of the model's
# forward covariance calculation, not lavaan's analytic delta matrix.
finite_difference_check <- function(f, label) {
  lv <- f$fit
  theta <- lavaan::coef(lv)
  forward <- function(x) {
    model <- getFromNamespace("lav_model_set_parameters", "lavaan")(lavmodel = lv@Model, x = x)
    getFromNamespace("lav_model_implied", "lavaan")(lavmodel = model)$cov[[1L]]
  }
  S <- lavaan::lavInspect(lv, "cov.ov")
  equal(forward(theta), S, 1e-10)
  D <- matrix(0, nrow(S)^2, length(theta))
  for (j in seq_along(theta)) {
    h <- 1e-5 * max(1, abs(theta[j]))
    plus <- minus <- theta
    plus[j] <- plus[j] + h; minus[j] <- minus[j] - h
    D[, j] <- as.vector((forward(plus) - forward(minus)) / (2 * h))
  }
  V <- lavaan::vcov(lv)
  numerical_se <- matrix(sqrt(pmax(rowSums((D %*% V) * D), 0)), nrow(S),
    dimnames = dimnames(S))
  detail <- ravel_sem_audit(wrap(f))$sem_implied_covariances
  observed <- subset(detail, space == "observed" & inference_method == "joint_parameter_delta")
  idx <- cbind(match(observed$lhs, rownames(S)), match(observed$rhs, colnames(S)))
  discrepancy <- max(abs(observed$se - numerical_se[idx]))
  check(paste(label, "implied observed covariance SEs agree with independent finite differences"), {
    stopifnot(nrow(observed) > 0L, discrepancy < 1e-7)
  })
  cat(label, "maximum implied covariance SE difference:", format(discrepancy, digits = 8), "\n")
  detail
}
continuous_detail <- finite_difference_check(continuous, "continuous")
check("exogenous latent covariance SEs match directly estimated covariance parameters", {
  latent <- subset(continuous_detail, space == "latent")
  pe <- lavaan::parameterEstimates(lav)
  lookup <- paste(pe$lhs, pe$rhs)
  index <- match(paste(latent$lhs, latent$rhs), lookup)
  missing <- is.na(index)
  index[missing] <- match(paste(latent$rhs[missing], latent$lhs[missing]), lookup)
  stopifnot(all(!is.na(index)), all(latent$status == "ok"))
  equal(latent$se, pe$se[index])
})
check("nonconverged sources cannot appear successful in any inferential output table", {
  failed <- continuous; failed$converged <- FALSE
  bad <- ravel_sem_audit(wrap(failed))
  for (nm in setdiff(names(bad), "sem_extraction_status"))
    stopifnot(all(bad[[nm]]$status == "not_converged"))
})
check("converged Heywood sources are explicitly labelled inadmissible", {
  S <- matrix(c(1, .7, .7, .7, 1, .1, .7, .1, 1), 3L,
    dimnames = list(c("x1", "x2", "x3"), c("x1", "x2", "x3")))
  heywood <- suppressWarnings(lavaan::cfa("g =~ x1 + x2 + x3", sample.cov = S,
    sample.nobs = 1000L, std.lv = TRUE))
  stopifnot(lavaan::lavInspect(heywood, "converged"),
    !suppressWarnings(lavaan::lavInspect(heywood, "post.check")))
  failed <- continuous; failed$fit <- heywood; failed$converged <- TRUE
  bad <- ravel_sem_audit(wrap(failed))
  for (nm in setdiff(names(bad), "sem_extraction_status"))
    stopifnot(all(bad[[nm]]$status == "inadmissible"))
})
check("missing source fits produce an explicit extraction status", {
  unavailable <- wrap(continuous)
  attr(unavailable, "fits")[[1]]$fit <- NULL
  bad <- ravel_sem_audit(unavailable)
  stopifnot(identical(bad$sem_extraction_status$status, "unavailable"))
})
check("standardization warnings are visible and cannot taint raw parameter SEs", {
  stage <- subset(audit$sem_extraction_status, stage == "standardized_parameters")
  stopifnot(nrow(stage) == 1L,
    all(audit$sem_standardized_parameters$status == stage$status),
    all(audit$sem_parameters$status == "ok"))
})

binary <- ravel_fit(dat, col_outcome = "Y_binary", outcome = "binary_liability",
  fit_constraints = "estimated_state", pi = .4)
binary_detail <- finite_difference_check(binary, "binary")
check("binary outcome variance absent from estimated moments has no invented SE", {
  yvar <- subset(binary_detail, space == "observed" & lhs == "Y" & rhs == "Y")
  stopifnot(nrow(yvar) == 1L, yvar$estimate == 1, is.na(yvar$se),
    yvar$inference_method == "unavailable", yvar$status == "se_unavailable")
})
check("exporting the audit preserves each standalone CSV table", {
  paths <- write_ravel_sem_audit(audit, tempfile("sem-audit-"))
  stopifnot(all(file.exists(paths)))
  for (nm in names(paths)) stopifnot(nrow(utils::read.csv(paths[[nm]])) == nrow(audit[[nm]]))
  unlink(paths)
})
utils::write.csv(do.call(rbind, checks), "dev/sem-audit-check/check-summary.csv", row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), "", "Source hashes:",
  capture.output(tools::md5sum(c("analysis/parental_deviation/ravel_sem_audit.R",
    "dev/run-sem-audit-checks.R")))), "dev/sem-audit-check/session.txt")
cat(length(checks), "SEM audit check groups passed.\n")
