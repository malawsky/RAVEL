# Validate observed-PGS and auxiliary-SEM reporting using retained synthetic fits.
# First run dev/check-paper-reporting.R --prepare-only if the synthetic caches are absent.
pkgload::load_all(".", quiet = TRUE)
RAVEL:::.ensure_lavaan_runtime()
options(ravel.workflow_functions_only = TRUE)
e <- new.env(parent = globalenv())
sys.source("analysis/parental_deviation/ravel_parental_deviation.R", e)
checks <- list()
check <- function(ok, label) {
  checks[[length(checks) + 1L]] <<- data.frame(check = label, passed = isTRUE(ok))
  if (!isTRUE(ok)) stop(label)
}
equal <- function(a, b, tol = 1e-7) isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))
for (name in c("child-sem", "parent-sem", "child-mixed", "child-mixed-no-covariates")) {
  result <- readRDS(file.path("/tmp/ravel-paper-reporting-cache", paste0(name, ".rds")))
  parent <- name == "parent-sem"
  records <- attr(result, if (parent) "fits" else "report_fit_records")
  for (key in names(records)) {
    fit <- records[[key]]$fit
    mixed <- inherits(fit, "ravel_mixed_results")
    settings <- records[[key]]$settings
    rows <- e$.ravel_report_pgs_structure(fit, settings)
    cov <- rows[rows$block == "pgs_input_covariance", ]
    cor <- rows[rows$block == "pgs_correlation", ]
    p <- length(unique(c(cov$lhs, cov$rhs)))
    check(nrow(cov) == p*(p+1)/2 && nrow(cor) == p*(p-1)/2, paste(name,key,"all_PGS_pairs"))
    check(all(is.finite(rows$estimate)) && all(is.finite(rows$se)), paste(name,key,"finite_estimates_SEs"))
    if (!mixed) {
      native <- lavaan::lavInspect(fit$fit,"cov.ov")
      check(equal(cov$estimate,native[cbind(cov$lhs,cov$rhs)]),paste(name,key,"native_covariances"))
      check(equal(cor$estimate,cov2cor(native)[cbind(cor$lhs,cor$rhs)]),paste(name,key,"native_correlations"))
      score_names <- unique(c(cov$lhs,cov$rhs))
      scales <- seq_along(score_names)*1e-5
      settings$score_scaling <- data.frame(canonical=score_names,column=paste0("custom",seq_along(scales)),scale=scales,center=0)
      raw <- e$.ravel_report_pgs_structure(fit,settings)
      factor <- scales[match(rows$lhs,score_names)]*scales[match(rows$rhs,score_names)]
      factor[rows$block=="pgs_correlation"] <- 1
      check(equal(raw$estimate,rows$estimate*factor,1e-10) && equal(raw$se,rows$se*factor,1e-10),paste(name,key,"raw_scale_and_SE_transformation"))
    } else {
      g <- RAVEL:::.ravel_mixed_genetic(2L,fit$parameters,fit$settings)
      labs <- c("mother","father","child","child_sib")
      labels <- unlist(lapply(1:2,function(k) paste0("PGS",k,"_",labs)))
      dimnames(g$P) <- list(labels,labels)
      check(equal(cor$estimate,cov2cor(g$P)[cbind(cor$lhs,cor$rhs)]),paste(name,key,"mixed_correlations"))
      pair <- c("PGS1_mother","PGS2_child_sib")
      f <- function(par) cov2cor(RAVEL:::.ravel_mixed_genetic(2L,par,fit$settings)$P)[1,8]
      grad <- numDeriv::grad(f,fit$parameters)
      expected <- sqrt(drop(grad %*% fit$vcov[names(fit$parameters), names(fit$parameters)] %*% grad))
      row <- cor[(cor$lhs==pair[1]&cor$rhs==pair[2]) | (cor$lhs==pair[2]&cor$rhs==pair[1]),]
      check(nrow(row)==1 && equal(row$se,expected,1e-6),paste(name,key,"independent_numerical_correlation_SE"))
    }
  }
  report <- e$ravel_reporting_tables(result,if(parent) "parental" else "child")
  aux <- report$fit_statistics[report$fit_statistics$engine=="lavaan_auxiliary",]
  mains <- records[vapply(records,function(z) identical(z$kind,"main"),logical(1))]
  check(length(unique(aux$source_fit))==length(mains),paste(name,"one_auxiliary_per_main"))
  for (key in names(mains)) {
    diagnostic <- mains[[key]]$diagnostic
    if (is.null(diagnostic)) diagnostic <- mains[[key]]$fit$one_factor
    lav <- diagnostic$auxiliary_fit$fit
    if (!inherits(lav,"lavaan")) next
    expected <- lavaan::fitMeasures(lav)
    rows <- aux[aux$source_fit==paste0(key,"|auxiliary_SEM"),]
    for (metric in c("aic","bic","rmsea","cfi","tli","srmr","chisq","df"))
      check(equal(rows$value[tolower(rows$statistic)==metric],unname(expected[metric])),paste(name,key,"auxiliary",metric))
    check(equal(rows$value[rows$statistic=="N"],sum(lavaan::lavInspect(lav,"nobs"))),paste(name,key,"auxiliary_N"))
  }
  check(!any(grepl("not a check|does not (test|validate|assess)|not a test",report$diagnostics$note,ignore.case=TRUE)),paste(name,"concise_diagnostic_notes"))
}
dir.create("dev/pgs-reporting-validation",showWarnings=FALSE)
write.csv(do.call(rbind,checks),"dev/pgs-reporting-validation/checks.csv",row.names=FALSE)
cat(length(checks),"checks passed\n")
