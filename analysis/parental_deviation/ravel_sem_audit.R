# Detailed source-SEM tables for the parental-phenotype workflow.
# Sourcing defines functions only. These tables do not change the primary results.

ravel_sem_audit <- function(result, level = .95) {
  records <- attr(result, "fits")
  if (!is.list(records) || !length(records))
    stop("result must retain the parental workflow's source fit records.")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be between zero and one.")
  tables <- list(sem_fit_statistics = list(), sem_parameters = list(),
    sem_covariances = list(), sem_standardized_parameters = list(),
    sem_implied_covariances = list(), sem_parameter_vcov = list(),
    sem_extraction_status = list())
  cutoff <- stats::qnorm((1 + level) / 2)
  for (key in names(records)) {
    record <- records[[key]]
    fit <- record$fit
    meta <- data.frame(source_fit = key,
      outcome = if (is.null(record$outcome)) attr(result, "settings")$outcome else record$outcome,
      siblings = if (is.null(record$parental_sibs)) NA_character_ else record$parental_sibs,
      model = if (is.null(record$model)) NA_character_ else record$model,
      assumption = record$kind, pi = record$pi,
      supplied_pi_se = if (is.null(record$pi_se)) 0 else record$pi_se,
      genetic_members = paste(record$genetic_roles, collapse = "+"),
      uncertainty_scope = "Source SEM conditional on supplied pi; external pi uncertainty is not included",
      stringsAsFactors = FALSE)
    add <- function(name, tab) {
      if (is.null(tab) || !nrow(tab)) return(invisible(NULL))
      tables[[name]][[length(tables[[name]]) + 1L]] <<-
        cbind(meta[rep(1L, nrow(tab)), , drop = FALSE], tab, row.names = NULL)
    }
    capture <- function(stage, expression) {
      warnings <- character()
      value <- tryCatch(withCallingHandlers(force(expression), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = identity)
      failed <- inherits(value, "error")
      status <- if (failed) "error" else if (length(warnings)) "warning" else "ok"
      note <- paste(unique(c(if (failed) conditionMessage(value), warnings)), collapse = "; ")
      add("sem_extraction_status", data.frame(stage = stage, status = status, note = note))
      list(value = if (failed) NULL else value, status = status, note = note)
    }
    if (is.null(fit) || !inherits(fit$fit, "lavaan")) {
      add("sem_extraction_status", data.frame(stage = "source_fit", status = "unavailable",
        note = if (is.null(record$note)) "No fitted lavaan object." else record$note))
      next
    }
    lav <- fit$fit
    source_status <- if (isTRUE(fit$converged)) "ok" else "not_converged"
    n <- capture("sample_size", sum(lavaan::lavInspect(lav, "nobs")))
    meta$N <- if (is.null(n$value)) NA_real_ else n$value
    check <- capture("admissibility", lavaan::lavInspect(lav, "post.check"))
    if (source_status == "ok" && !isTRUE(check$value))
      source_status <- if (is.null(check$value)) "admissibility_unavailable" else "inadmissible"
    measures <- capture("fit_statistics", lavaan::fitMeasures(lav))
    if (length(measures$value)) add("sem_fit_statistics", data.frame(
      statistic = names(measures$value), value = unname(measures$value),
      status = if (source_status != "ok") source_status else
        ifelse(is.finite(measures$value), measures$status, "unavailable"),
      note = measures$note, stringsAsFactors = FALSE))
    add("sem_fit_statistics", data.frame(statistic = c("source_converged", "source_admissible", "N"),
      value = c(as.numeric(isTRUE(fit$converged)),
        if (is.null(check$value)) NA_real_ else as.numeric(check$value), meta$N),
      status = source_status, note = "", stringsAsFactors = FALSE))

    # Unstandardized entries use the actual estimator's SEs (MLR/cluster robust
    # where requested), with fixed values explicitly identified by free index 0.
    pars <- capture("parameters", lavaan::parameterEstimates(lav, ci = TRUE, level = level))
    pt <- lavaan::parameterTable(lav)
    parameter_key <- function(x) paste(x$lhs, x$op, x$rhs, sep = "\r")
    if (!is.null(pars$value)) {
      p <- as.data.frame(pars$value)
      index <- match(parameter_key(p), parameter_key(pt))
      p$parameter_index <- index
      p$free_index <- pt$free[index]
      p$is_fixed <- p$free_index == 0L & !p$op %in% c(":=", "==", ">", "<")
      p$status <- if (source_status != "ok") source_status else pars$status
      p$note <- pars$note
      add("sem_parameters", p)
      covariance <- p[p$op == "~~", , drop = FALSE]
      latent <- lavaan::lavNames(lav, "lv")
      covariance$component <- ifelse(covariance$lhs %in% latent & covariance$rhs %in% latent,
        "latent", "observed residual/exogenous")
      covariance$quantity <- ifelse(covariance$lhs == covariance$rhs, "variance", "covariance")
      add("sem_covariances", covariance)
    } else p <- NULL

    standardized <- capture("standardized_parameters", lavaan::standardizedSolution(
      lav, type = "std.all", se = TRUE, ci = TRUE, level = level))
    if (!is.null(standardized$value)) {
      tab <- as.data.frame(standardized$value)
      tab$standardization <- "std.all"
      tab$status <- if (source_status != "ok") source_status else standardized$status
      tab$note <- standardized$note
      add("sem_standardized_parameters", tab)
    }

    variance <- capture("parameter_vcov", lavaan::vcov(lav))
    V <- variance$value
    if (is.matrix(V)) {
      ids <- which(lower.tri(V, diag = TRUE), arr.ind = TRUE)
      labels <- make.unique(names(lavaan::coef(lav)))
      add("sem_parameter_vcov", data.frame(parameter_i = ids[, 1], parameter_j = ids[, 2],
        label_i = labels[ids[, 1]], label_j = labels[ids[, 2]],
        covariance = V[ids], status = if (source_status != "ok") source_status else variance$status,
        note = "Sampling covariance of SEM parameter estimates, not a genetic covariance; lower triangle including diagonal",
        stringsAsFactors = FALSE))
    }

    # Delta rows are matched by their moment names, never an assumed vech order.
    # Validate them against lavaan's expected moments before reporting any SE.
    delta <- capture("implied_moment_delta", {
      D <- lavaan::lavInspect(lav, "delta")
      moments <- lavaan::lavInspect(lav, "wls.est")
      if (!is.matrix(D) || !is.matrix(V) || !is.numeric(moments) ||
          !identical(rownames(D), names(moments)) || ncol(D) != nrow(V) ||
          !identical(colnames(D), colnames(V)))
        stop("Moment names or free-parameter ordering do not match; implied-moment SEs are unavailable.")
      vv <- rowSums((D %*% V) * D)
      list(moments = moments, se = stats::setNames(ifelse(is.finite(vv) & vv >= -1e-10,
        sqrt(pmax(vv, 0)), NA_real_), names(moments)))
    })
    for (space in c("observed", "latent")) {
      implied <- capture(paste0("implied_", space, "_covariance"),
        lavaan::lavInspect(lav, if (space == "observed") "cov.ov" else "cov.lv"))
      S <- implied$value
      if (!is.matrix(S)) next
      ids <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE)
      lhs <- rownames(S)[ids[, 1]]; rhs <- colnames(S)[ids[, 2]]
      estimate <- unname(S[ids])
      se <- rep(NA_real_, length(estimate))
      method <- rep("unavailable", length(estimate))
      note <- rep("No validated derivative or direct covariance parameter matched this moment.", length(estimate))
      if (space == "observed" && !is.null(delta$value)) {
        d <- delta$value
        keys <- paste0(lhs, "~~", rhs)
        opposite <- paste0(rhs, "~~", lhs)
        keys[!keys %in% names(d$moments)] <- opposite[!keys %in% names(d$moments)]
        matched <- match(keys, names(d$moments))
        valid <- !is.na(matched) & abs(unname(d$moments[matched]) - estimate) <=
          1e-7 * pmax(1, abs(estimate))
        valid[is.na(valid)] <- FALSE
        se[valid] <- d$se[matched[valid]]
        method[valid] <- "joint_parameter_delta"
        note[valid] <- ""
      }
      if (space == "latent" && !is.null(p)) {
        pc <- p[p$op == "~~", , drop = FALSE]
        keys <- paste(lhs, rhs, sep = "\r")
        pos <- match(keys, paste(pc$lhs, pc$rhs, sep = "\r"))
        missing <- is.na(pos)
        pos[missing] <- match(keys[missing], paste(pc$rhs, pc$lhs, sep = "\r"))
        valid <- !is.na(pos) & abs(pc$est[pos] - estimate) <= 1e-7 * pmax(1, abs(estimate))
        valid[is.na(valid)] <- FALSE
        se[valid] <- pc$se[pos[valid]]
        method[valid] <- ifelse(pc$is_fixed[pos[valid]], "fixed_parameter", "direct_covariance_parameter")
        note[valid] <- ""
      }
      status <- ifelse(is.finite(se), "ok", "se_unavailable")
      dependencies <- if (space == "observed") list(variance, delta, implied, check) else
        list(pars, implied, check)
      warnings <- vapply(dependencies, function(x) x$note, "")
      warnings <- unique(warnings[nzchar(warnings)])
      if (any(vapply(dependencies, function(x) x$status != "ok", logical(1)))) {
        status[status == "ok"] <- "warning"
        note <- vapply(note, function(x) paste(c(x[nzchar(x)], warnings), collapse = "; "), "")
      }
      if (source_status != "ok") status[] <- source_status
      add("sem_implied_covariances", data.frame(space = space, lhs = lhs, rhs = rhs,
        quantity = ifelse(lhs == rhs, "variance", "covariance"), estimate = estimate,
        se = se, ci_low = estimate - cutoff * se, ci_high = estimate + cutoff * se,
        inference_method = method, status = status, note = note, stringsAsFactors = FALSE))
    }
  }
  # Different stages and estimators can supply extra columns. Preserve them all.
  bind <- function(parts) {
    if (!length(parts)) return(data.frame())
    columns <- unique(unlist(lapply(parts, names)))
    parts <- lapply(parts, function(x) {
      for (nm in setdiff(columns, names(x))) x[[nm]] <- NA
      x[columns]
    })
    out <- do.call(rbind, parts); rownames(out) <- NULL; out
  }
  lapply(tables, bind)
}

write_ravel_sem_audit <- function(audit, prefix) {
  paths <- stats::setNames(paste0(prefix, "_", names(audit), ".csv"), names(audit))
  for (nm in names(audit)) utils::write.csv(audit[[nm]], paths[[nm]], row.names = FALSE)
  invisible(paths)
}
