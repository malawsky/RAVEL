cov_test_cache <- new.env(parent = emptyenv())

cov_test_fixed_corr <- c(
  e1_mf = 0, e1_mc = 0.5, e1_fc = 0.5,
  e2_mf = 0, e2_mc = 0.5, e2_fc = 0.5
)

cov_test_fixture <- function(outcome = "binary_liability", correlated = TRUE,
                             with_covariates = TRUE, n = 8000L, seed = 71631L) {
  id <- paste(outcome, correlated, with_covariates, n, seed, sep = "_")
  key <- paste0("data_", id)
  if (exists(key, cov_test_cache, inherits = FALSE)) return(cov_test_cache[[key]])

  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit(if (had_seed) {
    assign(".Random.seed", old_seed, .GlobalEnv)
  } else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)

  roles <- c("child", "mother", "father")
  sigma_g <- matrix(c(1, .6, .6, .6, 1, .2, .6, .2, 1), 3, 3)
  dimnames(sigma_g) <- list(roles, roles)
  g <- matrix(rnorm(n * 3), n, 3) %*% chol(sigma_g)
  pc <- if (correlated) .5 * g[, 2] + .2 * g[, 3] + sqrt(.67) * rnorm(n) else rnorm(n)
  age <- .35 * pc + sqrt(1 - .35^2) * rnorm(n)
  sex <- rbinom(n, 1, .5)
  w <- cbind(PC1 = pc, age = age, sex = sex - .5)
  sigma_w <- matrix(c(1, .35, 0, .35, 1, 0, 0, 0, .25), 3, 3)
  cov_g_pc <- if (correlated) as.numeric(sigma_g %*% c(0, .5, .2)) else rep(0, 3)
  sigma_gw <- cbind(cov_g_pc, .35 * cov_g_pc, 0)
  b <- c(bD = .45, bm = .12, bf = 0)
  a <- if (with_covariates) c(.30, -.15, .20) else rep(0, 3)
  var_y <- as.numeric(t(b) %*% sigma_g %*% b + t(a) %*% sigma_w %*% a +
                        2 * t(b) %*% sigma_gw %*% a + .9^2)
  liability <- as.numeric(g %*% b + w %*% a + .9 * rnorm(n))
  y <- if (outcome == "continuous") liability else as.numeric(liability > .8 * sqrt(var_y))
  dat <- data.frame(FID = seq_len(n), Y = y, PC1 = pc, age = age, sex = sex)
  sigma_error <- matrix(c(1, .5, .5, .5, 1, 0, .5, 0, 1), 3, 3)
  reliability <- rbind(c(.70, .65, .60), c(.60, .75, .70))
  for (k in 1:2) {
    e <- matrix(rnorm(n * 3), n, 3) %*% chol(sigma_error)
    for (r in seq_along(roles)) {
      dat[[paste0("PGS", k, "_", roles[r])]] <-
        sqrt(reliability[k, r]) * g[, r] + sqrt(1 - reliability[k, r]) * e[, r]
    }
  }
  out <- list(id = id, dat = dat, outcome = outcome,
              covariates = if (with_covariates) c("PC1", "age", "sex") else character(),
              beta = if (outcome == "continuous") b else b / sqrt(var_y),
              var_y = var_y, sigma_g = sigma_g)
  cov_test_cache[[key]] <- out
  out
}

# This reference specifies the SEM independently of the package model builders.
cov_test_reference_model <- function(covariates = character(),
                                     crossperson = "free", fixed_corr = TRUE,
                                     fix_parent_paths = FALSE) {
  model <- "
    Gm =~ lm1*PGS1_mother + lm2*PGS2_mother
    Gf =~ lf1*PGS1_father + lf2*PGS2_father
    Gc =~ lc1*PGS1_child + lc2*PGS2_child
    PGS1_mother ~~ th1_m*PGS1_mother
    PGS1_father ~~ th1_f*PGS1_father
    PGS1_child ~~ th1_c*PGS1_child
    PGS2_mother ~~ th2_m*PGS2_mother
    PGS2_father ~~ th2_f*PGS2_father
    PGS2_child ~~ th2_c*PGS2_child
    Gm ~~ rhoMF*Gf
    Gc ~~ eta_m*Gm
    Gc ~~ eta_f*Gf
  "
  for (k in 1:2) {
    labels <- if (crossperson == "free") c("mf", "mc", "fc") else c("mf", "pc", "pc")
    model <- paste(model,
      sprintf("PGS%d_mother ~~ e%d_%s*PGS%d_father", k, k, labels[1], k),
      sprintf("PGS%d_mother ~~ e%d_%s*PGS%d_child", k, k, labels[2], k),
      sprintf("PGS%d_father ~~ e%d_%s*PGS%d_child", k, k, labels[3], k), sep = "\n")
    if (fixed_corr) {
      model <- paste(model,
        sprintf("e%d_mf == 0*sqrt(th%d_m*th%d_f)", k, k, k),
        sprintf("e%d_mc == 0.5*sqrt(th%d_m*th%d_c)", k, k, k),
        sprintf("e%d_fc == 0.5*sqrt(th%d_f*th%d_c)", k, k, k), sep = "\n")
    }
  }
  model <- paste(model, if (fix_parent_paths) "Y ~ bD*Gc + 0*Gm + 0*Gf" else
                   "Y ~ bD*Gc + bm*Gm + bf*Gf", sep = "\n")
  for (i in seq_along(covariates)) {
    model <- paste(model, sprintf("Y ~ c_%d*%s", i, covariates[i]), sep = "\n")
  }
  for (g in c("Gc", "Gm", "Gf")) {
    for (w in covariates) model <- paste(model, sprintf("%s ~~ %s", g, w), sep = "\n")
  }
  if (length(covariates) > 1L) {
    pairs <- combn(covariates, 2)
    for (j in seq_len(ncol(pairs))) {
      model <- paste(model, paste(pairs[, j], collapse = " ~~ "), sep = "\n")
    }
  }
  model
}

cov_test_fd_jacobian <- function(fun, x, step = 1e-5) {
  value <- fun(x)
  J <- matrix(0, length(value), length(x), dimnames = list(names(value), names(x)))
  for (j in seq_along(x)) {
    h <- step * max(1, abs(x[j]))
    plus <- minus <- x
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    J[, j] <- (fun(plus) - fun(minus)) / (2 * h)
  }
  J
}

# Independent first stage: raw-X probit and a numerical delta method on primitives.
cov_test_binary_reference <- function(data, model, step = 1e-5) {
  variables <- lavaan::lavNames(lavaan::lavaanify(model), "ov")
  predictors <- setdiff(variables, "Y")
  data <- data[stats::complete.cases(data[variables]), variables, drop = FALSE]
  X <- as.matrix(data[predictors])
  n <- nrow(X)
  p <- ncol(X)
  mu <- colMeans(X)
  centered <- sweep(X, 2L, mu)
  S <- crossprod(centered) / n
  design <- cbind("(Intercept)" = 1, X)
  probit <- stats::glm.fit(design, data$Y, family = stats::binomial("probit"))
  stopifnot(probit$converged, probit$rank == p + 1L,
            all(is.finite(probit$coefficients)))
  probability <- stats::pnorm(drop(design %*% probit$coefficients))
  density <- stats::dnorm(drop(design %*% probit$coefficients))
  score <- design * ((data$Y - probability) * density / (probability * (1 - probability)))
  information <- crossprod(design * (density / sqrt(probability * (1 - probability)))) / n
  IFbeta <- t(solve(information, t(score)))
  pairs <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE)
  covariance_names <- apply(pairs, 1L, function(pair) {
    paste0("cov:", predictors[pair[1]], ":", predictors[pair[2]])
  })
  primitive <- c(setNames(probit$coefficients, paste0("beta:", colnames(design))),
                 setNames(mu, paste0("mean:", predictors)),
                 setNames(S[pairs], covariance_names))
  IFcovariance <- vapply(seq_len(nrow(pairs)), function(j) {
    centered[, pairs[j, 1]] * centered[, pairs[j, 2]] - S[pairs[j, 1], pairs[j, 2]]
  }, numeric(n))
  primitive_if <- cbind(IFbeta, centered, IFcovariance)
  colnames(primitive_if) <- names(primitive)
  primitive_if <- sweep(primitive_if, 2L, colMeans(primitive_if))

  location_names <- vapply(variables, function(v) {
    if (v == "Y") "Y|t1" else paste0(v, "~1")
  }, character(1), USE.NAMES = FALSE)
  labels <- c(location_names, paste0(predictors, "~~", predictors))
  for (j in seq_len(length(variables) - 1L)) {
    for (i in seq.int(j + 1L, length(variables))) {
      labels <- c(labels, paste0(variables[j], "~~", variables[i]))
    }
  }
  decode <- function(value) {
    covariance <- matrix(0, p, p, dimnames = list(predictors, predictors))
    covariance[pairs] <- value[seq.int(2L * p + 2L, length(value))]
    covariance[upper.tri(covariance)] <- t(covariance)[upper.tri(covariance)]
    list(beta = value[seq_len(p + 1L)],
         mean = value[seq.int(p + 2L, 2L * p + 1L)], covariance = covariance)
  }
  map <- function(value) {
    primitive_value <- decode(value)
    delta <- primitive_value$beta[-1L]
    mean <- primitive_value$mean
    covariance <- primitive_value$covariance
    variance <- 1 + drop(crossprod(delta, covariance %*% delta))
    covariance_y <- drop(covariance %*% delta) / sqrt(variance)
    threshold <- -(primitive_value$beta[1L] + sum(mean * delta)) / sqrt(variance)
    out <- setNames(numeric(length(labels)), labels)
    for (j in seq_along(labels)) {
      label <- labels[j]
      if (label == "Y|t1") {
        out[j] <- threshold
      } else if (endsWith(label, "~1")) {
        out[j] <- mean[match(sub("~1$", "", label), predictors)]
      } else {
        pair <- strsplit(label, "~~", fixed = TRUE)[[1L]]
        out[j] <- if ("Y" %in% pair) covariance_y[match(setdiff(pair, "Y"), predictors)] else
          covariance[pair[1L], pair[2L]]
      }
    }
    out
  }
  moments <- map(primitive)
  J <- cov_test_fd_jacobian(map, primitive, step)
  moment_if <- primitive_if %*% t(J)
  moment_if <- sweep(moment_if, 2L, colMeans(moment_if))
  covariance <- diag(length(variables))
  dimnames(covariance) <- list(variables, variables)
  means <- setNames(numeric(length(variables)), variables)
  means[predictors] <- mu
  for (label in labels[grepl("~~", labels, fixed = TRUE)]) {
    pair <- strsplit(label, "~~", fixed = TRUE)[[1L]]
    covariance[pair[1L], pair[2L]] <- covariance[pair[2L], pair[1L]] <- moments[label]
  }
  th <- -means
  th["Y"] <- moments["Y|t1"]
  th_index <- integer(length(variables))
  th_index[match("Y", variables)] <- match("Y", variables)
  attr(th, "th.idx") <- setNames(th_index, location_names)
  working_variance <- vapply(labels, function(label) {
    if (label == "Y|t1") {
      tau <- moments[label]
      return(stats::pnorm(tau) * stats::pnorm(-tau) / stats::dnorm(tau)^2)
    }
    if (endsWith(label, "~1")) {
      v <- sub("~1$", "", label)
      return(covariance[v, v])
    }
    pair <- strsplit(label, "~~", fixed = TRUE)[[1L]]
    covariance[pair[1L], pair[1L]] * covariance[pair[2L], pair[2L]] +
      covariance[pair[1L], pair[2L]]^2
  }, numeric(1), USE.NAMES = FALSE)
  list(args = list(sample.cov = covariance, sample.mean = means, sample.th = th,
                   sample.nobs = n, NACOV = crossprod(moment_if) / n,
                   WLS.V = diag(1 / working_variance), sample.cov.rescale = FALSE),
       labels = labels, moments = moments, primitive = primitive, primitive_if = primitive_if,
       primitive_gamma = crossprod(primitive_if) / n, jacobian = J, moment_if = moment_if,
       map = map, decode = decode, pairs = pairs, predictors = predictors,
       probit = probit, information = information)
}

# Analytic differentials are checked against the numerical primitive Jacobian.
cov_test_binary_analytic_jacobian <- function(reference) {
  value <- reference$decode(reference$primitive)
  p <- length(reference$predictors)
  delta <- unname(value$beta[-1L])
  mu <- unname(value$mean)
  S <- value$covariance
  a <- drop(S %*% delta)
  v <- 1 + sum(delta * a)
  k <- value$beta[1L] + sum(mu * delta)
  J <- matrix(0, length(reference$labels), length(reference$primitive),
              dimnames = list(reference$labels, names(reference$primitive)))
  for (j in seq_along(reference$primitive)) {
    direction <- numeric(length(reference$primitive))
    direction[j] <- 1
    change <- reference$decode(direction)
    db <- change$beta[-1L]
    dm <- change$mean
    dS <- change$covariance
    dv <- drop(crossprod(delta, dS %*% delta)) + 2 * sum(a * db)
    dc <- drop(dS %*% delta + S %*% db) / sqrt(v) - a * dv / (2 * v^(3 / 2))
    dt <- -(change$beta[1L] + sum(mu * db) + sum(delta * dm)) / sqrt(v) +
      k * dv / (2 * v^(3 / 2))
    for (i in seq_along(reference$labels)) {
      label <- reference$labels[i]
      if (label == "Y|t1") {
        J[i, j] <- dt
      } else if (endsWith(label, "~1")) {
        J[i, j] <- dm[match(sub("~1$", "", label), reference$predictors)]
      } else {
        pair <- strsplit(label, "~~", fixed = TRUE)[[1L]]
        J[i, j] <- if ("Y" %in% pair) dc[match(setdiff(pair, "Y"), reference$predictors)] else
          dS[pair[1L], pair[2L]]
      }
    }
  }
  J
}

cov_test_reference_fit <- function(fixture, crossperson = "free", fixed_corr = TRUE,
                                   fix_parent_paths = FALSE, legacy_options = FALSE,
                                   start = NULL) {
  args <- list(model = cov_test_reference_model(fixture$covariates, crossperson,
                                               fixed_corr, fix_parent_paths),
               data = fixture$dat, std.lv = TRUE, fixed.x = FALSE,
               estimator = if (fixture$outcome == "continuous") "MLR" else "WLSMV",
               missing = if (fixture$outcome == "continuous") "fiml" else "listwise",
               ordered = if (fixture$outcome == "continuous") NULL else "Y",
               check.gradient = FALSE, warn = FALSE)
  if (!legacy_options) args$conditional.x <- FALSE
  if (fixture$outcome == "binary_liability" && length(fixture$covariates)) {
    args$data <- NULL
    args <- c(args, cov_test_binary_reference(fixture$dat, args$model)$args)
  }
  if (!is.null(start)) args$start <- start
  do.call(lavaan::sem, args)
}

cov_test_fit <- function(fixture, engine = "package", crossperson = "free",
                         fixed_corr = TRUE, fix_parent_paths = FALSE,
                         n_boot = 0L, h = 0, regime = "n_gen") {
  key <- paste("fit", fixture$id, engine, crossperson, fixed_corr,
               fix_parent_paths, n_boot, h, regime, sep = "_")
  if (exists(key, cov_test_cache, inherits = FALSE)) return(cov_test_cache[[key]])
  warnings <- character()
  out <- withCallingHandlers({
    if (engine == "reference") {
      fit <- cov_test_reference_fit(fixture, crossperson, fixed_corr, fix_parent_paths,
                                    legacy_options = !length(fixture$covariates))
      list(converged = isTRUE(lavaan::lavInspect(fit, "converged")), fit = fit)
    } else {
      RAVEL::ravel_fit(fixture$dat, outcome = fixture$outcome,
        fixed_effects = fixture$covariates, crossperson = crossperson,
        fixed_crossperson_corr = if (fixed_corr) cov_test_fixed_corr else NULL,
        fix_parent_paths = fix_parent_paths, h = h, am_correction = regime,
        am_n_generations = 1L, n_boot = n_boot, seed = 99173L)
    }
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  attr(out, "test_warnings") <- warnings
  cov_test_cache[[key]] <- out
  out
}

cov_test_converged <- function(result) {
  expect_true(isTRUE(result$converged), info = paste(c(result$error,
              attr(result, "test_warnings")), collapse = "\n"))
  isTRUE(result$converged)
}

cov_test_expect_reference <- function(result, reference, tolerance = 5e-5) {
  if (!cov_test_converged(result) || !cov_test_converged(reference)) return(invisible(NULL))
  got <- lavaan::coef(result$fit)
  expected <- lavaan::coef(reference$fit)
  expect_gt(length(got), 0L)
  expect_gt(length(expected), 0L)
  expect_setequal(names(got), names(expected))
  shared <- intersect(names(got), names(expected))
  expect_gt(length(shared), 0L)
  expect_true(all(is.finite(got)))
  expect_true(all(is.finite(expected)))
  expect_equal(got[shared], expected[shared], tolerance = tolerance)
  got_vcov <- lavaan::vcov(result$fit)
  reference_vcov <- lavaan::vcov(reference$fit)
  expect_true(is.matrix(got_vcov))
  expect_true(is.matrix(reference_vcov))
  if (!is.matrix(got_vcov) || !is.matrix(reference_vcov)) return(invisible(NULL))
  expect_equal(dim(got_vcov), rep(length(got), 2L))
  expect_equal(dim(reference_vcov), rep(length(expected), 2L))
  expect_true(all(is.finite(got_vcov)))
  expect_true(all(is.finite(reference_vcov)))
  expect_equal(got_vcov[shared, shared], reference_vcov[shared, shared], tolerance = tolerance)
  expect_true(isTRUE(lavaan::lavInspect(result$fit, "post.check")))
}

cov_test_genetic_fraction <- function(fit) {
  b <- lavaan::coef(fit)[c("bD", "bm", "bf")]
  sigma <- lavaan::lavInspect(fit, "cov.lv")[c("Gc", "Gm", "Gf"), c("Gc", "Gm", "Gf")]
  as.numeric(t(b) %*% sigma %*% b / lavaan::lavInspect(fit, "cov.ov")["Y", "Y"])
}
