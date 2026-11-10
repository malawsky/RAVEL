score_numerics_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    set.seed(9102)
    n <- 1000L
    G <- MASS::mvrnorm(n, rep(0, 3), matrix(c(1,.65,.65,.65,1,.3,.65,.3,1),3))
    K <- matrix(c(1,.5,.5,.5,1,0,.5,0,1),3)
    d <- data.frame(FID=seq_len(n), IID=paste0('i',seq_len(n)))
    for (j in 1:2) d[paste0('PGS',j,'_',c('child','mother','father'))] <-
      G + MASS::mvrnorm(n,rep(0,3),c(.4,1.3)[j]*K)
    d$Y <- drop(G %*% c(.32,.13,.08)) + rnorm(n)
    f <- ravel_fit(d,fit_constraints='equilibrium',pi=.4)
    cached <<- list(data=d,fit=f)
    cached
  }
})

test_that('converged starts are installed before any nonlinear SEM evaluation', {
  x <- score_numerics_fixture()
  model <- .ravel_score_signal_model(x$fit)
  pt <- lavaan::parameterTable(x$fit$fit)
  labels <- pt$label[pt$free>0 & nzchar(pt$label)]
  # Deliberately invalid explicit starts exercise the constructor failure window.
  # The converged reference must replace them before evaluating constraints.
  for (label in labels) model <- gsub(paste0('\\b',label,'\\*(?=[A-Z])'),
    paste0('start(0)*',label,'*'),model,perl=TRUE)
  dat <- as.data.frame(lavaan::lavInspect(x$fit$fit,'data'))
  before <- serialize(x$fit,NULL)
  ans <- .ravel_score_signal_capture(.ravel_score_signal_fit(model,dat,x$fit$fit,'MLR'))
  expect_null(ans$error)
  expect_true(.ravel_score_signal_validate_fit(ans$value)$valid)
  expect_equal(lavaan::fitMeasures(ans$value,'logl'),
    lavaan::fitMeasures(x$fit$fit,'logl'),tolerance=1e-6)
  expect_identical(serialize(x$fit,NULL),before)
})

test_that('a numerical refit error gets an audited same-model restart', {
  x <- score_numerics_fixture()
  real_fit <- .ravel_score_signal_fit
  calls <- 0L
  testthat::local_mocked_bindings(.ravel_score_signal_fit=function(...) {
    calls <<- calls+1L
    if (calls==1L) stop('NA/NaN/Inf in foreign function call (arg 1).')
    real_fit(...)
  },.package='RAVEL')
  z <- ravel_test_score_signal(x$fit)
  expect_true(z$status %in% c('ok','ok_with_warnings'),info=z$reason)
  expect_true(is.finite(z$p_value))
  retry <- z$diagnostics$refit_retries$H0
  if (!is.null(retry$numerical_error_restart)) retry <- retry$numerical_error_restart
  expect_true(retry$resolved)
  expect_match(retry$initial_error,'NA/NaN/Inf',fixed=TRUE)
  expect_true(retry$retry_validation$valid)
  expect_identical(retry$method,'same_model_data_derived_starts')
})

test_that('valid diagnostic p-values survive numerical warning metadata', {
  x <- score_numerics_fixture()
  real_checked <- .ravel_score_signal_fit_checked
  testthat::local_mocked_bindings(.ravel_score_signal_fit_checked=function(...) {
    ans <- real_checked(...)
    ans$warnings <- c(ans$warnings,'Synthetic warning for reporting regression test')
    ans
  },.package='RAVEL')
  z <- ravel_test_score_signal(x$fit)
  expect_identical(z$status,'ok_with_warnings',info=z$reason)
  expect_true(is.finite(z$statistic))
  expect_true(is.finite(z$p_value))
  expect_match(z$diagnostics$warnings$H0,'Synthetic warning')
  expect_output(print(z),'Statistic:')
})

test_that('auxiliary SEM survives a failed comparison and retains fit statistics', {
  x <- score_numerics_fixture()
  testthat::local_mocked_bindings(ravel_test_score_signal=function(...) {
    stop('Synthetic null-refit failure')
  },.package='RAVEL')
  z <- ravel_score_diagnostic(x$data,fit_constraints='equilibrium',pi=.4)
  expect_identical(z$status,'diagnostic_error')
  expect_s4_class(z$auxiliary_fit$fit,'lavaan')
  fm <- lavaan::fitMeasures(z$auxiliary_fit$fit,c('aic','bic','chisq','df','cfi','tli','rmsea','srmr'))
  expect_true(all(is.finite(fm)))
  expect_equal(z$selection$N_used,nrow(x$data))
})

test_that('degenerate selected variables fail explicitly without replacement', {
  d <- score_numerics_fixture()$data
  d$PGS1_child <- 1
  z <- ravel_score_diagnostic(d,fit_constraints='equilibrium',pi=.4)
  expect_identical(z$status,'invalid_data')
  expect_match(z$reason,'PGS1_child',fixed=TRUE)
  expect_equal(sort(z$selection$rows_used),seq_len(nrow(d)))
  d$PGS1_child <- rep(c(-1,1),length.out=nrow(d))*1e155
  z <- ravel_score_diagnostic(d,fit_constraints='equilibrium',pi=.4)
  expect_identical(z$status,'invalid_data')
  expect_match(z$reason,'finite positive variance',fixed=TRUE)
})

test_that('numerical restart cannot bypass model validation', {
  x <- score_numerics_fixture()
  real_fit <- .ravel_score_signal_fit
  calls <- 0L
  testthat::local_mocked_bindings(.ravel_score_signal_fit=function(...) {
    calls <<- calls+1L
    if (calls==1L) stop('NA/NaN/Inf in foreign function call (arg 1).')
    real_fit(...)
  },.ravel_score_signal_validate_fit=function(...) list(valid=FALSE,reason='Synthetic validation failure'),
    .package='RAVEL')
  z <- ravel_test_score_signal(x$fit)
  expect_identical(z$status,'fit_error')
  expect_true(is.na(z$p_value))
  expect_false(z$diagnostics$refit_retries$H0$resolved)
  expect_match(z$reason,'Synthetic validation failure',fixed=TRUE)
})
