# Independent scientific checks: the formulas below are derived from the
# rank-one assortment covariance recurrence, rather than the new fit builder.
adversarial_component_state <- function(rho_total, hidden_fraction, generation) {
  w <- c(sqrt(1 - hidden_fraction), sqrt(hidden_fraction))
  S <- diag(2L)
  advance <- function(S) {
    cT <- drop(S %*% w)
    spouse <- rho_total * tcrossprod(cT) / sum(w * cT)
    list(child = .5 * (S + spouse) + .5 * diag(2L), spouse = spouse)
  }
  if (generation > 1L) {
    for (j in seq_len(generation - 1L)) S <- advance(S)$child
  }
  next_state <- advance(S)
  spouse <- next_state$spouse
  child <- next_state$child
  pc <- .5 * (S + spouse)
  list(parent = S, child = child, spouse = spouse, parent_child = pc,
       rho_visible = spouse[1L, 1L] / S[1L, 1L],
       variance_ratio = child[1L, 1L] / S[1L, 1L],
       eta_visible = pc[1L, 1L] / sqrt(S[1L, 1L] * child[1L, 1L]))
}

adversarial_pi_packet <- function(beta = c(.3, .1, -.05), beta_cov = diag(1e-4, 3L)) {
  labels <- c("bD", "bm", "bf", "rhoMF", "eta_m", "eta_f")
  V <- matrix(0, 6L, 6L, dimnames = list(labels, labels))
  V[1:3, 1:3] <- beta_cov
  list(schema = "ravel_visible_transfer_v1", diagnostics = list(admissible = TRUE),
       theta = setNames(c(beta, 0, .5, .5), labels), vcov_core = V,
       settings = list(uncertainty = "lavaan", ci_level = .95),
       outcome = list(name = "adversarial_trait"),
       software = list(transform_fingerprint = RAVEL:::.report_transform_fingerprint()))
}

adversarial_exact_trios <- function(h = .5, rho = .4, generation = 5L,
                                   score_scales = c(.7, 1.5, .9, 2.3, .8, 1.2)) {
  state <- adversarial_component_state(rho, h, generation)
  # Raw component order is Gc,Hc,Gm,Hm,Gf,Hf. This six-dimensional
  # construction specifies inheritance before forming either observed PGS.
  K <- rbind(cbind(state$child, state$parent_child, state$parent_child),
             cbind(state$parent_child, state$parent, state$spouse),
             cbind(state$parent_child, state$spouse, state$parent))
  total_map <- matrix(0, 3L, 6L)
  for (j in 1:3) total_map[j, (2*j-1):(2*j)] <- c(sqrt(1-h), sqrt(h))
  visible <- c(1L, 3L, 5L)
  G <- K[visible, visible]
  GT <- K[visible, , drop=FALSE] %*% t(total_map)
  beta_raw <- c(bD=.3, bm=.08, bf=.04)
  gy <- drop(GT %*% beta_raw)
  E <- matrix(c(1,.5,.5,.5,1,0,.5,0,1),3L,3L)
  sigma <- rbind(cbind(G+.4*E,G), cbind(G,G+1.6*E))
  sigma <- rbind(cbind(sigma,c(gy,gy)),c(gy,gy,1))
  sigma <- sigma * tcrossprod(c(score_scales,1))
  vars <- c(paste0("PGS1_",c("child","mother","father")),
            paste0("PGS2_",c("child","mother","father")),"Y")
  dimnames(sigma) <- list(vars,vars)
  n <- 257L
  # Fixed orthogonal columns provide exactly the target ML covariance.
  x <- outer(seq_len(n),seq_len(7L),function(i,j) sin(i*j) + cos(i*(j+.5)))
  Q <- qr.Q(qr(scale(x,center=TRUE,scale=FALSE)))
  dat <- as.data.frame(sqrt(n)*Q %*% chol(sigma))
  names(dat) <- vars
  vg <- diag(G)
  Rg <- G/sqrt(outer(vg,vg))
  beta_visible <- drop(solve(Rg,gy/sqrt(vg)))
  names(beta_visible) <- names(beta_raw)
  list(dat=dat,state=state,beta_raw=beta_raw,beta_visible=beta_visible,
       beta_std=beta_raw*sqrt(diag(total_map %*% K %*% t(total_map))),
       reliability1=vg/(vg+.4),reliability2=vg/(vg+1.6))
}

# Independent Gaussian expected information. The fitted parameter of interest
# here is full-liability spouse r directly; no package correction derivative or
# lavaan information matrix is used to obtain the reference standard error.
adversarial_spouse_fisher <- function(h=.4,rho=.4,generation=3L,N=5000L,
                                      pi_se=.02,step=1e-5) {
  state_geometry <- function(r,h) {
    x <- adversarial_component_state(r,h,generation)
    K <- rbind(cbind(x$child,x$parent_child,x$parent_child),
               cbind(x$parent_child,x$parent,x$spouse),
               cbind(x$parent_child,x$spouse,x$parent))
    Tmap <- matrix(0,3L,6L)
    for(j in 1:3) Tmap[j,(2*j-1):(2*j)] <- c(sqrt(1-h),sqrt(h))
    G <- K[c(1L,3L,5L),c(1L,3L,5L)]
    list(G=G,GT=K[c(1L,3L,5L),,drop=FALSE]%*%t(Tmap))
  }
  geometry <- state_geometry(rho,h)
  G <- geometry$G
  R <- stats::cov2cor(G)
  beta <- drop(solve(R,geometry$GT %*% c(.3,.08,.04)/sqrt(diag(G))))
  lambda <- sqrt(G[2,2]/(G[2,2]+c(.4,1.6)))
  par <- c(rho,lambda,rep(0,6L),beta,1)
  covariance <- function(par,h) {
    G <- state_geometry(par[1],h)$G
    R <- stats::cov2cor(G)
    growth <- G[1,1]/G[2,2]
    scale <- exp(par[4:9])
    parent_lambda <- par[2:3]
    child_lambda <- sqrt(growth*parent_lambda^2/(1-parent_lambda^2+growth*parent_lambda^2))
    block1 <- diag(c(child_lambda[1],rep(parent_lambda[1],2L)))
    block2 <- diag(c(child_lambda[2],rep(parent_lambda[2],2L)))
    L <- rbind(block1,block2)*scale
    C <- L%*%R%*%t(L)
    E <- matrix(c(1,.5,.5,.5,1,0,.5,0,1),3L)
    for(j in 1:2) {
      ix <- (3*j-2):(3*j)
      residual <- sqrt(1-c(child_lambda[j],rep(parent_lambda[j],2L))^2)*scale[ix]
      C[ix,ix] <- C[ix,ix]+E*tcrossprod(residual)
    }
    score_y <- drop(L%*%R%*%par[10:12])
    rbind(cbind(C,score_y),c(score_y,par[13]))
  }
  K <- covariance(par,h)
  W <- solve(K)
  J <- lapply(seq_along(par),function(j) {
    plus <- minus <- par
    plus[j] <- plus[j]+step
    minus[j] <- minus[j]-step
    (covariance(plus,h)-covariance(minus,h))/(2*step)
  })
  H <- (covariance(par,h+step)-covariance(par,h-step))/(2*step)
  inner <- function(A,B) sum(diag(W%*%A%*%W%*%B))/2
  info <- outer(seq_along(par),seq_along(par),Vectorize(function(i,j) inner(J[[i]],J[[j]])))
  h_cross <- vapply(J,inner,numeric(1),B=H)
  sensitivity <- -drop(solve(info,h_cross))
  variance <- solve(info)/N
  total_variance <- variance[1,1]+sensitivity[1]^2*pi_se^2
  full_parent <- 1
  full_parent_derivative <- 0
  if(generation>1L) for(j in seq_len(generation-1L)) {
    full_parent_derivative <- .5*full_parent+.5*(1+rho)*full_parent_derivative
    full_parent <- .5*(1+rho)*full_parent+.5
  }
  raw_covariance_derivative <- full_parent+rho*full_parent_derivative
  names <- c(paste0("PGS1_",c("child","mother","father")),
             paste0("PGS2_",c("child","mother","father")),"Y")
  dimnames(K) <- list(names,names)
  list(covariance=K,information=info,condition_number=kappa(info,exact=TRUE),
       conditional_se=sqrt(variance[1,1]),rho_pi_sensitivity=sensitivity[1],
       total_se=sqrt(total_variance),pi_variance_fraction=sensitivity[1]^2*pi_se^2/total_variance,
       raw_covariance_conditional_se=raw_covariance_derivative*sqrt(variance[1,1]),
       raw_covariance_total_se=raw_covariance_derivative*sqrt(total_variance),
       raw_covariance_pi_sensitivity=raw_covariance_derivative*sensitivity[1],
       h=h,rho=rho,generation=generation,N=N,pi_se=pi_se,step=step)
}

adversarial_lavaan_spouse_check <- function(reference) {
  RAVEL:::.ensure_lavaan_runtime()
  h <- reference$h
  generation <- reference$generation
  labels <- c("bD","bm","bf","rhoMF","eta_m","eta_f")
  fixture <- adversarial_exact_trios(h=h,rho=reference$rho,generation=generation,
                                   score_scales=rep(1,6L))
  fixture$dat[,1:6] <- scale(fixture$dat[,1:6])
  fit_one <- function(h) {
    model <- RAVEL:::.build_model(fixed_crossperson=c(e1_mf=0,e2_mf=0),
      fixed_crossperson_corr=c(e1_mc=.5,e1_fc=.5,e2_mc=.5,e2_fc=.5))
    model <- paste(model,RAVEL:::.ravel_constraint_syntax("n_gen",h,generation),sep="\n")
    model <- RAVEL:::.ravel_constraint_starts(model,fixture$dat,"n_gen",h,generation)
    fit <- lavaan::sem(model,sample.cov=reference$covariance,sample.nobs=reference$N,
      sample.cov.rescale=FALSE,std.lv=TRUE,estimator="ML",se="standard",
      information="expected",warn=FALSE)
    if(!RAVEL:::.ravel_constraint_fit_valid(fit)) stop("Population moment fit failed")
    pe <- lavaan::parameterEstimates(fit)
    theta <- setNames(vapply(labels,function(label) pe$est[match(label,pe$label)],numeric(1)),labels)
    list(fit=fit,theta=theta)
  }
  center <- fit_one(h)
  step <- .001
  sensitivity <- (fit_one(h+step)$theta-fit_one(h-step)$theta)/(2*step)
  object <- structure(list(converged=TRUE,fit=center$fit,
    beta_vis=center$theta[1:3],latent_cov_vis=center$theta[4:6],
    am_correction="n_gen",am_n_generations=generation,
    constraint_spec=list(mode="n_gen",n_generations=generation,h=h,depends_on_pi=TRUE)),
    class="ravel_results")
  corrected <- ravel_correct_pi(object,h,h+c(-1,1)*qnorm(.975)*reference$pi_se,
    core_pi_sensitivity=sensitivity)
  index <- which(corrected$scale=="corrected_std" & corrected$parameter=="rho_corrected")
  raw_index <- which(corrected$scale=="corrected_raw" & corrected$parameter=="cov_corrected")
  list(corrected=corrected,conditional_se=sqrt(attr(corrected,"variance_components")$core[index,index]),
       total_se=corrected$se[index],rho_pi_sensitivity=attr(corrected,"jacobian")[index,"pi"],
       raw_covariance_conditional_se=sqrt(attr(corrected,"variance_components")$core[raw_index,raw_index]),
       raw_covariance_total_se=corrected$se[raw_index],
       raw_covariance_pi_sensitivity=attr(corrected,"jacobian")[raw_index,"pi"],
       rho_estimate=corrected$estimate[index])
}

test_that("the finite-generation state follows an independent two-component recurrence", {
  for (h in c(0, .2, .5, .8)) for (rho in c(-.7, -.2, 0, .2, .7, .95)) {
    for (generation in c(1L, 2L, 3L, 5L, 10L, 100L)) {
      ref <- adversarial_component_state(rho, h, generation)
      got <- RAVEL:::.am_state_components(rho, h, "n_gen", generation)
      expect_equal(got$vG, ref$parent[1, 1], tolerance = 5e-12)
      expect_equal(got$vGc, ref$child[1, 1], tolerance = 5e-12)
      expect_equal(got$rho_G, ref$rho_visible, tolerance = 5e-12)
      expect_equal(got$vGc / got$vG, ref$variance_ratio, tolerance = 5e-12)
      expect_equal(.5 * (1 + got$rho_G) * sqrt(got$vG / got$vGc),
                   ref$eta_visible, tolerance = 5e-12)
    }
  }
})

test_that("external pi variance gives the analytic random-mating correction variance", {
  pi <- .4
  pi_se <- .025
  beta <- c(.3, .1, -.05)
  Vbeta <- matrix(c(2e-4, 4e-5, -1e-5, 4e-5, 3e-4, 2e-5,
                    -1e-5, 2e-5, 1e-4), 3L, 3L)
  packet <- adversarial_pi_packet(beta, Vbeta)
  ci <- pi + c(-1, 1) * qnorm(.975) * pi_se
  got <- ravel_correct_pi(packet, pi, ci, am_correction = "equilibrium")
  expected <- Vbeta / (1 - pi) + tcrossprod(beta) * pi_se^2 / (4 * (1 - pi)^3)
  expect_equal(got$estimate[1:3], beta / sqrt(1 - pi), tolerance = 1e-10)
  expect_equal(unname(attr(got, "covariance")[1:3, 1:3]), expected, tolerance = 2e-8)
  expect_equal(unname(attr(got, "covariance")[5:7, 5:7]), expected, tolerance = 2e-8)
  expect_equal(got$pi_se, rep(pi_se, 8L), tolerance = 1e-12)

  # Matched working-scale SEs yield the same first-order raw-pi variance.
  logit_ci <- plogis(qlogis(pi) + c(-1, 1) * qnorm(.975) * pi_se / (pi * (1 - pi)))
  logit <- ravel_correct_pi(packet, pi, logit_ci, pi_ci_scale = "logit",
                            am_correction = "equilibrium")
  expect_equal(attr(logit, "covariance"), attr(got, "covariance"), tolerance = 1e-10)
})

test_that("cross covariance is propagated once on the conditional core scale", {
  pi <- .4
  pi_se <- .025
  beta <- c(.3, .1, -.05)
  packet <- adversarial_pi_packet(beta)
  cross <- setNames(c(1e-5, -2e-5, 3e-6, 0, 0, 0), names(packet$theta))
  got <- ravel_correct_pi(packet, pi, pi + c(-1, 1) * qnorm(.975) * pi_se,
                          core_pi_cov = cross, am_correction = "equilibrium")
  g <- beta / (2 * (1 - pi)^1.5)
  expected <- packet$vcov_core[1:3, 1:3] / (1 - pi) + tcrossprod(g) * pi_se^2 +
    (tcrossprod(cross[1:3], g) + tcrossprod(g, cross[1:3])) / sqrt(1 - pi)
  expect_equal(unname(attr(got, "covariance")[1:3, 1:3]), unname(expected), tolerance = 2e-8)
  expect_error(ravel_correct_pi(packet, pi, c(pi, pi), core_pi_cov = cross,
                                am_correction = "equilibrium"), "zero-variance pi")
})

test_that("pi uncertainty does not silently infer variance from incompatible intervals", {
  packet <- adversarial_pi_packet()
  expect_error(ravel_correct_pi(packet, .4, c(.35, .5)), "symmetric")
  expect_error(ravel_correct_pi(packet, .4, c(0, .6), pi_ci_scale = "logit"), "strictly")
  expect_error(ravel_correct_pi(packet, 0, c(0, .05)), "boundary")
  expect_error(ravel_correct_pi(packet, .4, c(.3, .5), am_n_generations = 1.5), "integer")
  zero <- ravel_correct_pi(packet, 0, c(0, 0), am_correction = "equilibrium")
  expect_equal(zero$estimate[1:3], unname(packet$theta[1:3]), tolerance = 1e-10)
  expect_equal(unname(attr(zero, "covariance")[1:3, 1:3]),
               unname(packet$vcov_core[1:3, 1:3]), tolerance = 1e-8)
})

test_that("constrained finite fits recover unequal score reliabilities and changing generations", {
  for (rho in c(-.5,.5)) {
    fixture <- adversarial_exact_trios(rho=rho,generation=3L)
    result <- suppressWarnings(ravel_fit(fixture$dat,h=.5,
      fit_constraints="n_gen",am_n_generations=3L))
    expect_true(isTRUE(result$converged))
    if (!isTRUE(result$converged)) next
    expect_equal(result$beta_vis,fixture$beta_visible,tolerance=2e-4)
    expect_equal(result$beta_total_assumed,fixture$beta_raw,tolerance=3e-4)
    expect_equal(result$beta_total_assumed_std,fixture$beta_std,tolerance=3e-4)
    expect_equal(unname(result$rho_vis),fixture$state$rho_visible,tolerance=2e-4)
    expect_equal(unname(result$latent_cov_vis[c("eta_m","eta_f")]),
                 rep(fixture$state$eta_visible,2L),tolerance=2e-4)
    # The score effects are deliberately in different raw units by role.
    # Equal raw loadings would therefore fail even though reliability is valid.
    rel <- setNames(result$reliability$reliability,
      paste(result$reliability$PGS,result$reliability$role,sep="_"))
    expect_equal(unname(rel[c("PGS1_child","PGS1_mother","PGS1_father")]),
                 unname(fixture$reliability1),tolerance=2e-4)
    expect_equal(unname(rel[c("PGS2_child","PGS2_mother","PGS2_father")]),
                 unname(fixture$reliability2),tolerance=2e-4)
  }
})

test_that("one-generation constrained measurement fitting is independent of hidden fraction", {
  fixture <- adversarial_exact_trios(h=.4,rho=.4,generation=1L)
  a <- suppressWarnings(ravel_fit(fixture$dat,h=.1,fit_constraints="n_gen",am_n_generations=1L))
  b <- suppressWarnings(ravel_fit(fixture$dat,h=.6,fit_constraints="n_gen",am_n_generations=1L))
  expect_true(isTRUE(a$converged))
  expect_true(isTRUE(b$converged))
  if (isTRUE(a$converged) && isTRUE(b$converged)) {
    expect_equal(a$beta_vis,b$beta_vis,tolerance=2e-5)
    expect_equal(a$latent_cov_vis,b$latent_cov_vis,tolerance=2e-5)
    expect_equal(a$reliability,b$reliability,tolerance=2e-5)
    expect_false(a$constraint_spec$depends_on_pi)
    expect_false(b$constraint_spec$depends_on_pi)
  }
})

test_that("finite pi propagation uses the full fitted-core chain derivative", {
  pi <- .4
  pi_se <- .02
  fixture <- adversarial_exact_trios(h=pi,rho=.5,generation=3L)
  packet <- adversarial_pi_packet()
  labels <- names(packet$theta)
  packet$theta <- setNames(c(fixture$beta_visible,fixture$state$rho_visible,
    fixture$state$eta_visible,fixture$state$eta_visible),labels)
  packet$settings$constraint_spec <- list(mode="n_gen",h=pi,n_generations=3L,
    depends_on_pi=TRUE)
  # A specified local fitted-core path isolates the chain rule from optimizer
  # noise. The conditional estimator covariance and pi covariance are distinct.
  sensitivity <- setNames(c(.04,-.01,.02,.015,-.008,-.008),labels)
  cross <- setNames(c(1e-5,-2e-6,0,0,0,0),labels)
  before <- packet
  result <- ravel_correct_pi(packet,pi,pi+c(-1,1)*qnorm(.975)*pi_se,
    am_n_generations=3L,core_pi_sensitivity=sensitivity,core_pi_cov=cross)
  expect_identical(packet,before)
  direct_path <- function(z) {
    core <- z[labels]+sensitivity*(z["pi"]-pi)
    estimate <- RAVEL:::.compute_scale_estimates(core[1:3],core["eta_m"],
      core["eta_f"],core["rhoMF"],z["pi"],"n_gen",3L)
    c(estimate$beta_total_assumed_std,estimate$rho_total_assumed,
      estimate$beta_total_assumed,estimate$cov_total_assumed)
  }
  x <- c(packet$theta,pi=pi)
  J <- matrix(0,8L,7L)
  for (j in seq_along(x)) {
    plus <- minus <- x
    plus[j] <- plus[j]+1e-5
    minus[j] <- minus[j]-1e-5
    J[,j] <- (direct_path(plus)-direct_path(minus))/(2e-5)
  }
  expect_equal(unname(attr(result,"jacobian")),J,tolerance=2e-6)
  V <- matrix(0,7L,7L)
  V[1:6,1:6] <- packet$vcov_core
  V[7,7] <- pi_se^2
  V[1:6,7] <- V[7,1:6] <- cross
  expect_equal(unname(attr(result,"covariance")),J %*% V %*% t(J),tolerance=2e-6)
  expect_gt(max(abs(attr(result,"jacobian")[,"pi"]-
                    attr(result,"partial_jacobian")[,"pi"])),.02)
  expect_error(ravel_correct_pi(packet,pi,pi+c(-.02,.02),am_n_generations=3L),
               "core_pi_sensitivity")
  expect_error(ravel_correct_pi(packet,pi+.01,pi+c(-.01,.03),am_n_generations=3L,
                 core_pi_sensitivity=sensitivity),"Refit")
})

test_that("finite spouse SE matches independent Gaussian expected information", {
  reference <- adversarial_spouse_fisher()
  actual <- adversarial_lavaan_spouse_check(reference)
  expect_equal(actual$rho_estimate,reference$rho,tolerance=3e-5)
  expect_equal(actual$conditional_se,reference$conditional_se,tolerance=3e-4)
  expect_equal(actual$rho_pi_sensitivity,reference$rho_pi_sensitivity,tolerance=3e-3)
  expect_equal(actual$total_se,reference$total_se,tolerance=5e-4)
  expect_equal(actual$raw_covariance_conditional_se,reference$raw_covariance_conditional_se,tolerance=3e-4)
  expect_equal(actual$raw_covariance_total_se,reference$raw_covariance_total_se,tolerance=5e-4)
  expect_equal(actual$raw_covariance_pi_sensitivity,reference$raw_covariance_pi_sensitivity,tolerance=3e-3)
})
