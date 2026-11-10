# Observed-PGS comparators; standalone base-R likelihood and FID score sandwich.
# Genetic score measurement error and hidden-variant corrections are not applied.
.ravel_pgs_delta <- function(fun, theta, V = NULL) {
  value <- fun(theta); J <- matrix(0, length(value), length(theta))
  for (j in seq_along(theta)) {
    h <- 1e-5 * max(1, abs(theta[j])); a <- b <- theta
    a[j] <- a[j] + h; b[j] <- b[j] - h
    J[, j] <- (fun(a) - fun(b)) / (2*h)
  }
  covariance <- if (is.null(V)) matrix(NA_real_, length(value), length(value)) else J %*% V %*% t(J)
  dimnames(covariance) <- list(names(value), names(value))
  list(value = value, vcov = covariance, jacobian = J)
}

.ravel_pgs_moments <- function(scores, iid, cluster, version) {
  take <- !duplicated(iid); x <- as.matrix(scores[take, , drop = FALSE]); id <- cluster[take]
  n <- nrow(x); q <- ncol(x); ix <- which(lower.tri(matrix(0,q,q),diag=TRUE), arr.ind=TRUE)
  z <- sweep(x,2L,colMeans(x),"-"); S <- crossprod(z)/n
  covariance <- S[ix]; correlation <- stats::cov2cor(S)[ix]
  influence <- vapply(seq_len(nrow(ix)),function(j) {
    a <- ix[j,1L]; b <- ix[j,2L]; z[,a]*z[,b]-S[a,b]
  },numeric(n))
  correlation_influence <- vapply(seq_len(nrow(ix)),function(j) {
    a <- ix[j,1L]; b <- ix[j,2L]
    influence[,j]/sqrt(S[a,a]*S[b,b])-correlation[j]/2*(
      (z[,a]^2-S[a,a])/S[a,a]+(z[,b]^2-S[b,b])/S[b,b])
  },numeric(n))
  grouped <- rowsum(cbind(influence,correlation_influence),id,reorder=FALSE)
  se <- if(nrow(grouped)>1L)sqrt(colSums(grouped^2))/n else rep(NA_real_,2*nrow(ix))
  estimate <- c(covariance,correlation)
  data.frame(regression=version, block=rep(c("observed_PGS_covariance","observed_PGS_correlation"),each=nrow(ix)),
    lhs=rep(colnames(x)[ix[,1L]],2L),rhs=rep(colnames(x)[ix[,2L]],2L),
    estimate=unname(estimate),se=se,ci_low=unname(estimate)-stats::qnorm(.975)*se,
    ci_high=unname(estimate)+stats::qnorm(.975)*se,
    scale=rep(c("Input PGS units (N-denominator covariance)","Correlation"),each=nrow(ix)),
    scope="Observed role scores; one row per IID, child-weighted parental moments",
    N=n,SE_clusters=length(unique(id)),stringsAsFactors=FALSE)
}

.ravel_pgs_groups <- function(X,y,family,iid,cluster,slope,family_random,individual_random) {
  members <- split(seq_along(y),family)
  patterns <- lapply(members,function(rows) {
    rows <- rows[order(iid[rows],rows)]; person <- match(iid[rows],unique(iid[rows]))
    list(rows=rows,person=person,slope=if(is.null(slope)) NULL else slope[rows])
  })
  keys <- vapply(patterns,function(z) paste(paste(z$person,collapse=","),
    if(is.null(z$slope)) "" else paste(format(z$slope,digits=17,scientific=TRUE,trim=TRUE),collapse=","),sep="|"),"")
  groups <- split(seq_along(patterns),factor(keys,levels=unique(keys)))
  lapply(groups,function(at) {
    z <- patterns[[at[1L]]]; rows <- do.call(rbind,lapply(patterns[at],`[[`,"rows"))
    m <- ncol(rows); nf <- nrow(rows); q <- ncol(X)
    design <- lapply(seq_len(m),function(j) X[rows[,j],,drop=FALSE])
    Y <- matrix(y[rows],nf,m)
    same <- outer(z$person,z$person,"==")*1
    B <- list(residual=diag(m))
    if(family_random) B$family <- matrix(1,m,m)
    if(individual_random) B$individual <- same
    if(!is.null(slope)) {
      B$slope <- same*outer(z$slope,z$slope)
      B$intercept_slope <- same*outer(z$slope,z$slope,"+")
    }
    list(Y=Y,X=design,n=nf,m=m,B=B,
      families=names(members)[at],clusters=cluster[rows[,1L]],
      XX=if(nf>=m*q) lapply(design,function(a) lapply(design,function(b) crossprod(a,b))) else NULL)
  })
}

.ravel_pgs_covariance <- function(par,B) {
  S <- exp(2*par[["residual_sd"]])*B$residual
  derivative <- list(residual_sd=2*S)
  for(nm in intersect(c("family","individual","slope"),names(B))) {
    name <- paste0(nm,"_sd"); component <- exp(2*par[[name]])*B[[nm]]
    S <- S+component; derivative[[name]] <- 2*component
  }
  if("intercept_slope"%in%names(B)) {
    a <- exp(par[["individual_sd"]]+par[["slope_sd"]]); rho <- tanh(par[["intercept_slope_cor"]])
    component <- a*rho*B$intercept_slope; S <- S+component
    derivative$individual_sd <- derivative$individual_sd+component
    derivative$slope_sd <- derivative$slope_sd+component
    derivative$intercept_slope_cor <- a*(1-rho^2)*B$intercept_slope
  }
  list(S=S,derivative=derivative[names(par)])
}

.ravel_pgs_evaluator <- function(groups) {
  q <- ncol(groups[[1L]]$X[[1L]]); cache <- new.env(parent=emptyenv())
  function(par,gradient=FALSE) {
    if(exists("par",cache,inherits=FALSE)&&identical(par,cache$par)&&
       (!gradient||!is.null(cache$value$gradient))) return(cache$value)
    value <- tryCatch({
      A <- matrix(0,q,q); b <- numeric(q); work <- vector("list",length(groups)); constant <- 0
      for(g in seq_along(groups)) {
        z <- groups[[g]]; c <- .ravel_pgs_covariance(par,z$B)
        ch <- chol(c$S); Q <- chol2inv(ch); U <- z$Y%*%Q
        for(i in seq_len(z$m)) {
          b <- b+drop(crossprod(z$X[[i]],U[,i]))
          for(j in seq_len(z$m)) if(Q[i,j]!=0) A <- A+Q[i,j]*
            if(is.null(z$XX)) crossprod(z$X[[i]],z$X[[j]]) else z$XX[[i]][[j]]
        }
        constant <- constant+z$n*(2*sum(log(diag(ch)))+z$m*log(2*pi))+sum(U*z$Y)
        work[[g]] <- c(list(Q=Q),c)
      }
      beta <- drop(solve(A,b)); value <- .5*(constant-sum(b*beta))
      grad <- stats::setNames(numeric(length(par)),names(par))
      for(g in seq_along(groups)) {
        z <- groups[[g]]; w <- work[[g]]
        mu <- do.call(cbind,lapply(z$X,function(x) drop(x%*%beta)))
        w$residual <- z$Y-mu; w$U <- w$residual%*%w$Q
        if(gradient) for(j in seq_along(par)) grad[j] <- grad[j]+.5*(
          z$n*sum(w$Q*w$derivative[[j]])-sum((w$U%*%w$derivative[[j]])*w$U))
        work[[g]] <- w
      }
      list(valid=TRUE,value=value,beta=beta,A=A,work=work,gradient=if(gradient) grad else NULL)
    },error=function(e) list(valid=FALSE,value=1e50,reason=conditionMessage(e),
      gradient=if(gradient) stats::setNames(rep(0,length(par)),names(par)) else NULL))
    cache$par <- par;cache$value <- value;value
  }
}

.ravel_pgs_continuous <- function(X,y,family,iid,cluster,slope,family_random,individual_random,control) {
  groups <- .ravel_pgs_groups(X,y,family,iid,cluster,slope,family_random,individual_random)
  evaluate <- .ravel_pgs_evaluator(groups)
  r <- stats::lm.fit(X,y)$residuals; variance <- mean(r^2)
  if(!is.finite(variance)||variance<=1e-12) stop("Comparator residual variance is not positive.")
  par <- c(residual_sd=log(sqrt(variance*.7)))
  if(family_random) par <- c(par,family_sd=log(sqrt(variance*.15)))
  if(individual_random) par <- c(par,individual_sd=log(sqrt(variance*.15)))
  if(!is.null(slope)) par <- c(par,slope_sd=log(sqrt(variance*.02)/stats::sd(slope)),intercept_slope_cor=0)
  lower <- stats::setNames(rep(log(1e-5),length(par)),names(par)); upper <- lower*0+log(1e3)
  if("intercept_slope_cor"%in%names(par)) { lower["intercept_slope_cor"] <- -5;upper["intercept_slope_cor"] <- 5 }
  if(length(par)==1L) { par[1L] <- log(sqrt(variance)); optimization <- list(convergence=0L,par=par,message="Closed-form ML") } else {
    opt <- utils::modifyList(list(iter.max=500L,eval.max=800L,rel.tol=1e-13,x.tol=1e-10),control)
    optimization <- stats::nlminb(par,function(x)evaluate(x)$value,
      function(x)evaluate(x,TRUE)$gradient,lower=lower,upper=upper,control=opt)
    par <- optimization$par
    # Resolve early objective-tolerance termination using the analytic gradient.
    for (refinement in seq_len(3L)) {
      current <- evaluate(par,TRUE)
      H <- tryCatch(stats::optimHess(par,function(x)evaluate(x)$value,
        function(x)evaluate(x,TRUE)$gradient),error=function(e)NULL)
      if(is.null(H)||any(!is.finite(H))||min(eigen(H,symmetric=TRUE,only.values=TRUE)$values)<=0)break
      step <- tryCatch(drop(solve(H,current$gradient)),error=function(e)NULL)
      if(is.null(step)||max(abs(step))>1)break
      candidate <- par-step
      if(any(candidate<=lower)||any(candidate>=upper))break
      refined <- evaluate(candidate,TRUE)
      if(!isTRUE(refined$valid)||refined$value>current$value+1e-10||
         max(abs(refined$gradient))>=max(abs(current$gradient)))break
      par <- candidate
      if(max(abs(refined$gradient))<1e-8)break
    }
    optimization$refined_par <- par
    optimization$refined_objective <- evaluate(par)$value
  }
  final <- evaluate(par,TRUE)
  if(!isTRUE(final$valid)) stop(final$reason)
  nf <- length(unique(family)); q <- ncol(X); p <- length(par)
  natural <- final$gradient[intersect(c("residual_sd","family_sd","individual_sd","slope_sd"),names(par))]
  natural <- natural/(2*exp(2*par[names(natural)]))
  near <- intersect(c("family_sd","individual_sd","slope_sd"),names(par))
  boundary <- near[exp(2*par[near])<1e-6]
  projected <- final$gradient
  projected[(par-lower<1e-6&projected>0)|(upper-par<1e-6&projected<0)] <- 0
  natural[intersect(boundary,names(natural)[natural>0])] <- 0
  if("intercept_slope_cor"%in%names(par)&&abs(tanh(par[["intercept_slope_cor"]]))>.999) boundary <- c(boundary,"intercept_slope_cor")
  converged <- all(is.finite(c(projected,natural)))&&max(abs(c(projected,natural)))/nf<1e-5
  status <- if(!converged) "not_converged" else if(length(boundary)) "variance_boundary" else "ok"
  note <- if(!converged) "Comparator gradient check failed." else if(length(boundary))
    paste("Variance boundary:",paste(boundary,collapse=", ")) else ""
  V <- NULL; scores <- NULL
  if(status=="ok") {
    inferred <- tryCatch({
      Hprofile <- stats::optimHess(par,function(x)evaluate(x)$value,
        function(x)evaluate(x,TRUE)$gradient,control=list(ndeps=rep(1e-5,p)))
      J <- matrix(0,q,p); scores <- matrix(0,nf,p+q)
      offset <- 0L; outer <- character(nf)
      for(g in seq_along(groups)) {
        z <- groups[[g]]; w <- final$work[[g]]; at <- offset+seq_len(z$n);offset <- offset+z$n
        outer[at] <- z$clusters
        for(j in seq_len(p)) scores[at,j] <- .5*(rowSums((w$U%*%w$derivative[[j]])*w$U)-sum(w$Q*w$derivative[[j]]))
        for(i in seq_len(z$m)) scores[at,p+seq_len(q)] <- scores[at,p+seq_len(q),drop=FALSE]+z$X[[i]]*w$U[,i]
      }
      for(j in seq_len(p)) {
        h <- 1e-5*max(1,abs(par[j]));a <- b <- par;a[j] <- a[j]+h;b[j] <- b[j]-h
        J[,j] <- (evaluate(a)$beta-evaluate(b)$beta)/(2*h)
      }
      H <- rbind(cbind(Hprofile+crossprod(J,final$A%*%J),-crossprod(J,final$A)),
        cbind(-final$A%*%J,final$A));H <- (H+t(H))/2
      e <- eigen(H,symmetric=TRUE,only.values=TRUE)$values
      if(min(e)<=max(e)*1e-9) stop("Comparator observed information is singular or indefinite.")
      bread <- solve(H); grouped <- rowsum(scores,outer,reorder=FALSE)
      if(nrow(grouped)<2L) stop("At least two SE clusters are required.")
      V <- bread%*%crossprod(grouped)%*%bread
      dimnames(V) <- list(c(names(par),colnames(X)),c(names(par),colnames(X)))
      list(vcov=(V+t(V))/2,cluster_scores=grouped,information=H)
    },error=identity)
    if(inherits(inferred,"error")) {status <- "inference_unavailable";note <- conditionMessage(inferred)} else V <- inferred$vcov
  }
  list(beta=stats::setNames(final$beta,colnames(X)),par=par,vcov=V,status=status,
    converged=converged,note=note,optimizer=optimization,
    fit_statistics=c(logLik=-final$value,n_parameters=p+q,AIC=2*final$value+2*(p+q)),
    method="Observed-PGS Gaussian mixed ML; FID-clustered CR0 joint score sandwich",
    likelihood_scope="Y | observed PGS, covariates",
    scope="Observed-PGS Gaussian regression")
}

.ravel_pgs_binary <- function(X,y,cluster) {
  messages <- character()
  fit <- withCallingHandlers(stats::glm.fit(X,y,family=stats::binomial("probit")),
    warning=function(w){messages <<- c(messages,conditionMessage(w));invokeRestart("muffleWarning")})
  beta <- fit$coefficients;eta <- drop(X%*%beta)
  positive <- exp(stats::dnorm(eta,log=TRUE)-stats::pnorm(eta,log.p=TRUE))
  negative <- exp(stats::dnorm(eta,log=TRUE)-stats::pnorm(-eta,log.p=TRUE))
  score <- ifelse(y==1,positive,-negative)
  curvature <- ifelse(y==1,positive*(eta+positive),negative*(-eta+negative))
  H <- crossprod(X,X*curvature); H <- (H+t(H))/2
  V <- NULL; note <- paste(unique(messages),collapse="; ")
  converged <- isTRUE(fit$converged)&&all(is.finite(beta))&&max(abs(colSums(X*score)))/nrow(X)<1e-5
  separation <- all((2*y-1)*eta > 0) || any(abs(beta)>30)
  status <- if(!converged) "not_converged" else if(separation) "separation" else "ok"
  if(separation) note <- paste(note,"Complete or numerical separation in the probit comparator.")
  if(status=="ok") {
    inference <- tryCatch({
      e <- eigen(H,symmetric=TRUE,only.values=TRUE)$values
      if(min(e)<=1e-9*max(e)) stop("Probit comparator information is singular or indefinite.")
      grouped <- rowsum(X*score,cluster,reorder=FALSE)
      if(nrow(grouped)<2L) stop("At least two SE clusters are required.")
      B <- solve(H);B%*%crossprod(grouped)%*%B
    },error=identity)
    if(inherits(inference,"error")) { status <- "inference_unavailable";note <- paste(note,conditionMessage(inference)) } else V <- inference
  }
  if(!is.null(V)) dimnames(V) <- list(colnames(X),colnames(X))
  list(beta=stats::setNames(beta,colnames(X)),par=numeric(),vcov=V,status=status,converged=converged,
    note=trimws(note),optimizer=fit,
    fit_statistics=c(logLik=-fit$deviance/2,n_parameters=fit$rank,AIC=fit$aic),
    method="Observed-PGS probit GLM; FID-clustered CR0 observed-score sandwich",
    likelihood_scope="Y | observed PGS, covariates",
    scope="Observed-PGS probit GLM")
}

.ravel_pgs_comparators <- function(dat,outcome,pgs_columns,covariates=NULL,
    family_id="M_ID",subject_id="IID",cluster="FID",family_random=FALSE,
    individual_random=FALSE,random_slope=NULL,interactions=NULL,interaction_center=NULL,
    binary_outcome=FALSE,control=list()) {
  roles <- unique(sub("[12]$","",names(pgs_columns)))
  if(!length(roles)||any(!c(outer(roles,1:2,paste0))%in%names(pgs_columns))) stop("Supply two named PGS columns for each role.")
  fixed <- unique(c(covariates,interactions,random_slope)); needed <- unique(c(outcome,unname(pgs_columns),fixed,family_id,subject_id,cluster))
  absent <- setdiff(needed,names(dat));if(length(absent))stop("Comparator columns missing: ",paste(absent,collapse=", "))
  for(id in c(family_id,subject_id,cluster)) if(anyNA(dat[[id]])||any(!nzchar(as.character(dat[[id]]))))stop("Comparator IDs must be observed: ",id)
  if(!is.null(random_slope)&&!individual_random)stop("A comparator random slope requires an individual intercept.")
  if(!all(vapply(dat[c(outcome,unname(pgs_columns),fixed)],is.numeric,logical(1))))stop("Comparator outcome, PGSs and covariates must be numeric.")
  numeric_columns <- unique(c(outcome,unname(pgs_columns),fixed))
  if(any(vapply(dat[numeric_columns],function(x)any(!is.na(x)&!is.finite(x)),logical(1))))stop("Comparator numeric inputs must be finite or NA.")
  keep <- stats::complete.cases(dat[numeric_columns]);d <- dat[keep,,drop=FALSE]
  if(nrow(d)<50L)stop("Fewer than 50 complete rows remain for PGS comparators.")
  y <- d[[outcome]];fam <- as.character(d[[family_id]]);iid <- as.character(d[[subject_id]]);fid <- as.character(d[[cluster]])
  map <- unique(data.frame(iid=iid,family=fam,cluster=fid))
  if(anyDuplicated(map$iid)||anyDuplicated(unique(map[c("family","cluster")])$family))stop("Comparator IID/family/FID nesting is inconsistent.")
  if(individual_random&&!any(table(iid)>1L)&&!binary_outcome)stop("An individual comparator effect requires repeated outcomes.")
  if(family_random&&individual_random&&!any(table(unique(map[c("iid","family")])$family)>1L)&&!binary_outcome)
    stop("Family and individual comparator effects cannot be separated with one child per family.")
  if(family_random&&!any(table(fam)>1L)&&!binary_outcome)stop("A family comparator effect needs some repeated family observations.")
  if(binary_outcome&&!setequal(unique(y),c(0,1)))stop("Binary comparator outcomes must contain 0 and 1.")
  centers <- stats::setNames(vapply(fixed,function(nm) mean(d[[nm]]),numeric(1)),fixed)
  if(length(interaction_center)) {
    if(is.null(names(interaction_center))||any(!is.finite(interaction_center))||any(!names(interaction_center)%in%fixed))stop("Comparator interaction centers must be named supplied covariates.")
    centers[names(interaction_center)] <- interaction_center
  }
  C <- lapply(fixed,function(nm)d[[nm]]-centers[[nm]]);names(C)<-fixed
  slope <- if(is.null(random_slope))NULL else C[[random_slope]]
  if(!is.null(slope)&&!binary_outcome&&!any(vapply(split(slope,iid),function(x)length(unique(x))>1L,logical(1))))
    stop("A comparator individual slope needs within-person age variation.")
  parameter <- c(child="beta_D",mother="beta_M",father="beta_F",mat_sib="beta_MS",pat_sib="beta_FS")
  if(any(!roles%in%names(parameter)))stop("Unsupported comparator PGS role.")
  rows <- fits <- scales <- moments <- list()
  for(k in 1:2) {
    version <- paste0("PGS",k); columns <- pgs_columns[paste0(roles,k)]
    raw <- d[unname(columns)];names(raw)<-roles; Z <- matrix(NA_real_,nrow(d),length(roles),dimnames=list(NULL,parameter[roles]))
    scaling <- list()
    for(j in seq_along(roles)) {
      id <- if(roles[j]=="child")iid else fam
      groups <- split(raw[[j]],id)
      if(any(vapply(groups,function(x)diff(range(x))>1e-10*max(abs(x)),logical(1))))
        stop("Comparator PGS changes within its person/family role: ",columns[j])
      unique_score <- raw[[j]][!duplicated(id)];center <- mean(unique_score);scale <- stats::sd(unique_score)
      if(!is.finite(scale)||scale<=0)stop("Comparator PGS has no variance: ",columns[j])
      Z[,j] <- (raw[[j]]-center)/scale
      scaling[[j]] <- data.frame(regression=version,role=roles[j],column=unname(columns[j]),center=center,scale=scale,
        N=length(unique_score),reference=if(roles[j]=="child")"Unique IID" else "One score per nuclear family",stringsAsFactors=FALSE)
    }
    X <- cbind(intercept=1,Z)
    if(length(C)) X <- cbind(X,stats::setNames(as.data.frame(C),paste0("fixed_",fixed)))
    X <- as.matrix(X)
    for(mod in interactions) X <- cbind(X,stats::setNames(as.data.frame(Z*C[[mod]]),paste0(colnames(Z),":",mod)))
    X <- as.matrix(X); xs <- apply(X,2,stats::sd);xs["intercept"]<-1
    if(any(!is.finite(xs))||any(xs<=0))stop("Comparator design contains a constant or redundant supplied column.")
    Xscaled <- sweep(X,2,xs,"/")
    if(qr(Xscaled,tol=1e-10,LAPACK=FALSE)$rank<ncol(Xscaled))stop("Comparator fixed-effect design is not full rank.")
    fit <- if(binary_outcome).ravel_pgs_binary(Xscaled,y,fid) else
      .ravel_pgs_continuous(Xscaled,y,fam,iid,fid,slope,family_random,individual_random,control)
    internal_beta <- fit$beta
    beta <- internal_beta/xs; names(beta)<-colnames(X)
    V <- if(is.null(fit$vcov))NULL else fit$vcov[names(beta),names(beta),drop=FALSE]/outer(xs,xs)
    fit$parameters <- beta;fit$coefficient_vcov<-V
    fit$internal_parameters<-c(fit$par,internal_beta)
    p <- length(fit$par)
    variance_value <- function(a) {
      out <- c(residual_variance=exp(2*a[["residual_sd"]]))
      for(nm in intersect(c("family_sd","individual_sd","slope_sd"),names(a)))out[sub("_sd$","_variance",nm)]<-exp(2*a[[nm]])
      if("intercept_slope_cor"%in%names(a))out<-c(out,intercept_slope_covariance=exp(a[["individual_sd"]]+a[["slope_sd"]])*tanh(a[["intercept_slope_cor"]]))
      out
    }
    if(p) {
      v <- .ravel_pgs_delta(variance_value,fit$par,if(is.null(fit$vcov))NULL else fit$vcov[seq_len(p),seq_len(p),drop=FALSE])
      se <- sqrt(pmax(0,diag(v$vcov)));fit$variance_vcov<-v$vcov
      fit$variance_components<-data.frame(parameter=names(v$value),estimate=unname(v$value),se=se,
        ci_low=unname(v$value)-stats::qnorm(.975)*se,ci_high=unname(v$value)+stats::qnorm(.975)*se,status=fit$status)
    } else fit$variance_components<-data.frame()
    genetic <- c(parameter[roles],unlist(lapply(interactions,function(m)paste0(parameter[roles],":",m))))
    contrasts <- diag(length(beta));dimnames(contrasts)<-list(names(beta),names(beta));contrasts<-contrasts[genetic,,drop=FALSE]
    for(suffix in c("",if(length(interactions))paste0(":",interactions))) {
      add <- function(name,pair,weights) {
        pair<-paste0(pair,suffix)
        if(all(pair%in%names(beta))) {
          a<-numeric(length(beta));names(a)<-names(beta);a[pair]<-weights
          contrasts<<-rbind(contrasts,stats::setNames(data.frame(t(a)),names(beta)))
          rownames(contrasts)[nrow(contrasts)]<<-paste0(name,suffix)
        }
      }
      add("beta_mid",c("beta_M","beta_F"),c(1,1));add("delta_M_mid",c("beta_M","beta_F"),c(1,-1))
      add("delta_P_mid",c("beta_F","beta_M"),c(1,-1));add("delta_M",c("beta_M","beta_MS"),c(1,-1));add("delta_F",c("beta_F","beta_FS"),c(1,-1))
    }
    contrasts<-as.matrix(contrasts);est<-drop(contrasts%*%beta)
    cv<-if(is.null(V))matrix(NA_real_,nrow(contrasts),nrow(contrasts))else contrasts%*%V%*%t(contrasts)
    se<-sqrt(pmax(0,diag(cv))); pv<-ifelse(is.finite(se)&se>0,2*stats::pnorm(-abs(est/se)),NA_real_)
    rows[[version]]<-data.frame(regression=version,scale="vanilla",parameter=names(est),estimate=unname(est),se=se,
      ci_low=unname(est)-stats::qnorm(.975)*se,ci_high=unname(est)+stats::qnorm(.975)*se,p_value=pv,
      N=nrow(d),status=fit$status,note=fit$note,method=fit$method,scope=fit$scope,
      outcome_scale=if(binary_outcome)"Probit latent response (residual SD = 1)" else "Input workflow outcome units",
      score_scale="Per 1 SD of observed role PGS",stringsAsFactors=FALSE)
    fit$N<-nrow(d);fit$n_subjects<-length(unique(iid));fit$n_families<-length(unique(fam));fit$n_clusters<-length(unique(fid))
    fit$rows_used<-which(keep);fit$centers<-centers;fit$design_scaling<-xs
    fit$family_random<-family_random&&!binary_outcome;fit$individual_random<-individual_random&&!binary_outcome
    fit$random_slope<-if(binary_outcome)NULL else random_slope
    fit$requested_random_effects<-list(family=family_random,individual=individual_random,slope=random_slope)
    fit$pgs_scaling<-do.call(rbind,scaling);fit$pgs_moments<-.ravel_pgs_moments(raw,iid,fid,version)
    fits[[version]]<-fit;scales[[version]]<-fit$pgs_scaling;moments[[version]]<-fit$pgs_moments
  }
  list(table=do.call(rbind,rows),fits=fits,pgs_scaling=do.call(rbind,scales),pgs_moments=do.call(rbind,moments))
}
