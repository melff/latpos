
latpos <- function(formula,data,subset,id,time,
                   unfold.method="Schoenemann",start=NULL,
                   sampler=mvnorm.sampler(),
                   ...){

  latent.dims <- all.vars(formula[c(1,3)])
  formula <- formula[1:2]
  formula <- update(formula,~.-1)
  cl <- match.call()
  m <- match.call(expand.dots=FALSE)
  mf <- m[c(1,match(c("formula","data","subset"),
          names(m),nomatch=0L))]
  mf[[1]] <- as.name("model.frame")
  mf$formula <- formula
  mf <- eval(mf,parent.frame())

  if(missing(id)) stop("argument 'id' is missing, with no default")
  if(missing(time)) stop("argument 'time' is missing, with no default")
  j <- eval(substitute(id),data,parent.frame())
  t <- eval(substitute(time),data,parent.frame())
  j <- match(j,unique(j))
  split(t,j) <- lapply(split(t,j),function(tj)match(tj,sort(unique(tj)))-1L)
  jt.order <- order(j,t)
  orig.order <- (1:nrow(mf))[jt.order]

  y <- t(model.matrix(formula,mf))[,jt.order]
  j <- j[jt.order]
  t <- t[jt.order]

  resp <- list()
  resp$n <- array(rep(colSums(y),each=nrow(y)),dim=dim(y))
  resp$y <- y/resp$n
  resp$y[resp$n==0] <- 0
  resp$j <- j
  resp$t0 <- t==0
  resp$t <- t
  s <- seq_along(t)
  resp$s0 <- ifelse(t==0,0,s-1)
  resp$s1 <- ifelse(t==0,0,s)

  if(!length(start)) start <- latpos.start.default(I=nrow(y),D=length(latent.dims),
                                                                              ...)

  start <- start[nzchar(names(start)) & !is.na(names(start))]
  
  start <- latpos.start(resp=resp,latent.dims=latent.dims,start=start,
                        unfold.method=unfold.method,...)
  
  fit <- latpos.fit(resp=resp,
                    start=start,
                    sampler=sampler,
                    control=latpos.control(...))
  fit$call <- cl
  fit$data <- data
  fit$jt.order <- jt.order
  fit$orig.order <- orig.order
  fit$start <- start

  fit$parm$covmat <- vcov.latpos(fit)
  structure(fit,
    class="latpos"
  )
}


latpos.start.default <- function(I,D,
                                 free.beta=FALSE,
                                 free.Sigma=c("diagonal","scale","full","none"),
                                 free.rho=FALSE,
                                 A,
                                 beta,
                                 Sigma,
                                 rho,
                                 zeta,
                                 Q.gamma,
                                 ...
                                 ){

  free.Sigma <- match.arg(free.Sigma)
  start <- list(
      free.beta=free.beta,
      free.Sigma=free.Sigma,
      free.rho=free.rho
    )

  if(!missing(beta))
    start$beta <- beta
  else if(!free.beta)
    start$beta <- numeric(length=D)
  
  ## Matrix enforcing restrictions on Gamma
  if(missing(Q.gamma)){

    if(D < 2){

      Q.gamma <- 1
    }
    else if(free.Sigma=="scale"){

      D1 <- D-1

      i <- 1:D1
      ij.1 <- i+D*(1:D1-1)
      ij.2 <- i+1+D*(2:D-1)

      Q.gamma <- matrix(0,nrow=D1,ncol=D^2)
      Q.gamma[cbind(i,ij.1)] <- 1
      Q.gamma[cbind(i,ij.2)] <- -1
      Q.gamma <- restrictor(Q.gamma)$reduction[,1,drop=FALSE]
    }
    else if(free.Sigma=="diagonal"){

      GammaPat <- matrix(0,D,D)

      Qutri.gamma <- diag(x=as.numeric(upper.tri(GammaPat)))
      Qutri.gamma <- Qutri.gamma[diag(Qutri.gamma)>0,,drop=FALSE]
      Qltri.gamma <- diag(x=as.numeric(lower.tri(GammaPat)))
      Qltri.gamma <- Qltri.gamma[diag(Qltri.gamma)>0,,drop=FALSE]

      Q.gamma <- rbind(Qutri.gamma,Qltri.gamma)
      Q.gamma <- restrictor(Q.gamma)$reduction
    }
    else if(free.Sigma=="full"){

      GammaPat <- matrix(0,D,D)
      Qltri.gamma <- diag(x=as.numeric(lower.tri(GammaPat)))
      Qltri.gamma <- Qltri.gamma[diag(Qltri.gamma)>0,,drop=FALSE]
      Q.gamma <- restrictor(Qltri.gamma)$reduction
    }
  }
  start$Q.gamma <- Q.gamma
  if(free.Sigma=="none"){

    if(!missing(Sigma)){

      start$Sigma <- Sigma
      start$Gamma <- chol(Sigma)
      start$Theta <- chol2inv(start$Gamma)
    }
    else {
      
      start$Sigma <- diag(nrow=D)
      start$Theta <- diag(nrow=D)
    }
  }
  else if(!missing(Sigma)){
    Gamma <- chol(Sigma)
    gamma <- crossprod(Q.gamma,as.vector(Gamma))
    Gamma[] <- Q.gamma %*% gamma
    start$gamma <- gamma
    start$Theta <- chol2inv(Gamma)
    start$Sigma <- crossprod(Gamma)
  }

  SigmaPat <- matrix(0,D,D)
  if(D<2)
    Q.sigma <- 1
  else if(free.Sigma=="full"){

    DD2 <- D*(D-1)/2
    i <- row(SigmaPat)[lower.tri(SigmaPat)]
    j <- col(SigmaPat)[lower.tri(SigmaPat)]
    ij <- i + (j-1)*D
    ji <- j + (i-1)*D
    ii <- 1:DD2

    Q.sigma <- matrix(0,nrow=DD2,ncol=D^2)
    Q.sigma[cbind(ii,ij)] <- 1
    Q.sigma[cbind(ii,ji)] <- -1
    Q.sigma <- restrictor(Q.sigma)$reduction

  } else if(free.Sigma=="diagonal"){

    Qutri.sigma <- diag(x=as.numeric(upper.tri(SigmaPat)))
    Qutri.sigma <- Qutri.sigma[diag(Qutri.sigma)>0,,drop=FALSE]
    Qltri.sigma <- diag(x=as.numeric(lower.tri(SigmaPat)))
    Qltri.sigma <- Qltri.sigma[diag(Qltri.sigma)>0,,drop=FALSE]

    Q.sigma <- rbind(Qutri.sigma,Qltri.sigma)
    Q.sigma <- restrictor(Q.sigma)$reduction

  } else { #if(free.Sigma)=="scale"

    D1 <- D-1

    i <- 1:D1
    ij.1 <- i+D*(1:D1-1)
    ij.2 <- i+1+D*(2:D-1)

    Q.sigma <- matrix(0,nrow=D1,ncol=D^2)
    Q.sigma[cbind(i,ij.1)] <- 1
    Q.sigma[cbind(i,ij.2)] <- -1
    Q.sigma <- restrictor(Q.sigma)$reduction[,1,drop=FALSE]
  }
  start$Q.sigma <- Q.sigma

  if(!missing(rho))
    start$rho <- rho
  else if(!start$free.rho)
    start$rho <- 1

  if(!missing(zeta))
    start$tau <- 1/zeta

  return(start)
}

latpos.control <- function(maxiter=200,
                            initial.size=101,
                            diff.Q.alpha=.05,
                            diff.Q.eps=1e-7,
                            diff.logLik.eps=1e-7,
                            abs.diff.psi.eps=0,
                            rel.diff.psi.eps=0,
                            max.size=Inf,
                            min.final.size=1000,
                            sparsity.eps=0,
                            ll.linesearch=FALSE,
                            Q.linesearch=FALSE,
                            ...){
  list(
    maxiter=maxiter,
    initial.size=initial.size,
    diff.Q.alpha=diff.Q.alpha,
    diff.Q.eps=diff.Q.eps,
    diff.logLik.eps=diff.logLik.eps,
    abs.diff.psi.eps=abs.diff.psi.eps,
    rel.diff.psi.eps=rel.diff.psi.eps,
    max.size=max.size,
    min.final.size=min.final.size,
    sparsity.eps=sparsity.eps,
    ll.linesearch=ll.linesearch,
    Q.linesearch=Q.linesearch
    )
}

latpos.start <- function(resp,latent.dims,start,unfold.method,restrictions=standard.restrictions,...){

  I <- nrow(resp$y)
  D <- length(latent.dims)

  y <- resp$y*resp$n
  start.diffs <- sqrt(-log(sweep(y+.5,2,colSums(y+.5),"/")))
  uf <- unfold(start.diffs,ndims=length(latent.dims),method=unfold.method,squared=TRUE)
  if("A" %in% names(start))
      A <- start$A
  else
      A <- uf$B
  B <- uf$A
  colnames(A) <- latent.dims
  colnames(B) <- latent.dims

  D <- ncol(A)
  I <- nrow(A)
  JT <- nrow(B)
  Tj1 <- tabulate(resp$j)
  Tj <- Tj1 - 1
  J <- length(Tj1)
  
  if(is.function(restrictions)) restrictions <- restrictions(A)
  else if(is.list(restrictions)) restrictions <- restrictions[c("C","d")]
  else stop("no support of restrictions with type",typeof(restrictions))
  ## Matrix enforcing the linear restrictions
  ## on A
  rest.C <- restrictions$C
  rest.d <- if(length(restrictions$d)) restrictions$d else numeric(nrow(rest.C))
  A.restrictor <- restrictor(C=rest.C,d=rest.d)
  Q.phi <- A.restrictor$reduction
  kappa.phi <- A.restrictor$offset

  transf <- rotate.to.restriction(X=A,C=restrictions$C,d=restrictions$d)
  A <- transf$transformed
  rot <- transf$rotation
  transl <- transf$translation

  phi <- crossprod(Q.phi,as.vector(A)-kappa.phi)
  A[] <- Q.phi%*%phi + kappa.phi
  B <- sweep(B%*%rot,2,transl,"-")
  
  start$A <- A
  start$phi <- phi
  start$Q.phi <- Q.phi
  start$kappa.phi <- kappa.phi

  Q.gamma <- start$Q.gamma

  U <- scale(B,scale=FALSE)
  beta <- attr(U,"scaled:center")
  if(!length(start$beta)) start$beta <- beta
  
  bS0 <- crossprod(U[resp$t0,,drop=FALSE])
  bS1 <- crossprod(U[resp$s0,,drop=FALSE])
  bS2 <- crossprod(U[resp$s0,,drop=FALSE],U[resp$s1,,drop=FALSE])
  bS3 <- crossprod(U[resp$s1,,drop=FALSE])

  if(length(start$Sigma)){
    Sigma <- start$Sigma
    Gamma <- start$Gamma
    mS0 <- bS0/sum(resp$t0)
    imS0 <- solve(mS0)
    L0 <- chol(imS0) %*% Gamma
    U <- U %*% L0
    }
  else Gamma <- chol(cov(U))

  if(length(start$rho)) rho <- start$rho
  else rho <- 1

  if(length(start$tau)) tau <- start$tau
  else tau <- 2

  bT1 <- sum(Tj1)
  bT <- sum(Tj)

  if(start$free.rho){
    par <- rho
    l.rho <- length(rho)
    i.rho <- 1:l.rho
  }
  else{
    par <- numeric(0)
    l.rho <- 0
    i.rho <- 0
  }
  if(start$free.Sigma!="none"){
      gamma <- crossprod(Q.gamma,as.vector(Gamma))
      par <- c(par,gamma)
      l.gamma <- length(gamma)
      i.gamma <- 1:l.gamma
  }
  else {

    l.gamma <- 0
    i.gamma <- 0
  }
  par <- c(par,tau)

  searchFun <- function(par){

      if(l.rho)
        rho <- par[i.rho]
      else
        rho <- 1
      if(l.gamma){
        gamma <- par[l.rho + i.gamma]
        Gamma[] <- Q.gamma %*% gamma
        Sigma <- crossprod(Gamma)
      }
      tau <- par[l.rho + l.gamma + 1]
      Theta <- chol2inv(Gamma)

      pen <- 0
      if(tau < 1) pen <- 1#pen + (tau - 1)^2
      if(rho < 0) pen <- 1#pen + rho^2
      if(rho > 1) pen <- 1#pen + (rho - 1)^2
      
      bV <- bS3 - 2*rho*bS2 + rho^2*bS1
      suppressWarnings(log.tau <- log(tau))

      ssq <- sum(diag((bS0+tau*bV)%*%Theta))
      logDet <- D*bT*log.tau + bT1*logdet(Theta)
      
      ll <- -ssq/2 + logDet/2
      #cat("\nll =",ll)
      -ll + pen*abs(ll)^2
  }
  opt.res <- optim(par,searchFun,method="BFGS")
  pseudo.ll <- -opt.res$value
  par <- opt.res$par
  if(l.rho)
    rho <- par[i.rho]
  else
    rho <- 1
  if(l.gamma){
    gamma <- par[l.rho + i.gamma]
    Gamma[] <- Q.gamma %*% gamma
    Sigma <- crossprod(Gamma)
  }
  tau <- par[l.rho + l.gamma + 1]


  if(!length(start$rho)) start$rho <- rho
  if(!length(start$tau)) start$tau <- tau
  if(!length(start$Sigma)) {
    start$gamma <- gamma
    start$Sigma <- Sigma
    start$Gamma <- Gamma
    start$Theta <- chol2inv(Gamma)
  }
  
  start$Utilde <- U
  

  start$Tj <- Tj
  start$Tj1 <- Tj1
  start$latent.dims <- latent.dims

  start
}

latpos.fit <- function(resp,start,
                            sampler=mvnorm.sampler(),
                            control=latpos.control()
          ){
## Fits a latent-positions martingale model with
## conditional multinomial distribution of policy objective emphases
## using Caffo, Jank & Jones' (2005) 'Ascend-based Monte Carlo expectation-maximization'
## JRSS B 67: 235-251

  maxiter           <- control$maxiter
  initial.size      <- control$initial.size
  diff.Q.alpha      <- control$diff.Q.alpha
  diff.Q.eps        <- control$diff.Q.eps
  diff.logLik.eps   <- control$diff.logLik.eps
  abs.diff.psi.eps  <- control$abs.diff.psi.eps
  rel.diff.psi.eps  <- control$rel.diff.psi.eps
  max.size          <- control$max.size
  min.final.size    <- control$min.final.size
  sparsity.eps      <- control$sparsity.eps
  ll.linesearch <- control$ll.linesearch
  Q.linesearch  <- control$Q.linesearch

  EM.maxiter <- control$EM.maxiter
  NR.threshold <- control$NR.threshold
  

  stopifnot(ncol(start$A)==ncol(start$U))

  latent.dims <- start$latent.dims
  ndims <- length(latent.dims)

  ## Deal with empirical zeroes

  if(sparsity.eps > 0){
#browser()
    y <- resp$y*resp$n
    nonnull <- which(colSums(resp$y)>0)
    y[,nonnull] <- y[,nonnull] + sparsity.eps
    n <- array(rep(colSums(y),each=nrow(resp$y)),dim=dim(y))
    resp$y[,nonnull] <- y[,nonnull]/n[,nonnull]
    resp$n <- n
    rm(y)
  }
  ##
  sample.size <- initial.size

  psi <- start$phi
  if(start$free.beta) psi <- c(psi,start$beta)
  if(start$free.Sigma!="none") psi <- c(psi,vech(start$Sigma))
  if(start$free.rho) psi <- c(psi,start$rho)
  psi <- c(psi,1/start$tau)

  trace <- list(Q = matrix(NA,nrow=maxiter,ncol=8),
                psi = matrix(NA,nrow=maxiter,ncol=length(psi))
            )

  parm <- start

  ### The Iterations ###################
  cat("\n\n\n****** Starting Monte Carlo EM algorithm ******")
  cat("\nGenerating Monte Carlo sample(s) - size",initial.size,"- for the first step")
  parm$Utilde <- latpos.utilde(resp=resp,parm=parm,maxiter=maxiter,verbose=FALSE)
  latent.data <- latpos.simul(resp=resp,parm=parm,
                        latent.data=list(sample.size=initial.size),
                        sampler=sampler)

  parm$logLik <- sum(latent.data$ll.j)

  use.NR <- FALSE

  for(iter in 1:maxiter){

    step.res <- latpos.MCEMstep(resp=resp,
                                  parm=parm,
                                  latent.data=latent.data,
                                  trace=trace,
                                  iter=iter,
                                  sampler=sampler,
                                  control=control)

    parm <- step.res$parm
    latent.data <- step.res$latent.data
    trace <- step.res$trace
    converged <- step.res$converged

    if(converged) break

  }

  GradInfo.Abeta <- latpos.GradInfo_Abeta(resp=resp,parm=parm,latent.data=latent.data)
  GradInfo.VarPar <- latpos.GradInfo_VarPar(resp=resp,parm=parm,latent.data=latent.data)

  Info <- as.matrix(bdiag(GradInfo.Abeta$Information,GradInfo.VarPar$Information))
  covmat <- solve(Info)
  Qmat <- as.matrix(bdiag(GradInfo.Abeta$restrictor,GradInfo.VarPar$restrictor))
  
  covmat <- Qmat %*% tcrossprod(covmat,Qmat)
  parm$covmat <- covmat

  list(
    resp=resp,parm=parm,
    latent.data=latent.data,
    trace=trace,
    control=control,
    sampler=sampler,
    converged=converged
    )
}





logdet <- function(x) determinant(x,logarithm=TRUE)$modulus




vech <- function(M){

  stopifnot(nrow(M)==ncol(M))
  M[lower.tri(M,diag=TRUE)]
}

"vech<-" <- function(M,value){

  stopifnot(nrow(M)==ncol(M))
  M[lower.tri(M,diag=TRUE)] <- value
  M
}

symmFill <- function(D,fill){

  M <- matrix(NA,D,D)
  M[lower.tri(M,diag=TRUE)] <- fill
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  M
}


print.latpos <- function(x,...){

  cat("\nSpatial model of latent positions")

  print.default(x$call)

  cat("\nPositions of objectives:\n")
  print.default(x$parm$A)
  cat("\nManifesto parameters:\n")
  cat("\nMean positions (beta):\n")
  print.default(x$parm$beta)
  cat("\nAutoregression coefficient (rho):\n")
  print.default(x$parm$rho)
  cat("\nSigma:\n")
  print.default(x$parm$Sigma)
  cat("\nzeta:\n")
  print.default(1/x$parm$tau)
  invisible(x)
}


vcov.latpos <- function(object,...) return(object$parm$covmat)
