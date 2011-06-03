
latpos <- function(formula,data,subset,id,time,
                   unfold.method="Schoenemann",start=NULL,
                   sampler=mvt.sampler(df=7*length(latent.dims)),
                   ...){

  latent.dims <- all.vars(formula[c(1,3)])
  formula <- formula[1:2]
  formula <- update(formula,~.-1)
  manifest <- all.vars(formula)
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


  if(length(start)){
    if(!is.list(start)) stop("'start' argument must be a list")

    start <- start[c("A","beta","Sigma","Sigma0","Sigma1","Gamma","tau")]
    sdf.args <- c(
                  list(
                    latent=latent.dims,manifest=manifest
                  ),
                  start,
                  list(...)
                )
    start <- do.call(latpos.start.default,sdf.args)
  }

  if(!length(start)) start <- latpos.start.default(latent=latent.dims,manifest=manifest,
                                                                              ...)

  start <- start[nzchar(names(start)) & !is.na(names(start))]
  
  start <- latpos.start(resp=resp,latent.dims=latent.dims,manifest=manifest,
                        start=start,
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

  structure(fit,
    class="latpos"
  )
}


latpos.start.default <- function(latent,manifest,
                                 free.beta=FALSE,
                                 free.Sigma=c("full","diagonal","scale","none"),
                                 free.Sigma01=c("full","scaled"),
                                 free.Gamma=c("full","none","scale","diagonal"),
                                 A,
                                 beta,
                                 Sigma0,
                                 Sigma1,
                                 Gamma,
                                 tau,
                                 ...
                                 ){

  D <- length(latent)
  I <- length(manifest)

  free.Sigma <- match.arg(free.Sigma)
  free.Sigma01 <- match.arg(free.Sigma01)
  free.Gamma <- match.arg(free.Gamma)
  start <- list(
      free.beta=free.beta,
      free.Sigma=free.Sigma,
      free.Sigma01=free.Sigma01,
      free.Gamma=free.Gamma
    )

  if(!missing(A)){

    if(is.matrix(A)){

      if(nrow(A)==I && ncol(A)==D)
        start$A <- A
      else {
        start$A <- matrix(0,I,D,dimnames=list(manifest,latent))
        start$A[rownames(A),colnames(A)] <- A
      }
    }
    else if(is.list(A)){
        start$A <- matrix(0,I,D,dimnames=list(manifest,latent))

        for(l in names(A)){
          A.l <- A[[l]]
          start$A[names(A.l),l] <- A.l
        }
    }
    else {
      stop("unsported type of starting values for A")
    }
  }

  if(!missing(beta))
    start$beta <- beta
  else if(!free.beta)
    start$beta <- numeric(length=D)
  
  ## Matrix enforcing restrictions on Lambda = chol(Sigma)
  D.sq <- D*D

  if(D < 2){

    Q.kappa <- 1
  }
  else if(free.Sigma=="scale"){

    D1 <- D-1

    i <- 1:D1
    ij.1 <- i+D*(1:D1-1)
    ij.2 <- i+1+D*(2:D-1)

    Q.kappa <- matrix(0,nrow=D1,ncol=D.sq)
    Q.kappa[cbind(i,ij.1)] <- 1
    Q.kappa[cbind(i,ij.2)] <- -1

    Pat <- matrix(0,D,D)

    Qutri <- diag(x=as.numeric(upper.tri(Pat)))
    Qutri <- Qutri[diag(Qutri)>0,,drop=FALSE]
    Qltri <- diag(x=as.numeric(lower.tri(Pat)))
    Qltri <- Qltri[diag(Qltri)>0,,drop=FALSE]

    Q.kappa <- rbind(Q.kappa,Qutri,Qltri)

  }
  else if(free.Sigma=="diagonal"){

    Pat <- matrix(0,D,D)

    Qutri <- diag(x=as.numeric(upper.tri(Pat)))
    Qutri <- Qutri[diag(Qutri)>0,,drop=FALSE]
    Qltri <- diag(x=as.numeric(lower.tri(Pat)))
    Qltri <- Qltri[diag(Qltri)>0,,drop=FALSE]

    Q.kappa <- rbind(Qutri,Qltri)
  }
  else if(free.Sigma=="full"){

    Pat <- matrix(0,D,D)
    Qltri <- diag(x=as.numeric(lower.tri(Pat)))
    Q.kappa <- Qltri[diag(Qltri)>0,,drop=FALSE]
  }
  if(free.Sigma01=="full"){

    Q.kappa <- as.matrix(bdiag(Q.kappa,Q.kappa))
  }
  else {

    Q.01 <- diag(nrow=D.sq)
    Q.01 <- cbind(Q.01,-Q.01)

    Q.kappa <- rbind(
                    cbind(Q.kappa,matrix(0,nrow=nrow(Q.kappa),ncol=ncol(Q.kappa))),
                    Q.01)
  }
  Q.kappa <- restrictor(Q.kappa)$reduction
  Q.kappa0 <- Q.kappa[1:D.sq,,drop=FALSE]
  Q.kappa1 <- Q.kappa[D.sq + 1:D.sq,,drop=FALSE]
  
  start$Q.kappa <- Q.kappa
  start$Q.kappa0 <- Q.kappa0
  start$Q.kappa1 <- Q.kappa1

  if(!missing(Sigma0))
    start$Sigma0 <- Sigma0
  else if(free.Sigma=="none")
    start$Sigma0 <- diag(nrow=D)
    
  if(!missing(Sigma1))
    start$Sigma1 <- Sigma1
  else if(free.Sigma=="none")
    start$Sigma1 <- diag(nrow=D)
  else if(free.Sigma01 == "scaled")
    start$Sigma1 <- start$Sigma0
  

  if(!missing(tau))
    start$tau <- tau

  ## Matrix enforcing restrictions on Gamma
  if(D < 2){

    Q.rho <- 1
  }
  else if(free.Gamma=="scale"){

    D1 <- D-1

    i <- 1:D1
    ij.1 <- i+D*(1:D1-1)
    ij.2 <- i+1+D*(2:D-1)

    Q.rho <- matrix(0,nrow=D1,ncol=D.sq)
    Q.rho[cbind(i,ij.1)] <- 1
    Q.rho[cbind(i,ij.2)] <- -1

    Pat <- matrix(0,D,D)

    Qutri <- diag(x=as.numeric(upper.tri(Pat)))
    Qutri <- Qutri[diag(Qutri)>0,,drop=FALSE]
    Qltri <- diag(x=as.numeric(lower.tri(Pat)))
    Qltri <- Qltri[diag(Qltri)>0,,drop=FALSE]

    Q.rho <- rbind(Q.rho,Qutri,Qltri)

  }
  else if(free.Gamma=="diagonal"){

    Pat <- matrix(0,D,D)

    Qutri <- diag(x=as.numeric(upper.tri(Pat)))
    Qutri <- Qutri[diag(Qutri)>0,,drop=FALSE]
    Qltri <- diag(x=as.numeric(lower.tri(Pat)))
    Qltri <- Qltri[diag(Qltri)>0,,drop=FALSE]

    Q.rho <- rbind(Qutri,Qltri)
  }
  else if(free.Gamma=="full"){

    start$Q.rho <- diag(nrow=D*D)
  }

  if(free.Gamma == "diagonal"){

    Q.rho <- restrictor(Q.rho)$reduction
    start$Q.rho <- Q.rho
  }

  if(!missing(Gamma))
    start$Gamma <- Gamma
  else 
    start$Gamma <- diag(nrow=D)

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
                            force.increase=FALSE,
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
    force.increase=force.increase,
    ll.linesearch=ll.linesearch,
    Q.linesearch=Q.linesearch
    )
}

latpos.start <- function(resp,latent.dims,manifest,start,unfold.method,restrictions=standard.restrictions,...){

  I <- length(manifest)
  D <- length(latent.dims)

  y <- resp$y*resp$n
  start.diffs <- sqrt(-log(sweep(y+.5,2,colSums(y+.5),"/")))
  uf <- unfold(start.diffs,ndims=length(latent.dims),method=unfold.method,squared=TRUE)
  if("A" %in% names(start)){
      A <- start$A
      trans <- ProcrustesTrans(uf$B,A,translate=TRUE)
      A. <- sweep(uf$B%*%trans$A,2,trans$b,"-")
      B <- sweep(uf$A%*%trans$A,2,trans$b,"-")
  }
  else {
      A <- uf$B
      B <- uf$A
  }
  colnames(A) <- latent.dims
  rownames(A) <- manifest
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
  A.restrictor <- tryCatch(restrictor(C=rest.C,d=rest.d),
                           error=function(e) stop("Your model does not seem to be identified. Reduce the number of latent dimensions or set more restrictions.",
                                                  call.=FALSE))
  Q.phi <- A.restrictor$reduction
  r.phi <- A.restrictor$offset

  transf <- rotate.to.restriction(X=A,C=restrictions$C,d=restrictions$d)
  dimnA <- dimnames(A)
  A <- transf$transformed
  rot <- transf$rotation
  transl <- transf$translation
  dimnames(A) <- dimnA

  phi <- crossprod(Q.phi,as.vector(A)-r.phi)
  A[] <- Q.phi%*%phi + r.phi
  B <- sweep(B%*%rot,2,transl,"-")
   
  start$A <- A
  start$phi <- phi
  start$Q.phi <- Q.phi
  start$r.phi <- r.phi

  U <- scale(B,scale=FALSE)
  beta <- attr(U,"scaled:center")
  if(!length(start$beta)) start$beta <- beta

  Gamma <- start$Gamma

  ## Enforce restrictions on Gamma

  if(start$free.Gamma=="none") l.rho <- 0
  else {

    Gamma <- start$Gamma
    Q.rho <- start$Q.rho
    rho <- crossprod(Q.rho,as.vector(Gamma))
    l.rho <- length(rho)
  }

  if(length(start$Sigma0))
    Sigma0 <- start$Sigma0
  else 
    Sigma0 <- cov(U[resp$t0,,drop=FALSE])
  if(length(start$Sigma1))
    Sigma1 <- start$Sigma1
  else {
    if(start$free.Sigma01=="scaled")
      Sigma1 <- Sigma0
    else 
      Sigma1 <- cov(U[resp$s1,,drop=FALSE]-tcrossprod(U[resp$s0,,drop=FALSE],Gamma))
  }


  if(length(start$tau)) tau <- start$tau
  else if(start$free.Sigma01=="scaled"){
    tau <- solve(cov(U[resp$s1,,drop=FALSE]-tcrossprod(U[resp$s0,,drop=FALSE],Gamma)),Sigma0)
    tau <- sqrt(mean(diag(as.matrix(tau))))
  }
  else tau <- 1

  start$Sigma0 <- Sigma0
  start$Sigma1 <- Sigma1

  start$Gamma <- Gamma

  bT1 <- sum(Tj1)
  bT <- sum(Tj)

  bS00 <- crossprod(U[resp$t0,,drop=FALSE])
  bS11 <- crossprod(U[resp$s0,,drop=FALSE])
  bS12 <- crossprod(U[resp$s0,,drop=FALSE],U[resp$s1,,drop=FALSE])
  bS21 <- t(bS12)
  bS22 <- crossprod(U[resp$s1,,drop=FALSE])

  start$Tj <- Tj
  start$Tj1 <- Tj1
  start$latent.dims <- latent.dims

  start <- varParMax(start,S22=bS22,S12=bS12,S11=bS11,S00=bS00)

  start$Utilde <- U


  start
}

latpos.fit <- function(resp,start,
                            sampler=mvnorm.sampler(),
                            control=latpos.control()
          ){
## Fits a latent-positions state-space model with
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

  if(start$free.Sigma!="none"){

    Lambda0 <- chol(start$Sigma0)
    Lambda1 <- chol(start$Sigma1)
    Q.kappa0 <- start$Q.kappa0
    Q.kappa1 <- start$Q.kappa1

    kappa <- crossprod(Q.kappa0,as.vector(Lambda0)) + crossprod(Q.kappa1,as.vector(Lambda1))
  }

  if(start$free.Gamma!="none"){

    Gamma <- start$Gamma
    Q.rho <- start$Q.rho
    rho <- crossprod(Q.rho,as.vector(Gamma))
  }

  psi <- parm2psi(start)

  trace <- list(Q = matrix(NA,nrow=maxiter,ncol=8),
                psi = matrix(NA,nrow=maxiter,ncol=length(psi))
            )

  parm <- start

  print(parm$A)

  ### The Iterations ###################
  cat("\n\n\n****** Starting Monte Carlo EM algorithm ******")
  cat("\nGenerating Monte Carlo sample(s) - size",initial.size,"- for the first step")
  parm$Utilde <- latpos.utilde(resp=resp,parm=parm,maxiter=maxiter,verbose=FALSE)
  latent.data <- latpos.simul(resp=resp,parm=parm,
                        latent.data=list(sample.size=initial.size),
                        sampler=sampler)

  parm$logLik <- sum(latent.data$ll.j)
  parm$sample.size <- latent.data$sample.size

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
    cat(paste("\n",format(Sys.time(),usetz=TRUE),"\n",sep=""))
    cat("A:\n")
    print(parm$A)
    cat("Sigma0:\n")
    print(parm$Sigma0)
    cat("Sigma1:\n")
    print(parm$Sigma1/parm$tau^2)
    cat("Gamma:\n")
    print(parm$Gamma)
    if(converged) break

  }
  parm$zeta <- 1/parm$tau

  GradInfo.Abeta <- latpos.GradInfo_Abeta(resp=resp,parm=parm,latent.data=latent.data)
  GradInfo.VarPar <- latpos.GradInfo_VarPar(resp=resp,parm=parm,latent.data=latent.data)

  parm$Information <- list(Abeta=GradInfo.Abeta$Information,
                           VarPar=GradInfo.VarPar$Information)

  parm$Info.cpl <- list(Abeta=GradInfo.Abeta$Info.cpl,
                        VarPar=GradInfo.VarPar$Info.cpl)

  parm$Info.miss <- list(Abeta=GradInfo.Abeta$Info.miss,
                        VarPar=GradInfo.VarPar$Info.miss)

  parm$Info.restr <- list(Abeta=GradInfo.Abeta$restrictor,
                        VarPar=GradInfo.VarPar$restrictor)

  parm$gradient <- list(Abeta=GradInfo.Abeta$gradient,
                        VarPar=GradInfo.VarPar$gradient)
  
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

  cat("\nSpatial model of latent positions\n")

  print.default(x$call)

  cat("\nPositions of objectives:\n")
  print.default(x$parm$A)
  cat("\nManifesto parameters:\n")
  cat("\nMean positions (beta):\n")
  print.default(x$parm$beta)
  cat("\nAutoregression coefficients (Gamma):\n")
  print.default(x$parm$Gamma)
  cat("\nStarting position variance:\n")
  print.default(x$parm$Sigma0)
  cat("\nPosition change variance:\n")
  if(x$parm$free.Sigma01=="scaled")
    print.default(x$parm$Sigma1/x$parm$tau^2)
  else
    print.default(x$parm$Sigma1)
  invisible(x)
}



