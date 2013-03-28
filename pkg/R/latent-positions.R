
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

  parm.names <- c("A","beta","Sigma","Sigma0","Sigma1","Gamma","tau")

  dots <- list(...)
  if(any(parm.names %in% names(dots))){
    if(!length(start))
      start <- dots[parm.names]
    else
      start <- c(start,dots[setdiff(parm.names,names(start))])
  }


  if(length(start)){
    if(!is.list(start)) stop("'start' argument must be a list")

    start <- start[parm.names]
    sdf.args <- c(
                  list(
                    latent=latent.dims,manifest=manifest
                  ),
                  start,
                  dots
                )
    start <- do.call(latpos.start.default,sdf.args)
  }
  else
    start <- latpos.start.default(latent=latent.dims,manifest=manifest,...)

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
  fit$sampler <- sampler

  structure(fit,
    class="latpos"
  )
}


latpos.start.default <- function(latent,manifest,
                                 free.beta=TRUE,
                                 free.Gamma=TRUE,
                                 A,
                                 beta,
                                 Gamma,
                                 ...
                                 ){

  D <- length(latent)
  I <- length(manifest)

  start <- list(
      free.beta=free.beta,
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
      stop("unsupported type of starting values for A")
    }
  }

  if(!missing(beta))
    start$beta <- beta
  else if(!free.beta)
    start$beta <- numeric(length=D)
  

  ## Matrix enforcing restrictions on Gamma -- so
  ## far none implemented ...
  if(D < 2){

    start$Q.rho <- as.matrix(1)
  }
  else {

    start$Q.rho <- diag(D*D)

  }


  if(!missing(Gamma))
    start$Gamma <- Gamma
  else 
    start$Gamma <- diag(nrow=D)

  return(start)
}

latpos.control <- function(maxiter=200,
                            initial.size=101,
                            Lambda.alpha=.05,
                            Lambda.eps=1e-7,
                            diff.logLik.eps=1e-7,
                            abs.diff.psi.eps=0,
                            rel.diff.psi.eps=0,
                            max.size=Inf,
                            min.final.size=1000,
                            force.increase=TRUE,
                            Q.linesearch=TRUE,
                            ...){
  list(
    maxiter=maxiter,
    initial.size=initial.size,
    Lambda.alpha=Lambda.alpha,
    Lambda.eps=Lambda.eps,
    diff.logLik.eps=diff.logLik.eps,
    abs.diff.psi.eps=abs.diff.psi.eps,
    rel.diff.psi.eps=rel.diff.psi.eps,
    max.size=max.size,
    min.final.size=min.final.size,
    force.increase=force.increase,
    Q.linesearch=Q.linesearch
    )
}

latpos.start <- function(resp,latent.dims,manifest,start,unfold.method,restrictions=standard.restrictions,maxiter,...){

  I <- length(manifest)
  D <- length(latent.dims)

  y <- resp$y*resp$n
  start.diffs <- sqrt(-log(sweep(y+.5,2,colSums(y+.5),"/")))
  uf <- unfold(start.diffs,ndims=start$ndims,method=unfold.method,squared=TRUE)
  start$uf <- uf
  if("A" %in% names(start)){
      A <- start$A
      trans <- ProcrustesTrans(uf$B,A,translate=TRUE)
      A <- sweep(uf$B%*%trans$A,2,trans$b,"-")
      B <- sweep(uf$A%*%trans$A,2,trans$b,"-")
  }
  else {
      A <- uf$B[,1:D,drop=FALSE]
      B <- uf$A[,1:D,drop=FALSE]
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
  else stop("no support for restrictions of type",typeof(restrictions))

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

  beta <- colMeans(B)
  
  if(!length(start$beta)) start$beta <- beta

  Gamma <- start$Gamma

  ## Enforce restrictions on Gamma

  if(!start$free.Gamma) l.rho <- 0
  else {

    Gamma <- start$Gamma
    Q.rho <- start$Q.rho
    rho <- crossprod(Q.rho,as.vector(Gamma))
    l.rho <- length(rho)
  }

  start$Gamma <- Gamma

  ## Enforce symmetry restrictions on Sigma0 and Sigma1

  if(D>1){

    Q.kappa <- restr.to.symm(D)
    start$Q.kappa0 <- start$Q.kappa1 <- Q.kappa$reduction
  }
  else {

    start$Q.kappa0 <- start$Q.kappa1 <- as.matrix(1)
  }
  
  ## Get starting values
  
  bT1 <- sum(Tj1)
  bT <- sum(Tj)

  start$Tj <- Tj
  start$Tj1 <- Tj1
  start$latent.dims <- latent.dims

  bb0 <- colSums(B[resp$t0,,drop=FALSE])
  bs1 <- colSums(B[resp$s0,,drop=FALSE])
  bs2 <- colSums(B[resp$s1,,drop=FALSE])

  bS00 <- crossprod(B[resp$t0,,drop=FALSE])
  bS11 <- crossprod(B[resp$s0,,drop=FALSE])
  bS12 <- crossprod(B[resp$s0,,drop=FALSE],B[resp$s1,,drop=FALSE])
  bS22 <- crossprod(B[resp$s1,,drop=FALSE])

  start <- LVdistMax(start,
                      b0=bb0,s1=bs1,s2=bs2,
                      S22=bS22,S12=bS12,S11=bS11,S00=bS00,maxiter=100,verbose=TRUE,starting=TRUE)
  cat("\n")

  start$Btilde <- list(B=B)

  start
}

latpos.fit <- function(resp,start,
                            sampler=mvt.sampler(df=27),
                            control=latpos.control()
          ){
## Fits a latent-positions state-space model with
## conditional multinomial distribution of policy objective emphases
## using Caffo, Jank & Jones' (2005) 'Ascend-based Monte Carlo expectation-maximization'
## JRSS B 67: 235-251

  maxiter           <- control$maxiter
  initial.size      <- control$initial.size
  Lambda.alpha      <- control$Lambda.alpha
  Lambda.eps        <- control$Lambda.eps
  diff.logLik.eps   <- control$diff.logLik.eps
  abs.diff.psi.eps  <- control$abs.diff.psi.eps
  rel.diff.psi.eps  <- control$rel.diff.psi.eps
  max.size          <- control$max.size
  min.final.size    <- control$min.final.size
  Q.linesearch  <- control$Q.linesearch


  stopifnot(ncol(start$A)==ncol(start$U))

  latent.dims <- start$latent.dims
  ndims <- length(latent.dims)

  ##
  sample.size <- initial.size

  psi <- parm2psi(start)

  trace <- list(Q = matrix(NA,nrow=maxiter,ncol=8),
                psi = matrix(NA,nrow=maxiter,ncol=length(psi))
            )

  parm <- start

  print(parm$A)

  ### The Iterations ###################
  cat("\n\n\n****** Starting Monte Carlo EM algorithm ******")
  cat("\nGenerating Monte Carlo sample(s) - size",initial.size,"- for the first step")
  parm$Btilde <- latpos.Btilde(resp=resp,parm=parm,maxiter=maxiter,verbose=FALSE)
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
    cat("beta:\n")
    print(parm$beta)
    cat("Sigma0:\n")
    print(parm$Sigma0)
    cat("Sigma1:\n")
    print(parm$Sigma1)
    cat("Gamma:\n")
    print(parm$Gamma)
    if(converged) break

  }

  CplInfo.phi <- latpos.CplInfo_phi(resp=resp,parm=parm,latent.data=latent.data)
  CplInfo.LVdist <- latpos.CplInfo_LVdist(resp=resp,parm=parm,latent.data=latent.data)
  MissInfo <- latpos.missinfo(resp=resp,parm=parm,latent.data=latent.data)

  
  Info.cpl <- bdiag(CplInfo.phi$Information,
                       CplInfo.LVdist$Information)
  Info.miss <- MissInfo$var.gradient

  Info.restr <- bdiag(CplInfo.phi$restrictor,
                      CplInfo.LVdist$restrictor)

  Info.obs <- Info.cpl - Info.miss

  parm$covmat <- Info.restr%*%tcrossprod(solve(Info.obs),Info.restr)

  parm$Information <- list( complete.data=Info.cpl,
                            missing=Info.miss,
                            observed.data=Info.obs,
                            restrictor=Info.restr
                            )
  
  list(
    resp=resp,parm=parm,
    sample.size=latent.data$sample.size,
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
  cat("\nStarting position variance (Sigma0):\n")
  print.default(x$parm$Sigma0)
  cat("\nPosition change variance (Sigma1):\n")
  print.default(x$parm$Sigma1)
  invisible(x)
}



