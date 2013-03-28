latpos.eval.parms <- function(y,n,j,t,parm,B,compute,weights){

  ndim <- length(parm$latent.dims)
  #nmiss.weights <- !missing(weights)
  if(missing(weights)) weights <- rep(1,nrow(B))
  else weights <- as.vector(weights)

  A <- parm$A
  U <- sweep(B,2,parm$beta,"-")

  Gamma <- parm$Gamma
  Sigma0 <- parm$Sigma0
  Sigma1 <- parm$Sigma1
  Lambda0 <- chol(Sigma0)
  Lambda1 <- chol(Sigma1)
  Theta0 <- chol2inv(Lambda0)
  Theta1 <- chol2inv(Lambda1)

  ssq <- numeric(nrow(U))

  if(any(c("deviance","logLik.j") %in% compute)){

    s <- 1:nrow(B)
    s0 <- ifelse(t==0,0,s-1)
    s1 <- ifelse(t==0,0,s)
    t0 <- which(t==0)
    tnon0 <- t > 0

    w0 <- weights[t0]
    wnon0 <- weights[tnon0]
    bT0 <- sum(w0)
    bT <- sum(wnon0)
    bT1 <- bT + bT0

    U0 <- U[t0,,drop=FALSE]
    diffU <- U[s1,,drop=FALSE] - tcrossprod(U[s0,,drop=FALSE],parm$Gamma)

    ssq <- 0
    ssq[t0] <- rowSums(w0*U0*(U0%*%Theta0))
    ssq[tnon0] <- rowSums(wnon0*diffU*(diffU%*%Theta1))

    logdet.Theta0 <- -2*sum(log(diag(Lambda0)))
    logdet.Theta1 <- -2*sum(log(diag(Lambda1)))

    logDet <- 0
    logDet[t0] <- w0*logdet.Theta0
    logDet[tnon0] <- wnon0*logdet.Theta1

  }

  if(any(c("deviance","p","logLik.j","XWX","XWy","g.j") %in% compute)){

    p <- latpos_p(A,B)
    ll <- ll_p(p,y,n,weights)
  }

  res <- list()
  if("deviance" %in% compute){

    dev.resid <- ll_p(y,y,n,weights) - ll
    dev.resid[y==0] <- 0
    res$deviance <- 2*sum(dev.resid) + sum(ssq) - sum(logDet)
  }
  if("p" %in% compute)
    res$p <- p

  if("logLik.j" %in% compute){

    log.constpart <- (lfactorial(n[1,]) - colSums(lfactorial(y*n)))*as.vector(weights)

    ll <- colSums(ll) + log.constpart - ssq/2 + logDet/2
    res$logLik.j <- drop(rowsum(ll,j))

  }

  if(any(c("XWX","XWy","g.j") %in% compute)){

    X <- d.eta.d.phi(A,B,parm$Q.phi)

    jtk <- rep(1:ncol(y),each=nrow(y))

  }

  if(any(c("XWX","XWy") %in% compute)){
    res$XWX <- latpos_XWX(X,p,n,weights)
  }
  if("XWy" %in% compute || "g.j" %in% compute)
    r <- latpos_resid(p,y,n,weights)
  
  if("XWy" %in% compute){

    Xr <- crossprod(X,r)
    phi <- parm$phi
    res$XWy <- res$XWX%*%phi + Xr
  }
  if("g.j" %in% compute){

    ji <- rep(j,each=nrow(y))
    res$g.j <- rowsum(X*r,ji)
  }

  res
}


latpos.integ.ll <- function(resp,parm,latent.data,compute){

    w.sim <- latent.data$w.sim
    B.sim <- latent.data$B.sim

    chunk.size <- getOption("latpos.chunk.size")
    batch.size <- chunk.size%/%(4*prod(dim(B.sim)[c(1,3)]))
    batch.size <- 2*(batch.size%/%2)
    m <- latent.data$sample.size %/% batch.size
    r <- latent.data$sample.size %% batch.size

    dev <- 0
    XWX <- 0
    XWy <- 0

    compute.XWX <- "XWX" %in% compute
    compute.XWy <- "XWy" %in% compute

    I <- nrow(resp$y)
    JT <- ncol(resp$y)
    JTK <- JT*batch.size
    D <- length(parm$latent.dims)

    j. <- resp$j

    if(m>0){

      y <- rep(resp$y,batch.size)
      n <- rep(resp$n,batch.size)
      dim(y) <- dim(n) <- c(I,JTK)
      j <- rep(resp$j,batch.size)
      t <- rep(resp$t,batch.size)

      kk <- 1:batch.size
      for(k in 1:m){

        B <- array(B.sim[,kk,,drop=FALSE],c(JTK,D))
        w <- w.sim[j.,kk,drop=FALSE]
        res <- latpos.eval.parms(y=y,n=n,j=j,t=t,
                            parm=parm,B=B,weights=w,
                            compute=compute)
        dev <- dev + res$deviance
        if(compute.XWX)
          XWX <- XWX + res$XWX
        if(compute.XWy)
          XWy <- XWy + res$XWy
        kk <- kk + batch.size
      }
    }
    if(r>0){

      kk <- m*batch.size + 1:r
      y <- rep(resp$y,r)
      n <- rep(resp$n,r)
      dim(y) <- dim(n) <- c(I,JT*r)
      j <- rep(resp$j,r)
      t <- rep(resp$t,r)

      B <- array(B.sim[,kk,,drop=FALSE],c(JT*r,D))
      w <- w.sim[j.,kk,drop=FALSE]
      repl <- rep(kk,each=JT)
      res <- latpos.eval.parms(y=y,n=n,j=j,t=t,
                          parm=parm,B=B,weights=w,
                          compute=compute)
      dev <- dev + res$deviance
      if(compute.XWX)
        XWX <- XWX + res$XWX
      if(compute.XWy)
        XWy <- XWy + res$XWy
    }

    res <- list()
    res$deviance <- dev
    if(compute.XWX)
      res$XWX <- XWX
    if(compute.XWy)
      res$XWy <- XWy

    res
}
