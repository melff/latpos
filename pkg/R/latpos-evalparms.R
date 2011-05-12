latpos.eval.parms <- function(y,n,j,t,parm,U,replications,compute,weights){

  ndim <- length(parm$latent.dims)
  #nmiss.weights <- !missing(weights)
  if(missing(weights)) weights <- rep(1,nrow(U))
  else weights <- as.vector(weights)

  s <- 1:nrow(U)
  s0 <- ifelse(t==0,0,s-1)
  s1 <- ifelse(t==0,0,s)
  t0 <- which(t==0)
  tnon0 <- t > 0

  A <- parm$A
  B <- sweep(U,2,parm$beta,"+")


  ssq <- numeric(nrow(U))

  if(any(c("deviance","logLik.j") %in% compute)){

    w0 <- weights[t0]
    wnon0 <- weights[tnon0]

    U0 <- U[t0,,drop=FALSE]
    diffU <- U[s1,,drop=FALSE] - parm$rho*U[s0,,drop=FALSE]

    ssq <- 0
    ssq[t0] <- rowSums(w0*U0*(U0%*%parm$Theta))
    ssq[tnon0] <- parm$tau*rowSums(wnon0*diffU*(diffU%*%parm$Theta))

    logdet.Theta <- logdet(parm$Theta)
    log.tau <- log(parm$tau)

    logDet <- 0
    logDet[t0] <- w0*logdet.Theta
    logDet[tnon0] <- ndim*wnon0*log.tau + wnon0*logdet.Theta

    #browser()

  }

  if(any(c("deviance","p","logLik.j","XWX","XWy","G.j") %in% compute)){

#     eta <- eta(A,B)
#     exp.eta <- exp(eta)
#     p <- sweep(exp.eta,2,colSums(exp.eta),"/")
#     p.oldstyle <- p
    p <- latpos_p(A,B)
# browser()
    #I <- nrow(A)
    #wn <- as.vector(rep(weights,each=I)*n)
    #ll.oldstyle <- wn*y*log(p)
    ll <- ll_p(p,y,n,weights)
  }

  res <- list()
  if("deviance" %in% compute){

    dev.resid <- ll_p(y,y,n,weights) - ll
    dev.resid[y==0] <- 0
  #if(any(weights!=1)) browser()
    res$deviance <- 2*sum(dev.resid) + sum(ssq) - sum(logDet)
  }
  if("p" %in% compute)
    res$p <- p

  if("logLik.j" %in% compute){

    log.constpart <- (lfactorial(n[1,]) - colSums(lfactorial(y*n)))*as.vector(weights)

    #if(nmiss.weights)browser()
    ll <- colSums(ll) + log.constpart - ssq/2 + logDet/2
    res$logLik.j <- drop(rowsum(ll,j))

  }

  if(any(c("XWX","XWy","G.j") %in% compute)){

#     X <- d.eta.d.alpha(A,B)%*%parm$Q.phi
#     if(parm$free.beta){
# 
#       X.beta <- d.eta.d.beta(A,U)
#       X <- cbind(X,X.beta)
#     }
#     X.oldstyle <- X

    if(parm$free.beta){
      X <- d.eta.d.phibeta(A,B,parm$Q.phi)
    }
    else
      X <- d.eta.d.phi(A,B,parm$Q.phi)

    jtk <- rep(1:ncol(y),each=nrow(y))

# browser()
  }

  if(any(c("XWX","XWy") %in% compute)){
    #wX <- X - rowsum(as.vector(p)*X,jtk)[jtk,]
    #res$XWX <- crossprod(wX,(wn*as.vector(p))*wX)
    #XWX.oldstyle <- crossprod(wX,(wn*as.vector(p))*wX)
    res$XWX <- latpos_XWX(X,p,n,weights)
    #browser()
  }
  if("XWy" %in% compute){

    #Xr.oldstyle <- crossprod(X,wn*as.vector(y-p))
    Xr <- crossprod(X,latpos_resid(p,y,n,weights))
    psi <- parm$phi
    if(parm$free.beta)
      psi <- c(psi,parm$beta)

    res$XWy <- res$XWX%*%psi + Xr
  }
  if("G.j" %in% compute){

    I <- nrow(y)
    p <- as.vector(p)
    y <- as.vector(y)
    n <- as.vector(n)
    w <- as.vector(rep(weights,each=I))
    j <- rep(j,each=I)
    repl <- rep(replications,each=I)
    R <- X*n*(y-p)
    resid.j <- lapply(split(1:nrow(R),j),function(jj){
                    R.j <- R[jj,,drop=FALSE]
                    w.j <- w[jj]
                    rep.j <- repl[jj]
                    AbetaRes.j(R.j,w.j,rep.j)
                  })

    res$g.j <- sapply(resid.j,"[[",i="R")
    res$G.j <- Sapply(resid.j,"[[",i="RR")
    res$wgwg.j<- Sapply(resid.j,"[[",i="wRwR")
    res$wwg.j <- sapply(resid.j,"[[",i="wwR")
    res$ww.j <- sapply(resid.j,"[[",i="ww")
   
  }

  res
}


latpos.integ.ll <- function(resp,parm,latent.data,compute){

    w.sim <- latent.data$w.sim
    U.sim <- latent.data$U.sim

    chunk.size <- getOption("latpos.chunk.size")
    batch.size <- chunk.size%/%(4*prod(dim(U.sim)[c(1,3)]))
    batch.size <- 2*(batch.size%/%2)
    m <- latent.data$sample.size %/% batch.size
    r <- latent.data$sample.size %% batch.size

    dev <- 0
    XWX <- 0
    XWy <- 0
    G.j <- 0
    g.j <- 0
    wgwg.j <- 0
    wwg.j <- 0
    ww.j <- 0

    compute.XWX <- "XWX" %in% compute
    compute.XWy <- "XWy" %in% compute
    compute.G.j <- "G.j" %in% compute

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

        U <- array(U.sim[,kk,,drop=FALSE],c(JTK,D))
        w <- w.sim[j.,kk,drop=FALSE]
        repl <- rep(kk,each=JT)
        res <- latpos.eval.parms(y=y,n=n,j=j,t=t,
                            parm=parm,U=U,weights=w,
                            replications=repl,
                            compute=compute)
        dev <- dev + res$deviance
        if(compute.XWX)
          XWX <- XWX + res$XWX
        if(compute.XWy)
          XWy <- XWy + res$XWy
        if(compute.G.j){
          G.j <- G.j + res$G.j
          g.j <- g.j + res$g.j
          wgwg.j <- wgwg.j + res$wgwg.j
          wwg.j <- wwg.j + res$wwg.j
          ww.j <- ww.j + res$ww.j
          }
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

      U <- array(U.sim[,kk,,drop=FALSE],c(JT*r,D))
      w <- w.sim[j.,kk,drop=FALSE]
      repl <- rep(kk,each=JT)
      res <- latpos.eval.parms(y=y,n=n,j=j,t=t,
                          parm=parm,U=U,weights=w,
                          replications=repl,
                          compute=compute)
      dev <- dev + res$deviance
      if(compute.XWX)
        XWX <- XWX + res$XWX
      if(compute.XWy)
        XWy <- XWy + res$XWy
      if(compute.G.j){
        G.j <- G.j + res$G.j
        g.j <- g.j + res$g.j
        wgwg.j <- wgwg.j + res$wgwg.j
        wwg.j <- wwg.j + res$wwg.j
        ww.j <- ww.j + res$ww.j
        }
    }

    res <- list()
    res$deviance <- dev
    if(compute.XWX)
      res$XWX <- XWX
    if(compute.XWy)
      res$XWy <- XWy
    if(compute.G.j){
      res$G.j <- G.j
      res$g.j <- g.j
      res$wgwg.j <- wgwg.j
      res$wwg.j <- wwg.j
      res$ww.j <- ww.j
    }

    res
}
