latpos.VarParMstep <- function(resp,parm,latent.data,maxiter){


    w.sim <- latent.data$w.sim
    U.sim <- latent.data$U.sim
    sample.size <- latent.data$sample.size

    chunk.size <- getOption("latpos.chunk.size")
    batch.size <- chunk.size%/%(4*prod(dim(U.sim)[c(1,3)]))
    batch.size <- 2*(batch.size%/%2)
    m <- sample.size %/% batch.size
    r <- sample.size %% batch.size

    dev <- 0

    I <- nrow(resp$y)
    JT <- ncol(resp$y)
    JTK <- JT*batch.size
    D <- length(parm$latent.dims)

    j. <- resp$j

    bS0 <- 0
    bS1 <- 0
    bS2 <- 0
    bS3 <- 0
    bT0 <- 0
    bT <- 0

    if(m>0){

      kk <- 1:batch.size
      for(k in 1:m){

        t0 <- rep(resp$t0,batch.size)
        s0 <- rep(resp$s0,batch.size)
        s1 <- rep(resp$s1,batch.size)

        U <- array(U.sim[,kk,,drop=FALSE],c(JTK,D))
        w <- as.vector(w.sim[j.,kk,drop=FALSE])
        w0 <- w[t0]
        w <- w[s0]

        bS0 <- bS0 + crossprod(U[t0,,drop=FALSE],w0*U[t0,,drop=FALSE])
        bS1 <- bS1 + crossprod(U[s0,,drop=FALSE], w*U[s0,,drop=FALSE])
        bS2 <- bS2 + crossprod(U[s0,,drop=FALSE], w*U[s1,,drop=FALSE])
        bS3 <- bS3 + crossprod(U[s1,,drop=FALSE], w*U[s1,,drop=FALSE])

        bT0 <- bT0 + sum(w0)
        bT <- bT + sum(w)
        kk <- kk + batch.size
      }
    }
    if(r>0){

      kk <- m*batch.size + 1:r

      t0 <- rep(resp$t0,r)
      s0 <- rep(resp$s0,r)
      s1 <- rep(resp$s1,r)

      U <- array(U.sim[,kk,,drop=FALSE],c(JT*r,D))
      w <- as.vector(w.sim[j.,kk,drop=FALSE])
      w0 <- w[t0]
      w <- w[s0]

      bS0 <- bS0 + crossprod(U[t0,,drop=FALSE],w0*U[t0,,drop=FALSE])
      bS1 <- bS1 + crossprod(U[s0,,drop=FALSE], w*U[s0,,drop=FALSE])
      bS2 <- bS2 + crossprod(U[s0,,drop=FALSE], w*U[s1,,drop=FALSE])
      bS3 <- bS3 + crossprod(U[s1,,drop=FALSE], w*U[s1,,drop=FALSE])

      bT0 <- bT0 + sum(w0)
      bT <- bT + sum(w)
    }
    bS2 <- (bS2 + t(bS2))/2

    Gamma <- parm$Gamma
    Sigma <- crossprod(Gamma)
    rho <- parm$rho
    tau <- parm$tau

    #cat("\nvech(Sigma) =",vech(Sigma),"rho =",rho,"zeta =",1/tau)

    free.rho <- parm$free.rho
    free.Sigma <- parm$free.Sigma
    Q.gamma <- parm$Q.gamma

    bT1 <- bT + bT0

    if(free.rho){
      par <- rho
      lower.par <- 0
      upper.par <- 1
      l.rho <- length(rho)
      i.rho <- 1:l.rho
    }
    else{
      par <- numeric(0)
      lower.par <- numeric(0)
      upper.par <- numeric(0)
      l.rho <- 0
      i.rho <- 0
    }
    if(free.Sigma!="none"){
        gamma <- crossprod(Q.gamma,as.vector(Gamma))
        par <- c(par,gamma)
        lower.par <- c(lower.par,rep(-Inf,length(gamma)))
        upper.par <- c(upper.par,rep(Inf,length(gamma)))
        l.gamma <- length(gamma)
        i.gamma <- 1:l.gamma
    }
    else {

      l.gamma <- 0
      i.gamma <- 0
    }
    par <- c(par,tau)
    lower.par <- c(lower.par,1)
    upper.par <- c(upper.par,Inf)

    searchFun <- function(par){

        if(l.rho)
          rho <- par[i.rho]
        else
          rho <- 1
        if(l.gamma){
          gamma <- par[l.rho + i.gamma]
          Gamma[] <- Q.gamma %*% gamma
          Sigma <- crossprod(Gamma)
          Theta <- chol2inv(Gamma)
        }
        tau <- par[l.rho + l.gamma + 1]

        #pen <- 0
        #if(tau < 1) pen <- 1#pen + (tau - 1)^2
        #if(rho < 0) pen <- 1#pen + rho^2
        #if(rho > 1) pen <- 1#pen + (rho - 1)^2
        #if(any(diag(Sigma)<0)) pen <- 1#pen + sum(abs(diag(Sigma)[diag(Sigma)<0]))^2

        log.tau <- log(tau)
        bV <- bS3 - 2*rho*bS2 + rho^2*bS1

        ssq <- sum(diag((bS0+tau*bV)%*%Theta))
        logDet <- D*bT*log.tau + bT1*logdet(Theta)

        ll <- -ssq/2 + logDet/2
        #cat("\nvech(Sigma) =",vech(Sigma),"rho =",rho,"zeta =",1/tau)
        #cat(" ll =",ll)
        -ll
    }
    cat("\n")
    opt.res <- nlminb(start=par,objective=searchFun,grad=NULL,lower=lower.par,upper=upper.par,control=list(trace=10))
    ll <- -opt.res$objective
    par <- opt.res$par
    if(l.rho)
      rho <- par[i.rho]
    else
      rho <- 1
    if(l.gamma){
      gamma <- par[l.rho + i.gamma]
      Gamma[] <- Q.gamma %*% gamma
      Sigma <- crossprod(Gamma)
      Theta <- chol2inv(Gamma)
    }
    tau <- par[l.rho + l.gamma + 1]

    #cat("\nvech(Sigma) =",vech(Sigma),"rho =",rho,"zeta =",1/tau)

    parm$rho <- rho
    parm$gamma <- gamma
    parm$Gamma <- Gamma
    parm$Theta <- Theta
    parm$Sigma <- Sigma
    parm$tau <- tau

    parm
}



