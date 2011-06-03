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
    J <- length(resp$Tj)
    JT <- ncol(resp$y)
    JTK <- JT*batch.size
    D <- length(parm$latent.dims)

    j. <- resp$j

    bS00 <- 0
    bS11 <- 0
    bS12 <- 0
    bS22 <- 0
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

        bS00 <- bS00 + crossprod(U[t0,,drop=FALSE],w0*U[t0,,drop=FALSE])
        bS11 <- bS11 + crossprod(U[s0,,drop=FALSE], w*U[s0,,drop=FALSE])
        bS12 <- bS12 + crossprod(U[s0,,drop=FALSE], w*U[s1,,drop=FALSE])
        bS22 <- bS22 + crossprod(U[s1,,drop=FALSE], w*U[s1,,drop=FALSE])

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

      bS00 <- bS00 + crossprod(U[t0,,drop=FALSE],w0*U[t0,,drop=FALSE])
      bS11 <- bS11 + crossprod(U[s0,,drop=FALSE], w*U[s0,,drop=FALSE])
      bS12 <- bS12 + crossprod(U[s0,,drop=FALSE], w*U[s1,,drop=FALSE])
      bS22 <- bS22 + crossprod(U[s1,,drop=FALSE], w*U[s1,,drop=FALSE])

      bT0 <- bT0 + sum(w0)
      bT <- bT + sum(w)
    }
    bS21 <- t(bS12)
    
    bT1 <- bT + bT0

    parm <- varParMax(parm,S22=bS22,S12=bS12,S11=bS11,S00=bS00)

    parm
}



