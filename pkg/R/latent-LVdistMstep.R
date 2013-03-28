latpos.LVdistMstep <- function(resp,parm,latent.data,maxiter){


    w.sim <- latent.data$w.sim
    B.sim <- latent.data$B.sim
    sample.size <- latent.data$sample.size

    chunk.size <- getOption("latpos.chunk.size")
    batch.size <- chunk.size%/%(4*prod(dim(B.sim)[c(1,3)]))
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

    bb0 <- 0
    bs1 <- 0
    bs2 <- 0
    bS00 <- 0
    bS11 <- 0
    bS12 <- 0
    bS22 <- 0

    if(m>0){

      kk <- 1:batch.size
      for(k in 1:m){

        t0 <- rep(resp$t0,batch.size)
        s0 <- rep(resp$s0,batch.size)
        s1 <- rep(resp$s1,batch.size)

        B <- array(B.sim[,kk,,drop=FALSE],c(JTK,D))
        w <- as.vector(w.sim[j.,kk,drop=FALSE])

        w0 <- w[t0]
        w <- w[s0]

        B0 <- B[t0,,drop=FALSE]
        B1 <- B[s0,,drop=FALSE]
        B2 <- B[s1,,drop=FALSE]
        
        bb0 <- bb0 + colSums(w0*B0)
        bs1 <- bs1 + colSums(w*B1)
        bs2 <- bs2 + colSums(w*B2)
        
        bS00 <- bS00 + crossprod(B0,w0*B0)
        bS11 <- bS11 + crossprod(B1, w*B1)
        bS12 <- bS12 + crossprod(B1, w*B2)
        bS22 <- bS22 + crossprod(B2, w*B2)

        kk <- kk + batch.size
      }
    }
    if(r>0){

      kk <- m*batch.size + 1:r

      t0 <- rep(resp$t0,r)
      s0 <- rep(resp$s0,r)
      s1 <- rep(resp$s1,r)

      B <- array(B.sim[,kk,,drop=FALSE],c(JT*r,D))
      w <- as.vector(w.sim[j.,kk,drop=FALSE])

      w0 <- w[t0]
      w <- w[s0]

      B0 <- B[t0,,drop=FALSE]
      B1 <- B[s0,,drop=FALSE]
      B2 <- B[s1,,drop=FALSE]

      bb0 <- bb0 + colSums(w0*B0)
      bs1 <- bs1 + colSums(w*B1)
      bs2 <- bs2 + colSums(w*B2)

      bS00 <- bS00 + crossprod(B0,w0*B0)
      bS11 <- bS11 + crossprod(B1, w*B1)
      bS12 <- bS12 + crossprod(B1, w*B2)
      bS22 <- bS22 + crossprod(B2, w*B2)

    }
    bS21 <- t(bS12)

    parm <- LVdistMax(parm,b0=bb0,s1=bs1,s2=bs2,
                S22=bS22,S12=bS12,S11=bS11,S00=bS00,maxiter=maxiter,verbose=TRUE)

    parm
}



