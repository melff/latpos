U1U1 <- function(resp,parm,latent.data){

  B.sim <- latent.data$B.sim
  w.sim <- latent.data$w.sim

  sample.size <- latent.data$sample.size

  chunk.size <- getOption("latpos.chunk.size")
  batch.size <- chunk.size%/%(4*prod(dim(B.sim)[c(1,3)]))
  batch.size <- 2*(batch.size%/%2)
  m <- sample.size %/% batch.size
  r <- sample.size %% batch.size

  Tj1 <- parm$Tj1
  Tj <- parm$Tj

  I <- nrow(resp$y)
  J <- length(Tj)
  JT <- ncol(resp$y)
  JTK <- JT*batch.size
  JK <- J*batch.size
  D <- length(parm$latent.dims)

  beta <- parm$beta
  U1U1.j <- 0

  if(m>0){

    t0 <- rep(resp$t0,batch.size)
    s0 <- rep(resp$s0,batch.size)
    s1 <- rep(resp$s1,batch.size)
    jr <- rep(resp$j,batch.size) + rep(J*(1:batch.size-1),each=JT)
    j.r <- rep(1:J,batch.size)
    Tj.r <- rep(Tj,batch.size)

    s0. <- s0
    s0 <- s0. + JT*rep(1:batch.size-1,each=JT)
    s0[t0] <- 0

    jr1 <- jr[s0]
    
    skel <- matrix(0,nrow=JK,D*D)

    kk <- 1:batch.size

    for(k in 1:m){

      B <- array(B.sim[,kk,,drop=FALSE],c(JTK,D))
      w <- as.vector(w.sim[,kk,drop=FALSE])

      B1 <- B[s0,,drop=FALSE]

      U1 <- sweep(B1,2,beta,"-")

      U1U1.jr <- skel
      U1U1.jr[Tj.r>0,] <- flat.gcrossprod(U1,U1,jr1)

      U1U1.j <- U1U1.j + rowsum(w*U1U1.jr,j.r)

      kk <- kk + batch.size
    }
  }
  if(r>0){

    t0 <- rep(resp$t0,r)
    s0 <- rep(resp$s0,r)
    s1 <- rep(resp$s1,r)
    jr <- rep(resp$j,r) + rep(J*(1:r-1),each=JT)
    j.r <- rep(1:J,r)

    Tj.r <- rep(Tj,r)

    s0. <- s0
    s0 <- s0. + JT*rep(1:r-1,each=JT)
    s0[t0] <- 0

    jr1 <- jr[s0]
    
    skel <- matrix(0,nrow=J*r,D*D)

    kk <- m*batch.size + 1:r

    B <- array(B.sim[,kk,,drop=FALSE],c(JT*r,D))
    w <- as.vector(w.sim[,kk,drop=FALSE])

    B1 <- B[s0,,drop=FALSE]

    U1 <- sweep(B1,2,beta,"-")

    U1U1.jr <- skel
    U1U1.jr[Tj.r>0,] <- flat.gcrossprod(U1,U1,jr1)

    U1U1.j <- U1U1.j + rowsum(w*U1U1.jr,j.r)

  }

  U1U1 <- colSums(U1U1.j)
  dim(U1U1) <- c(D,D)

  U1U1
}




latpos.CplInfo_LVdist <- function(resp,parm,latent.data){

  free.beta <- parm$free.beta
  free.Gamma <- parm$free.Gamma

  Sigma0 <- parm$Sigma0
  Sigma1 <- parm$Sigma1
  beta   <- parm$beta
  Gamma  <- parm$Gamma

  Lambda0 <- chol(Sigma0)
  Lambda1 <- chol(Sigma1)

  Theta0 <- chol2inv(Lambda0)
  Theta1 <- chol2inv(Lambda1)

  Tj <-  parm$Tj
  Tj1 <- parm$Tj1
  D <- ncol(Sigma0)
  D.sq <- D*D

  J <- length(Tj1)
  bT1 <- sum(Tj1)
  bT <- sum(Tj)

  U1U1 <- U1U1(resp=resp,parm=parm,latent.data=latent.data)

  IGamma <- diag(D) - Gamma
  Info.beta <- J*Theta0+bT*crossprod(IGamma,Theta1%*%IGamma)

  Info.Gamma <- U1U1 %x% Theta1

  Theta0.x.Theta0 <- Theta0 %x% Theta0
  Theta1.x.Theta1 <- Theta1 %x% Theta1

  Info.Sigma0 <-  J/2*Theta0.x.Theta0
  Info.Sigma1 <- bT/2*Theta1.x.Theta1

  Q.kappa0 <- parm$Q.kappa0
  Q.kappa1 <- parm$Q.kappa1
  Q.rho    <- parm$Q.rho

  Info.cpl.unrest <- list()
  Qmat <- list()
  Info.cpl.unrest <- c(Info.cpl.unrest,list(Info.beta))

  if(free.beta)
    Qmat <- c(Qmat,list(diag(length(beta))))
  else
    Qmat <- c(Qmat,list(matrix(0,nrow=length(beta),ncol=0)))

  Info.cpl.unrest <- c(Info.cpl.unrest,list(Info.Gamma))
  if(free.Gamma)
    Qmat <- c(Qmat,list(Q.rho))
  else
    Qmat <- c(Qmat,list(matrix(0,nrow=length(Gamma),ncol=0)))

  Info.cpl.unrest <- c(Info.cpl.unrest,list(Info.Sigma0,Info.Sigma1))
  Info.cpl.unrest <- as.matrix(bdiag(Info.cpl.unrest))

  Qmat <- c(Qmat,list(Q.kappa0,Q.kappa1))
  Qmat <- as.matrix(bdiag(Qmat))

  Info.cpl <- crossprod(Qmat,Info.cpl.unrest%*%Qmat)

  list(
    Information=as.matrix(Info.cpl),
    restrictor=Qmat)
}
