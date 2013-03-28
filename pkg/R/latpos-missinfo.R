latpos.missinfo <- function(resp,parm,latent.data){

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

  j. <- resp$j

  Q.kappa0 <- parm$Q.kappa0
  Q.kappa1 <- parm$Q.kappa1
  
  beta <- parm$beta
  Gamma <- parm$Gamma
  Sigma0 <- parm$Sigma0
  Sigma1 <- parm$Sigma1
  Theta0 <- solve(Sigma0)
  Theta1 <- solve(Sigma1)
  Theta0.x.Theta0 <- Theta0 %x% Theta0
  Theta1.x.Theta1 <- Theta1 %x% Theta1
  I.x.Theta1 <- diag(D) %x% Theta1

  I.Gamma <- diag(D)-Gamma

  g.j <- 0
  G.j <- 0

  if(m>0){

    y <- rep(resp$y,batch.size)
    n <- rep(resp$n,batch.size)
    dim(y) <- dim(n) <- c(I,JTK)

    tr <- rep(resp$t,batch.size)
    t0 <- rep(resp$t0,batch.size)
    s0 <- rep(resp$s0,batch.size)
    s1 <- rep(resp$s1,batch.size)
    jr <- rep(resp$j,batch.size) + rep(J*(1:batch.size-1),each=JT)
    j.r <- rep(1:J,batch.size)
    Tj.r <- rep(Tj,batch.size)

    s0. <- s0
    s0 <- s0. + JT*rep(1:batch.size-1,each=JT)
    s0[t0] <- 0
    s1. <- s1
    s1 <- s1. + JT*rep(1:batch.size-1,each=JT)
    s1[t0] <- 0

    jr0 <- jr[t0]
    jr1 <- jr[s0]
    jr2 <- jr[s1]

    g.jr.skel <- matrix(0,nrow=JK,D*D)

    kk <- 1:batch.size

    for(k in 1:m){

      B <- array(B.sim[,kk,,drop=FALSE],c(JTK,D))
      w <- w.sim[,kk,drop=FALSE]

      g.phi.jr <- latpos.eval.parms(y,n,j=jr,t=tr,parm=parm,B=B,weights=as.vector(w[j.,]),
                                  compute="g.j")$g.j
      
      B0 <- B[t0,,drop=FALSE]
      B1 <- B[s0,,drop=FALSE]
      B2 <- B[s1,,drop=FALSE]

      bb0.jr <- B0
      bs1.jr <- rowsum(B1,jr1)
      bs2.jr <- rowsum(B2,jr1)

      g.beta.jr <- -t(matrix(Theta0%*%beta,D,JK))
      g.beta.jr <- g.beta.jr - t(matrix(crossprod(I.Gamma,Theta0%*%I.Gamma)%*%beta,D,JK))*rep(Tj,batch.size)
      g.beta.jr <- g.beta.jr + bb0.jr%*%Theta0
      g.beta.jr[Tj.r>0] <- g.beta.jr[Tj.r>0] + (bs2.jr - tcrossprod(bs1.jr,Gamma)%*%Theta1%*%I.Gamma)

      U0 <- sweep(B0,2,beta,"-")
      U1 <- sweep(B1,2,beta,"-")
      U2 <- sweep(B2,2,beta,"-")

      U2.GammaU1 <- U2-tcrossprod(U1,Gamma)

      R.Gamma <- flat.gcrossprod(U2.GammaU1,U1,jr1)
      g.Gamma.jr <- g.jr.skel
      g.Gamma.jr[Tj.r>0,] <- R.Gamma %*% I.x.Theta1

      R.Sigma0 <- sweep(flat.outerprod(U0,U0),2,as.vector(Sigma0),"-")
      g.Sigma0.jr <- 0.5*R.Sigma0 %*% Theta0.x.Theta0
      R.Sigma1 <- sweep(flat.outerprod(U2.GammaU1,U2.GammaU1),2,as.vector(Sigma1),"-")
      g.Sigma1.jr <- g.jr.skel
      g.Sigma1.jr[Tj.r>0,] <- 0.5*rowsum(R.Sigma1,jr1) %*% Theta1.x.Theta1

      g.jr <- cbind(g.phi.jr,g.beta.jr,g.Gamma.jr,g.Sigma0.jr%*%Q.kappa0,g.Sigma1.jr%*%Q.kappa1)

      w <- as.vector(w)
      g.j <- g.j + rowsum(w*g.jr,j.r)
      G.j <- G.j + gcrossprod(g.jr,w*g.jr,j.r)

      kk <- kk + batch.size
    }
  }
  if(r>0){

    kk <- m*batch.size + 1:r
    y <- rep(resp$y,r)
    n <- rep(resp$n,r)
    dim(y) <- dim(n) <- c(I,JT*r)

    tr <- rep(resp$t,r)
    t0 <- rep(resp$t0,r)
    s0 <- rep(resp$s0,r)
    s1 <- rep(resp$s1,r)
    jr <- rep(resp$j,r) + rep(J*(1:r-1),each=JT)
    j.r <- rep(1:J,r)

    Tj.r <- rep(Tj,r)

    s0. <- s0
    s0 <- s0. + JT*rep(1:r-1,each=JT)
    s0[t0] <- 0
    s1. <- s1
    s1 <- s1. + JT*rep(1:r-1,each=JT)
    s1[t0] <- 0

    jr0 <- jr[t0]
    jr1 <- jr[s0]
    jr2 <- jr[s1]

    g.jr.skel <- matrix(0,nrow=J*r,D*D)

    B <- array(B.sim[,kk,,drop=FALSE],c(JT*r,D))
    w <- w.sim[,kk,drop=FALSE]

    g.phi.jr <- latpos.eval.parms(y,n,j=jr,t=tr,parm=parm,B=B,weights=as.vector(w[j.,]),
                                compute="g.j")$g.j

    B0 <- B[t0,,drop=FALSE]
    B1 <- B[s0,,drop=FALSE]
    B2 <- B[s1,,drop=FALSE]

    bb0.jr <- B0
    bs1.jr <- rowsum(B1,jr1)
    bs2.jr <- rowsum(B2,jr1)

    g.beta.jr <- -t(matrix(Theta0%*%beta,D,J*r))
    g.beta.jr <- g.beta.jr - t(matrix(crossprod(I.Gamma,Theta0%*%I.Gamma)%*%beta,D,J*r))*rep(Tj,r)
    g.beta.jr <- g.beta.jr + bb0.jr%*%Theta0
    g.beta.jr[Tj.r>0] <- g.beta.jr[Tj.r>0] + (bs2.jr - tcrossprod(bs1.jr,Gamma)%*%Theta1%*%I.Gamma)

    U0 <- sweep(B0,2,beta,"-")
    U1 <- sweep(B1,2,beta,"-")
    U2 <- sweep(B2,2,beta,"-")

    U2.GammaU1 <- U2-tcrossprod(U1,Gamma)

    R.Gamma <- flat.gcrossprod(U2.GammaU1,U1,jr1)
    g.Gamma.jr <- g.jr.skel
    g.Gamma.jr[Tj.r>0,] <- R.Gamma %*% I.x.Theta1

    R.Sigma0 <- sweep(flat.outerprod(U0,U0),2,as.vector(Sigma0),"-")
    g.Sigma0.jr <- 0.5*R.Sigma0 %*% Theta0.x.Theta0
    R.Sigma1 <- sweep(flat.outerprod(U2.GammaU1,U2.GammaU1),2,as.vector(Sigma1),"-")
    g.Sigma1.jr <- g.jr.skel
    g.Sigma1.jr[Tj.r>0,] <- 0.5*rowsum(R.Sigma1,jr1) %*% Theta1.x.Theta1

    g.jr <- cbind(g.phi.jr,g.beta.jr,g.Gamma.jr,g.Sigma0.jr%*%Q.kappa0,g.Sigma1.jr%*%Q.kappa1)

    w <- as.vector(w)
    g.j <- g.j + rowsum(w*g.jr,j.r)
    G.j <- G.j + gcrossprod(g.jr,w*g.jr,j.r)

  }

  g <- colSums(g.j)
  var.gradient <- colSums(G.j) - crossprod(g.j)

  list(gradient=g,var.gradient=var.gradient)
}
