CrossProd <- function(X,Y=X){

    stopifnot(nrow(X)==nrow(Y))
    Y <- array(Y,dim=c(nrow(Y),ncol(Y),ncol(X)))
    Y <- aperm(Y,c(1,3,2))
    Y*as.vector(X)
}

VarPar_GradVar <- function(resp,parm,latent.data){

  U <- latent.data$U.sim
  w <- latent.data$w.sim
  j <- resp$j
  t <- resp$t
  t0 <- resp$t0
  s0 <- resp$s0
  s1 <- resp$s1
  t0[t0] <- which(t0)

  Sigma0 <- parm$Sigma0
  Sigma1 <- parm$Sigma1
  Theta0 <- solve(Sigma0)
  Theta1 <- solve(Sigma1)
  Gamma <- parm$Gamma
  tau <- parm$tau

  Tj1 <- parm$Tj1
  Tj <- Tj1-1

  J <- max(j)
  D <- dim(U)[3]
  nsim <- dim(U)[2]

  npars <- length(Theta0)+length(Theta1)+length(Gamma)+1

  G <- array(0,c(J,npars,npars))
  g <- array(0,c(J,npars))
  S11 <- array(0,c(J,D,D))

  for(jj in 1:J){

    t0.j <- t0[j==jj]
    s0.j <- s0[j==jj]
    s1.j <- s1[j==jj]
    T.j <- Tj[jj]
    T1.j <- Tj1[jj]
    
    U.j0 <- U[t0.j,,,drop=FALSE]
    U.js0 <- U[s0.j,,,drop=FALSE]
    U.js1 <- U[s1.j,,,drop=FALSE]

    dim(U.j0) <- dim(U.j0)[2:3]
    dim(U.js0) <- c(prod(dim(U.js0)[1:2]),dim(U.js0)[3])
    dim(U.js1) <- c(prod(dim(U.js1)[1:2]),dim(U.js1)[3])
    diff.U.js <- U.js1 - tcrossprod(U.js0,Gamma)

    S00j <- CrossProd(U.j0)
    
    V22j <- CrossProd(diff.U.js)
    R12j <- CrossProd(U.js0,diff.U.js)

    dim(V22j) <- c(T.j,nsim,D,D)
    dim(R12j) <- c(T.j,nsim,D,D)
    V22j <- colSums(V22j)
    R12j <- colSums(R12j)

    D.sq <- D*D

    grad.Theta0.j <- sweep(-S00j,c(2,3),Sigma0,"+")
    grad.Theta1.j <- sweep(-tau^2*V22j,c(2,3),T.j*Sigma1,"+")
    grad.Gamma.j <- tau^2*matrix(R12j,nrow=prod(dim(R12j)[1:2]),ncol=D)%*%Theta1
    grad.tau.j <- D*T.j/tau - tau*(matrix(V22j,nrow=dim(V22j)[1],ncol=D.sq)%*%as.vector(Theta1))

    dim(grad.Theta0.j) <- c(nsim,D.sq)
    dim(grad.Theta1.j) <- c(nsim,D.sq)
    dim(grad.Gamma.j) <- c(nsim,D.sq)
    
    grad.j <- cbind(grad.Theta0.j,grad.Theta1.j,grad.Gamma.j,grad.tau.j)

    w.j <- w[jj,]

    G[jj,,] <- crossprod(grad.j,w.j*grad.j)
    g[jj,] <- crossprod(w.j,grad.j)

    S11[jj,,] <- crossprod(U.js0,w.j*U.js0)
  }

  var.gradient <- colSums(G)-crossprod(g)
  S11 <- colSums(S11)
  
  list(gradient=g,var.gradient=var.gradient,S11=S11)
}




latpos.GradInfo_VarPar <- function(resp,parm,latent.data){

  Sigma0 <- parm$Sigma0
  Sigma1 <- parm$Sigma1
  Gamma  <- parm$Gamma

  Lambda0 <- chol(Sigma0)
  Lambda1 <- chol(Sigma1)

  Theta0 <- chol2inv(Lambda0)
  Theta1 <- chol2inv(Lambda1)

  tau <- parm$tau
  
  free.Sigma <- parm$free.Sigma
  free.rho <-   parm$free.rho

  
  Tj <-  parm$Tj
  Tj1 <- parm$Tj1
  D <- ncol(Sigma0)
  D.sq <- D*D

  bT1 <- sum(Tj1)
  bT <- sum(Tj)

  Theta0.x.Theta0 <- Theta0 %x% Theta0
  Theta1.x.Theta1 <- Theta1 %x% Theta1

  Info.Sigma0 <- Theta0.x.Theta0
  Info.Sigma1 <- Theta1.x.Theta1

  d.theta0.d.sigma0 <- Sigma0 %x% Sigma0
  d.theta1.d.sigma1 <- Sigma1 %x% Sigma1

  d.sigma0.d.lambda0 <- array(0,c(D,D,D,D))
  d.sigma1.d.lambda1 <- array(0,c(D,D,D,D))
  defg <- quick.grid(d=1:D,e=1:D,f=1:D,g=1:D)
  d <- defg[,1]
  e <- defg[,2]
  f <- defg[,3]
  g <- defg[,4]
  delta.de <- as.numeric(d==e)
  delta.ge <- as.numeric(g==e)
  dg <- defg[,c(1,4)]
  df <- defg[,c(1,3)]
  d.sigma0.d.lambda0[defg] <- delta.de*Lambda0[dg]+Lambda0[df]*delta.ge
  d.sigma1.d.lambda1[defg] <- delta.de*Lambda1[dg]+Lambda1[df]*delta.ge
  dim(d.sigma0.d.lambda0) <- c(D.sq,D.sq)
  dim(d.sigma1.d.lambda1) <- c(D.sq,D.sq)

  vp <- VarPar_GradVar(resp=resp,parm=parm,latent.data=latent.data)
  gradient <- vp$gradient
  var.gradient <- vp$var.gradient
  S11 <- vp$S11


  Info.Gamma <- S11 %x% (tau^2*Theta0)

  Info.tau <- D*bT/tau^2

  Info.Lambda0 <- d.sigma0.d.lambda0%*%tcrossprod(Info.Sigma0,d.sigma0.d.lambda0)
  Info.Lambda1 <- d.sigma1.d.lambda1%*%tcrossprod(Info.Sigma1,d.sigma1.d.lambda1)


  Info.Lambda <- as.matrix(bdiag(Info.Lambda0,Info.Lambda1))

  Q.kappa <- parm$Q.kappa
  Info.kappa <- crossprod(Q.kappa,Info.Lambda%*%Q.kappa)

  Info.Theta1.tau <- bT/tau*as.vector(Sigma1)
  Info.Sigma1.tau <- d.theta1.d.sigma1%*%Info.Theta1.tau
  Info.Lambda1.tau <- d.sigma1.d.lambda1%*%Info.Sigma1.tau
  Info.Lambda.tau <- c(rep(0,D.sq),Info.Lambda1.tau)
  Q.kappa0 <- parm$Q.kappa0
  Q.kappa1 <- parm$Q.kappa1
  Info.kappa.tau <- c(crossprod(Q.kappa,Info.Lambda.tau))

  if(length(parm$Q.rho)){

    Q.rho <- parm$Q.rho
    Info.rho <- crossprod(Q.rho,Info.Gamma%*%Q.rho)
    Info.rho.tau <- rep(0,ncol(Q.rho))
  }
  else {

    Info.rho <- Info.rho.tau <- Q.rho <- matrix(0,0,0)
  }

  Info.cpl <- as.matrix(bdiag(Info.kappa,Info.rho))

  if(parm$free.Sigma01=="scaled"){

    Info.cpl <- cbind(Info.cpl,c(Info.kappa.tau,Info.rho.tau))
    Info.cpl <- rbind(Info.cpl,c(Info.kappa.tau,Info.rho.tau,Info.tau))
  }

  Info.miss <- var.gradient

  if(parm$free.Gamma=="none"){
      ii.Gamma <- 2*D.sq + 1:D.sq
      Info.miss <- Info.miss[-ii.Gamma,-ii.Gamma]
      Qmat <- Q.kappa
  }
  else
      Qmat <- as.matrix(bdiag(Q.kappa,Q.rho))

  if(parm$free.Sigma01=="scaled"){

    Qmat <- as.matrix(bdiag(Qmat,1))
  }
  else {

    np <- ncol(Info.miss)
    Info.miss <- Info.miss[-np,-np]
  }
  Info.miss <- crossprod(Qmat,Info.miss%*%Qmat)

  Info.obs <- Info.cpl - Info.miss

  covmat <- solve(Info.obs)
  covmat <- Qmat%*%tcrossprod(covmat,Qmat)

  d.parm.d.parmred <- bdiag(d.sigma0.d.lambda0,d.sigma1.d.lambda1)
  if(parm$free.Gamma!="none"){

    d.parm.d.parmred <- bdiag(d.parm.d.parmred,diag(nrow=D.sq))
  }
  if(parm$free.Sigma01=="scaled"){

    d.parm.d.parmred <- bdiag(d.parm.d.parmred,1)
  }
  d.parm.d.parmred <- as.matrix(d.parm.d.parmred)
  covmat <- d.parm.d.parmred %*% tcrossprod(covmat,d.parm.d.parmred)
  
  list(
    covmat=covmat,
    gradient=as.matrix(gradient),
    Information=as.matrix(Info.obs),
    Info.cpl=as.matrix(Info.cpl),
    Info.miss=as.matrix(Info.miss),
    var.gradient=as.matrix(var.gradient),
    restrictor=Qmat)
}