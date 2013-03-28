LVdistMax <- function(parm,b0,s1,s2,S22,S12,S11,S00,maxiter=100,verbose=FALSE,starting=FALSE){

  b0 <- as.vector(b0)
  s1 <- as.vector(s1)
  s2 <- as.vector(s2)

  S22 <- as.matrix(S22)
  S12 <- as.matrix(S12)
  S11 <- as.matrix(S11)
  S00 <- as.matrix(S00)

  S21 <- t(S12)
  Tj <- parm$Tj
  bT <- sum(Tj)
  J <- length(Tj)
  D <- nrow(S00)
  latent.dims <- parm$latent.dims

  Q.rho <- parm$Q.rho

  free.beta <- parm$free.beta
  free.Gamma <- parm$free.Gamma

  l.rho <- ncol(Q.rho)

  start <- numeric(0)

  Gamma <- parm$Gamma
  beta <- parm$beta
  
  if(free.Gamma){

    ii.rho <- 1:l.rho

    
    rho <- crossprod(Q.rho,as.vector(Gamma))
    start <- c(start,rho)
  } 
  if(free.beta){

    l.beta <- length(beta)
    ii.beta <- l.rho + 1:l.beta
    start <- c(start,beta)
  }
  else{

    l.beta <- 0
    ii.beta <- 0
  }

  Imat <- diag(D)

  if(free.Gamma || free.beta){
  
    obFun <- function(par){

      if(free.Gamma)
        Gamma <- matrix(Q.rho%*%par[ii.rho],D,D)
      else
        Gamma <- parm$Gamma

      if(free.beta)
        beta <- par[ii.beta]
      else
        beta <- parm$beta

      betab0 <- tcrossprod(beta,b0)
      Sigma0 <- (S00 - betab0 - t(betab0))/J + tcrossprod(beta)

      IGbeta <- (Imat - Gamma)%*%beta
      s2Gs1 <- s2 - Gamma%*%s1

      IGbeta.s2Gs1 <- tcrossprod(IGbeta,s2Gs1)
      Gamma.S12 <- Gamma%*%S12

      V2211 <- S22 - Gamma.S12 - t(Gamma.S12) + Gamma%*%tcrossprod(S11,Gamma)
      Sigma1 <- (V2211 - IGbeta.s2Gs1 - t(IGbeta.s2Gs1))/bT + tcrossprod(IGbeta)

      Lambda0 <- chol(Sigma0)
      Lambda1 <- chol(Sigma1)

      logdet.Sigma0 <- 2*sum(log(diag(Lambda0)))
      logdet.Sigma1 <- 2*sum(log(diag(Lambda1)))
      structure(J*logdet.Sigma0 + bT*logdet.Sigma1,
          Sigma0=Sigma0,Sigma1=Sigma1)

    }

    cat("\n")
    ores <- nlminb(start=start,objective=obFun,control=list(trace=1))
    par <- ores$par

    if(free.Gamma){

        Gamma <- matrix(Q.rho%*%par[ii.rho],D,D)
        parm$Gamma <- Gamma
    }
    if(free.beta){

        beta <- par[ii.beta]
        parm$beta <- beta
    }

  }

  betab0 <- tcrossprod(beta,b0)
  Sigma0 <- (S00 - betab0 - t(betab0))/J + tcrossprod(beta)

  IGbeta <- (Imat - Gamma)%*%beta
  s2Gs1 <- s2 - Gamma%*%s1

  IGbeta.s2Gs1 <- tcrossprod(IGbeta,s2Gs1)
  Gamma.S12 <- Gamma%*%S12

  V2211 <- S22 - Gamma.S12 - t(Gamma.S12) + Gamma%*%tcrossprod(S11,Gamma)
  Sigma1 <- (V2211 - IGbeta.s2Gs1 - t(IGbeta.s2Gs1))/bT + tcrossprod(IGbeta)

  parm$Sigma0 <- Sigma0
  parm$Sigma1 <- Sigma1

  parm
}