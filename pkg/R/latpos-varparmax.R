varParMax <- function(parm,S22,S12,S11,S00){

  S21 <- t(S12)
  bT <- sum(parm$Tj)
  J <- length(parm$Tj)
  D <- nrow(S00)
  latent.dims <- parm$latent.dims
  
  if(parm$free.Sigma=="full" && parm$free.Sigma01=="full"
      && parm$free.Gamma %in% c("none","full")){

    if(parm$free.Gamma=="none"){

      V22 <- S22 - S21 - S12 + S11
    }
    else {

      Gamma <- t(solve(S11,S12))
      dimnames(Gamma) <- list(latent.dims,latent.dims)
      parm$Gamma <- Gamma
      Gamma.S12 <- Gamma%*%S12
      V22 <- S22 - Gamma.S12 - t(Gamma.S12) + Gamma%*%tcrossprod(S11,Gamma)
    }

    Sigma0 <- S00/J
    Sigma1 <- V22/bT

    dimnames(Sigma0) <- list(latent.dims,latent.dims)
    dimnames(Sigma1) <- list(latent.dims,latent.dims)
    parm$Sigma1 <- Sigma1
    parm$Sigma0 <- Sigma0

    parm$tau <- 1
  }
  else {

    if(length(parm$Gamma))
      Gamma <- parm$Gamma
    else {

      if(parm$free.Gamma=="none")
        Gamma <- diag(nrow=D)
      else
        Gamma <- t(solve(S11,S12))
    }

    if(parm$free.Gamma=="none") l.rho <- 0
    else {

      Q.rho <- parm$Q.rho
      rho <- crossprod(Q.rho,as.vector(Gamma))
      l.rho <- length(rho)
    }

    Sigma0 <- parm$Sigma0
    Sigma1 <- parm$Sigma1
    if(!length(Sigma0)) Sigma0 <- S00/J
    if(!length(Sigma1)) {

      Gamma.S12 <- Gamma%*%S12
      V22 <- S22 - Gamma.S12 - t(Gamma.S12) + Gamma%*%tcrossprod(S11,Gamma)
      Sigma1 <- V22/bT
    }

    Lambda0 <- chol(Sigma0)
    Lambda1 <- chol(Sigma1)

    Q.kappa0 <- parm$Q.kappa0
    Q.kappa1 <- parm$Q.kappa1

    kappa <- crossprod(Q.kappa0,as.vector(Lambda0)) + crossprod(Q.kappa1,as.vector(Lambda1))
    l.kappa <- length(kappa)

    tau <- parm$tau
    if(!length(tau)) tau <- 1

    if(parm$free.Sigma01=="scaled")
      l.tau <- 1
    else
      l.tau <- 0

    par <- kappa
    i.kappa <- 1:l.kappa
    if(l.rho) {

      par <- c(par,rho)
      i.rho <- l.kappa + 1:l.rho
    }
    if(l.tau) {

      par <- c(par,tau)
      i.tau <- l.kappa + l.rho + 1
    }

    searchFun <- function(par){

      kappa <- par[i.kappa]
      Lambda0[] <- Q.kappa0%*%kappa
      Lambda1[] <- Q.kappa1%*%kappa

      Theta0 <- chol2inv(Lambda0)
      Theta1 <- chol2inv(Lambda1)

      if(l.rho){

        rho <- par[i.rho]
        Gamma[] <- Q.rho%*%rho
      }
      if(l.tau)
        tau <- par[i.tau]

      logDet <- 2*suppressWarnings(D*bT*log(tau) - J*sum(log(diag(Lambda0))) - bT*sum(log(diag(Lambda1))))

      R <- S22 - tcrossprod(S21,Gamma) - crossprod(Gamma,S12) + Gamma%*%tcrossprod(S11,Gamma)
      ssq <- tau^2*tr(R%*%Theta1) + tr(S00%*%Theta0)

      ssq - logDet
    }

    opt.res <- optim(par,searchFun,method="BFGS",control=list(trace=1,REPORT=1))
    pseudo.ll <- -opt.res$value
    par <- opt.res$par

    kappa <- par[i.kappa]
    Lambda0[] <- Q.kappa0%*%kappa
    Lambda1[] <- Q.kappa1%*%kappa

    Sigma0 <- crossprod(Lambda0)
    Sigma1 <- crossprod(Lambda1)

    dimnames(Sigma0) <- list(latent.dims,latent.dims)
    dimnames(Sigma1) <- list(latent.dims,latent.dims)

    parm$Sigma0 <- Sigma0
    parm$Sigma1 <- Sigma1

    if(l.rho){

      rho <- par[i.rho]
      dimnames(parm$Gamma) <- list(latent.dims,latent.dims)
      parm$Gamma[] <- Q.rho%*%rho
    }
    if(l.tau){
      tau <- par[i.tau]
      parm$tau <- tau
    }
    else
      parm$tau <- 1
  }

  parm 
}