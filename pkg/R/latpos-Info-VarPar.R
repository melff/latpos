VarParResids <- function(resp,parm,latent.data){

  U.sim <- latent.data$U.sim
  w.sim <- latent.data$w.sim
  j <- resp$j
  t <- resp$t
  t0 <- resp$t0
  s0 <- resp$s0
  s1 <- resp$s1

  Sigma <- parm$Sigma
  Theta <- parm$Theta
  rho <- parm$rho
  tau <- parm$tau

  Tj1 <- parm$Tj1

  free.rho <- parm$free.rho

  t0[t0] <- which(t0)
  u.j <- unique(j)
  jt <- 1:nrow(U.sim)

  jt <- split(jt,j)

  res <- lapply(u.j,function(j){
                    jt <- jt[[j]]
                    t0 <- t0[jt]
                    s0 <- s0[jt]
                    s1 <- s1[jt]
                    t0 <- t0[t0>0]
                    s0 <- s0[s0>0]
                    s1 <- s1[s1>0]
                    VarParRes.j(
                        Sigma=Sigma,
                        Theta=Theta,
                        rho=rho,
                        tau=tau,
                        Tj1=Tj1[j],
                        w=w.sim[j,],
                        U.t0=U.sim[t0,,,drop=FALSE],
                        U.s0=U.sim[s0,,,drop=FALSE],
                        U.s1=U.sim[s1,,,drop=FALSE],
                        free.rho=free.rho
                      )})

}


VarParRes.j <- function(Sigma,Theta,rho,tau,Tj1,w,U.t0,U.s0,U.s1,free.rho){

  Tj <- Tj1-1

  S0 <- array(U.t0,dim(U.t0)[c(1,2,3,3)])
  S0 <- S0*aperm(S0,c(1,2,4,3))
  S0 <- colSums(S0)

  S1 <- array(U.s0,dim(U.s0)[c(1,2,3,3)])
  S2 <- array(U.s1,dim(U.s1)[c(1,2,3,3)])

  S3 <- colSums(S2*aperm(S2,c(1,2,4,3)))
  S2 <- colSums(S1*aperm(S2,c(1,2,4,3)))
  S1 <- colSums(S1*aperm(S1,c(1,2,4,3)))

  S0 <- aperm(S0,c(2,3,1))
  S1 <- aperm(S1,c(2,3,1))
  S2 <- (aperm(S2,c(2,3,1))+aperm(S2,c(3,2,1)))/2
  S3 <- aperm(S3,c(2,3,1))

  nsim <- dim(S1)[3]
  ndim <- dim(S1)[1]
  ndim2 <- ndim*ndim

  V <- S3 - 2*rho*S2 + rho^2*S1

  R.Theta <- (as.vector(Tj1*Sigma) - as.vector(S0 + tau*V))/2
  dim(R.Theta) <- c(ndim^2,nsim)
  dim(V) <- c(ndim2,nsim)
  dim(S2) <- c(ndim2,nsim)
  dim(S1) <- c(ndim2,nsim)
  vecTheta <- as.vector(Theta)

  R.Theta <- t(R.Theta)
  S2 <- t(S2)
  S1 <- t(S1)
  V <- t(V)

  R.rho <- tau*((S2-rho*S1)%*%vecTheta)
  R.tau <- (ndim*Tj/tau-V%*%vecTheta)/2

  wR.Theta <- w*R.Theta
  wR.rho <- w*R.rho
  wR.tau <- w*R.tau

  RR.Theta <- crossprod(R.Theta,wR.Theta)
  RR.Theta.rho <- crossprod(R.Theta,wR.rho)
  RR.Theta.tau <- crossprod(R.Theta,wR.tau)
  RR.rho.tau <- crossprod(R.rho,wR.tau)
  RR.rho <- crossprod(R.rho,wR.rho)
  RR.tau <- crossprod(R.tau,wR.tau)

  wRwR.Theta <- crossprod(wR.Theta,wR.Theta)
  wRwR.Theta.rho <- crossprod(wR.Theta,wR.rho)
  wRwR.Theta.tau <- crossprod(wR.Theta,wR.tau)
  wRwR.rho.tau <- crossprod(wR.rho,wR.tau)
  wRwR.rho <- crossprod(wR.rho,wR.rho)
  wRwR.tau <- crossprod(wR.tau,wR.tau)

  R.Theta <- colSums(wR.Theta)
  R.rho <- colSums(wR.rho)
  R.tau <- colSums(wR.tau)

  wwR.Theta <- drop(crossprod(w,wR.Theta))
  wwR.rho <- drop(crossprod(w,wR.rho))
  wwR.tau <- drop(crossprod(w,wR.tau))

  wwR.R.Theta <- tcrossprod(wwR.Theta,R.Theta)
  wwR.R.Theta.rho <- tcrossprod(wwR.Theta,R.rho)
  wwR.R.Theta.tau <- tcrossprod(wwR.Theta,R.tau)
  wwR.R.rho.Theta <- tcrossprod(wwR.rho,R.Theta)
  wwR.R.rho <- tcrossprod(wwR.rho,R.rho)
  wwR.R.rho.tau <- tcrossprod(wwR.rho,R.tau)
  wwR.R.tau.Theta <- tcrossprod(wwR.tau,R.Theta)
  wwR.R.tau.rho <- tcrossprod(wwR.tau,R.rho)
  wwR.R.tau <- tcrossprod(wwR.tau,R.tau)

  ww <- sum(w^2)

  ww.RR.Theta <- ww*tcrossprod(R.Theta,R.Theta)
  ww.RR.Theta.rho <- ww*tcrossprod(R.Theta,R.rho)
  ww.RR.Theta.tau <- ww*tcrossprod(R.Theta,R.tau)
  ww.RR.rho <- ww*tcrossprod(R.rho,R.rho)
  ww.RR.rho.tau <- ww*tcrossprod(R.rho,R.tau)
  ww.RR.tau <- ww*tcrossprod(R.tau,R.tau)


  S1 <- crossprod(S1,w)
  dim(S1) <- c(ndim,ndim)
  S2 <- crossprod(S2,w)
  dim(S2) <- c(ndim,ndim)
  V <- crossprod(V,w)
  dim(V) <- c(ndim,ndim)

  if(free.rho){

    R <- c(R.Theta, R.rho, R.tau)

    RR <- rbind(cbind(RR.Theta,RR.Theta.rho,RR.Theta.tau),
                cbind(t(RR.Theta.rho),RR.rho,RR.rho.tau),
                cbind(t(RR.Theta.tau),RR.rho.tau,RR.tau))

    wRwR <- rbind(cbind(wRwR.Theta,wRwR.Theta.rho,wRwR.Theta.tau),
                cbind(t(wRwR.Theta.rho),wRwR.rho,wRwR.rho.tau),
                cbind(t(wRwR.Theta.tau),wRwR.rho.tau,wRwR.tau))

    wwR.R <- rbind(cbind(wwR.R.Theta,wwR.R.Theta.rho,wwR.R.Theta.tau),
                cbind(wwR.R.rho.Theta,wwR.R.rho,wwR.R.rho.tau),
                cbind(wwR.R.tau.Theta,wwR.R.tau.rho,wwR.R.tau))

    ww.RR <- rbind(cbind(ww.RR.Theta,ww.RR.Theta.rho,ww.RR.Theta.tau),
                cbind(t(ww.RR.Theta.rho),ww.RR.rho,ww.RR.rho.tau),
                cbind(t(ww.RR.Theta.tau),ww.RR.rho.tau,ww.RR.tau))

  }
  else {

    R <- c(R.Theta, R.tau)

    RR <- rbind(cbind(RR.Theta,RR.Theta.tau),
                cbind(t(RR.Theta.tau),RR.tau))

    wRwR <- rbind(cbind(wRwR.Theta,wRwR.Theta.tau),
                cbind(t(wRwR.Theta.tau),wRwR.tau))

    wwR.R <- rbind(cbind(wwR.R.Theta,wwR.R.Theta.tau),
                cbind(wwR.R.tau.Theta,wwR.R.tau))

    ww.RR <- rbind(cbind(ww.RR.Theta,ww.RR.Theta.tau),
                cbind(t(ww.RR.Theta.tau),ww.RR.tau))

  }

  var.R <- wRwR - wwR.R - t(wwR.R) + ww.RR

  list(R=R,RR=RR,S1=S1,S2=S2,V=V,var.R=var.R)
}

latpos.GradInfo_VarPar <- function(resp,parm,latent.data){

  Sigma <- parm$Sigma
  Theta <- parm$Theta

  tau <- parm$tau
  rho <- parm$rho

  free.Sigma <- parm$free.Sigma
  free.rho <-   parm$free.rho

  zeta <- 1/tau
  
  Tj <-  parm$Tj
  Tj1 <- parm$Tj1
  D <- ncol(Sigma)

  Theta.x.Theta <- Theta %x% Theta
  d.theta.d.sigma <- -Theta.x.Theta
  d.tau.d.zeta <- -tau^2
  
  Info.cpl.Theta <- sum(Tj1)/2*(Sigma %x% Sigma)
  Info.cpl.Sigma <- sum(Tj1)/2*Theta.x.Theta

  Info.cpl.tau <- D*sum(Tj)/2*zeta^2
  Info.cpl.zeta <- D*sum(Tj)/2*tau^2

  resids <- VarParResids(resp=resp,parm=parm,latent.data=latent.data)

  R <- sapply(resids,"[[", i="R")
  RR <- Sapply(resids,"[[", i="RR")

  S1 <- sapply(resids,"[[", i="S1")
  S1 <- if(is.matrix(S1))rowSums(S1) else sum(S1)

  S2 <- sapply(resids,"[[", i="S2")
  S2 <- if(is.matrix(S2))rowSums(S2) else sum(S2)

  V <- sapply(resids,"[[", i="V")
  V <- if(is.matrix(V))rowSums(V) else sum(V)

  var.R <- Sapply(resids,"[[", i="var.R")
  
  S2rhoS1 <- S2-rho*S1

  Info.cpl.rho <- crossprod(S1,as.vector(Theta))
  Info.cpl.Theta.rho <- tau/2 * as.vector(S2rhoS1)
  Info.cpl.Theta.tau <- as.vector(V)/2
  Info.cpl.rho.tau <- -sum(S2rhoS1*Theta)

  Info.cpl.Sigma.rho <- d.theta.d.sigma%*%Info.cpl.Theta.rho
  Info.cpl.Sigma.zeta <- d.theta.d.sigma%*%Info.cpl.Theta.tau*d.tau.d.zeta
  Info.cpl.rho.zeta <- Info.cpl.rho.tau*d.tau.d.zeta

  RR <- rowSums(RR,dims=2) - tcrossprod(R)
  R <- rowSums(R)
  var.R <- rowSums(var.R,dims=2)

  if(free.rho){
    Info.cpl <- rbind(
                  cbind(Info.cpl.Sigma,Info.cpl.Sigma.rho,Info.cpl.Sigma.zeta),
                  cbind(t(Info.cpl.Sigma.rho),Info.cpl.rho,Info.cpl.rho.zeta),
                  cbind(t(Info.cpl.Sigma.zeta),Info.cpl.rho.zeta,Info.cpl.zeta)
                  )
    Trans <- bdiag(d.theta.d.sigma,1, d.tau.d.zeta)
  }
  else{
    Info.cpl <- rbind(
                  cbind(Info.cpl.Sigma,Info.cpl.Sigma.zeta),
                  cbind(t(Info.cpl.Sigma.zeta),Info.cpl.zeta)
                  )
    Trans <- bdiag(d.theta.d.sigma,d.tau.d.zeta)
  }

  Info.miss <- Trans%*%tcrossprod(RR,Trans)
  gradient <- Trans%*%R
  var.gradient <- Trans%*%tcrossprod(var.R,Trans)

  Info.obs <- Info.cpl - Info.miss

  Q.sigma <- parm$Q.sigma
  if(free.rho)
    Qmat <- bdiag(Q.sigma,1,1)
  else
    Qmat <- bdiag(Q.sigma,1)

  Info.obs <- crossprod(Qmat,Info.obs%*%Qmat)
  gradient <- crossprod(Qmat,gradient)
  var.gradient <- crossprod(Qmat,var.gradient%*%Qmat)

  list(
    gradient=as.vector(gradient),
    Information=as.matrix(Info.obs),
    var.gradient=as.matrix(var.gradient),
    restrictor=as.matrix(Qmat))
}