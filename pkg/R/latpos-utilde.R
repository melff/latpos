latpos.utilde <- function(resp,parm,maxiter=100,verbose=FALSE){

  A <- parm$A
  y <- resp$y
  n <- resp$n
  t <- resp$t

  I <- nrow(y)
  J <- ncol(y)
  D <- ncol(A)
  IJ <- I*J

  n <- as.vector(n)

  P. <- Matrix(0,nrow=J,ncol=IJ)
  ij <- 1:IJ
  j <- rep(1:J,each=I)

  U <- parm$Utilde
  if(is.list(U))
    U <- U$U
  vecU <- as.vector(t(U))
  
  repl <- rep(1,length(y))

  res <- latpos.eval.parms(y=y,n=resp$n,j=resp$j,t=t,parm=parm,U=U,
                replications=repl,
                compute=c("deviance","p"))
  p <- res$p
  dev <- res$deviance

  Omega <- OmegaMat(Sigma0=parm$Sigma0,
                    Sigma1=parm$Sigma1,
                    Gamma=parm$Gamma,
                    tau=parm$tau,
                    ndim=length(parm$latent.dims),Tj=parm$Tj)

  for(iter in 1:maxiter){

    last.dev <- dev
    if(verbose) cat("\nIntegrand iteration",iter)

    p <- as.vector(p)
    Z <- d.eta.d.u(A=A,U=U)

    P.[cbind(j,ij)] <- sqrt(n)*p
    ZWZ <- symmpart(crossprod(Z,n*p*Z) - crossprod(P.%*%Z))
    ZWZOmega <- ZWZ + Omega

    Zr <- crossprod(Z,n*as.vector(y-p))
    #ZWy <- ZWZ %*% vecU + Zr
    #vecU <- as.vector(solve(ZWZOmega,ZWy))
    vecU <- vecU + as.vector(solve(ZWZOmega,Zr-Omega%*%vecU))

    last.U <- U
    U <- t(structure(vecU,dim=c(D,J)))

    res <- latpos.eval.parms(y=y,n=resp$n,j=resp$j,t=t,parm=parm,U=U,
                replications=repl,
                compute=c("deviance","p"))
    p <- res$p
    dev <- res$deviance

    if(!is.finite(dev)){

      trouble <- apply(!is.finite(p),2,any)
      if(!all(trouble)){

        U[trouble,] <- last.U[trouble,]
        res <- latpos.eval.parms(y=y,n=resp$n,j=resp$j,t=t,parm=parm,U=U,
                      replications=repl,
                      compute=c("deviance","p"))
        p <- res$p
        dev <- res$deviance
      }
      for(iiter in 1:maxiter){

        if(is.finite(dev)) break
        if(verbose) cat("\nStep halved")

        U <- (U + last.U)/2
        res <- latpos.eval.parms(y=y,n=resp$n,j=resp$j,t=t,parm=parm,U=U,
                      replications=repl,
                      compute=c("deviance","p"))
        p <- res$p
        dev <- res$deviance
      }

    }
    if(verbose) cat(" - deviance:",dev)

    crit <- abs(dev-last.dev)/abs(.1+last.dev)
    if(verbose && is.finite(crit)) cat(" criterion",crit)
    if(is.finite(crit) && crit < 1e-7){

      if(verbose) cat("\nConverged\n")
      break
    }

    if(is.finite(crit) && dev > last.dev){

      if(verbose) cat("\nCannot decrease deviance, backing up\n")
      U <- last.U
      res <- latpos.eval.parms(y=y,n=resp$n,j=resp$j,t=t,parm=parm,U=U,
                              replications=repl,
                              compute=c("deviance","p"))
      p <- res$p
      dev <- res$deviance
      p <- as.vector(p)
      Z <- d.eta.d.u(A=A,U=U)

      P.[cbind(j,ij)] <- sqrt(n)*p
      ZWZ <- symmpart(crossprod(Z,n*p*Z) - crossprod(P.%*%Z))
      ZWZOmega <- ZWZ + Omega
      break
    }

  }

  list(
    iK2 = t(solve(chol(symmpart(ZWZOmega)))),
    U   = U
  )
}
