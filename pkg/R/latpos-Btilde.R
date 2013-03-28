latpos.Btilde <- function(resp,parm,maxiter=100,verbose=FALSE){

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

  B <- parm$Btilde$B
  
  res <- latpos.eval.parms(y=y,n=resp$n,j=resp$j,t=t,parm=parm,B=B,
                compute=c("deviance","p"))
  p <- res$p
  dev <- res$deviance

  Omega <- OmegaMat(Sigma0=parm$Sigma0,
                    Sigma1=parm$Sigma1,
                    Gamma=parm$Gamma,
                    ndim=length(parm$latent.dims),Tj=parm$Tj)

  for(iter in 1:maxiter){

    last.dev <- dev
    if(verbose) cat("\nIntegrand iteration",iter)

    p <- as.vector(p)
    Z <- d.eta.d.b(A=A,B=B)

    P.[cbind(j,ij)] <- sqrt(n)*p
    ZWZ <- symmpart(crossprod(Z,n*p*Z) - crossprod(P.%*%Z))
    ZWZOmega <- ZWZ + Omega

    Zr <- crossprod(Z,n*as.vector(y-p))

    U <- sweep(B,2,parm$beta,"-")

    vecB <- as.vector(t(B))
    vecU <- as.vector(t(U))
    
    vecB <- vecB + as.vector(solve(ZWZOmega,Zr-Omega%*%vecU))

    last.B <- B
    B <- t(structure(vecB,dim=c(D,J)))
    
    res <- latpos.eval.parms(y=y,n=resp$n,j=resp$j,t=t,parm=parm,B=B,
                compute=c("deviance","p"))
    p <- res$p
    dev <- res$deviance

    if(!is.finite(dev)){

      trouble <- apply(!is.finite(p),2,any)
      if(!all(trouble)){

        B[trouble,] <- last.B[trouble,]
        res <- latpos.eval.parms(y=y,n=resp$n,j=resp$j,t=t,parm=parm,B=B,
                      compute=c("deviance","p"))
        p <- res$p
        dev <- res$deviance
      }
      for(iiter in 1:maxiter){

        if(is.finite(dev)) break
        if(verbose) cat("\nStep halved")

        B <- (B + last.B)/2
        res <- latpos.eval.parms(y=y,n=resp$n,j=resp$j,t=t,parm=parm,B=B,
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
      B <- last.B
      res <- latpos.eval.parms(y=y,n=resp$n,j=resp$j,t=t,parm=parm,B=B,
                              compute=c("deviance","p"))
      p <- res$p
      dev <- res$deviance
      p <- as.vector(p)
      Z <- d.eta.d.b(A=A,B=B)

      P.[cbind(j,ij)] <- sqrt(n)*p
      ZWZ <- symmpart(crossprod(Z,n*p*Z) - crossprod(P.%*%Z))
      ZWZOmega <- ZWZ + Omega
      break
    }

  }

  list(
    iK2 = t(solve(chol(symmpart(ZWZOmega)))),
    B   = B
  )
}
