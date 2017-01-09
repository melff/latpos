latpos.LaplaceEM <- function(resp,start,
                            control=latpos.control()
          ){

    maxiter           <- control$maxiter
    eps        <- control$eps

    stopifnot(ncol(start$A)==ncol(start$U))

    latent.dims <- start$latent.dims
    ndims <- length(latent.dims)

    parm <- start
    
    ## cat(paste("\n",format(Sys.time(),usetz=TRUE),"\n",sep=""))
    ## cat("A:\n")
    ## print(parm$A)
    ## cat("Sigma0:\n")
    ## print(parm$Sigma0)
    ## cat("Sigma1:\n")
    ## print(parm$Sigma1)
    ## cat("Gamma:\n")
    ## print(parm$Gamma)

    
    Q.phi <- parm$Q.phi
    r.phi <- parm$r.phi

    y <- resp$y
    n <- resp$n
    weights <- rep(1,ncol(n))

    j <- resp$j
    t1 <- resp$t1
    t2 <- resp$t2

    Tj <- table(j)
    TT <- sum(Tj)
    TT1 <- TT-1
    
    P. <- Matrix(0,nrow=length(y),ncol=ncol(y))
    ik <- 1:length(y)
    k <- as.vector(col(y))

    
  ### The Iterations ###################
    A <- parm$A
    phi <- parm$phi
    Sigma0 <- parm$Sigma0
    Sigma1 <- parm$Sigma1
    Gamma <- parm$Gamma
    
    issq0 <- diag(x=1/sqrt(diag(Sigma0)),nrow=nrow(Sigma0))
    Sigma0 <- issq0 %*% Sigma0 %*% issq0
    Sigma1 <- issq0 %*% Sigma1 %*% issq0

    Theta0 <- solve(Sigma0)
    Theta1 <- solve(Sigma1)
    
    ll_y <- ll_p(y,y,n,weights)
    ll_y[y==0] <- 0
    ll_y <- sum(ll_y)

    ll <- LaplLogLik(resp,parm)
    
    cat("\n****** Starting Laplace EM algorithm ******")
    cat("\n\tApprox. log-likelihood:",ll)

    converged <- FALSE
    for(iter in 1:maxiter){
        
        cat("\nIteration",iter)
        Btilde <- latpos.Btilde(resp=resp,parm=parm,maxiter=maxiter,verbose=TRUE)
        parm$Btilde <- Btilde
        B0 <- Btilde$B0
        B1 <- Btilde$B1
        B <- B0[j,] + B1

        p <- latpos_p(A,B)
        logf <- ll_p(p,y,n,weights)
        logf <- sum(logf)  
        
        dev <- 2*(ll_y - logf)

        iconverged <- FALSE

        last.A <- A
        for(iiter in 1:maxiter){

            cat("\n\tInner iteration",iiter)
            P.[cbind(ik,k)] <- sqrt(n)*p
            W <- Diagonal(x=as.vector(n*p)) - tcrossprod(P.)
            ## iW <- Diagonal(x=as.vector(1/(n*p)))
            r <- as.vector(n*(y-p))

            X <- d.eta.d.phi(A,B,Q.phi)
            Xr <- crossprod(X,r)
            XWX <- crossprod(X,W%*%X)
            XWy <- XWX%*%phi + Xr

            last.phi <- phi
            phi <- solve(XWX,XWy)
            A[] <- Q.phi%*%phi + r.phi

            p <- latpos_p(A,B)
            logf <- ll_p(p,y,n,weights)
            logf <- sum(logf)  

            last.dev <- dev
            dev <- 2*(ll_y - logf)

            cat(" - deviance:",dev)
            while(!is.finite(dev) || dev > last.dev){

                cat("\n\t\tStep halved")
                phi <- (phi + last.phi)/2
                A[] <- Q.phi%*%phi + r.phi

                p <- latpos_p(A,B)
                logf <- ll_p(p,y,n,weights)
                logf <- sum(logf)  
                
                logf <- sum(logf)  
                dev <- 2*(ll_y - logf)

                cat(" new deviance",dev)
                crit <- abs(dev-last.dev)/abs(.1+last.dev)
                if(crit < eps) {
                    cat("\n")
                    break
                }
            }
            crit <- abs(dev-last.dev)/abs(.1+last.dev)
            cat(" criterion",crit)
            if(crit < eps) {
                cat("\n\tConverged")
                break
            }
        }

        S00 <- crossprod(B0)
        S11 <- crossprod(B1[t1,])
        S12 <- crossprod(B1[t1,],B1[t2,])
        S22 <- crossprod(B1[t2,])

        offdiag <- row(S00)!=col(S00)
        S00[offdiag] <- 0
        S11[offdiag] <- 0
        S12[offdiag] <- 0
        S22[offdiag] <- 0

        last.Gamma <- Gamma
        last.Sigma1 <- Sigma1
        Gamma <- solve(S11,S12)
        Sigma1 <-(S22 - Gamma%*%S11%*%t(Gamma))/TT1

        parm$A <- A
        parm$phi <- phi
        parm$Sigma0 <- Sigma0
        parm$Sigma1 <- Sigma1
        parm$Gamma <- Gamma
        
        ## cat(paste("\n",format(Sys.time(),usetz=TRUE),"\n",sep=""))
        ## cat("A:\n")
        ## print(parm$A)
        ## cat("Sigma0:\n")
        ## print(parm$Sigma0)
        ## cat("Sigma1:\n")
        ## print(parm$Sigma1)
        ## cat("Gamma:\n")
        ## print(parm$Gamma)

        last.ll <- ll 
        ll <- LaplLogLik(resp,parm)
        cat("\n\tApprox. log-likelihood:",ll)

        while(ll < last.ll){

            cat("\n\tStep halved")
            A <- (A + last.A)/2
            phi <- crossprod(Q.phi,as.vector(A)-r.phi)
            Sigma1 <- (Sigma1 + last.Sigma1)/2
            Gamma <- (Gamma + last.Gamma)/2

            parm$A <- A
            parm$phi <- phi
            parm$Sigma1 <- Sigma1
            parm$Gamma <- Gamma

            ll <- LaplLogLik(resp,parm)
            cat(" approx. log-likelihood:",ll)
            
            crit <- abs(ll - last.ll)/(.1 + abs(last.ll))
            if(crit < eps) break
            
        }

        
        crit <- abs(ll - last.ll)/(.1 + abs(last.ll))
        cat(" criterion:",crit)
        if(crit < eps){
            cat("\nConverged")
            converged <- TRUE
            break
        }       
    } # end for

    cat("\n")
    Info.phi <- LapInfo.phi(resp,parm)
    Info.VPar <- LapInfo.VPar(resp,parm)
    
    covmat.A <- Q.phi%*%solve(Info.phi)%*%t(Q.phi)
    covmat.sigma1 <- solve(Info.VPar$sigma1)
    covmat.gamma <- solve(Info.VPar$gamma)

  list(
      resp=resp,parm=parm,
      logLik=ll,
    covmat=list(A=covmat.A,
                sigma1=covmat.sigma1,
                gamma=covmat.gamma),
    control=control,
    converged=converged
    )
}

LaplLogLik <- function(resp,parm){

    res <- 0

    j <- resp$j
    t1 <- resp$t1
    t2 <- resp$t2
    
    Tj <- table(j)
    TT <- sum(Tj)
    TT1 <- TT - 1
    J <- length(Tj)
    
    y <- resp$y
    n <- resp$n
    weights <- resp$weights
    if(!length(weights))
        weights <- rep(1,ncol(y))
    
    A <- parm$A
    B0 <- parm$Btilde$B0
    B1 <- parm$Btilde$B1

    B <- B0[j,] + B1
    
    p <- latpos_p(A,B)
    logf <- ll_p(p,y,n,weights)
    logf <- sum(logf)  
    
    Theta0 <- solve(parm$Sigma0)
    Theta1 <- solve(parm$Sigma1)
    Gamma <- parm$Gamma
    
    logDetTheta0 <- logdet(Theta0)
    logDetTheta1 <- logdet(Theta1)

    S00 <- crossprod(B0)
    S11 <- crossprod(B1[t1,])
    S12 <- crossprod(B1[t1,],B1[t2,])
    S22 <- crossprod(B1[t2,])
    GS12 <- Gamma%*%S12
    R1212 <- S22 + Gamma%*%S11%*%t(Gamma) - (GS12 +t(GS12))
    
    logg <- (J*logDetTheta0 + sum(Theta0*S00)
        + TT1*logDetTheta1
        - sum(Theta1*R1212))/2

    ik <- 1:length(y)
    k <- as.vector(col(y))

    P. <- Matrix(0,nrow=length(y),ncol=ncol(y))
    P.[cbind(ik,k)] <- sqrt(n)*p
    W <- Diagonal(x=as.vector(n*p)) - tcrossprod(P.)

    D <- ncol(B)
    D2 <- D*D
    M01 <- Matrix(0,nrow=length(B),ncol=D*J)
    jd <- rep((j-1)*D,each=D) + rep(1:D,length(j))
    k <- 1:nrow(M01)
    M01[cbind(k,jd)] <- 1
    
    Z1 <- d.eta.d.b(A=A,B=B)
    Z0 <- Z1 %*% M01
    Z <- cbind(Z0,Z1)
    ZWZ <- crossprod(Z,W%*%Z)

    Omega0 <- Diagonal(n=J) %x% Theta0
    Omega1 <- lapply(Tj,Omega.1j,Theta1,Gamma)
    Omega1 <- do.call(bdiag,Omega1)
    
    Omega <- bdiag(Omega0,Omega1)
    ZWZOmega <- ZWZ + Omega
   
    
    logdetZWZOmega <- Matrix::determinant(ZWZOmega)$modulus

    res <- logf + logg - logdetZWZOmega
    
    return(res)
}

LapInfo.phi <- function(resp,parm){

    y <- resp$y
    n <- resp$n

    j <- resp$j
    Tj <- table(j)
    J <- length(Tj)

    A <- parm$A
    Btilde <- parm$Btilde
    B0 <- Btilde$B0
    B1 <- Btilde$B1
    B <- B0[j,] + B1
    Q.phi <- parm$Q.phi
    
    p <- latpos_p(A,B)
    P. <- Matrix(0,nrow=length(y),ncol=ncol(y))
    ik <- 1:length(y)
    k <- as.vector(col(y))

    P.[cbind(ik,k)] <- sqrt(n)*p
    W <- Diagonal(x=as.vector(n*p)) - tcrossprod(P.)

    X <- d.eta.d.phi(A,B,Q.phi)  
    XWX <- crossprod(X,W%*%X)

    D <- ncol(B)
    D2 <- D*D
    M01 <- Matrix(0,nrow=length(B),ncol=D*J)
    jd <- rep((j-1)*D,each=D) + rep(1:D,length(j))
    k <- 1:nrow(M01)
    M01[cbind(k,jd)] <- 1
    
    Z1 <- d.eta.d.b(A=A,B=B)
    Z0 <- Z1 %*% M01
    Z <- cbind(Z0,Z1)
    ZWZ <- crossprod(Z,W%*%Z)

    XWZ <- crossprod(X,W%*%Z)
    ZWX <- Matrix::t(XWZ)

    Theta0 <- solve(parm$Sigma0)
    Theta1 <- solve(parm$Sigma1)
    Gamma <- parm$Gamma
    
    Omega0 <- Diagonal(n=J) %x% Theta0
    Omega1 <- lapply(Tj,Omega.1j,Theta1,Gamma)
    Omega1 <- do.call(bdiag,Omega1)
    
    Omega <- bdiag(Omega0,Omega1)
    ZWZOmega <- ZWZ + Omega
    
    XiVX <- XWX - XWZ%*%solve(ZWZOmega,ZWX)
    return(XiVX)
}

LapInfo.VPar <- function(resp,parm){

    y <- resp$y
    n <- resp$n

    j <- resp$j
    Tj <- table(j)
    J <- length(Tj)
    t1 <- resp$t1
    t2 <- resp$t2

    TT <- sum(Tj)
    TT1 <- TT - 1
    
    A <- parm$A
    Btilde <- parm$Btilde
    B0 <- Btilde$B0
    B1 <- Btilde$B1
    B <- B0[j,] + B1
    Q.phi <- parm$Q.phi
    
    p <- latpos_p(A,B)
    P. <- Matrix(0,nrow=length(y),ncol=ncol(y))
    ik <- 1:length(y)
    k <- as.vector(col(y))
    
    P.[cbind(ik,k)] <- sqrt(n)*p
    W <- Diagonal(x=as.vector(n*p)) - tcrossprod(P.)

    D <- ncol(B)
    D2 <- D*D
    M01 <- Matrix(0,nrow=length(B),ncol=D*J)
    jd <- rep((j-1)*D,each=D) + rep(1:D,length(j))
    k <- 1:nrow(M01)
    M01[cbind(k,jd)] <- 1
    
    Z1 <- d.eta.d.b(A=A,B=B)
    Z0 <- Z1 %*% M01
    Z <- cbind(Z0,Z1)
    ZWZ <- crossprod(Z,W%*%Z)

    Sigma0 <- parm$Sigma0
    Sigma1 <- parm$Sigma1
    Theta0 <- solve(Sigma0)
    Theta1 <- solve(Sigma1)
    Gamma <- parm$Gamma
    
    Omega0 <- Diagonal(n=J) %x% Theta0
    Omega1 <- lapply(Tj,Omega.1j,Theta1,Gamma)
    Omega1 <- do.call(bdiag,Omega1)
    
    Omega <- bdiag(Omega0,Omega1)
    ZWZOmega <- ZWZ + Omega
    
    iZWZOmega <- solve(ZWZOmega)

    jj0 <- 1:(D*J)
    jj1 <- D*J + 1:length(B)

    jj. <- cbind(rep(rep(1:D,D),J),
                 rep(rep(1:D,each=D),J)) + rep((1:J-1)*D,each=D2)

    tt. <- cbind(rep(rep(1:D,D),TT),rep(rep(1:D,each=D),TT)) + rep((1:TT-1)*D,each=D2)

    Z.sigma <- matrix(0,nrow=D,ncol=D2)
    Z.sigma[cbind(1:D,(1:D)^2)] <- 1
    Z.gamma <- Z.sigma
    
    K00 <- iZWZOmega[jj0,jj0]
    K11 <- iZWZOmega[jj1,jj1]

    K00jj <- K00[jj.]
    dim(K00jj) <- c(D,D,J)
    sK00 <- rowSums(K00jj,dims=2)
    KK00 <- kronprod_sum(aperm(K00jj,c(3,1,2)))

    ThetaTheta00 <- Theta0%x%Theta0
    
    Info.sigma0 <- tcrossprod(Z.sigma%*%ThetaTheta00%*%(J*Sigma0%x%Sigma0-KK00)%*%ThetaTheta00,Z.sigma)/2

    S11 <- crossprod(B1[t1,])
    
    K11tt <- K11[tt.]
    dim(K11tt) <- c(D,D,TT)
    sK11 <- rowSums(K11tt[,,t1,drop=FALSE],dims=2)
    
    Info.gamma <- tcrossprod(Z.gamma%*%(TT*Theta1%x%(S11+sK11)),Z.gamma)

    KK11 <- kronprod_sum(aperm(K11tt[,,t1,drop=FALSE],c(3,1,2)))
    KK22 <- kronprod_sum(aperm(K11tt[,,t2,drop=FALSE],c(3,1,2)))
    
    ThetaTheta11 <- Theta1%x%Theta1
    GG <- Gamma%x%Gamma
    Info.sigma1 <- tcrossprod(
        Z.sigma%*%
        ThetaTheta11%*%(TT1*Sigma1%x%Sigma1-KK22 - GG%*%KK11%*%t(GG))%*%ThetaTheta11,
        Z.sigma)/2
    
    list(sigma0=Info.sigma0,
         sigma1=Info.sigma1,
         gamma=Info.gamma)
}

kronprod_sum <- function(x){
    k <- dim(x)[1]
    n <- dim(x)[2]
    m <- dim(x)[3]
    dim(x) <- c(k,n*m)
    xx <- crossprod(x)
    dim(xx) <- c(n,m,n,m)
    xx <- aperm(xx,c(1,3,2,4))
    dim(xx) <- c(n*n,m*m)
    xx
}

