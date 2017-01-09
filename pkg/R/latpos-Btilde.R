latpos.Btilde <- function(resp,parm,maxiter=100,verbose=FALSE){

    A <- parm$A
    y <- resp$y
    n <- resp$n
    j <- resp$j
    t <- resp$t
    t1 <- resp$t1
    t2 <- resp$t2

    D <- ncol(A)
    
    Tj <- table(j)
    J <- length(Tj)
    TT <- sum(Tj)
    TT1 <- TT-J

    jj0 <- 1:(D*J)
    jj1 <- D*J + 1:(D*TT)

    Q0 <- sum0Mat(J) 
    Q1 <- lapply(Tj,sum0Mat)
    Q1 <- bdiag(Q1)
    Q01 <- bdiag(Q0,Q1) %x% Diagonal(D)
    
    Theta0 <- solve(parm$Sigma0)
    Theta1 <- solve(parm$Sigma1)
    Gamma <- parm$Gamma
    
    B0 <- parm$Btilde$B0
    B1 <- parm$Btilde$B1
    
    B <- B0[j,] + B1
    b <- c(t(B0),t(B1))

    D2 <- D*D
    M01 <- Matrix(0,nrow=length(B),ncol=D*J)
    jd <- rep((j-1)*D,each=D) + rep(1:D,length(j))
    k <- 1:nrow(M01)
    M01[cbind(k,jd)] <- 1

    Omega0 <- Diagonal(n=J) %x% Theta0
    Omega1 <- lapply(Tj,Omega.1j,Theta1,Gamma)
    Omega1 <- do.call(bdiag,Omega1)
    Omega <- bdiag(Omega0,Omega1)
    
    weights <- rep(1,ncol(y))
    ll_y <- ll_p(y,y,n,weights)
    ll_y[y==0] <- 0

    p <- latpos_p(A,B)
    ll <- ll_p(p,y,n,weights)
    dev <- sum(ll_y-ll)

    P. <- Matrix(0,nrow=length(y),ncol=ncol(y))
    ik <- 1:length(y)
    k <- as.vector(col(y))
    P.[cbind(ik,k)] <- sqrt(n)*p
    W <- Diagonal(x=as.vector(n*p)) - tcrossprod(P.)

    
    for(iter in 1:maxiter){

        last.dev <- dev

        Z1 <- d.eta.d.b(A=A,B=B)
        Z0 <- Z1 %*% M01
        Z <- cbind(Z0,Z1)

        ZWZ <- crossprod(Z,W%*%Z)
        
        ZWZOmega <- ZWZ + Omega

        Zr <- crossprod(Z,as.vector(n*(y-p)))

        last.b <- b

        ZWy <- ZWZ%*%last.b + Zr

        b <- Q01%*%solve(crossprod(Q01,ZWZOmega%*%Q01),crossprod(Q01,ZWy))
        b0 <- b[jj0]
        b1 <- b[jj1]
        B0 <- t(matrix(b0,ncol=nrow(B0),nrow=ncol(B0)))
        B1 <- t(matrix(b1,ncol=nrow(B1),nrow=ncol(B1)))

        B <- B0[j,] + B1

        p <- latpos_p(A,B)
        ll <- ll_p(p,y,n,weights)
        dev <- sum(ll_y-ll)
        
        #plot(B1[,1],type="l")
        #Sys.sleep(0.1)
        
        #cat("\nIteration:",iter," Deviance:",dev)
        while(!is.finite(dev) || dev > last.dev){

            b <- (b + last.b)/2
            b0 <- b[jj0]
            b1 <- b[jj1]
            B0 <- t(matrix(b0,ncol=nrow(B0),nrow=ncol(B0)))
            B1 <- t(matrix(b1,ncol=nrow(B1),nrow=ncol(B1)))
            
            B <- B0[j,] + B1

            p <- latpos_p(A,B)
            ll <- ll_p(p,y,n,weights)
            dev <- sum(ll_y-ll)
            ## cat("\n\t Stepsize halved - new deviance:",dev)
            if(is.finite(dev)){
                if(dev < last.dev) break
                crit <- abs(dev-last.dev)/abs(.1+last.dev)
                ## cat(" criterion:",crit)
                if(crit < 1e-7) break
            }
        }
        
        crit <- abs(dev-last.dev)/abs(.1+last.dev)
        #cat(" criterion:",crit)
        if(crit < 1e-7) break
    }

    res <- list(B0=B0,
                B1=B1)
    
    return(res)
}





sum0Mat <- function(n){
    x <- Diagonal(n)
    x[n,] <- -1
    x[,-n]
}
