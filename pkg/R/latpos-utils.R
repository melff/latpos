
latpos_p <- function(A,B) .Call("latpos_p",A,B)
ll_p <- function(p,y,n,weights) .Call("ll_p",p,y,n,weights)
latpos_resid <- function(p,y,n,weights) .Call("latpos_resid",p,y,n,weights)

d.eta.d.phi <- function(A,B,Q) .Call("d_eta_d_phi",A,B,Q)

latpos_XWX <- function(X,p,n,weights) .Call("latpos_XWX",X,p,n,weights)

eta <- function(A,B){

    D <- ncol(A)
    I <- nrow(A)
    J <- nrow(B)


    ijd <- quick.grid(i=1:I,j=1:J,d=1:D)

    i <- ijd[,1]
    j <- ijd[,2]
    d <- ijd[,3]
    ij <- i+I*(j-1)

    AB <- matrix(0,nrow=c(I*J),ncol=D)
    AB[cbind(ij,d)] <- A[cbind(i,d)]*B[cbind(j,d)]
    eta <- rowSums(AB)
    eta <- eta-rep(rowSums(A*A/2),J)
    dim(eta) <- c(I,J)
    eta
}

d.eta.d.alpha <- function(A,B){

    D <- ncol(A)
    I <- nrow(A)
    J <- nrow(B)

    res <- matrix(0,nrow=c(I*J),ncol=c(I*D))

    ijd <- quick.grid(i=1:I,j=1:J,d=1:D)

    i <- ijd[,1]
    j <- ijd[,2]
    d <- ijd[,3]
    ij <- i+I*(j-1)
    id <- i+I*(d-1)

    res[cbind(ij,id)] <- B[cbind(j,d)] - A[cbind(i,d)]
    res
}

d.eta.d.b <- function(A,B){


    D <- ncol(A)
    I <- nrow(A)
    J <- nrow(B)

    res <- Matrix(0,nrow=c(I*J),ncol=c(J*D))

    ijd <- quick.grid(i=1:I,j=1:J,d=1:D)

    i <- ijd[,1]
    j <- ijd[,2]
    d <- ijd[,3]
    ij <- i+I*(j-1)
    dj <- d+D*(j-1)
    
    res[cbind(ij,dj)] <- A[cbind(i,d)]  - B[cbind(j,d)]
    res
}



quick.grid <- function(...){

    my.args <- list(...)

    for(i in seq(along=my.args)){

        x <- my.args[[i]]
        if(i==1) {

            res <- as.matrix(x)
        }
        else {

            II <- nrow(res)
            JJ <- length(x)
            ii <- 1:II
            jj <- 1:JJ

            res <- cbind(res[rep(ii,JJ),],rep(x,each=II))
        }
    }

    colnames(res) <- names(my.args)
    res
}

QR <- function(M){
    ## standardised QR decomposition:
    ## diagonal elements of R are always
    ## positive, only "significant"
    ## columns in Q are returned

    qrM <- qr(M)
    rnk <- qrM$rank
    Q <- qr.Q(qrM)[,1:rnk,drop=FALSE]
    R <- qr.R(qrM)[1:rnk,,drop=FALSE]

    sgndR <- sign(diag(R))
    sgndR <- diag(x=sgndR,nrow=length(sgndR))
    Q <- Q%*%sgndR
    R <- sgndR%*%R

    list(Q=Q, R=R)
}




split.matrix <- function(x, f, drop = FALSE, along=c("rows","columns"), ...){

   along <- match.arg(along)
   i <- switch(along,
       rows=split(seq_len(nrow(x)),f,drop,...),
       columns=split(seq_len(ncol(x)),f,drop,...))
   switch(along,
      rows=lapply(i,function(i)x[i,,drop=FALSE]),
      columns=lapply(i,function(i)x[,i,drop=FALSE]))
}

dblocks <- function(x,f,drop=FALSE){

  stopifnot(nrow(x)==ncol(x))
  i <- split(seq_len(nrow(x)),f,drop)
  lapply(i,function(i)x[i,i,drop=FALSE])
}

tr <- function(x) sum(diag(as.matrix(x)))
logSymmDet <- function(M) 2*sum(log(diag(chol(M))))
solveSymm <- function(M) chol2inv(chol(M))

parm2psi <- function(parm){

    psi <- parm$phi

    Lambda0 <- chol(parm$Sigma0)
    Lambda1 <- chol(parm$Sigma1)
    kappa <- c(
                Lambda0[upper.tri(Lambda0,diag=TRUE)],
                Lambda1[upper.tri(Lambda0,diag=TRUE)]
              )
    
    psi <- c(psi,kappa)
        
    if(parm$free.Gamma) {
      rho <- as.vector(parm$Gamma)
      psi <- c(psi,rho)
    }

    psi
}

gcrossprod <- function(X,Y,f){

  r <- nrow(X)
  stopifnot(nrow(Y)==r && length(f)==r)
  
  m <- ncol(X)
  n <- ncol(Y)

  i <- rep(1:m,n)
  j <- rep(1:n,each=m)

  R <- rowsum(X[,i,drop=FALSE]*Y[,j,drop=FALSE],f)
  dim(R) <- c(dim(R)[1],m,n)
  R
}

flat.gcrossprod <- function(X,Y,f){

  r <- nrow(X)
  stopifnot(nrow(Y)==r && length(f)==r)

  m <- ncol(X)
  n <- ncol(Y)

  i <- rep(1:m,n)
  j <- rep(1:n,each=m)

  rowsum(X[,i,drop=FALSE]*Y[,j,drop=FALSE],f)
}


flat.outerprod <- function(X,Y){

  r <- nrow(X)
  stopifnot(nrow(Y)==r)

  m <- ncol(X)
  n <- ncol(Y)

  i <- rep(1:m,n)
  j <- rep(1:n,each=m)

  X[,i,drop=FALSE]*Y[,j,drop=FALSE]
}

