
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

d.eta.d.beta <- function(A,U){


    D <- ncol(A)
    I <- nrow(A)
    J <- nrow(U)

    res <- matrix(0,nrow=c(I*J),ncol=D)

    ijd <- quick.grid(i=1:I,j=1:J,d=1:D)

    i <- ijd[,1]
    j <- ijd[,2]
    d <- ijd[,3]
    ij <- i+I*(j-1)

    res[cbind(ij,d)] <- A[cbind(i,d)]
    res
}



d.eta.d.u <- function(A,U){


    D <- ncol(A)
    I <- nrow(A)
    J <- nrow(U)

    res <- Matrix(0,nrow=c(I*J),ncol=c(J*D))

    ijd <- quick.grid(i=1:I,j=1:J,d=1:D)

    i <- ijd[,1]
    j <- ijd[,2]
    d <- ijd[,3]
    ij <- i+I*(j-1)
    dj <- d+D*(j-1)

    res[cbind(ij,dj)] <- A[cbind(i,d)]
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



restrictor <- function(C,d=numeric(m),sign=7){

  ## Create a matrix that maps a small set of
  ## linearly unrestriced parameters
  ## to a larger set of linearly restricted
  ## parameters

  if(!is.matrix(C)) C <- t(as.vector(C))

  m <- nrow(C)
  n <- ncol(C)
  Cm <- t(solve(tcrossprod(C),C))
  #Cm <- MASS::ginv(C) ## too imprecise
  
  M <- diag(n) - Cm %*% C
  QRM <- QR(M)
  list(
    reduction=round(QRM$Q,sign),
    projection=M,
    offset=Cm%*%d,
    C=C,
    ginv.C=Cm,
    d=d
    )
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

rotate.to.restriction <- function(X,C,d,maxiter=100,eps=1e-7,verbose=TRUE){

  restr <- restrictor(C=C,d=d)
  k <- restr$offset
  M <- restr$projection

  Xrest <- array(k + M%*%as.vector(X),dim=dim(X))
  rot <- ProcrustesTrans(X,target=Xrest,translate=TRUE)
  A <- rot$A
  b <- rot$b
  Xtrans <- sweep(X%*%A,2,b,"-")

  ssq <- crossprod(C%*%as.vector(Xtrans)-d)
  if(verbose)cat("\nInitial ssq:",ssq)

  for(iter in 1:maxiter){

    Xrest <- array(k + M%*%as.vector(Xtrans),dim=dim(X))
    rot <- ProcrustesTrans(X,target=Xrest,translate=TRUE)
    last.A <- A
    last.b <- b
    A <- rot$A
    b <- rot$b
    Xtrans <- sweep(X%*%A,2,b,"-")

    crit <- max(max(abs(A - last.A)),max(abs(b-last.b)))
    ssq <- crossprod(C%*%as.vector(Xtrans)-d)
    if(verbose)cat("\nIteration",iter,"ssq:",ssq,"crit:",crit)
    
    if(crit <eps) {
      if(verbose)cat("\nConverged\n")
      break
    }
  }

  Xrest <- array(k + M%*%as.vector(Xtrans),dim=dim(X))
  list(
    transformed = Xtrans,
    restricted = Xrest,
    rotation = A,
    translation = b
    )
}

