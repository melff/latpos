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



standard.restrictions <- function(A){

  upper.tri.A <- upper.tri(A)

  C.tri <- diag(x=as.numeric(upper.tri.A))
  C.tri <- C.tri[diag(C.tri)>0,]

  C.sum <- t(diag(nrow=ncol(A)) %x% rep(1,nrow(A)))
  C <- rbind(C.tri,C.sum)

  list(C=C,d=numeric(nrow(C)))
}
