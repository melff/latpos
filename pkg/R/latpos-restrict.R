rotate.to.restriction <- function(X,C,d,maxiter=100,eps=1e-7,verbose=FALSE){

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
  C.tri <- C.tri[diag(C.tri)>0,,drop=FALSE]

  C.sum <- t(diag(nrow=ncol(A)) %x% rep(1,nrow(A)))
  C <- rbind(C.tri,C.sum)

  list(C=C,d=numeric(nrow(C)))
}

set.parms.free <- function(...){

  pat <- list(...)

  function(A) set.parms.free2(pat,A)
}

set.parms.free2 <- function(pat,A){

  Apat <- array(1,dim=dim(A),dimnames=dimnames(A))

  for(i in seq_along(pat)){
    latent.dim <- names(pat)[i]
    free.indicators <- as.character(pat[[i]])
    if(!all(free.indicators %in% rownames(Apat))) stop("undefined indicator")
    Arows <- na.omit(match(free.indicators,rownames(Apat)))
    Acol <- match(latent.dim,colnames(Apat))
    if(!length(Acol)) stop("undefined latent dimension")
    if(length(Arows) && length(Acol)){
      Apat[Arows,Acol] <- 0
      }
  }

  C <- diag(x=as.vector(Apat),nrow=length(as.vector(Apat)))
  use <- diag(C)>0
  C <- C[use,,drop=FALSE]

  C.sum <- t(diag(nrow=ncol(A)) %x% rep(1,nrow(A)))
  C <- rbind(C,C.sum)
  
  d <- numeric(nrow(C))
  list(C=C,d=d)
}

restr.to.symm <- function(D){

  if(D < 2) stop("argument 'D' must be >=2")
  
  Q <- matrix(0,D^2,D^2)
  d.e <- quick.grid(d=1:D,e=1:D)
  d.e <- d.e[d.e[,1]>d.e[,2],,drop=FALSE]

  d <- d.e[,1]
  e <- d.e[,2]

  de <- d + D*(e-1)
  ed <- e + D*(d-1)

  Q[cbind(de,de)] <- 1
  Q[cbind(de,ed)] <- -1

  Q <- Q[diag(Q)>0,,drop=FALSE]
  restrictor(Q)
}