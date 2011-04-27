ProcrustesTrans <- function(X,target,translate=FALSE,dilate=FALSE) {
        C <- crossprod(target,X)
        S <- svd(C)
        A <- S$v %*% t(S$u)
        if(dilate) A <- A * sum(S$d) / sum(diag(crossprod(t(X))))
        if(translate){
            centX <- apply(X,2,mean)
            centY <- apply(target,2,mean)
            b <- centY - crossprod(A,centX)
        } else
        b <- numeric(ncol(X))
        list(A=A,b=b)
}


procrustes.latpos <- function(x,target,...){

    X <- x$parm$A

    if(is.list(target)){

        nms <- unique(unlist(lapply(target,names)))
        if(length(nms)) nc <- length(nms)
        else nc <- max(unlist(lapply(target,length)))

        nr <- length(target)
        tmp <- matrix(0,nrow=nr,ncol=nc)

        if(length(lnms <- names(target)))
          rownames(tmp) <- lnms

        if(length(nms)){

          colnames(tmp) <- nms
          for(i in seq_along(target)){

            ti <- target[[i]]
            tmp[i,names(ti)] <- ti
          }
        }
        else {

          for(i in seq_along(target)){

            ti <- target[[i]]
            tmp[i,seq_along(ti)] <- ti
          }
        }
        target <- tmp
    }


    if(ncol(target)>ncol(X)) stop("to many targets specified")

    if(length(colnames(target))) {

        rotdims <- match(colnames(target),colnames(X))
        if(any(is.na(rotdims))) stop("undefined columns selected")
        }
    else rotdims <- seq_len(ncol(target))

    tmp <- matrix(0,nrow=nrow(X),ncol=ncol(target))

    if(length(rownames(target))){

      ii <- match(rownames(target),rownames(X))
      if(any(is.na(ii))) stop("undefined rows selected")
    }
    else {

      ii <- seq_len(nrow(target))
    }

    tmp[ii,] <- target
    target <- tmp


    transf <- ProcrustesTrans(X[,rotdims,drop=FALSE],target,...)

    Rot <- diag(nrow=ncol(X))
    trnsl <- numeric(ncol(X))
    Rot[rotdims,rotdims] <- transf$A
    trnsl <- trnsl[rotdims] <- transf$b

    x$parm.orig <- x$parm
    x$parm.orig$U.sim <- NULL
    x$parm.orig$w.sim <- NULL
    
    x$parm <- x$parm
    x$parm$U.sim <- NULL
    x$parm$w.sim <- NULL

    x$parm$A <- sweep(x$parm$A%*%Rot,2,-trnsl)
    x$parm$beta <- as.vector(x$parm$beta%*%Rot) - trnsl
    x$parm$Utilde$U <- x$parm$Utilde$U%*%Rot
    x$parm$Utilde <- latpos.utilde(parm=x$parm,resp=x$resp,verbose=FALSE)

    x$parm$Sigma <- crossprod(Rot,x$parm$Sigma%*%Rot)

    if(!inherits(x,"latposProcrustes"))
      class(x) <- c("latposProcrustes",class(x))

    attr(x,"rotation") <- Rot
    attr(x,"translation") <- trnsl

    x$parm.transf$covmat <- vcov.latposProcrustes(x)

    return(x)
}


print.latposProcrustes <- function(x,...){

  cat("\nProcrustes-rotated spatial model of latent positions")

  print.default(x$call)

  cat("\nPositions of objectives:\n")
  print.default(x$parm.transf$A)
  cat("\nManifesto parameters:\n")
  cat("\nMean positions (beta):\n")
  print.default(x$parm.transf$beta)
  cat("\nAutoregression coefficient (rho):\n")
  print.default(x$parm.transf$rho)
  cat("\nSigma:\n")
  print.default(x$parm.transf$Sigma)
  cat("\nzeta:\n")
  print.default(1/x$parm.transf$tau)
  invisible(x)
}

vcov.latposProcrustes <- function(object,...){

  free.beta <- object$parm$free.beta
  free.Sigma <- object$parm$free.Sigma
  free.rho <- object$parm$free.rho
  Rotmat <- attr(object,"rotation")

  if(free.beta)
    ii.Abeta <- length(object$parm$A) + length(object$parm$beta)
  else
    ii.Abeta <- length(object$parm$A)
  ii.Abeta <- 1:ii.Abeta

  if(free.Sigma=="diagonal")
    ii.varPar <- ncol(object$parm$Sigma)
  else if (free.Sigma=="scale")
    ii.varPar <- 1
  else
    ii.varPar <- length(object$parm$Sigma)

  if(free.rho)
    ii.varPar <- ii.varPar + 2
  else
    ii.varPar <- ii.varPar + 1
  ii.varPar <- length(ii.Abeta) + 1:ii.varPar

  vcov.Abeta <- object$parm$covmat[ii.Abeta,ii.Abeta]
  vcov.varPar <- object$parm$covmat[ii.varPar,ii.varPar]

  Trans.A <- t(Rotmat) %x% diag(nrow=nrow(object$parm$A))
  if(free.beta)
    Trans.Abeta <- bdiag(Trans.A,diag(nrow=length(object$parm$beta)))
  else
    Trans.Abeta <- Trans.A
  vcov.Abeta <- tcrossprod(Trans.Abeta%*%vcov.Abeta,Trans.Abeta)

  Trans.sigma <- t(Rotmat) %x% diag(nrow=nrow(object$parm$Sigma))

  SigmaPat <- matrix(0,D,D)

  if(free.Sigma=="diagonal"){

    Qutri.sigma <- diag(x=as.numeric(upper.tri(SigmaPat)))
    Qutri.sigma <- Qutri.sigma[diag(Qutri.sigma)>0,,drop=FALSE]
    Qltri.sigma <- diag(x=as.numeric(lower.tri(SigmaPat)))
    Qltri.sigma <- Qltri.sigma[diag(Qltri.sigma)>0,,drop=FALSE]

    Q.sigma <- rbind(Qutri.sigma,Qltri.sigma)
    Q.sigma <- restrictor(Q.sigma)

  } else if(free.Sigma=="scale"){

    D1 <- D-1

    i <- 1:D1
    ij.1 <- i+D*(1:D1-1)
    ij.2 <- i+1+D*(2:D-1)

    Q.sigma <- matrix(0,nrow=D1,ncol=D^2)
    Q.sigma[cbind(i,ij.1)] <- 1
    Q.sigma[cbind(i,ij.2)] <- -1
    Q.sigma <- restrictor(Q.sigma)[,1,drop=FALSE]
  }

  if(free.Sigma!="full")
    Trans.sigma <- Trans.sigma %*% Q.sigma

  if(free.rho)
    Trans.varPar <- bdiag(Trans.sigma,1,1)
  else
    Trans.varPar <- bdiag(Trans.sigma,1)
  vcov.varPar <- tcrossprod(Trans.varPar%*%vcov.varPar,Trans.varPar)

  vcov <- as.matrix(bdiag(vcov.Abeta,vcov.varPar))
  return(vcov)
}

