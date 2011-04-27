OmegaMat1 <- function(Theta,tau,Dmat,ndim,Tj1){

   Thetaj <- c(list(Theta),rep(list(tau*Theta),Tj))
   Thetaj <- bdiag(Thetaj)
   crossprod(Dmat,Thetaj %*% Dmat)
}



Dmat1 <- function(Tj1,D=1,rho=1){

   Dmat <- .Dmat1(Tj1,rho=rho)
   bdiag(Dmat) %x% Diagonal(n=D)
}




OmegaMat <- function(Theta,tau,rho,ndim,Tj){

   Dmat <- Dmat(Tj+1,D=ndim,rho=rho)
   Thetaj <-unlist(lapply(Tj,
              function(Tj)
                c(list(Theta),rep(list(tau*Theta),Tj))),
                recursive=FALSE)
   Thetaj <- bdiag(Thetaj)
   crossprod(Dmat,Thetaj %*% Dmat)
}


Dmat <- function(Tj,D=1,rho=1){

   Dmat <- lapply(Tj,.Dmat1,rho=rho)
   bdiag(Dmat) %x% Diagonal(n=D)
}

.Dmat1 <- function(Tj1,rho=1){

      Dmat <- Diagonal(n=Tj1)
      if(Tj1>1){
            i <- 2:nrow(Dmat)
            j <- i-1
            Dmat[cbind(i,j)] <- -rho
      }
      Dmat
 }

LagMat <- function(Tj,D=1){

   LagMat <- lapply(Tj,.LagMat1)
   bdiag(LagMat) %x% Diagonal(n=D)
}

.LagMat1 <- function(Tj1){

      Diagonal(n=Tj1)[-Tj1,,drop=FALSE]
}

Drop0Mat <- function(Tj,D=1){

   Drop0Mat <- lapply(Tj,.Drop0Mat1)
   bdiag(Drop0Mat) %x% Diagonal(n=D)
}


.Drop0Mat1 <- function(Tj1){

      Diagonal(n=Tj1)[-1,,drop=FALSE]
}
