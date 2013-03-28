OmegaMat <- function(Sigma0,Sigma1,Gamma,ndim,Tj){

   Dmat <- Dmat(Tj+1,D=ndim,Gamma=Gamma)
   Theta0 <- solve(Sigma0)
   Theta1 <- solve(Sigma1)
   Xi.j <- unlist(lapply(Tj,
              function(Tj)
                c(list(Theta0),rep(list(Theta1),Tj))),
                recursive=FALSE)
   Xi.j <- bdiag(Xi.j)
   crossprod(Dmat,Xi.j %*% Dmat)
}


Dmat <- function(Tj,D=1,Gamma=diag(nrow=D)){

   Dmat <- lapply(Tj,Dmat1,D=D,Gamma=Gamma)
   bdiag(Dmat)
}

Dmat1 <- function(Tj1,D,Gamma){

      Dmat <- Diagonal(n=Tj1*D)
      if(Tj1>1){
        Tj <- Tj1-1
        tmp <- rep(list(Gamma),Tj)
        tmp <- bdiag(tmp)
        ii <- seq_len(Tj*D)
        Dmat[D+ii,ii] <- Dmat[D+ii,ii] - tmp
      }
      Dmat
 }


