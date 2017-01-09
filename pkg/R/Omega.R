Omega.1j <- function(Tj,Theta1,Gamma){

    d <- ncol(Theta1)
    I <- Diagonal(Tj)
    Theta.1j <- I[-Tj,-Tj,drop=FALSE] %x% Theta1
    Delta.j <- I[-1,,drop=FALSE] %x% Diagonal(d) - I[-Tj,,drop=FALSE]%x%Gamma
    Omega.1j <- crossprod(Delta.j,Theta.1j%*%Delta.j)
    return(Omega.1j)
}

dbind <- function(...){

    args <- list(...)
    ncol <- sapply(args,ncol)
    ncol <- sum(ncol)
    nrow <- sapply(args,nrow)
    nrow <- sum(nrow)

    y <- matrix(0,nrow,ncol)
    
    r <- 0
    c <- 0

    for(x in args){

        i <- r + 1:nrow(x)
        j <- c + 1:ncol(x)
        y[i,j] <- x

        r <- r + nrow(x)
        c <- c + ncol(x)
    }
    return(y)
}




