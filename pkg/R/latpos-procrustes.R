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


