latpos.start.default <- function(latent,manifest,
                                 free.Gamma=TRUE,
                                 A,
                                 Gamma,
                                 ...
                                 ){

    D <- length(latent)
    I <- length(manifest)

    start <- list(
        free.Gamma=free.Gamma
    )

    if(!missing(A)){

        if(is.matrix(A)){

            if(nrow(A)==I && ncol(A)==D)
                start$A <- A
            else {
                start$A <- matrix(0,I,D,dimnames=list(manifest,latent))
                start$A[rownames(A),colnames(A)] <- A
            }
        }
        else if(is.list(A)){
            start$A <- matrix(0,I,D,dimnames=list(manifest,latent))

            for(l in names(A)){
                A.l <- A[[l]]
                start$A[names(A.l),l] <- A.l
            }
        }
        else {
            stop("unsupported type of starting values for A")
        }
    }


    if(!missing(Gamma))
        start$Gamma <- Gamma
    else 
        start$Gamma <- matrix(0,D,D)

    return(start)
}

latpos.start <- function(resp,latent.dims,manifest,start,unfold.method,restrictions=standard.restrictions,maxiter,...){

    I <- length(manifest)
    D <- length(latent.dims)

    y <- resp$y*resp$n
    start.diffs <- sqrt(-log(sweep(y+.5,2,colSums(y+.5),"/")))
    uf <- unfold(start.diffs,ndims=D,method=unfold.method,squared=TRUE)
    start$uf <- uf
    if("A" %in% names(start)){
        A <- start$A
        trans <- ProcrustesTrans(uf$B,A,translate=TRUE)
        A <- sweep(uf$B%*%trans$A,2,trans$b,"-")
        B <- sweep(uf$A%*%trans$A,2,trans$b,"-")
    }
    else {
        A <- uf$B[,1:D,drop=FALSE]
        B <- uf$A[,1:D,drop=FALSE]
    }

    A <- scale(A)
    B <- scale(B)
    colnames(A) <- latent.dims
    rownames(A) <- manifest
    colnames(B) <- latent.dims
    
    D <- ncol(A)
    I <- nrow(A)
    JT <- nrow(B)
    Tj<- tabulate(resp$j)
    J <- length(Tj)

    if(is.function(restrictions)) restrictions <- restrictions(A)
    else if(is.list(restrictions)) restrictions <- restrictions[c("C","d")]
    else stop("no support for restrictions of type",typeof(restrictions))

    ## Matrix enforcing the linear restrictions
    ## on A
    rest.C <- restrictions$C
    rest.d <- if(length(restrictions$d)) restrictions$d else numeric(nrow(rest.C))
    A.restrictor <- tryCatch(restrictor(C=rest.C,d=rest.d),
                             error=function(e) stop("Your model does not seem to be identified. Reduce the number of latent dimensions or set more restrictions.",
                                                    call.=FALSE))
    Q.phi <- A.restrictor$reduction
    r.phi <- A.restrictor$offset

    transf <- rotate.to.restriction(X=A,C=restrictions$C,d=restrictions$d)
    dimnA <- dimnames(A)
    A <- transf$transformed
    rot <- transf$rotation
    transl <- transf$translation
    dimnames(A) <- dimnA

    phi <- crossprod(Q.phi,as.vector(A)-r.phi)
    A[] <- Q.phi%*%phi + r.phi
    B <- sweep(B%*%rot,2,transl,"-")
    
    start$A <- A
    start$phi <- phi
    start$Q.phi <- Q.phi
    start$r.phi <- r.phi
    
    ## Enforce symmetry restrictions on Sigma0 and Sigma1

    if(D>1){

        Q.kappa <- restr.to.symm(D)
        start$Q.kappa0 <- start$Q.kappa1 <- Q.kappa$reduction
    }
    else {

        start$Q.kappa0 <- start$Q.kappa1 <- as.matrix(1)
    }
    
    ## Get starting values
    
    start$Tj <- Tj
    start$latent.dims <- latent.dims

    BB <- split.matrix(B,resp$j)
    B0 <- lapply(BB,colMeans)
    B1 <- Map(qcentr,BB,B0)
 
    B0 <- do.call(rbind,B0)
    B1 <- do.call(rbind,B1)
    start$Btilde <- list(B0=B0,B1=B1)   
    
    S00 <- crossprod(B0)
    S11 <- crossprod(B1[resp$s1,,drop=FALSE])
    S12 <- crossprod(B1[resp$s1,,drop=FALSE],B1[resp$s2,,drop=FALSE])
    S22 <- crossprod(B1[resp$s2,,drop=FALSE])
    
    offdiag <- row(S00)!=col(S00)
    S00[offdiag] <- 0
    S11[offdiag] <- 0
    S12[offdiag] <- 0
    S22[offdiag] <- 0

    Gamma <- t(solve(S11,S12))
    R22 <- S22 - Gamma%*%S12

    Sigma0 <- S00/J
    Sigma1 <- R22/sum(Tj-1)
    
    issq0 <- diag(x=1/sqrt(diag(Sigma0)),nrow=nrow(Sigma0))
    Sigma0 <- issq0 %*% Sigma0 %*% issq0
    Sigma1 <- issq0 %*% Sigma1 %*% issq0
    
    start$Sigma0 <- Sigma0#cov2cor(Sigma0)
    start$Sigma1 <- Sigma1
    start$Gamma <- Gamma
    
    start
}

qcentr <- function(X,mean){
    t(t(X)-c(mean))
}

quncentr <- function(X,mean){
    t(t(X)+c(mean))
}
