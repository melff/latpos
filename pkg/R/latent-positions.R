
latpos <- function(formula,data,subset,id,time,
                   unfold.method="Schoenemann",start=NULL,
                   ...){

  latent.dims <- all.vars(formula[c(1,3)])
  formula <- formula[1:2]
  formula <- update(formula,~.-1)
  manifest <- all.vars(formula)
  cl <- match.call()
  m <- match.call(expand.dots=FALSE)
  mf <- m[c(1,match(c("formula","data","subset"),
          names(m),nomatch=0L))]
  mf[[1]] <- as.name("model.frame")
  mf$formula <- formula
  mf <- eval(mf,parent.frame())

  if(missing(id)) stop("argument 'id' is missing, with no default")
  if(missing(time)) stop("argument 'time' is missing, with no default")
  j <- eval(substitute(id),data,parent.frame())
  t <- eval(substitute(time),data,parent.frame())
  j <- match(j,unique(j))
  split(t,j) <- lapply(split(t,j),function(tj)match(tj,sort(unique(tj))))
  jt.order <- order(j,t)
  orig.order <- (1:nrow(mf))[jt.order]

  y <- t(model.matrix(formula,mf))[,jt.order]
  j <- j[jt.order]
  t <- t[jt.order]

  resp <- list()
  resp$n <- array(rep(colSums(y),each=nrow(y)),dim=dim(y))
  resp$y <- y/resp$n
  resp$y[resp$n==0] <- 0
  resp$j <- j
  resp$t0 <- t==1
  resp$t <- t
  resp$t1 <- ifelse(t==1,0,t-1)
  resp$t2 <- ifelse(t==1,0,t)
  s <- seq_along(t)
  resp$s <- s  
  resp$s1 <- ifelse(t==1,0,s-1)
  resp$s2 <- ifelse(t==1,0,s)

  parm.names <- c("A","Sigma","Sigma0","Sigma1","Gamma","tau")

  dots <- list(...)
  if(any(parm.names %in% names(dots))){
    if(!length(start))
      start <- dots[parm.names]
    else
      start <- c(start,dots[setdiff(parm.names,names(start))])
  }

    
  if(length(start)){
    if(!is.list(start)) stop("'start' argument must be a list")

    start <- start[parm.names]
    sdf.args <- c(
                  list(
                    latent=latent.dims,manifest=manifest
                  ),
                  start,
                  dots
                )
    start <- do.call(latpos.start.default,sdf.args)
  }
  else
    start <- latpos.start.default(latent=latent.dims,manifest=manifest,...)

  start <- start[nzchar(names(start)) & !is.na(names(start))]
  
  start <- latpos.start(resp=resp,latent.dims=latent.dims,manifest=manifest,
                        start=start,
                        unfold.method=unfold.method,...)
  
  fit <- latpos.fit(resp=resp,
                    start=start,
                    control=latpos.control(...))
    
  fit$call <- cl
  fit$data <- data
  fit$jt.order <- jt.order
  fit$orig.order <- orig.order
  fit$start <- start

  structure(fit,
    class="latpos"
  )
}



latpos.control <- function(maxiter=200,
                            eps=1e-7,
                            ...){
  list(
    maxiter=maxiter,
    eps=eps
    )
}


latpos.fit <- function(resp,start,
                       control=latpos.control(),
                       ...){

    parm <- start

    fit <- latpos.LaplaceEM(resp,start,control)
        
    return(fit)
}





logdet <- function(x) determinant(x,logarithm=TRUE)$modulus




vech <- function(M){

  stopifnot(nrow(M)==ncol(M))
  M[lower.tri(M,diag=TRUE)]
}

"vech<-" <- function(M,value){

  stopifnot(nrow(M)==ncol(M))
  M[lower.tri(M,diag=TRUE)] <- value
  M
}

symmFill <- function(D,fill){

  M <- matrix(NA,D,D)
  M[lower.tri(M,diag=TRUE)] <- fill
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  M
}


print.latpos <- function(x,...){

  cat("\nSpatial model of latent positions\n")

  print.default(x$call)

  cat("\nLocations of criteria:\n")
  print.default(x$parm$A)
  cat("\nLatent positions distribution:\n")
  cat("\nAutoregression coefficients (Gamma):\n")
  print.default(x$parm$Gamma)
  cat("\nActor average position variance (Sigma0):\n")
  print.default(x$parm$Sigma0)
  cat("\nPosition change variance (Sigma1):\n")
  print.default(x$parm$Sigma1)
  invisible(x)
}

logLik.latpos <- function(object,...){
    ll <- object$logLik
    parm <- object$parm
    df <- with(parm, length(phi) + ncol(Sigma1) + ncol(Gamma))

    structure(ll,
              df=df,
              class="logLik")
}
