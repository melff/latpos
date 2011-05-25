latpos.impute <- function(resp,parm,sampler,sample.size){

  j <- resp$j
  u.j <- unique(j)
  J <- length(u.j)

  Utilde <- parm$Utilde

  JT <- ncol(resp$y)
  ndim <- length(parm$latent.dims)

  sample.size <- 2*(sample.size%/%2+sample.size%%2)
  chunk.size <- getOption("latpos.chunk.size")

  U.sim <- array(0,dim=c(JT,sample.size,ndim))

  jD <- rep(j,each=ndim)
  ll.j <- numeric(J)

  sampler$reset()
  for(jj in 1:J){

    j. <- which(j==jj)
    jD. <- which(jD==jj)

    y.j <- resp$y[,j.,drop=FALSE]
    n.j <- resp$n[,j.,drop=FALSE]
    j.j <- j[j.]
    t.j <- resp$t[j.]

    Utilde.j <- Utilde$U[j.,,drop=FALSE]
    iK2.j <- Utilde$iK2[jD.,jD.,drop=FALSE]

    batch.size <- chunk.size%/%(length(Utilde.j)*4)
    batch.size <- 2*(batch.size%/%2)

    kept <- 0

    log.thresh.j <- 0
    cat("\nUnit",jj,"\n")
    repeat{

        sim.tmp <- simul.U.imp(y=y.j,n=n.j,j=j.j,t=t.j,parm=parm,
                            Utilde=Utilde.j,
                            iK2=iK2.j,
                            size=batch.size,
                            sampler=sampler)
        log.w.tmp <- sim.tmp$log.w
        ll.tmp <- sim.tmp$ll.cpl
        Utmp <- sim.tmp$U

        if(log.thresh.j==0) log.thresh.j <- max(log.w.tmp)+1
        else if(max(log.w.tmp)>log.thresh.j){

          log.thresh.j <- max(log.w.tmp)+1
          kept <- 0
          cat("  Restarting\n")
          next
        }

        keep <- log.w.tmp - log.thresh.j > log(runif(n=batch.size))

        Utmp <- Utmp[,keep,,drop=FALSE]
        ll.tmp <- ll.tmp[keep]

        n_keep <- sum(keep)
        if(n_keep > 0){

          if(kept + n_keep > sample.size) {

            n_keep <- sample.size - kept
            keep.tmp <- 1:n_keep
            Utmp <- Utmp[,keep.tmp,,drop=FALSE]
            ll.tmp <- ll.tmp[keep.tmp]

          }

          kk <- kept + 1:n_keep
          U.sim[j.,kk,] <- Utmp

          ll.j[j.] <- ll.j[j.] + sum(ll.tmp)/sample.size
        }
        kept <- kept + n_keep
        cat(" ",n_keep,"of",batch.size,"random vectors kept -",kept,"in total\n")

        if(kept>=sample.size) break
    }

  }

  list(U.sim=U.sim, sample.size=sample.size)
}




predict.latpos <- function(object, newdata = NULL, id=NULL, time=NULL,
                            type=c("posterior modes","posterior means","multiple imputation"),
                            se.fit=FALSE, interval=c("none","normal","percentile"), level=0.95,
                            sample.size = object$sample.size,
                            batch.size=object$sample.size,
                            sampler=mvnorm.sampler(),
                            control=latpos.control()){

  type <- match.arg(type)
  interval <- match.arg(interval)

  formula <- as.formula(object$call$formula)
  latent.dims <- all.vars(formula[c(1,3)])

  parm <- object$parm

  maxiter <- control$maxiter

  if(length(newdata)){

    ndims <- length(latent.dims)
    formula <- formula[1:2]
    formula <- update(formula,~.-1)

    m <- match.call(expand.dots=FALSE)
    mf <- m[c(1,match(c("formula","data","subset"),
            names(m),nomatch=0L))]
    mf[[1]] <- as.name("model.frame")
    mf$formula <- formula
    mf$data <- newdata
    mf <- eval(mf,parent.frame())

    y <- t(model.matrix(formula,mf))
    if(missing(id)) id <- object$call$id
    if(missing(time)) time <- object$call$time
    j <- eval(substitute(id),newdata,parent.frame())
    t <- eval(substitute(time),newdata,parent.frame())

    j <- match(j,unique(j))
    split(t,j) <- lapply(split(t,j),function(tj)match(tj,sort(unique(tj)))-1L)
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
    resp$t0 <- t==0
    resp$t <- t
    s <- seq_along(t)
    resp$s0 <- ifelse(t==0,0,s-1)
    resp$s1 <- ifelse(t==0,0,s)

    parm$Utilde <- latpos.utilde(resp=resp,parm=parm,maxiter=maxiter)
  }
  else {

    orig.order <- object$orig.order
    resp <- object$resp
    if(!is.list(parm$Utilde))
      parm$Utilde <- latpos.utilde(resp=resp,parm=parm,maxiter=maxiter)
  }

  beta <- parm$beta
  if(type=="posterior modes"){

    Utilde <- parm$Utilde
    iK2 <- Utilde$iK2
    Utilde <- Utilde$U

    if(length(parm$vcov))
      vcov <- parm$vcov
    else
      vcov <- vcov(object)

    if(parm$free.beta){

      l.beta <- length(beta)
      i.beta <- length(parm$A) + 1:l.beta
      vcov.beta <- vcov[i.beta,i.beta]
      se.beta <- sqrt(diag(vcov.beta))
    }
    else
      se.beta <- 0*beta

    B <- sweep(Utilde,2,beta,"+")
    colnames(B) <- latent.dims

    if(se.fit || interval!="none"){

      var.Utilde <- array(diag(crossprod(iK2)),dim(Utilde))
      if(object$fix.beta)
        se.B <- sqrt(var.Utilde)
      else
        se.B <- sqrt(sweep(var.Utilde,2,se.beta^2,"+"))
      colnames(se.B) <- latent.dims
    }

    if(!(se.fit || interval !="none")) return(B[orig.order,])
    else if(se.fit) return(list(fit=B[orig.order,],se.fit=se.B[orig.order,]))
    else {

      res <- array(B,c(dim(B),3))

      zval <- qnorm((1 - level)/2)
      res[,,2] <- B + zval*se.B
      res[,,3] <- B - zval*se.B
      dimnames(res) <- list(dimnames(B)[[1]],latent.dims,c("fit","lwr","upr"))
      return(res[orig.order,,])
    }
  }
  else {

    Usim <- latpos.impute(resp=resp,parm=parm,sampler=sampler,sample.size)
    Usim <- Usim$U.sim

    if(type=="multiple imputation"){

      B <- beta + aperm(Usim,c(3,1,2))
      return(aperm(B,c(2,1,3))[orig.order,,])
    }
    else{ # type=="posterior means"

      Usim <- aperm(Usim,c(2,1,3))
      means.Usim <- colMeans(Usim)
      B <- t(beta+t(means.Usim))
      colnames(B) <- latent.dims

      if(se.fit || interval=="normal"){

            se.B <- t(sqrt(colMeans(Usim^2) - means.Usim^2))
            colnames(se.B) <- latent.dims
      }
      if(se.fit){

        return(list(fit=B[orig.order,],se.fit=se.B[orig.order,]))
      }
      if(interval=="normal") {

        res <- array(B,c(dim(B),3))

        zval <- qnorm((1 - level)/2)
        res[,,2] <- B + zval*se.B
        res[,,3] <- B - zval*se.B
        dimnames(res) <- list(dimnames(B)[[1]],latent.dims,c("fit","lwr","upr"))
        return(res[orig.order,,])
      }
      else if(interval == "percentile"){

        res <- array(B,c(dim(B),3))
        prob <- (1 - level)/2
        prob <- c(prob,1-prob)
        B.lowup <- apply(Usim,c(2,3),quantile,probs=prob)
        B.lowup <- beta+aperm(B.lowup,c(3,2,1))
        B.lowup <- aperm(B.lowup,c(2,1,3))
        res[,,2:3] <- B.lowup
        dimnames(res) <- list(dimnames(B)[[1]],latent.dims,c("fit","lwr","upr"))
        return(res[orig.order,,])
      }
      return(B[orig.order,])
    }
  }
}

