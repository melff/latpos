latpos.predict_simul <- function(resp,parm,sampler,sample.size,
                                 mc.cores = getOption("mc.cores", 2L),
                                 verbose=FALSE,
                                 very.silent=FALSE){

  j <- resp$j
  u.j <- unique(j)
  J <- length(u.j)

  Btilde <- parm$Btilde

  JT <- ncol(resp$y)
  ndim <- length(parm$latent.dims)

  sample.size <- 2*(sample.size%/%2+sample.size%%2)
  chunk.size <- getOption("latpos.chunk.size")

  B.sim <- array(0,dim=c(JT,sample.size,ndim))

  jD <- rep(j,each=ndim)
  ll.j <- numeric(J)

  sampler$reset()
  do.unit <- function(jj){

    j. <- which(j==jj)
    jD. <- which(jD==jj)

    y.j <- resp$y[,j.,drop=FALSE]
    n.j <- resp$n[,j.,drop=FALSE]
    j.j <- j[j.]
    t.j <- resp$t[j.]

    Btilde.j <- Btilde$B[j.,,drop=FALSE]
    iK2.j <- Btilde$iK2[jD.,jD.,drop=FALSE]

    batch.size <- chunk.size%/%(length(Btilde.j)*4)
    batch.size <- 2*(batch.size%/%2)

    kept <- 0

    log.thresh.j <- 0
    
    B.sim.j <- B.sim[j.,,,drop=FALSE]
    ll.j <- ll.j[j.]
    
    #cat("\nUnit",jj,"\n")
    repeat{

        sim.tmp <- simul.B.imp(y=y.j,n=n.j,j=j.j,t=t.j,parm=parm,
                            Btilde=Btilde.j,
                            iK2=iK2.j,
                            size=batch.size,
                            sampler=sampler)
        log.w.tmp <- sim.tmp$log.w
        ll.tmp <- sim.tmp$ll.cpl
        Btmp <- sim.tmp$B

        if(log.thresh.j==0) log.thresh.j <- max(log.w.tmp)+1
        else if(max(log.w.tmp)>log.thresh.j){

          log.thresh.j <- max(log.w.tmp)+1
          kept <- 0
          if(verbose) cat("  Restarting in Unit",jj,"\n")
          next
        }

        keep <- log.w.tmp - log.thresh.j > log(runif(n=batch.size))

        Btmp <- Btmp[,keep,,drop=FALSE]
        ll.tmp <- ll.tmp[keep]

        n_keep <- sum(keep)
        if(n_keep > 0){

          if(kept + n_keep > sample.size) {

            n_keep <- sample.size - kept
            keep.tmp <- 1:n_keep
            Btmp <- Btmp[,keep.tmp,,drop=FALSE]
            ll.tmp <- ll.tmp[keep.tmp]

          }

          kk <- kept + 1:n_keep
          B.sim.j[,kk,] <- Btmp

          ll.j <- ll.j + sum(ll.tmp)/sample.size
        }
        kept <- kept + n_keep
        if(verbose)
         cat(" ",n_keep,"of",batch.size,"random vectors kept -",kept,"of",sample.size,"for Unit",jj,"\n")

        if(kept>=sample.size){ 
          if(!very.silent) cat("Done with Unit",jj,"\n")
          break
        }
    }
    list(B.sim=B.sim.j,ll=ll.j)
  }
  require(parallel)
  res <- mclapply(1:J,do.unit,mc.cores=mc.cores)
  for(jj in 1:J){
    j. <- which(j==jj)
    B.sim[j.,,] <- res[[jj]]$B.sim[,,]
    ll.j[j.] <- res[[jj]]$ll
  }

  list(B.sim=B.sim, sample.size=sample.size)
}




predict.latpos <- function(object, newdata = NULL, id=NULL, time=NULL,
                            type=c("posterior modes","posterior means","simulate"),
                            se.fit=FALSE, interval=c("none","normal","percentile"), level=0.95,
                            sample.size = object$sample.size,
                            sampler=object$sampler,
                            maxiter=100,
                            ...){

  type <- match.arg(type)
  interval <- match.arg(interval)

  formula <- as.formula(object$call$formula)
  latent.dims <- all.vars(formula[c(1,3)])

  parm <- object$parm

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

    parm$Btilde <- latpos.Btilde(resp=resp,parm=parm,maxiter=maxiter)
  }
  else {

    orig.order <- object$orig.order
    resp <- object$resp
    if(!is.list(parm$Utilde))
      parm$Btilde <- latpos.Btilde(resp=resp,parm=parm,maxiter=maxiter)
  }

  beta <- parm$beta
  if(type=="posterior modes"){

    Btilde <- parm$Btilde
    iK2 <- Btilde$iK2
    B <- Btilde$B

    colnames(B) <- latent.dims

    if(se.fit || interval!="none"){

      var.B <- array(diag(crossprod(iK2)),dim(B))
      se.B <- sqrt(var.B)
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

    Bsim <- latpos.predict_simul(resp=resp,parm=parm,sampler=sampler,sample.size,...)
    Bsim <- Bsim$B.sim

    if(type=="simulate"){

      return(aperm(Bsim,c(1,3,2))[orig.order,,])
    }
    else{ # type=="posterior means"

      Bsim <- aperm(Bsim,c(2,1,3))
      B <- colMeans(Bsim)
      colnames(B) <- latent.dims

      if(se.fit || interval=="normal"){

            se.B <- t(sqrt(colMeans(Bsim^2) - B^2))
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
        B.lowup <- apply(Bsim,c(2,3),quantile,probs=prob)
        
        res[,,2:3] <- aperm(B.lowup,c(2,3,1))
        dimnames(res) <- list(dimnames(B)[[1]],latent.dims,c("fit","lwr","upr"))
        return(res[orig.order,,])
      }
      return(B[orig.order,])
    }
  }
}

