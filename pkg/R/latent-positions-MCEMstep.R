latpos.MCEMstep <- function(resp,parm,
                            latent.data,
                            trace,
                            iter,
                            sampler,
                            control){

  maxiter           <- control$maxiter
  initial.size      <- control$initial.size
  diff.Q.alpha      <- control$diff.Q.alpha
  diff.Q.eps        <- control$diff.Q.eps
  diff.logLik.eps   <- control$diff.logLik.eps
  abs.diff.psi.eps  <- control$abs.diff.psi.eps
  rel.diff.psi.eps  <- control$rel.diff.psi.eps
  max.size          <- control$max.size
  min.final.size    <- control$min.final.size
  sparsity.eps      <- control$sparsity.eps
  ll.linesearch <- control$ll.linesearch
  Q.linesearch  <- control$Q.linesearch

  last.parm <- parm

  last.psi <- last.parm$phi
  if(last.parm$free.beta) last.psi <- c(last.psi,last.parm$beta)
  if(last.parm$free.Sigma!="none") last.psi <- c(last.psi,vech(last.parm$Sigma))
  if(last.parm$free.rho) last.psi <- c(last.psi,last.parm$rho)
  last.psi <- c(last.psi,1/last.parm$tau)

  if(!latent.data$sample.size) latent.data$sample.size <- initial.size
  sample.size <- latent.data$sample.size
  last.sample.size <- sample.size
  sample.size.start <- sample.size
  
  
  last.logLik <- last.parm$logLik

  diffQsign <- FALSE
  converged <- FALSE
  used.sample.size <- latent.data$sample.size

  cat("\n\n**** Monte Carlo EM Iteration",iter,"****")
  while(!diffQsign && !converged) {

    used.sample.size <- latent.data$sample.size
    maybe.converged <- FALSE
    ## The M-step: Location parameters
    
    parm <- latpos.AbetaMstep(resp=resp,parm=parm,latent.data=latent.data,maxiter=maxiter)

    ## The M-step: Variance parameters
    
    parm <- latpos.VarParMstep(resp=resp,parm=parm,latent.data=latent.data,maxiter=maxiter)

    ## Checking the amount of increase in the Q function
    
    Lambda.Q.res <- latpos.Lambda.Q(resp=resp,latent.data=latent.data,parm=parm,last.parm=last.parm)

    Q.psi <- Lambda.Q.res$Q
    diff.Q.psi <- Lambda.Q.res$diff.Q
    ASE.diff.Q.psi <- Lambda.Q.res$ASE.diff.Q
    cat("\nQ improvement:",diff.Q.psi)

    psi <- parm$phi
    if(parm$free.beta) psi <- c(psi,parm$beta)
    if(parm$free.Sigma!="none") psi <- c(psi,vech(parm$Sigma))
    if(parm$free.rho) psi <- c(psi,parm$rho)
    psi <- c(psi,1/parm$tau)

    zval.diff.Q.psi <- diff.Q.psi/ASE.diff.Q.psi
    crit.diff.Q.psi <- abs(diff.Q.psi + qnorm(diff.Q.alpha)*ASE.diff.Q.psi)/(.1+abs(Q.psi))
    diff.psi <- psi - last.psi
    psi.max.diff <- max(abs(diff.psi))
    max.diff.id <- which(abs(diff.psi)==psi.max.diff)[1]
    psi.crit <- max(abs(psi-last.psi)/(.0001+last.psi))

    if(Q.linesearch &&
      abs(zval.diff.Q.psi) > qnorm(1-diff.Q.alpha) &&
      crit.diff.Q.psi > diff.Q.eps &&
      psi.crit > rel.diff.psi.eps &&
      abs(psi.max.diff) > abs.diff.psi.eps &&
      diff.Q.psi < 0 &&
      sample.size.start == latent.data$sample.size ){

      cat(" -- need to conduct line search")

      searchFun <- function(lambda){

        cat("\nlambda =",lambda)
        mlambda <- 1-lambda
        parm$A <- lambda*parm$A + mlambda*last.parm$A
        parm$phi <- lambda*parm$phi + mlambda*last.parm$phi
        parm$beta <- lambda*parm$beta + mlambda*last.parm$beta
        parm$Gamma <- lambda*parm$Gamma + mlambda*last.parm$Gamma
        parm$Sigma <- crossprod(parm$Gamma)
        parm$Theta <- chol2inv(parm$Gamma)
        parm$rho <- lambda*parm$rho + mlambda*last.parm$rho
        parm$tau <- lambda*parm$tau + mlambda*last.parm$tau
        Lambda.Q.res <- latpos.Lambda.Q(resp=resp,latent.data=latent.data,parm=parm,last.parm=last.parm)
        diff.Q.psi <- Lambda.Q.res$diff.Q
        cat(" Q improvement:",diff.Q.psi)
        structure(-diff.Q.psi,Lambda.Q.res=Lambda.Q.res,parm=parm)
      }

      opt.res <- optimise(searchFun,interval=c(0,1))
      parm <- attr(opt.res$objective,"parm")
      Lambda.Q.res <- attr(opt.res$objective,"Lambda.Q.res")
      Q.psi <- Lambda.Q.res$Q
      diff.Q.psi <- Lambda.Q.res$diff.Q
      ASE.diff.Q.psi <- Lambda.Q.res$ASE.diff.Q

      zval.diff.Q.psi <- diff.Q.psi/ASE.diff.Q.psi
      crit.diff.Q.psi <- abs(diff.Q.psi + qnorm(diff.Q.alpha)*ASE.diff.Q.psi)/(.1+abs(Q.psi))
      diff.psi <- psi - last.psi
      psi.max.diff <- max(abs(diff.psi))
      max.diff.id <- which(abs(diff.psi)==psi.max.diff)[1]
      psi.crit <- max(abs(diff.psi)/(.0001+last.psi))

    }

    cat("\nLast value of psi:  ",last.psi)
    cat("\nCurrent value of psi:",psi)
    cat("\nSignificance criterion:",zval.diff.Q.psi)
    cat("\nQ difference criterion:",crit.diff.Q.psi)
    cat("\nMax abs change in psi:",psi.max.diff)
    cat("\nMax rel change in psi:",psi.crit)

    diffQsign <- abs(zval.diff.Q.psi) > qnorm(1-diff.Q.alpha)
    if(crit.diff.Q.psi < diff.Q.eps && crit.diff.Q.psi > 0) maybe.converged <- TRUE
    if(psi.crit < rel.diff.psi.eps)  maybe.converged <- TRUE
    if(abs(psi.max.diff) < abs.diff.psi.eps)  maybe.converged <- TRUE

    if(diff.Q.psi<0){

      cat("\nCannot improve Q-function, stepping back")
      parm <- last.parm
      psi <- last.psi
      maybe.converged <- TRUE
    }
    else if(!diffQsign) {

      parm <- last.parm
      cat("\nSample size increase from",latent.data$sample.size)
      latent.data$sample.size <- ceiling(latent.data$sample.size*1.5)
      cat(" to",latent.data$sample.size)
      maybe.converged <- FALSE
    }

    ## First find the maximum of the integrand - for
    ## optimal importance sampling
    Utilde <- latpos.utilde(resp=resp,parm=parm,maxiter=maxiter,verbose=FALSE)
    parm$Utilde <- Utilde

    cat("\nGenerating Monte Carlo sample for the next step - size",latent.data$sample.size)
    ## The Monte-Carlo E-step: Simulating from
    ## the posterior distribution of the latent data

    gc();gc()
    latent.data <- latpos.simul(resp=resp,parm=parm,
                           latent.data=latent.data,
                           sampler=sampler)
    sample.size <- latent.data$sample.size
    parm$logLik <- logLik <- sum(latent.data$ll.j)

    diff.logLik <- logLik - last.logLik
    crit.logLik <- abs(diff.logLik)/abs(last.logLik)

    cat("\nCurrent log-likelihood:",logLik)
    cat(" - increase:",diff.logLik)
    cat(" - relative increase:",crit.logLik*sign(diff.logLik))

    
  }

  trace$Q[iter,1] <- diff.Q.psi
  trace$Q[iter,2] <- ASE.diff.Q.psi
  trace$Q[iter,3] <- used.sample.size
  trace$Q[iter,4] <- psi.crit
  trace$Q[iter,5] <- psi.max.diff
  trace$Q[iter,6] <- diff.psi[max.diff.id]
  trace$Q[iter,7] <- psi[max.diff.id]
  trace$Q[iter,8] <- logLik

  trace$psi[iter,] <- psi

  if(interactive()){
    #first.plot.iter <- 1
    if(iter < 11) first.plot.iter <- 1
    else first.plot.iter <- iter - 10
    plot.ii <- first.plot.iter:iter

    plot.trace(plot.ii,trace$Q[plot.ii,,drop=FALSE],trace$psi[plot.ii,,drop=FALSE])
  }

  if(maybe.converged){

    if(latent.data$sample.size >= min.final.size)
        converged<-TRUE
    else{

      cat("\nMCEM algorithm may have converged, setting sample size to",min.final.size,"to be sure")
      latent.data$sample.size <- min.final.size
    }
  }
  else{

    sample.size.new <- ceiling(sample.size.start*
                              (2*qnorm(1-abs(diff.Q.alpha))/zval.diff.Q.psi)^2
                            )
    if(sample.size.new > max.size){
      cat("\nMaximum sample size reached")
      sample.size.new <- max.size
      }

    if(sample.size.new > latent.data$sample.size){
      cat("\nNew sample size",sample.size.new)
      latent.data$sample.size <- sample.size.new
      }
  }

  if(converged){

    cat("\nMCEM algorithm has converged\n")
  }


  if(converged){
  trace$Q <- trace$Q[1:iter,,drop=FALSE]
  trace$psi <- trace$psi[1:iter,,drop=FALSE]

  colnames(trace$Q) <- c("diff Q", "ASE diff Q","Sample size",
                         "Max rel diff in psi","Max abs diff in psi",
                         "Max diff psi","Psi with max diff","Log-likelihood")
  }
  
  list(
    parm=parm,
    latent.data=latent.data,
    trace=trace,
    diff.logLik = diff.logLik,
    diff.psi = diff.psi,
    psi.crit = psi.crit,
    psi.max.diff = psi.max.diff,
    diff.Q.psi = diff.Q.psi,
    converged=converged
    )
}

latpos.Lambda.Q <- function(resp,latent.data,parm,last.parm){

    w.sim <- latent.data$w.sim
    U.sim <- latent.data$U.sim

    chunk.size <- getOption("latpos.chunk.size")
    batch.size <- chunk.size%/%(4*prod(dim(U.sim)[c(1,3)]))
    batch.size <- 2*(batch.size%/%2-1)
    m <- latent.data$sample.size %/% batch.size
    r <- latent.data$sample.size %% batch.size

    I <- nrow(resp$y)
    J <- length(unique(resp$j))
    JT <- length(resp$j)
    JTK <- JT*batch.size
    sample.size <- dim(U.sim)[2]
    D <- dim(U.sim)[3]

    j. <- resp$j

    Q <- 0
    sum.wLambda <- 0
    sum.sq.wLambda <- 0
    sum.w.wLambda <- 0
    sum.sq.w <- 0

    if(m>0){

      y <- rep(resp$y,batch.size)
      n <- rep(resp$n,batch.size)
      dim(y) <- dim(n) <- c(I,JTK)
      j <- rep(resp$j,batch.size) + rep(JT*(1:batch.size-1),each=JT)
      t <- rep(resp$t,batch.size)

      kk <- 1:batch.size
      for(k in 1:m){
        #cat("#")
        U <- array(U.sim[,kk,,drop=FALSE],c(JTK,D))
        w.kk <- w.sim[,kk,drop=FALSE]
        w <- w.kk[j.,]
        res <- latpos.eval.parms(y=y,n=n,j=j,t=t,
                                parm=parm,U=U,weights=w,
                                compute="logLik.j")
        last.res <- latpos.eval.parms(y=y,n=n,j=j,t=t,
                                parm=last.parm,U=U,weights=w,
                                compute="logLik.j")

        Q <- Q + sum(res$logLik.j)
        wLambda.kk <- res$logLik.j - last.res$logLik.j
        dim(wLambda.kk) <- dim(w.kk)
        sum.wLambda <- sum.wLambda + rowSums(wLambda.kk)
        sum.sq.wLambda <- sum.sq.wLambda + rowSums(wLambda.kk^2)
        sum.w.wLambda <- sum.w.wLambda + rowSums(w.kk*wLambda.kk)
        sum.sq.w <- sum.sq.w + rowSums(w.kk^2)

        kk <- kk + batch.size
      }
    }
    if(r>0){
      #cat("#")
        y <- rep(resp$y,r)
        n <- rep(resp$n,r)
        dim(y) <- dim(n) <- c(I,JT*r)
        j <- rep(resp$j,r) + rep(JT*(1:r-1),each=JT)
        t <- rep(resp$t,r)

        kk <- m*batch.size + 1:r
        U <- array(U.sim[,kk,,drop=FALSE],c(JT*r,D))
        w.kk <- w.sim[,kk,drop=FALSE]
        w <- w.kk[j.,]
        res <- latpos.eval.parms(y=y,n=n,j=j,t=t,
                                parm=parm,U=U,weights=w,
                                compute="logLik.j")
        last.res <- latpos.eval.parms(y=y,n=n,j=j,t=t,
                                parm=last.parm,U=U,weights=w,
                                compute="logLik.j")

        Q <- Q + sum(res$logLik.j)
        wLambda.kk <- res$logLik.j - last.res$logLik.j
        dim(wLambda.kk) <- dim(w.kk)
        sum.wLambda <- sum.wLambda + rowSums(wLambda.kk)
        sum.sq.wLambda <- sum.sq.wLambda + rowSums(wLambda.kk^2)
        sum.w.wLambda <- sum.w.wLambda + rowSums(w.kk*wLambda.kk)
        sum.sq.w <- sum.sq.w + rowSums(w.kk^2)
    }


    diff.Q <- sum(sum.wLambda)
    aVar.diff.Q <- sum(sum.sq.wLambda-2*sum.wLambda*sum.w.wLambda+sum.wLambda^2*sum.sq.w)

    ASE.diff.Q <- sqrt(aVar.diff.Q)

    list(Q=Q,
      diff.Q=diff.Q,
      ASE.diff.Q=ASE.diff.Q)
}
