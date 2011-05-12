latpos.simul <- function(resp,parm,latent.data,sampler){

  j <- resp$j
  u.j <- unique(j)
  J <- length(u.j)

  JT <- ncol(resp$y)
  ndim <- length(parm$latent.dims)

  sample.size <- latent.data$sample.size
  sample.size <- 2*(sample.size%/%2+sample.size%%2)
  chunk.size <- getOption("latpos.chunk.size")

  Utilde <- parm$Utilde

  if(length(latent.data$U.sim)){

    last.size <- dim(latent.data$U.sim)[2]
    if(sample.size>last.size){

      U.sim <- array(0,dim=c(JT,sample.size,ndim))
      w.sim <- array(0,dim=c(J,sample.size))
    }
    else{

      U.sim <- latent.data$U.sim
      w.sim <- latent.data$w.sim
    }
  }
  else {

    U.sim <- array(0,dim=c(JT,sample.size,ndim))
    w.sim <- array(0,dim=c(J,sample.size))
  }

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
    m <- sample.size %/% batch.size
    r <- sample.size %% batch.size

    max.log.w.tmp <- -Inf

    if(m > 0){

      kk <- 1:batch.size
      for(k in 1:m){
        #cat(".")
        Utmp <- simul.U.imp(y=y.j,n=n.j,j=j.j,t=t.j,parm=parm,
                            Utilde=Utilde.j,
                            iK2=iK2.j,
                            size=batch.size,
                            sampler=sampler)

        U.sim[j.,kk,] <- Utmp$U
        w.sim[jj,kk] <- Utmp$log.w
        max.log.w.tmp <- max(max.log.w.tmp,Utmp$log.w)
        kk <- kk + batch.size
      }
    }

    if(r > 0){

      #cat(".")
      kk <- m*batch.size + 1:r
      Utmp <- simul.U.imp(y=y.j,n=n.j,j=j.j,t=t.j,parm=parm,
                          Utilde=Utilde.j,
                          iK2=iK2.j,
                          size=r,
                          sampler=sampler)

      U.sim[j.,kk,] <- Utmp$U
      max.log.w.tmp <- max(max.log.w.tmp,Utmp$log.w)
      w.sim[jj,kk] <- Utmp$log.w
    }

    batch.size <- chunk.size%/%4
    batch.size <- 2*(batch.size%/%2)
    m <- sample.size %/% batch.size
    r <- sample.size %% batch.size

    #ll.sim[jj,] <- Utmp$ll


    sum.w.tmp <- 0
    if(m > 0){

      kk <- 1:batch.size
      for(k in 1:m){

        w.sim[jj,kk] <- exp(w.sim[jj,kk] - max.log.w.tmp)
        sum.w.tmp <- sum.w.tmp + sum(w.sim[jj,kk])
        kk <- kk + batch.size
      }
    }
    if(r > 0){

      kk <- m*batch.size + 1:r
      w.sim[jj,kk] <- exp(w.sim[jj,kk] - max.log.w.tmp)
      sum.w.tmp <- sum.w.tmp + sum(w.sim[jj,kk])
    }
    mean.w.tmp <- sum.w.tmp/sample.size
    ll.j[jj] <- log(mean.w.tmp)+max.log.w.tmp
    if(m > 0){

      kk <- 1:batch.size
      for(k in 1:m){
        w.sim[jj,kk] <- w.sim[jj,kk]/sum.w.tmp
        kk <- kk + batch.size
      }
    }
    if(r > 0){

      kk <- m*batch.size + 1:r
      w.sim[jj,kk] <- w.sim[jj,kk]/sum.w.tmp
    }
  }

  list(U.sim=U.sim, ll.j=ll.j, w.sim=w.sim, sample.size=sample.size)
}


simul.U.imp <- function(y,n,j,t,parm,Utilde,iK2,size,sampler){

  D <- ncol(parm$A)
  I <- nrow(parm$A)
  Tj <- ncol(y)

  utilde <- as.vector(t(Utilde))
  TjD <- length(utilde)

  U <- sampler$sample(size=size,ndim=TjD)
  log.f.U <- sampler$log.density(U)

  U <- sweep(as.matrix(U%*%iK2),2,utilde,"+")

  dim(U) <- c(size,D,Tj)
  U <- aperm(U,c(3,1,2))

  y <- as.vector(y)
  y <- rep(y,size)
  dim(y) <- c(I,Tj*size)
  n <- rep(n,size)
  dim(n) <- dim(y)
  j <- rep(1:size,each=Tj)
  t <- rep(t,size)
  repl <- rep(1,length(j))

  dimU <- dim(U)
  dim(U) <- c(prod(dimU[1:2]),dimU[3])

  ll <- latpos.eval.parms(y=y,n=n,j=j,t=t,parm=parm,U=U,
                    replications=repl,
                    compute=c("logLik.j")) # Each replication is a "group"
  ll <- ll$logLik.j

  dim(U) <- dimU
  log.iw <- -log.f.U + determinant(iK2)$modulus

  list(
    U = U,
    ll.cpl = ll,
    log.iw = log.iw,
    log.w = ll + log.iw
  )
}

