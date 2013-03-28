latpos.simul <- function(resp,parm,latent.data,sampler){

  j <- resp$j
  u.j <- unique(j)
  J <- length(u.j)

  JT <- ncol(resp$y)
  ndim <- length(parm$latent.dims)

  sample.size <- latent.data$sample.size
  sample.size <- 2*(sample.size%/%2+sample.size%%2)
  chunk.size <- getOption("latpos.chunk.size")

  Btilde <- parm$Btilde

  if(length(latent.data$B.sim)){

    last.size <- dim(latent.data$B.sim)[2]
    if(sample.size>last.size){

      B.sim <- array(0,dim=c(JT,sample.size,ndim))
      w.sim <- array(0,dim=c(J,sample.size))
    }
    else{

      B.sim <- latent.data$B.sim
      w.sim <- latent.data$w.sim
    }
  }
  else {

    B.sim <- array(0,dim=c(JT,sample.size,ndim))
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

    Btilde.j <- Btilde$B[j.,,drop=FALSE]
    iK2.j <- Btilde$iK2[jD.,jD.,drop=FALSE]

    batch.size <- chunk.size%/%(length(Btilde.j)*4)
    batch.size <- 2*(batch.size%/%2)
    m <- sample.size %/% batch.size
    r <- sample.size %% batch.size

    max.log.w.tmp <- -Inf

    if(m > 0){

      kk <- 1:batch.size
      for(k in 1:m){
        #cat(".")
        Btmp <- simul.B.imp(y=y.j,n=n.j,j=j.j,t=t.j,parm=parm,
                            Btilde=Btilde.j,
                            iK2=iK2.j,
                            size=batch.size,
                            sampler=sampler)

        B.sim[j.,kk,] <- Btmp$B
        w.sim[jj,kk] <- Btmp$log.w
        max.log.w.tmp <- max(max.log.w.tmp,Btmp$log.w)
        kk <- kk + batch.size
      }
    }

    if(r > 0){

      #cat(".")
      kk <- m*batch.size + 1:r
      Btmp <- simul.B.imp(y=y.j,n=n.j,j=j.j,t=t.j,parm=parm,
                          Btilde=Btilde.j,
                          iK2=iK2.j,
                          size=r,
                          sampler=sampler)

      B.sim[j.,kk,] <- Btmp$B
      max.log.w.tmp <- max(max.log.w.tmp,Btmp$log.w)
      w.sim[jj,kk] <- Btmp$log.w
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

  list(B.sim=B.sim, ll.j=ll.j, w.sim=w.sim, sample.size=sample.size)
}


simul.B.imp <- function(y,n,j,t,parm,Btilde,iK2,size,sampler){

  D <- ncol(parm$A)
  I <- nrow(parm$A)
  Tj <- ncol(y)

  btilde <- as.vector(t(Btilde))
  TjD <- length(btilde)

  B <- sampler$sample(size=size,ndim=TjD)
  log.f.B <- sampler$log.density(B)

  B <- sweep(as.matrix(B%*%iK2),2,btilde,"+")

  dim(B) <- c(size,D,Tj)
  B <- aperm(B,c(3,1,2))

  y <- as.vector(y)
  y <- rep(y,size)
  dim(y) <- c(I,Tj*size)
  n <- rep(n,size)
  dim(n) <- dim(y)
  j <- rep(1:size,each=Tj)
  t <- rep(t,size)

  dimB <- dim(B)
  dim(B) <- c(prod(dimB[1:2]),dimB[3])

  ll <- latpos.eval.parms(y=y,n=n,j=j,t=t,parm=parm,B=B,
                    compute=c("logLik.j")) # Each replication is a "group"
  ll <- ll$logLik.j

  dim(B) <- dimB
  log.iw <- -log.f.B + determinant(iK2)$modulus

  list(
    B = B,
    ll.cpl = ll,
    log.iw = log.iw,
    log.w = ll + log.iw
  )
}

