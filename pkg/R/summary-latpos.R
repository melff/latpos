vcov_Abeta <- function(parm,use.eigen=FALSE,...){

  Info.Abeta <- parm$Information$Abeta
  restr.Abeta <- parm$Info.restr$Abeta

  cov.Abeta <- solve(Info.Abeta)
  if(any(diag(cov.Abeta)<0) && use.eigen){

    eigen.Abeta <- eigen(Info.Abeta)
    evals <- eigen.Abeta$values
    evecs <- eigen.Abeta$vectors

    evals <- 1/evals
    evals[evals < 0] <- 0

    cov.Abeta <- tcrossprod(evecs%*%diag(evals),evecs)

  }
  tcrossprod(restr.Abeta%*%cov.Abeta,restr.Abeta)
}

vcov_VarPar <- function(parm,use.eigen=FALSE,...){

  Info.VarPar <- parm$Information$VarPar
  restr.VarPar <- parm$Info.restr$VarPar
  Info.VarPar <- crossprod(restr.VarPar,Info.VarPar%*%restr.VarPar)

  cov.VarPar <- solve(Info.VarPar)

  if(any(diag(cov.VarPar)<0) && use.eigen){

    eigen.VarPar <- eigen(Info.VarPar)
    evals <- eigen.VarPar$values
    evecs <- eigen.VarPar$vectors

    evals <- 1/evals
    evals[evals < 0] <- 0

    cov.VarPar <- tcrossprod(evecs%*%diag(evals),evecs)

  }
  tcrossprod(restr.VarPar%*%cov.VarPar,restr.VarPar)
}



vcov.latpos <- function(object,use.eigen=FALSE,...) {

  parm <- object$parm

  cov.Abeta <- vcov_Abeta(parm,use.eigen=use.eigen,...)
  cov.VarPar <- vcov_VarPar(parm,use.eigen=use.eigen,...)

  covmat <- as.matrix(bdiag(cov.Abeta,cov.VarPar))

  return(covmat)
}

summary.latpos <- function(object,...){

  parm <- object$parm

  A <- parm$A
  beta <- parm$beta
  Sigma <- parm$Sigma
  rho <- parm$rho

  free.beta <- parm$free.beta
  free.rho <- parm$free.rho
  free.Sigma <- parm$free.Sigma

  zeta <- if(length(parm$zeta)) parm$zeta else 1/parm$tau

  l.A <- length(A)
  l.beta <- length(beta)
  l.rho <- length(rho)
  l.Sigma <- length(Sigma)

  i.A <- 1:l.A
  i.beta <- if(free.beta) l.A + 1:l.beta else integer(0)

  se.Abeta <- sqrt(diag(vcov_Abeta(parm,...)))
  se.A <- se.Abeta[i.A]
  se.beta <- if(free.beta) se.Abeta[i.beta] else 0*beta

  i.Sigma <- 1:l.Sigma
  i.rho <- if(free.rho) l.Sigma + 1 else integer(0)
  i.zeta <- if(free.rho) l.Sigma + 2 else l.Sigma + 1

  vcovVarPar <- vcov_VarPar(parm,...)
  se.VarPar <- sqrt(diag(vcovVarPar))
  se.Sigma <- se.VarPar[i.Sigma]
  se.rho <- if(free.rho) se.VarPar[i.rho] else 0
  se.zeta <- se.VarPar[i.zeta]

  Sigma1 <- zeta*Sigma
  i.Sigma.zeta <- c(i.Sigma,i.zeta)
  vcovSigmazeta <- vcovVarPar[i.Sigma.zeta,i.Sigma.zeta]
  d.Sigma1.d.sigma_zeta <- cbind(zeta*diag(nrow=length(Sigma)),as.vector(Sigma))
  vcov_Sigma1 <- d.Sigma1.d.sigma_zeta %*% tcrossprod(vcovSigmazeta,d.Sigma1.d.sigma_zeta)
  se.Sigma1 <- sqrt(diag(vcov_Sigma1))

  zval.A <- as.vector(A)/se.A
  pval.A <- 2*pnorm(abs(zval.A),lower.tail=FALSE)
  zval.A[as.vector(A)==0] <- 0
  pval.A[zval.A==0] <- 1

  zval.beta <- if(free.beta) beta/se.beta else 0*beta
  pval.beta <- if(free.beta) rep(1,length(beta)) else 2*pnorm(abs(zval.beta),lower.tail=FALSE)

  if(!length(colnames(A))) colnames(A) <- paste("Dim",seq_len(ncol(A)),sep=".")
  A.tab <- cbind(as.vector(A),se.A,zval.A,pval.A)
  dim(A.tab) <- c(dim(A),4)
  dimnames(A.tab) <- c(dimnames(A),list(c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))

  beta.tab <- cbind(beta,se.beta,zval.beta,pval.beta)
  if(!length(names(beta))) names(beta) <- colnames(A)
  dim(beta.tab) <- c(length(beta),4)
  dimnames(beta.tab) <- list(names(beta),c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

  Sigma.tab <- cbind(as.vector(Sigma),se.Sigma)
  dim(Sigma.tab) <- c(dim(Sigma),2)
  if(!length(dimnames(Sigma))){
      dimnames(Sigma) <- list(colnames(A),colnames(A))
  }
  dimnames(Sigma.tab) <- c(dimnames(Sigma),list(c("Estimate", "Std. Error")))

  rho.tab <- structure(c(rho,se.rho),names=c("Estimate", "Std. Error"))
  zeta.tab <- structure(c(zeta,se.rho),names=c("Estimate", "Std. Error"))

  Sigma1.tab <- cbind(as.vector(Sigma1),se.Sigma1)
  dim(Sigma1.tab) <- c(dim(Sigma1),2)
  if(!length(dimnames(Sigma1))){
      dimnames(Sigma1) <- list(colnames(A),colnames(A))
  }
  dimnames(Sigma1.tab) <- c(dimnames(Sigma1),list(c("Estimate", "Std. Error")))

  n <- object$resp$n
  y <- object$resp$y
  n.y <- n*y
  n.log.y <- n.y*log(y)
  n.log.y[y==0] <- 0
  n <- n[1,]
  constpart <- lfactorial(n) - colSums(lfactorial(n.y))
  ll <- object$parm$logLik - sum(constpart)
  
  deviance <- 2*(sum(n.log.y) - ll)

  total.obs <- ncol(y)
  total.units <- length(object$parm$Tj)
  total.counts <- sum(n)

  structure(
    c(object,
      list(
        deviance = deviance,
        total.obs = total.obs,
        total.units = total.units,
        total.counts = total.counts,
        tabs=list(
          A = A.tab,
          beta = beta.tab,
          Sigma = Sigma.tab,
          rho = rho.tab,
          zeta = zeta.tab,
          zeta.Sigma = Sigma1.tab
          )
        )),
     class=c("summary.latpos","latpos")
   )
}


printOld.summary.latpos <- function(x,...){

  cat("\nSpatial model of latent positions\n")

  print.default(x$call)

  cat("\nPositions of objectives:\n")
  print.default(x$tabs$A)
  cat("\nManifesto parameters:\n")
  cat("\nMean positions (beta):\n")
  print.default(x$tabs$beta)
  cat("\nSigma:\n")
  print.default(x$tabs$Sigma)
  cat("\nAutoregression coefficient (rho):\n")
  print.default(x$tabs$rho)
  cat("\nzeta:\n")
  print.default(x$tabs$zeta)
  cat("\nzeta x Sigma:\n")
  print.default(x$tabs$zeta.Sigma)
  invisible(x)
}

formatEST <- function(x) formatC(x,digits=getOption("digits"),width=-1,format="f")
formatSE <- function(x) paste("(",formatC(x,digits=getOption("digits"),width=-1,format="f"),")",sep="")
formatSumStat <- function(x) formatC(x,digits=1,width=-1,format="f")

print.summary.latpos <- function(x,...){

  cat("\nSpatial model of latent positions\n\n")

  cat("Call:\n  ")
  print.default(x$call)
  cat("\n")

  A.tab <- x$tabs$A[,,1:2,drop=FALSE]
  rown  <- dimnames(A.tab)[[1]]
  dimn  <- dimnames(A.tab)[[2]]
  statn <- dimnames(A.tab)[[3]]
  tmp <- A.tab

  ndim <- length(dimn)

  for(i in 1:ndim){
    A.tab[,i,1] <- formatEST(tmp[,i,1])
    A.tab[,i,2] <- formatSE(tmp[,i,2])
    A.tab[tmp[,i,2]==0,i,2] <- ""
    A.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
    }
  A.tab <- aperm(A.tab,c(1,3,2))

  tmp <- array(A.tab,c(dim(A.tab)[1],prod(dim(A.tab)[2:3])))
  tmp <- cbind(rown,tmp)
  tmp <- rbind(c("",rep(statn,ndim)),
               tmp)
  A.tab <- rbind(c("",as.vector(rbind(dimn,rep("",ndim)))),
                tmp)
  dimnames(A.tab) <- NULL

  Sigma.tab <- x$tabs$Sigma
  zSigma.tab <- x$tabs$zeta.Sigma

  if(x$parm$free.Sigma=="full"){
  
    tmp <- Sigma.tab
    for(i in 1:ndim){
      Sigma.tab[,i,1] <- formatEST(tmp[,i,1])
      Sigma.tab[,i,2] <- formatSE(tmp[,i,2])
      Sigma.tab[tmp[,i,2]==0,i,2] <- ""
      Sigma.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
      }
    Sigma.tab <- aperm(Sigma.tab,c(1,3,2))

    Sigma.tab <- array(Sigma.tab,c(dim(Sigma.tab)[1],prod(dim(Sigma.tab)[2:3])))
    Sigma.tab <- cbind(dimn,Sigma.tab)
    Sigma.tab <- rbind(c("Between units",rep("",ndim*2)),Sigma.tab)

    tmp <- zSigma.tab
    for(i in 1:ndim){
      zSigma.tab[,i,1] <- formatEST(tmp[,i,1])
      zSigma.tab[,i,2] <- formatSE(tmp[,i,2])
      zSigma.tab[tmp[,i,2]==0,i,2] <- ""
      zSigma.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
      }
    zSigma.tab <- aperm(zSigma.tab,c(1,3,2))

    zSigma.tab <- array(zSigma.tab,c(dim(zSigma.tab)[1],prod(dim(zSigma.tab)[2:3])))
    zSigma.tab <- cbind(dimn,zSigma.tab)
    zSigma.tab <- rbind(c("Between timepoints",rep("",ndim*2)),zSigma.tab)

    posPar.tab <- rbind(c("(Co-)Variance",rep("",ndim*2)),
                        Sigma.tab,zSigma.tab)

  }
  else{

    iii <- quick.grid(1:ndim,1:2)[,c(1,1,2)]

    Sigma.tab <- Sigma.tab[iii]
    zSigma.tab <- zSigma.tab[iii]
    dim(Sigma.tab) <- c(ndim,2)
    dim(zSigma.tab) <- c(ndim,2)
    
    tmp <- Sigma.tab
    for(i in 1:ndim){
      Sigma.tab[i,1] <- formatEST(tmp[i,1])
      Sigma.tab[i,2] <- formatSE(tmp[i,2])
      }
    
    tmp <- zSigma.tab
    for(i in 1:ndim){
      zSigma.tab[i,1] <- formatEST(tmp[i,1])
      zSigma.tab[i,2] <- formatSE(tmp[i,2])
      }

    posPar.tab <- rbind(c("Variance",rep("",ndim*2)),
                        c("Between units",t(Sigma.tab)),
                        c("Between timepoints",t(zSigma.tab))
                        )

  }

  if(x$parm$free.rho){

    rho.tab <- x$tabs$rho[1:2]
    tmp <- rho.tab
    rho.tab[1] <- formatEST(tmp[1])
    rho.tab[2] <- formatSE(tmp[2])
    tmp <- t(rep(rho.tab,ndim))

    rho.tab <- cbind("Autoregression slope",tmp)
    dimnames(rho.tab) <- NULL

    posPar.tab <- rbind(rho.tab,posPar.tab)
  }


  if(x$parm$free.beta){

    beta.tab <- x$tabs$beta[,1:2]
    tmp <- beta.tab
    for(i in 1:ndim){
      beta.tab[i,1] <- formatEST(tmp[i,1])
      beta.tab[i,2] <- formatSE(tmp[i,2])
      }
    tmp <- t(as.matrix(as.vector(t(beta.tab))))

    beta.tab <- cbind("Means",tmp)
    dimnames(beta.tab) <- NULL

    posPar.tab <- rbind(beta.tab,posPar.tab)
  }

  tmp <- rbind(c("",rep(statn,ndim)),
               posPar.tab)
  posPar.tab <- rbind(c("",as.vector(rbind(dimn,rep("",ndim)))),
                tmp)
  dimnames(posPar.tab) <- NULL


  nrows.A.tab <- nrow(A.tab)
  nrows.posPar.tab <- nrow(posPar.tab)

  i.A.tab <- 1:nrows.A.tab
  i.posPar.tab <- nrows.A.tab + 1:nrows.posPar.tab

  tab <- rbind(A.tab,posPar.tab)
  tab <- apply(tab,2,format,justify="right")
  tab <- apply(tab,1,paste,collapse="  ")

  tab <- c(
    "Positions of objectives:",
    "",
    paste("  ",tab[i.A.tab]),
    "",
    "Positions of manifestos:",
    "",
    paste("  ",tab[i.posPar.tab])
    )

  tab <- paste(tab,"\n",collapse="")
  cat(tab)

  sumstats <- matrix(c(
                    "Likelihood:",formatSumStat(x$parm$logLik),
                    "Deviance:",formatSumStat(x$deviance)),
                    byrow=TRUE,
                    nrow=2,ncol=2)

  #sumstats <- apply(sumstats,2,format,justify="right")
  #sumstats <- apply(sumstats,1,paste,collapse="  ")

  #sumstats <- c("Summary statistics:","",paste("  ",sumstats))
  #sumstats <- paste(sumstats,"\n",collapse="")

  #cat("\n")
  #cat(sumstats)

  totals <- matrix(c(
                    "N. of units:",x$total.units,
                    "N. of observations:",x$total.obs,
                    "N. of counts:",x$total.counts
                    ),byrow=TRUE,ncol=2)

  #totals <- apply(totals,2,format,justify="right")
  #totals <- apply(totals,1,paste,collapse="  ")

  #totals <- c("Totals:","",paste("  ",totals))
  #totals <- paste(totals,"\n",collapse="")
  #cat("\n")
  #cat(totals)

  summaries <- rbind(sumstats,totals)
  summaries <- apply(summaries,2,format,justify="right")
  summaries <- apply(summaries,1,paste,collapse="  ")

  sumstats <- paste("  ",summaries[1:2])
  totals <- paste("  ",summaries[2+1:3])

  summaries <- c("Summary statistics:","",
                 sumstats,
                 "Totals:","",
                 totals
                )
  summaries <- paste(summaries,"\n",collapse="")

  cat("\n")
  cat(summaries)
  
  invisible(x)
}

toLatex.summary.latpos <- function(x,...){

  out <- "%Spatial model of latent positions\n%\n%"


  A.tab <- x$tabs$A[,,1:2,drop=FALSE]
  rown  <- dimnames(A.tab)[[1]]
  dimn  <- dimnames(A.tab)[[2]]
  statn <- dimnames(A.tab)[[3]]
  tmp <- A.tab

  ndim <- length(dimn)

  for(i in 1:ndim){
    A.tab[,i,1] <- formatEST(tmp[,i,1])
    A.tab[,i,2] <- formatSE(tmp[,i,2])
    A.tab[tmp[,i,2]==0,i,2] <- ""
    A.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
    }
  A.tab <- aperm(A.tab,c(1,3,2))

  tmp <- array(A.tab,c(dim(A.tab)[1],prod(dim(A.tab)[2:3])))
  tmp <- cbind(rown,tmp)
  tmp <- rbind(c(#"",
                  "Positions of objectives:",
                  rep(
                      paste("\\multicolumn{1}{c}{",statn,"}",sep=""),
                      ndim)),
               tmp)
#
  A.tab <- if(ndim>1) rbind(c("",as.vector(rbind(dimn,rep("",ndim)))),
                            tmp)
           else tmp
  dimnames(A.tab) <- NULL

  Sigma.tab <- x$tabs$Sigma
  zSigma.tab <- x$tabs$zeta.Sigma

  if(x$parm$free.Sigma=="full"){

    tmp <- Sigma.tab
    for(i in 1:ndim){
      Sigma.tab[,i,1] <- formatEST(tmp[,i,1])
      Sigma.tab[,i,2] <- formatSE(tmp[,i,2])
      Sigma.tab[tmp[,i,2]==0,i,2] <- ""
      Sigma.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
      }
    Sigma.tab <- aperm(Sigma.tab,c(1,3,2))

    Sigma.tab <- array(Sigma.tab,c(dim(Sigma.tab)[1],prod(dim(Sigma.tab)[2:3])))
    Sigma.tab <- cbind(dimn,Sigma.tab)
    Sigma.tab <- rbind(c("Between units",rep("",ndim*2)),Sigma.tab)

    tmp <- zSigma.tab
    for(i in 1:ndim){
      zSigma.tab[,i,1] <- formatEST(tmp[,i,1])
      zSigma.tab[,i,2] <- formatSE(tmp[,i,2])
      zSigma.tab[tmp[,i,2]==0,i,2] <- ""
      zSigma.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
      }
    zSigma.tab <- aperm(zSigma.tab,c(1,3,2))

    zSigma.tab <- array(zSigma.tab,c(dim(zSigma.tab)[1],prod(dim(zSigma.tab)[2:3])))
    zSigma.tab <- cbind(dimn,zSigma.tab)
    zSigma.tab <- rbind(c("Between timepoints",rep("",ndim*2)),zSigma.tab)

    posPar.tab <- rbind(c("(Co-)Variance",rep("",ndim*2)),
                        Sigma.tab,zSigma.tab)

  }
  else{

    iii <- quick.grid(1:ndim,1:2)[,c(1,1,2)]

    Sigma.tab <- Sigma.tab[iii]
    zSigma.tab <- zSigma.tab[iii]
    dim(Sigma.tab) <- c(ndim,2)
    dim(zSigma.tab) <- c(ndim,2)

    tmp <- Sigma.tab
    for(i in 1:ndim){
      Sigma.tab[i,1] <- formatEST(tmp[i,1])
      Sigma.tab[i,2] <- formatSE(tmp[i,2])
      }

    tmp <- zSigma.tab
    for(i in 1:ndim){
      zSigma.tab[i,1] <- formatEST(tmp[i,1])
      zSigma.tab[i,2] <- formatSE(tmp[i,2])
      }

    posPar.tab <- rbind(c("Variance",rep("",ndim*2)),
                        c("Between units",t(Sigma.tab)),
                        c("Between timepoints",t(zSigma.tab))
                        )

  }

  if(x$parm$free.rho){

    rho.tab <- x$tabs$rho[1:2]
    tmp <- rho.tab
    rho.tab[1] <- formatEST(tmp[1])
    rho.tab[2] <- formatSE(tmp[2])
    tmp <- t(rep(rho.tab,ndim))

    rho.tab <- cbind("Autoregression slope",tmp)
    dimnames(rho.tab) <- NULL

    posPar.tab <- rbind(rho.tab,posPar.tab)
  }


  if(x$parm$free.beta){

    beta.tab <- x$tabs$beta[,1:2]
    tmp <- beta.tab
    for(i in 1:ndim){
      beta.tab[i,1] <- formatEST(tmp[i,1])
      beta.tab[i,2] <- formatSE(tmp[i,2])
      }
    tmp <- t(as.matrix(as.vector(t(beta.tab))))

    beta.tab <- cbind("Means",tmp)
    dimnames(beta.tab) <- NULL

    posPar.tab <- rbind(beta.tab,posPar.tab)
  }

  tmp <- rbind(c(#"",
                  "Positions of manifestos:",
                      rep(
                      paste("\\multicolumn{1}{c}{",statn,"}",sep=""),
                      ndim)),
               posPar.tab)
  posPar.tab <- tmp #if(ndim>1) rbind(c("",as.vector(rbind(dimn,rep("",ndim)))),
                #                 tmp)
                #else tmp
  dimnames(posPar.tab) <- NULL


  nrows.A.tab <- nrow(A.tab)
  nrows.posPar.tab <- nrow(posPar.tab)

  i.A.tab <- 1:nrows.A.tab
  i.posPar.tab <- nrows.A.tab + 1:nrows.posPar.tab

  tab <- rbind(A.tab,posPar.tab)
  tab <- apply(tab,2,format,justify="right")
  tab <- apply(tab,1,paste,collapse=" & ")

  A.tab <- tab[i.A.tab]
  posPar.tab <- tab[i.posPar.tab]

  A.tab <- paste(A.tab,"\\\\",sep="")
  A.tab <- if(ndim > 1) c(A.tab[1:2],"\\midrule",A.tab[-(1:2)])
           else  c(A.tab[1],"\\midrule",A.tab[-1])

  posPar.tab <- paste(posPar.tab,"\\\\",sep="")
  posPar.tab <- if(ndim > 1) c(posPar.tab[1:2],"\\midrule",posPar.tab[-(1:2)])
                else c(posPar.tab[1],"\\midrule",posPar.tab[-1])

  nc <- 2*length(dimn)
  nc1 <- nc+1

  head <- paste("\\begin{tabular}{",paste(c("l",rep("D{.}{.}{-1}",nc)),collapse=""),"}",sep="")
  tail <- paste("\\end{tabular}")

  tab <- c(
    head,
    "\\toprule",
    #paste("\\multicolumn{",nc,"}{l}{Positions of objectives:}",sep=""),
    #"Positions of objectives:",
    A.tab,
    "\\midrule",
    #paste("\\multicolumn{",nc,"}{l}{Positions of manifestos:}",sep=""),
    #"Positions of manifestos:",
    posPar.tab#,
    #"\\bottomrule",
    #tail
    )

  if(ndim>1)
    for(dn in dimn){

      pat <- paste(dn,"&")
      subst <- paste("\\multicolumn{2}{c}{",dn,"}",sep="")
      tab <- gsub(pat,subst,tab,fixed=TRUE)
    }

  out <- c(out,tab)

  sumstats <- matrix(c(
                    "Likelihood:",formatSumStat(x$parm$logLik),
                    "Deviance:",formatSumStat(x$deviance)),
                    byrow=TRUE,
                    nrow=2,ncol=2)


  totals <- matrix(c(
                    "N. of units:",x$total.units,
                    "N. of observations:",x$total.obs,
                    "N. of counts:",x$total.counts
                    ),byrow=TRUE,ncol=2)


  summaries <- rbind(sumstats,totals)
  summaries <- apply(summaries,2,format,justify="right")
  summaries <- apply(summaries,1,paste,collapse=" & ")

  sumstats <- summaries[1:2]
  totals <- summaries[2+1:3]

  summaries <- c(#"\\par",
                 #"\\begin{tabular}{lD{.}{.}{1}}",
                 #"\\toprule",
                 "\\midrule",
                 "\\multicolumn{2}{l}{Summary statistics:}\\\\",
                 paste(sumstats,"\\\\",sep=""),
                 "\\midrule",
                 "\\multicolumn{2}{l}{Totals:}\\\\",
                 paste(totals,"\\\\",sep=""),
                 "\\bottomrule",
                 "\\end{tabular}"
                )

  out <- c(out,summaries)

  structure(out,class="Latex")
}

relabel.summary.latpos <- function(x,...,gsub=FALSE,fixed=TRUE,warn=FALSE){

  tabs <- x$tabs

  tabs$A <- colrename(tabs$A,...,warn=warn)
  tabs$A <- rowrename(tabs$A,...,warn=warn)

  tabs$beta <- rowrename(tabs$beta,...,warn=warn)

  tabs$Sigma <- colrename(tabs$Sigma,...,warn=warn)
  tabs$Sigma <- rowrename(tabs$Sigma,...,warn=warn)

  tabs$zeta.Sigma <- colrename(tabs$zeta.Sigma,...,warn=warn)
  tabs$zeta.Sigma <- rowrename(tabs$zeta.Sigma,...,warn=warn)

  x$tabs <- tabs

  return(x)
}
