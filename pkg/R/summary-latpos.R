vcov_A <- function(parm,use.eigen=FALSE,...){

  Info.A <- parm$Information$A
  restr.A <- parm$Info.restr$A

  cov.A <- solve(Info.A)
  if(any(diag(cov.A)<0) && use.eigen){

    eigen.A <- eigen(Info.A)
    evals <- eigen.A$values
    evecs <- eigen.A$vectors

    evals <- 1/evals
    evals[evals < 0] <- 0

    cov.A <- tcrossprod(evecs%*%diag(evals),evecs)

  }
  tcrossprod(restr.A%*%cov.A,restr.A)
}

vcov_LVdist <- function(parm,use.eigen=FALSE,...){

  Info.LVdist <- parm$Information$LVdist
  restr.LVdist <- parm$Info.restr$LVdist

  cov.LVdist <- solve(Info.LVdist)

  if(any(diag(cov.LVdist)<0) && use.eigen){

    eigen.LVdist <- eigen(Info.LVdist)
    evals <- eigen.LVdist$values
    evecs <- eigen.LVdist$vectors

    evals <- 1/evals
    evals[evals < 0] <- 0

    cov.LVdist <- tcrossprod(evecs%*%diag(evals),evecs)

  }
  tcrossprod(restr.LVdist%*%cov.LVdist,restr.LVdist)
}



vcov.latpos <- function(object,use.eigen=FALSE,...) {

  return(object$parm$covmat)
}

summary.latpos <- function(object,...){

  parm <- object$parm

  A <- parm$A
  beta <- parm$beta
  Sigma0 <- parm$Sigma0
  Sigma1 <- parm$Sigma1
  Gamma <- parm$Gamma
  Sigma1 <- Sigma1
  
  free.beta <- parm$free.beta
  free.Gamma <- parm$free.Gamma
  free.Sigma0 <- parm$free.Sigma0
  free.Sigma1 <- parm$free.Sigma1

  l.A <- length(A)
  l.beta <- length(beta)

  i.A <- 1:l.A

  se <- sqrt(diag(parm$covmat))
  
  se.A <- se[i.A]

  zval.A <- as.vector(A)/se.A
  pval.A <- 2*pnorm(abs(zval.A),lower.tail=FALSE)
  zval.A[as.vector(A)==0] <- 0
  pval.A[zval.A==0] <- 1


  l.Sigma0 <- length(Sigma0)
  l.Sigma1 <- length(Sigma1)
  l.Gamma <- length(Gamma)

  i.beta <- if(free.beta) l.A + 1:l.beta else integer(0)
  i.Gamma <- l.A + l.beta + 1:l.Gamma
  i.Sigma0 <- l.A + l.beta + l.Gamma + 1:l.Sigma0
  i.Sigma1 <- l.A + l.beta + l.Gamma + l.Sigma0 + 1:l.Sigma1
  
  se.beta <- if(free.beta) se[i.beta] else 0*beta
  zval.beta <- if(free.beta) beta/se.beta else 0*beta
  pval.beta <- if(free.beta) rep(1,length(beta)) else 2*pnorm(abs(zval.beta),lower.tail=FALSE)
  
  se.Sigma0 <- se[i.Sigma0]
  se.Sigma1 <- se[i.Sigma1]
  if(!free.Gamma)
    se.Gamma <- matrix(0,0,0)
  else
    se.Gamma <- se[i.Gamma]

  if(!length(colnames(A))) colnames(A) <- paste("Dim",seq_len(ncol(A)),sep=".")
  A.tab <- cbind(as.vector(A),se.A,zval.A,pval.A)
  dim(A.tab) <- c(dim(A),4)
  dimnames(A.tab) <- c(dimnames(A),list(c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))

  beta.tab <- cbind(beta,se.beta,zval.beta,pval.beta)
  if(!length(names(beta))) names(beta) <- colnames(A)
  dim(beta.tab) <- c(length(beta),4)
  dimnames(beta.tab) <- list(names(beta),c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

  Sigma0.tab <- cbind(as.vector(Sigma0),se.Sigma0)
  dim(Sigma0.tab) <- c(dim(Sigma0),2)
  if(!length(dimnames(Sigma0))){
      dimnames(Sigma0) <- list(colnames(A),colnames(A))
  }
  dimnames(Sigma0.tab) <- c(dimnames(Sigma0),list(c("Estimate", "Std. Error")))

  Sigma1.tab <- cbind(as.vector(Sigma1),se.Sigma1)
  dim(Sigma1.tab) <- c(dim(Sigma1),2)
  if(!length(dimnames(Sigma1))){
      dimnames(Sigma1) <- list(colnames(A),colnames(A))
  }
  dimnames(Sigma1.tab) <- c(dimnames(Sigma1),list(c("Estimate", "Std. Error")))

  if(length(se.Gamma)){
    Gamma.tab <- cbind(as.vector(Gamma),se.Gamma)
    dim(Gamma.tab) <- c(dim(Gamma),2)
    if(!length(dimnames(Gamma))){
        dimnames(Gamma) <- list(colnames(A),colnames(A))
    }
    dimnames(Gamma.tab) <- c(dimnames(Gamma),list(c("Estimate", "Std. Error")))
  }
  else{
    Gamma.tab <- NULL
  }

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
          Sigma0 = Sigma0.tab,
          Sigma1 = Sigma1.tab,
          Gamma = Gamma.tab
          )
        )),
     class=c("summary.latpos","latpos")
   )
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

  Sigma0.tab <- x$tabs$Sigma0
  Sigma1.tab <- x$tabs$Sigma1
  Gamma.tab <- x$tabs$Gamma


  tmp <- Sigma0.tab
  for(i in 1:ndim){
    Sigma0.tab[,i,1] <- formatEST(tmp[,i,1])
    Sigma0.tab[,i,2] <- formatSE(tmp[,i,2])
    Sigma0.tab[tmp[,i,2]==0,i,2] <- ""
    Sigma0.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
    }
  Sigma0.tab <- aperm(Sigma0.tab,c(1,3,2))

  Sigma0.tab <- array(Sigma0.tab,c(dim(Sigma0.tab)[1],prod(dim(Sigma0.tab)[2:3])))
  Sigma0.tab <- cbind(dimn,Sigma0.tab)
  Sigma0.tab <- rbind(c("Between units",rep("",ndim*2)),Sigma0.tab)

  tmp <- Sigma1.tab
  for(i in 1:ndim){
    Sigma1.tab[,i,1] <- formatEST(tmp[,i,1])
    Sigma1.tab[,i,2] <- formatSE(tmp[,i,2])
    Sigma1.tab[tmp[,i,2]==0,i,2] <- ""
    Sigma1.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
    }
  Sigma1.tab <- aperm(Sigma1.tab,c(1,3,2))

  Sigma1.tab <- array(Sigma1.tab,c(dim(Sigma1.tab)[1],prod(dim(Sigma1.tab)[2:3])))
  Sigma1.tab <- cbind(dimn,Sigma1.tab)
  Sigma1.tab <- rbind(c("Between timepoints",rep("",ndim*2)),Sigma1.tab)

  posPar.tab <- rbind(c("(Co-)Variance",rep("",ndim*2)),
                      Sigma0.tab,Sigma1.tab)


  if(x$parm$free.Gamma){

    tmp <- Gamma.tab
    for(i in 1:ndim){
      Gamma.tab[,i,1] <- formatEST(tmp[,i,1])
      Gamma.tab[,i,2] <- formatSE(tmp[,i,2])
      Gamma.tab[tmp[,i,2]==0,i,2] <- ""
      Gamma.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
      }
    Gamma.tab <- aperm(Gamma.tab,c(1,3,2))

    Gamma.tab <- array(Gamma.tab,c(dim(Gamma.tab)[1],prod(dim(Gamma.tab)[2:3])))
    Gamma.tab <- cbind(dimn,Gamma.tab)
    Gamma.tab <- rbind(c("Autoregression slope",rep("",ndim*2)),Gamma.tab)

    posPar.tab <- rbind(Gamma.tab,posPar.tab)

  }

  if(x$parm$free.beta){

    beta.tab <- x$tabs$beta[,1:2,drop=FALSE]
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

  totals <- matrix(c(
                    "N. of units:",x$total.units,
                    "N. of observations:",x$total.obs,
                    "N. of counts:",x$total.counts
                    ),byrow=TRUE,ncol=2)

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

toLatex.summary.latpos <- function(object,...){

  out <- "%Spatial model of latent positions\n%\n%"


  A.tab <- object$tabs$A[,,1:2,drop=FALSE]
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

  Sigma0.tab <- object$tabs$Sigma0
  Sigma1.tab <- object$tabs$Sigma1
  Gamma.tab <- object$tabs$Gamma

  tmp <- Sigma0.tab
  for(i in 1:ndim){
    Sigma0.tab[,i,1] <- formatEST(tmp[,i,1])
    Sigma0.tab[,i,2] <- formatSE(tmp[,i,2])
    Sigma0.tab[tmp[,i,2]==0,i,2] <- ""
    Sigma0.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
    }
  Sigma0.tab <- aperm(Sigma0.tab,c(1,3,2))

  Sigma0.tab <- array(Sigma0.tab,c(dim(Sigma0.tab)[1],prod(dim(Sigma0.tab)[2:3])))
  Sigma0.tab <- cbind(dimn,Sigma0.tab)
  Sigma0.tab <- rbind(c("Between units",rep("",ndim*2)),Sigma0.tab)

  tmp <- Sigma1.tab
  for(i in 1:ndim){
    Sigma1.tab[,i,1] <- formatEST(tmp[,i,1])
    Sigma1.tab[,i,2] <- formatSE(tmp[,i,2])
    Sigma1.tab[tmp[,i,2]==0,i,2] <- ""
    Sigma1.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
    }
  Sigma1.tab <- aperm(Sigma1.tab,c(1,3,2))

  Sigma1.tab <- array(Sigma1.tab,c(dim(Sigma1.tab)[1],prod(dim(Sigma1.tab)[2:3])))
  Sigma1.tab <- cbind(dimn,Sigma1.tab)
  Sigma1.tab <- rbind(c("Between timepoints",rep("",ndim*2)),Sigma1.tab)

  posPar.tab <- rbind(c("(Co-)Variance",rep("",ndim*2)),
                      Sigma0.tab,Sigma1.tab)


  if(object$parm$free.Gamma){

    tmp <- Gamma.tab
    for(i in 1:ndim){
      Gamma.tab[,i,1] <- formatEST(tmp[,i,1])
      Gamma.tab[,i,2] <- formatSE(tmp[,i,2])
      Gamma.tab[tmp[,i,2]==0,i,2] <- ""
      Gamma.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
      }
    Gamma.tab <- aperm(Gamma.tab,c(1,3,2))

    Gamma.tab <- array(Gamma.tab,c(dim(Gamma.tab)[1],prod(dim(Gamma.tab)[2:3])))
    Gamma.tab <- cbind(dimn,Gamma.tab)
    Gamma.tab <- rbind(c("Autoregression slope",rep("",ndim*2)),Gamma.tab)

    posPar.tab <- rbind(Gamma.tab,posPar.tab)

  }

  if(object$parm$free.beta){

    beta.tab <- object$tabs$beta[,1:2,drop=FALSE]
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

  if(ndim>1)
    for(dn in dimn){

      pat <- paste(dn,"&")
      subst <- paste("\\multicolumn{2}{c}{",dn,"}",sep="")
      A.tab <- gsub(pat,subst,A.tab,fixed=TRUE)
    }
    
  posPar.tab <- paste(posPar.tab,"\\\\",sep="")
  posPar.tab <- c(posPar.tab[1],"\\midrule",posPar.tab[-1])

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

  out <- c(out,tab)

  sumstats <- matrix(c(
                    "Likelihood:",formatSumStat(object$parm$logLik),
                    "Deviance:",formatSumStat(object$deviance)),
                    byrow=TRUE,
                    nrow=2,ncol=2)


  totals <- matrix(c(
                    "Number of units:",object$total.units,
                    "Number of observations:",object$total.obs,
                    "Sum of counts:",object$total.counts
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

  tabs$Sigma0 <- colrename(tabs$Sigma0,...,warn=warn)
  tabs$Sigma0 <- rowrename(tabs$Sigma0,...,warn=warn)

  tabs$Sigma1 <- colrename(tabs$Sigma1,...,warn=warn)
  tabs$Sigma1 <- rowrename(tabs$Sigma1,...,warn=warn)

  x$tabs <- tabs

  return(x)
}
