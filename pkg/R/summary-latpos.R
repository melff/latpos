vcov.latpos <- function(object,use.eigen=FALSE,...) {

  return(object$parm$covmat)
}

summary.latpos <- function(object,...){

  parm <- object$parm
    
  A <- parm$A
  Sigma0 <- parm$Sigma0
  Sigma1 <- parm$Sigma1
  Gamma <- parm$Gamma
  
  se.A <- sqrt(diag(object$covmat$A))

  zval.A <- as.vector(A)/se.A
  pval.A <- 2*pnorm(abs(zval.A),lower.tail=FALSE)
  zval.A[as.vector(A)==0] <- 0
  pval.A[zval.A==0] <- 1

  
  se.sigma1 <- sqrt(diag(object$covmat$sigma1))
  se.gamma <- sqrt(diag(object$covmat$gamma))

  if(!length(colnames(A))) colnames(A) <- paste("Dim",seq_len(ncol(A)),sep=".")
  A.tab <- cbind(as.vector(A),se.A,zval.A,pval.A)
  dim(A.tab) <- c(dim(A),4)
  dimnames(A.tab) <- c(dimnames(A),list(c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))


  sigma0 <- diag(Sigma0)
  Sigma0.tab <- cbind(sigma0,0)
  dim(Sigma0.tab) <- c(1,length(sigma0),2)
  dimnames(Sigma0.tab) <- list("Sigma_0",colnames(A),c("Estimate", "Std. Error"))
    
  sigma1 <- diag(Sigma1)
  Sigma1.tab <- cbind(sigma1,se.sigma1)
  dim(Sigma1.tab) <- c(1,length(sigma1),2)
  dimnames(Sigma1.tab) <- list("Sigma_1",colnames(A),c("Estimate", "Std. Error"))

  gamma <- diag(Gamma)
  zval.gamma <- gamma/se.gamma
  pval.gamma <- 2*pnorm(abs(zval.gamma),lower.tail=FALSE)
  Gamma.tab <- cbind(gamma,se.gamma,zval.gamma,pval.gamma)

  dim(Gamma.tab) <- c(1,length(gamma),4)
  dimnames(Gamma.tab) <- list("Autoregression",colnames(A),c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

  n <- object$resp$n
  y <- object$resp$y
  n <- n[1,]

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
    
  Sigma1.tab <- x$tabs$Sigma1
  Gamma.tab <- x$tabs$Gamma[,,1:2,drop=FALSE]

  tmp <- Sigma1.tab
  for(i in 1:ndim){
    Sigma1.tab[,i,1] <- formatEST(tmp[,i,1])
    Sigma1.tab[,i,2] <- formatSE(tmp[,i,2])
    Sigma1.tab[tmp[,i,2]==0,i,2] <- ""
    Sigma1.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
    }
  Sigma1.tab <- aperm(Sigma1.tab,c(1,3,2))
    
  tmp <- array(Sigma1.tab,c(dim(Sigma1.tab)[1],prod(dim(Sigma1.tab)[2:3])))
  Sigma1.tab <- cbind("Between timepoints",tmp)
  dimnames(Sigma1.tab) <- NULL

  Sigma0.tab <- cbind("Between units",t(rep(c("1",""),ndim)))
  dimnames(Sigma0.tab) <- NULL

  tmp <- c("",rep(statn,ndim))
  posPar.hdr <- rbind(c("",as.vector(rbind(dimn,rep("",ndim)))),
                tmp)
    
    posPar.tab <- rbind(
        posPar.hdr,
        c("Variance",rep("",ncol(Sigma0.tab)-1)),
        Sigma0.tab,
        Sigma1.tab)


    tmp <- Gamma.tab
    for(i in 1:ndim){
      Gamma.tab[,i,1] <- formatEST(tmp[,i,1])
      Gamma.tab[,i,2] <- formatSE(tmp[,i,2])
      Gamma.tab[tmp[,i,2]==0,i,2] <- ""
      Gamma.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
      }
    Gamma.tab <- aperm(Gamma.tab,c(1,3,2))

    tmp <- array(Gamma.tab,c(dim(Gamma.tab)[1],prod(dim(Gamma.tab)[2:3])))
    Gamma.tab <- cbind("Autoregress. coef.",tmp)
    
    posPar.tab <- rbind(posPar.tab,
                        rep("",ncol(Gamma.tab)),
                        Gamma.tab)

    
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
                    "Likelihood:",formatSumStat(x$logLik)),
                    byrow=TRUE,
                    nrow=1,ncol=2)

  totals <- matrix(c(
                    "N. of units:",x$total.units,
                    "N. of observations:",x$total.obs,
                    "N. of counts:",x$total.counts
                    ),byrow=TRUE,ncol=2)

  summaries <- rbind(sumstats,totals)
  summaries <- apply(summaries,2,format,justify="right")
  summaries <- apply(summaries,1,paste,collapse="  ")
    
  summaries <- c("Summary statistics:\n",summaries)
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
  tmp <- rbind(c( "Positions of objectives:",
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
  Gamma.tab <- object$tabs$Gamma[,,1:2,drop=FALSE]

  tmp <- Sigma0.tab
  for(i in 1:ndim){
    Sigma0.tab[,i,1] <- formatEST(tmp[,i,1])
    Sigma0.tab[,i,2] <- formatSE(tmp[,i,2])
    Sigma0.tab[tmp[,i,2]==0,i,2] <- ""
    Sigma0.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
    }
  Sigma0.tab <- aperm(Sigma0.tab,c(1,3,2))

  Sigma0.tab <- array(Sigma0.tab,c(dim(Sigma0.tab)[1],prod(dim(Sigma0.tab)[2:3])))
    Sigma0.tab <- cbind("Between units",Sigma0.tab)

  tmp <- Sigma1.tab
  for(i in 1:ndim){
    Sigma1.tab[,i,1] <- formatEST(tmp[,i,1])
    Sigma1.tab[,i,2] <- formatSE(tmp[,i,2])
    Sigma1.tab[tmp[,i,2]==0,i,2] <- ""
    Sigma1.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
    }
  Sigma1.tab <- aperm(Sigma1.tab,c(1,3,2))

  Sigma1.tab <- array(Sigma1.tab,c(dim(Sigma1.tab)[1],prod(dim(Sigma1.tab)[2:3])))
  Sigma1.tab <- cbind("Between timepoints",Sigma1.tab)
    
  posPar.tab <- rbind(c("Variance",rep("",ndim*2)),
                      Sigma0.tab,Sigma1.tab)


    tmp <- Gamma.tab
    for(i in 1:ndim){
      Gamma.tab[,i,1] <- formatEST(tmp[,i,1])
      Gamma.tab[,i,2] <- formatSE(tmp[,i,2])
      Gamma.tab[tmp[,i,2]==0,i,2] <- ""
      Gamma.tab[tmp[,i,1]==0 & tmp[,i,2]==0,i,1] <- ""
      }
    Gamma.tab <- aperm(Gamma.tab,c(1,3,2))

    tmp <- array(Gamma.tab,c(dim(Gamma.tab)[1],prod(dim(Gamma.tab)[2:3])))
    Gamma.tab <- cbind("Autoregress. coef.",tmp)
    
    posPar.tab <- rbind(posPar.tab,
                        rep("",ncol(Gamma.tab)),
                        Gamma.tab)
    
  tmp <- rbind(c("Positions of manifestos:",
                      rep(c("",""),ndim)),
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
                    "Likelihood:",formatSumStat(object$logLik)),
                    byrow=TRUE,
                    nrow=1,ncol=2)


  totals <- matrix(c(
                    "Number of units:",object$total.units,
                    "Number of observations:",object$total.obs,
                    "Sum of counts:",object$total.counts
                    ),byrow=TRUE,ncol=2)


  summaries <- rbind(sumstats,totals)
  summaries <- apply(summaries,2,format,justify="right")
  summaries <- apply(summaries,1,paste,collapse=" & ")

  summaries <- c(#"\\par",
                 #"\\begin{tabular}{lD{.}{.}{1}}",
                 #"\\toprule",
                 "\\midrule",
                 "\\multicolumn{2}{l}{Summary statistics:}\\\\",
                 paste(summaries,"\\\\",sep=""),
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

  tabs$Sigma0 <- colrename(tabs$Sigma0,...,warn=warn)
  tabs$Sigma0 <- rowrename(tabs$Sigma0,...,warn=warn)

  tabs$Sigma1 <- colrename(tabs$Sigma1,...,warn=warn)
  tabs$Sigma1 <- rowrename(tabs$Sigma1,...,warn=warn)

  x$tabs <- tabs

  return(x)
}
