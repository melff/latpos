latpos.CplInfo_A <- function(resp,parm,latent.data){

  res <- latpos.integ.ll(resp=resp,
                         parm=parm,
                         latent.data=latent.data,
                         compute=c("XWX"))

  Q.phi <- parm$Q.phi
  Information <- Q.phi%*%tcrossprod(res$XWX,Q.phi)

  list(
    Information=as.matrix(Information),
    restrictor=as.matrix(Q.phi)
    )
}


latpos.GradInfo_A <- function(resp,parm,latent.data){

  res <- latpos.integ.ll(resp=resp,
                         parm=parm,
                         latent.data=latent.data,
                         compute=c("XWX","G.j"))
  
  g.j <- res$g.j
  G.j <- res$G.j

  XWX <- res$XWX

  gg.j <- tcrossprod(g.j)
  G.j <- aperm(G.j,c(3,1,2))
  G.j <- colSums(G.j)

  gradient <- g.j
  Info.miss <- G.j - gg.j
  Info <- XWX - Info.miss

  Q.phi <- parm$Q.phi

  list(
    gradient=gradient,
    Information=as.matrix(Info),
    Info.cpl=as.matrix(XWX),
    Info.miss=as.matrix(Info.miss),
    restrictor=as.matrix(Q.phi)
    )
}


ARes.j <- function(R,w,repl){

  wR <- w*R
  R <- rowsum(R,repl)
  wR <- rowsum(wR,repl)
  RR <- crossprod(R,wR)
  R <- colSums(wR)

  res <- list(
    R = R,
    RR = RR
    )

  res
}