latpos.GradInfo_Abeta <- function(resp,parm,latent.data){

  res <- latpos.integ.ll(resp=resp,
                         parm=parm,
                         latent.data=latent.data,
                         compute=c("XWX","G.j"))
  
  g.j <- res$g.j
  G.j <- res$G.j
  wgwg.j <- res$wgwg.j
  wwg.j <- res$wwg.j
  ww.j <- res$ww.j

  XWX <- res$XWX

  gg.j <- tcrossprod(g.j)
  G.j <- aperm(G.j,c(3,1,2))
  ww.G.j <- colSums(ww.j*G.j)
  G.j <- colSums(G.j)
  
  wgwg.j <- rowSums(wgwg.j,dims=2)
  wwg.g.j <- tcrossprod(wwg.j,g.j)

  gradient <- rowSums(g.j)
  Info.miss <- G.j - gg.j
  Info <- XWX - Info.miss

  var.g <- wgwg.j - wwg.g.j - t(wwg.g.j) + ww.G.j

  Q.phi <- parm$Q.phi
  beta <- parm$beta
  free.beta <- parm$free.beta

  if(free.beta)
    Qmat <- bdiag(Q.phi,diag(nrow=length(beta)))
  else
    Qmat <- Q.phi

  list(
    gradient=as.vector(gradient),
    Information=as.matrix(Info),
    Info.cpl=as.matrix(XWX),
    Info.miss=as.matrix(Info.miss),
    var.gradient=as.matrix(var.g),
    restrictor=as.matrix(Qmat)
    )
}


AbetaRes.j <- function(R,w,repl){

  wR <- w*R
  R <- rowsum(R,repl)
  wR <- rowsum(wR,repl)
  w <- tapply(w,repl,"[",i=1)
  RR <- crossprod(R,wR)
  wRwR <- crossprod(wR)
  wwR <- drop(crossprod(w,wR))
  R <- colSums(wR)
  #wwR.R <- tcrossprod(wwR,R)
  ww <- sum(w^2)

  res <- list(
    R = R,
    RR = RR,
    #var.g.j = wRwR - wwR.R - t(wwR.R) + ww.RR
    wRwR = wRwR,
    wwR = wwR,
    ww = ww
    )

  res
}