latpos.CplInfo_phi <- function(resp,parm,latent.data){

  res <- latpos.integ.ll(resp=resp,
                         parm=parm,
                         latent.data=latent.data,
                         compute=c("XWX"))

  Q.phi <- parm$Q.phi

  list(
    Information=as.matrix(res$XWX),
    restrictor=as.matrix(Q.phi)
    )
}

