latpos.AMstep <- function(resp,parm,latent.data,maxiter){

    A <- parm$A
    beta <- parm$beta
    Q.phi <- parm$Q.phi
    r.phi <- parm$r.phi

    JT <- ncol(resp$y)
    I <- nrow(resp$y)


    l.phi <- ncol(Q.phi)

    res <- latpos.integ.ll(resp=resp,parm=parm,latent.data=latent.data,compute=c("deviance","XWX","XWy"))
    dev <- res$deviance

    cat("\nM-step initial deviance:",dev)

    for(iter in 1:maxiter){

      cat("\nM-step Iteration for alpha",iter)

      last.parm <- parm
      XWX <- res$XWX
      XWy <- res$XWy
# browser()
      phi <- solve(XWX,XWy)
      #cat(" psi =",psi)
      A[] <- Q.phi %*% phi + r.phi
      parm$phi <- phi
      parm$A <- A

      res <- latpos.integ.ll(resp=resp,parm=parm,latent.data=latent.data,compute=c("deviance","XWX","XWy"))
      last.dev <- dev
      dev <- res$deviance
# browser()
      if(!is.finite(dev)){

        for(iiter in 1:maxiter){

          parm$A <- (last.parm$A + parm$A)/2
          parm$phi <- (last.parm$phi + parm$phi)/2

          res <- latpos.integ.ll(resp=resp,parm=parm,latent.data=latent.data,compute=c("deviance","XWX","XWy"))
          dev <- res$deviance
          cat("\nStep halved")
          if(is.finite(dev)){

            cat(" deviance:",dev)
            break
          }
          psi <- parm$phi
          if(parm$free.beta)
            psi <- c(psi,parm$beta)
          cat(" psi =",psi)
        }

        if(!is.finite(dev)) stop("could not find parameters for finite deviance")

      }
      if(dev > last.dev){

        cat("\nCould not decrease deviance - stepping back")
        parm <- last.parm
        break
      }

      crit <- abs(dev-last.dev)/abs(.1+last.dev)
      cat(" - deviance:",dev)
      cat(" criterion:",crit)
      if(is.finite(crit) && crit < 1e-7){

        cat("\nConverged")
        break
      }

    }

    return(parm)
}

