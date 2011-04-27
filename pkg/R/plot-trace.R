plot.trace <- function(x,y,z,plot.Q=TRUE){

    layout(matrix(1:6,nrow=2,ncol=3))
    if(plot.Q){

      plot(x,y[,1],type="n",
          xlab="Iteration",
          ylab="diff Q",
          ylim=range(c(y[,1]+qnorm(0.025)*y[,2],y[,1]+qnorm(0.975)*y[,2]))
          )
      segments(x0=x,
              x1=x,
              y0=y[,1]+qnorm(0.025)*y[,2],
              y1=y[,1]+qnorm(0.975)*y[,2]
              )
      lines(x,y[,1],type="b")
      #abline(h=Q.eps,col="red")
    } else {

      plot(1:2,1:2,type="n",axes=FALSE,ann=FALSE)
    }
    
    plot(x,y[,3],type="b",
         xlab="Iteration",
         ylab="Simulation sample size")
    plot(x,y[,4],type="b",
         xlab="Iteration",
         ylab="max rel diff psi")
    #abline(h=rel.psi.eps,col="red")

    plot(x,y[,5],type="b",
         xlab="Iteration",
         ylab="max abs diff psi")
    #abline(h=psi.eps,col="red")

    plot(x,y[,8],type="b",
         xlab="Iteration",
         ylab="Log-likelihood")

    plot(x,z[,1],type="l",
         xlab="Iteration",
         ylab="psi",
         ylim=range(z))
    for(i in 2:ncol(z))
         lines(x,z[,i],type="l",lty=i)


    layout(1)
}
