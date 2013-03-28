mvt.sampler <- function(df,fix.seed=TRUE){

  check.size <- function(size) (size%/%2)*2 + 1
  if(fix.seed){

    if(!exists(".Random.seed",globalenv())) runif(1)
    save.seed <- get(".Random.seed",globalenv())

  }

  reset <- function(){

      if(fix.seed){

        ## Reset random seed for smoothness
        ## (This is important especially for early small samples)
        set.seed(save.seed)
        message("\nResetting seed to ",save.seed[1]," ...")
      }
  }

  sample <- function(size,ndim){

      #size1 <- size-1
      if(size%%2) stop("need even sample size for anti-thetic variates")
      U <- matrix(NA,nrow=size,ncol=ndim)
      #size1 <- size1 %/% 2
      size1 <- size %/% 2
      ii <- (1:size1)*2
      #U[1,] <- 0
      U[ii,] <- rmvt(size1,sigma=diag(ndim),df=df)
      U[ii-1,] <- -U[ii,]
      U
  }

  log.density <- function(U) dmvt(U,df=df,log=TRUE)

  list(
    df=df,
    check.size=check.size,
    sample=sample,
    reset=reset,
    log.density=log.density
  )
}


mvnorm.sampler <- function(fix.seed=TRUE){

  ## alpha > 0 can be used to exclude rare events with
  ## with heavy importance weights that may drive up
  ## the Monte Carlo variance - that may be dubious but anyway ...

  check.size <- function(size) (size%/%2)*2 
  if(fix.seed){

    if(!exists(".Random.seed",globalenv())) runif(1)
    save.seed <- get(".Random.seed",globalenv())

  }

  reset <- function(){

      if(fix.seed){

        ## Reset random seed for smoothness
        ## (This is important especially for early small samples)
        set.seed(save.seed)
        message("\nResetting seed to",save.seed)
      }
  }

  sample <- function(size,ndim){

      #size1 <- size-1
      if(size%%2) stop("need even sample size for anti-thetic variates")
      U <- matrix(NA,nrow=size,ncol=ndim)
      #size1 <- size1 %/% 2
      size1 <- size %/% 2
      ii <- (1:size1)*2
      #U[1,] <- 0
      U[ii,] <- qnorm(runif(n=size1*ndim,min=0,max=1))
      U[ii-1,] <- -U[ii,]
      U
  }

  log.density <- function(U) rowSums(dnorm(U,log=TRUE))

  list(
    df=df,
    check.size=check.size,
    sample=sample,
    reset=reset,
    log.density=log.density
  )
}
