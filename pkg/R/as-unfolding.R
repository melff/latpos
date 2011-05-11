as.unfolding <- function(object,...) UseMethod("as.unfolding")
as.unfolding.latpos <- function(object,
        prediction=c("posterior modes","posterior means","multiple imputation"),
        biplot_type=c("text","density"),
        sample.size = object$latent.data$sample.size,
        batch.size=object$latent.data$sample.size,
        ...){

  if(is.character(prediction)){

    prediction <- match.arg(prediction)
    prediction <- predict(object,type=prediction,sample.size=sample.size,batch.size=batch.size)
  }

  if(is.list(prediction)){

    if(is.null(prediction$fit)) stop("prediction does not have a 'fit' element")
    B <- prediction$fit 
  }
  else if(is.array(prediction)){

    if(length(dim(prediction))==3){ # prediction type was "multiple imputation"

      B <- matrix(aperm(prediction,c(3,1,2)),ncol=ncol(prediction))
    }
    else if(is.matrix(prediction)){

      B <- prediction
    }
    else stop("prediction has too many dimensions")
  }
  else {
     B <- as.matrix(prediction)
  }

  res <- structure(list(
    A=object$parm$A,
    B=B
    ),
    class="unfolding"
    )
  attr(res,"plot_discrete") <- c(TRUE,FALSE)
  attr(res,"biplot_type") <- c("text","density")
  attr(res,"procrustes_use") <- "A"
  res
}