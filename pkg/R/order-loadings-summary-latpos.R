OrderLoadings <- function(x,...) UseMethod("OrderLoadings")
OrderLoadings.summary.latpos <- function(x,...){

  criterion <- x$tabs$A[,,1,drop=FALSE]
  criterion <- as.data.frame(criterion)

  ii <- do.call(order,rev(criterion))

  x$tabs$A <- x$tabs$A[ii,,,drop=FALSE]
  x
}
