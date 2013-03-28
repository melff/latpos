\name{latpos_restrictions}
\alias{restrictor}
\alias{set.parms.free}
\title{Specify and Implement Restrictions on the Location Parameters}
\description{
  \code{restrictor} creates a linear transformation that maps a small set of
  linearly unrestriced parameters to a larger set of linearly restricted
  parameters. 

  \code{set.parms.free} specifies a call to \code{restrictor} in which
  only the specified parameters are allowed to differ from zero.
}
\usage{

  restrictor(C,d=numeric(m),sign=7)

  set.parms.free(...)
}
\arguments{
  \item{C}{
    The matrix involved in the linear restriction \eqn{C\alpha=d}{C*alpha=d}.
  }
  \item{d}{
    The vector involved in the linear restriction  \eqn{C\alpha=d}{C*alpha=d}.
  }
  \item{sign}{
    The number of significant digits to use for rounding to compensate
    finite machine precision in computing the QR decomposition.
  }
  \item{\dots}{
    several character vector arguments. Each character vector corresponds
    to one of the axes of the latent space, and each character string
    in a vector corresponds to the name of a policy objective that can
    obtain coordinate values different from zero.
  }
}
\value{
  The values of these functions are for internal use only. 

  If \eqn{C\alpha=d}{C*alpha=d} then \eqn{\alpha=Q\phi+r}{alpha=Q*phi+r}.
  The function \code{restrictor} returns a list with the
    elements "reduction" (which equals \eqn{Q}{Q}) and "offset" (which equals \eqn{r}{r})

  The function \code{set.parms.free} returns a function that generates
  arguments \eqn{C}{C} and \eqn{d}{d} with which the \code{restrictor}
  is called inside of the function \code{\link{latpos}}.
}
