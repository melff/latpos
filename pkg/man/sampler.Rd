\name{sampler}
\alias{mvt.sampler}
\alias{mvnorm.sampler}
\title{Control the Generation of Random Numbers for Sampling from the Proposal Distribution}
\description{
  The functions \code{mvt.sampler} and \code{mvnorm.sampler} can
  be used to fine-tune the construction and sampling from the proposal distribution
  in the importance sampling algorithm used in the MCEM algorithm to fit
  \code{\link{latpos}} models.
}
\usage{
  mvt.sampler(df,fix.seed=TRUE) # Multivariate Student's t-distribution
  mvnorm.sampler(fix.seed=TRUE) # Multivariate normal distribution
}
\arguments{
  \item{df}{a number, the degrees of freedom of the multivariate Student's t-distribution.}
  \item{fix.seed}{a logical value; if TRUE the seed of the random number generated
    is reset in each MCEM iteration to make sure that approximate
    marginal log-likelihood surface is continuous while the simulation sample
    size is constant.}
}
\value{
  Both functions return a list that contains a function to
  generate random numbers and to compute the log-density function of the
  distribution from which random numbers are generated.
}
