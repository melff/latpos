\name{latpos_tunefit}
\alias{latpos.control}
\alias{latpos.start}
\title{Starting Values and Fine Tuning of the MCEM algorithm}
\description{

  Function \code{latpos.start} is used to construct "good" starting values,
  while function \code{latpos.control} provides settings for 
  the numerical aspects of the MCEM algorithm, with reasonable defaults.
}
\usage{
  latpos.start(resp,latent.dims,manifest,start,
               unfold.method,restrictions=standard.restrictions,
               maxiter,...)
  latpos.control(maxiter=200,initial.size=101,
                            Lambda.alpha=.05,
                            Lambda.eps=1e-7,
                            diff.logLik.eps=1e-7,
                            abs.diff.psi.eps=0,
                            rel.diff.psi.eps=0,
                            max.size=Inf,
                            min.final.size=1000,
                            force.increase=TRUE,
                            Q.linesearch=TRUE,
                            ...)
}
\arguments{
  \item{resp}{an internal representation of the observed data.}
  \item{latent.dims}{a character vector with the names of the axes of
    the latent space.}
  \item{manifest}{a character vector with the names of the observed
    variables, i.e. emphasis counts of policy objectives.}
  \item{start}{an optional list with starting values for the model parameters}
  \item{unfold.method}{the unfolding method to be used to generate reasonable
    starting values.}
  \item{restrictions}{an object representing restrictions on the model
    parameters, see \code{\link{restrictor}}.}
  \item{maxiter}{the maximum number of iterations to use,
    in \code{latpos.start} to get initial values for the posterior modes,
    in \code{latpos.control} to set the maximum number of MCEM iterations.}
  \item{initial.size}{a positive number, the simulation sample size
    to use in the first MCEM iteration.}
  \item{Lambda.alpha}{a "significance level" for the increase of the
    Q-function. If the increase is not "statistically significant" at this level,
    the sample size is automatically increased.}
  \item{Lambda.eps}{a non-negative number as convergence critierion.
    If the increase of the Q-function is smaller than this value,
    convergence of the MCEM is declared.}
  \item{diff.logLik.eps}{a non-negative number as convergence critierion.
    If the increase of the marginal log-likelihood is smaller than this value,
    convergence of the MCEM is declared.}
  \item{abs.diff.psi.eps}{a non-negative number as an alternative convergence critierion.
    if the absolute change of the model parameters is smaller than this value,
    convergence of the MCEM is declared.}
  \item{rel.diff.psi.eps}{a non-negative number as an alternative convergence critierion.
    if the absolute change of the model parameters is smaller than this value,
    convergence of the MCEM is declared.}
  \item{max.size}{a positive number, the maximum simulation sample size to be used.}
  \item{min.final.size}{a positive number, the minimal simulation sample size to be
    used in the final iterations of the MCEM algorithm.}
  \item{force.increase}{logical; if TRUE and the likelihood or the Q-function cannot be increased
    then conduct a line search for the optimal step size.}
  \item{Q.linesearch}{logical; if TRUE,
    force.increase==TRUE and the likelihood or the Q-function cannot be increased
    then conduct a line search for the optimal step size;
    if FALSE, but force.increase==TRUE and the likelihood or the Q-function cannot be increased
    then step back to the values of the previous iteration.}
  \item{\dots}{other arguments, ignored.}
}
