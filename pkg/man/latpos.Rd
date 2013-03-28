\name{latpos}
\alias{latpos}
\alias{print.latpos}
\alias{summary.latpos}
\alias{print.summary.latpos}
\alias{relabel.summary.latpos}
\alias{predict.latpos}
\title{Dynamic State-Space Models of Coded Political Texts}
\description{
 The function \code{latpos} allows to specify and estimate
 a dynamic state-space model of political texts and returns
 an object from which the (latent) positions, which actors take
 in their texts, can be predicted.
}
\usage{

  latpos(formula,data,subset,id,time,
         unfold.method="Schoenemann",start=NULL,
         sampler=mvt.sampler(df=7*length(latent.dims)),
         ...)

  \method{predict}{latpos}(object, newdata = NULL, id=NULL, time=NULL,
          type=c("posterior modes","posterior means","simulate"),
          se.fit=FALSE, interval=c("none","normal","percentile"), level=0.95,
          sample.size = object$sample.size,
          sampler=object$sampler,
          maxiter=100,...)
}
\arguments{
  \item{formula}{a formula describing the model.}
  \item{data}{an optional data frame that contains data to which the model is fitted.}
  \item{subset}{a logical vector; an optional condition that defines a subset of the
                optiona data frame to which the analysis is to be restricted.}
  \item{id}{name of a variable that identifies the actors.}
  \item{time}{name of a variable that identifies the occasions on which actors
    publish their political texts.}
  \item{unfold.method}{a text string, name of the method which is used in the unfolding
    procedure to obtain starting values.}
  \item{start}{an optional list giving starting values.}
  \item{sampler}{an object that generates random numbers to be used in
    the MCEM algorithm.}
  \item{\dots}{further optional arguments,
    that are passed to \code{\link{latpos.control}} and \code{\link{latpos.start}}
    by the function \code{latpos} and ignored by the \code{predict} method.
  }
  \item{object}{an object of class "latpos".}
  \item{newdata}{an optional data frame to predict actors positions from, based
    on the fitted model contained in \code{object}.}
  \item{type}{a character string identifying the type of predictions about
    the latent positions, either "posterior modes","posterior means", or "simulate".
    If the type is "posterior means" or "posterior modes", only one prediction
    is generated per coded text. If the type is "simulate",
    then \code{sample.size} predictions are generated from the
    empirical Bayes posterior distribution of the positions.
    }
  \item{se.fit}{logical, should standard errors be provided with the posterior means?}
  \item{interval}{a character vector. If "none", then no prediction intervals are
    returned, if "normal" than prediction intervals are contstructed based
      on normal quantiles, if "percentile", then predicion intervals are
      constructed based on the empirical percentiles of samples from the posterior.}
  \item{level}{a number between zero and one, the nominal coverage level of the
    prediction intervals.}
  \item{sample.size}{a positive number, the sample size for the simulated posterior means,
    or the number of imputed values generated from the posterior.
    }
  \item{maxiter}{a positive number, the maximum number of iterations used to compute
    posterior modes.}
}
\value{
  \code{latpos} returns an object of class "latpos" that contains parameter
  estimates etc.

  \code{predict.latpos} returns a vector or a matrix depending on the arguments
  provided to the function.
}
\examples{
\dontrun{

load(file="manifesto-counts.RData")
options(latpos.chunk.size=4*1024*1024) # determines the size of
# data chunks used in internal computations.

econ.latpos <- latpos(
   nationalize+controlecon+econplan+marketregul+incentives+econorthodox+freeenterp~econlr,
   data=manifesto.counts,
   id=party,
   time=year*100 + month,
   free.beta=TRUE,
   initial.size=50
 )
 
summary(econ.latpos)


nonmat.latpos <- latpos(
    laworder+tradmoralpos+tradmoralneg+natwaylifepos+
    democracy+freedmhumrights~authlib+tradmod,
    data=manifesto.counts,
    id=party,
    time=year*100 + month,
    initial.size=50,
    restrictions=set.parms.free(
      authlib=c("natwaylifepos","laworder","democracy","freedmhumrights"),
      tradmod=c("tradmoralpos","tradmoralneg","democracy","freedmhumrights")
    ),
    start = list(
          A=list(
              authlib=c(natwaylifepos=1.4,laworder=.7,democracy=-1),
              tradmod=c(tradmoralpos=1,tradmoralneg=-1)
            )
        )
)
summary(nonmat.latpos)
}
}