% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SL.flexsurv.R
\name{predict.SL.flexsurv}
\alias{predict.SL.flexsurv}
\title{Prediction method for accelerated failure time models}
\usage{
\method{predict}{SL.flexsurv}(object, newdata, ...)
}
\arguments{
\item{object}{An object of class \code{SL.smearglm}}

\item{newdata}{A \code{data.frame} to generate predictions for}

\item{...}{Other arguments (unused)}
}
\value{
A vector of numeric predictions
}
\description{
The numerical integration step can
be unstable and we have tried to build in checks to guard against this. In particular, we
first try to integrate with upper limit = Inf, but if that fails move to 1e8,
which sometimes is able to provide a sane answer when upper limit = Inf fails.
The function keeps trying smaller and smaller values, but will not go smaller than 1e6.
In that case it returns a random number between 0 and \code{object$randomUpper}.
}
