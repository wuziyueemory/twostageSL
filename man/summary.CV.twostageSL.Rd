% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.CV.twostageSL.R
\name{summary.CV.twostageSL}
\alias{summary.CV.twostageSL}
\alias{print.summary.CV.twostageSL}
\title{Summary Function For Cross-Validated Two stage Super Learner}
\usage{
\method{summary}{CV.twostageSL}(object, obsWeights = NULL, ...)

\method{print}{summary.CV.twostageSL}(x, digits = max(2, getOption("digits") - 2), ...)
}
\arguments{
\item{object}{An object of class \code{CV.twostageSL}, the result of a call to \code{CV.twostageSL}.}

\item{obsWeights}{Optional vector for observation weights.}

\item{...}{additional arguments …}

\item{x}{An object of class \code{summary.CV.twostageSL}, the result of a call to \code{summary.CV.twostageSL}.}

\item{digits}{The number of significant digits to use when printing.}
}
\value{
\code{summary.CV.twostageSL} returns a list with components
\item{call}{The function call from \code{CV.twostageSL}.}
\item{method}{Describes the loss function used. Currently default is CC.LS.scale.}
\item{V}{Number of folds.}
\item{Risk.SL}{Risk estimate for the two stage super learner.}
\item{Risk.dSL}{Risk estimate for the discrete super learner (the cross-validation selector).}
\item{Risk.library}{A matrix with the risk estimates for each algorithm in the library.}
\item{Table}{A table with the mean risk estimate and standard deviation across the folds for the two stage super learner and all algorithms in the library.}
}
\description{
Summary method for the \code{CV.twostageSL} function
}
\details{
Summary method for \code{CV.twostageSL}. Calculates the V-fold cross-validated estimate of either the mean squared error or the -2*log(L) depending on the loss function used.
}
\seealso{
\code{\link{CV.twostageSL}}
}
\author{
Ziyue Wu
}
