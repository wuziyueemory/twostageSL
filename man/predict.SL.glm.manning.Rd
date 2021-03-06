% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SL.glm.manning.R
\name{predict.SL.glm.manning}
\alias{predict.SL.glm.manning}
\title{Prediction method for \code{SL.glm.manning}}
\usage{
\method{predict}{SL.glm.manning}(object, newdata, ...)
}
\arguments{
\item{object}{An object of class \code{SL.glm.manning}}

\item{newdata}{A \code{data.frame} to generate predictions for}

\item{...}{Other arguments (unused)}
}
\value{
A numeric vector of predictions
}
\description{
Prediction method for \code{SL.glm.manning}
}
\examples{
# load cost data
data(cost_data)
# fit manning model
fit_glm.manning <- SL.glm.manning(Y = cost_data$totalcost, X = cost_data[, c("female", "race")],
                              newX = cost_data[, c("female", "race")])
# get back predictions
pred_glm.manning <- predict(fit_glm.manning$fit, newdata = cost_data[,c("female", "race")])
}
