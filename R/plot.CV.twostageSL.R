#' Graphical Display Of The V-Fold CV Risk Estimates
#'
#' The function plots the V-fold cross-validated risk estimates for the two stage super learner, the standard super learner, the discrete super learner and each algorithm in the library. By default the estimates will be sorted and include an asymptotic 95% confidence interval.
#'
#' @param x The output from \code{CV.twostageSL}.
#' @param package Either "ggplot2" or "lattice". The package selected must be available.
#' @param constant A numeric value. The confidence interval is defined as p +/- constant * se, where p is the point estimate and se is the standard error. The default is the quantile of the standard normal corresponding to a 95\% CI.
#' @param sort Logical. Should the rows in the plot be sorted from the smallest to the largest point estimate. If FALSE, then the order is two stage super learner, discrete super learner, then the estimators in \code{library.2stage} and \code{library.1stage} and the last is standard super learner.
#' @param ... Additional arguments for \code{summary.CV.twostageSL}
#'
#' @details see \code{summary.CV.twostageSL} for details on how the estimates are computed
#' @return Returns the plot (either a ggplot2 object (class \code{ggplot}) or a lattice object (class \code{trellis}))
#' @export
#' @author Ziyue Wu
#' @seealso \code{\link{summary.CV.twostageSL}} and \code{\link{CV.twostageSL}}
#'
plot.CV.twostageSL <- function(x, package = "ggplot2", constant = qnorm(0.975), sort = TRUE, ...) {
  sumx <- summary(x, ...)
  # if(sort) sumx$Table$Algorithm <- stats:::reorder.default(sumx$Table$Algorithm, -sumx$Table$Ave)\
  if(sort) sumx$Table$Algorithm <- reorder(sumx$Table$Algorithm, -sumx$Table$Ave)
  Mean <- sumx$Table$Ave
  se <- sumx$Table$se
  Lower <- Mean - constant*se
  Upper <- Mean + constant*se
  # d <- data.frame(Y = Mean, X = sumx$Table$Algorithm, Lower = Lower, Upper = Upper)
  assign("d", data.frame(Y = Mean, X = sumx$Table$Algorithm, Lower = Lower, Upper = Upper))

  if(package == "lattice") {
    .SL.require("lattice")
    p <- lattice::dotplot(X ~ Y, data = d, xlim = c(min(d$Lower) - 0.02, max(d$Upper) + 0.02), xlab = "V-fold CV Risk Estimate", ylab = "Method", panel = function(x, y){
      lattice::panel.xyplot(x, y, pch = 16, cex = 1)
      lattice::panel.segments(d$Lower, y, d$Upper, y, lty = 1)
    })
  }
  if(package == "ggplot2") {
    .SL.require("ggplot2")
    p <- ggplot2::ggplot(d, ggplot2::aes_string(x = "X", y = "Y", ymin = "Lower", ymax = "Upper")) + ggplot2::geom_pointrange() + ggplot2::coord_flip() + ggplot2::ylab("V-fold CV Risk Estimate") + ggplot2::xlab("Method")
  }
  return(p)
}
