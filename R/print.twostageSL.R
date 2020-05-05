#' @export
print.twostageSL <- function(x, ...) {
  cat("\nCall: ", deparse(x$call, width.cutoff = .9*getOption("width")), "\n\n", fill = getOption("width"))
  print(cbind(Risk = x$cvRisk, Coef = x$coef))
}
#' @export
coef.twostageSL <- function(object, ...) {
  object$coef
}
#' @export
print.CV.twostageSL <- function(x, ...) {
  cat("\nCall: ", deparse(x$call, width.cutoff = .9*getOption("width")), "\n\n", fill = getOption("width"))
  cat("Cross-validated predictions from the two-stage Super Learner:  SL.predict \n\nCross-validated predictions from the discrete super learner (cross-validation selector):  discreteSL.predict \n\nWhich library algorithm was the discrete super learner:  whichDiscreteSL \n\nCross-validated prediction for all algorithms in the library and standard super learner:  library.predict\n")
}
#' @export
coef.CV.twostageSL <- function(object, ...) {
  object$coef
}

