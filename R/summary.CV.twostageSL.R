#' Summary Function For Cross-Validated Two stage Super Learner
#'
#' Summary method for the \code{CV.twostageSL} function
#' @rdname summary.CV.twostageSL
#' @param object An object of class \code{CV.twostageSL}, the result of a call to \code{CV.twostageSL}.
#' @param x An object of class \code{summary.CV.twostageSL}, the result of a call to \code{summary.CV.twostageSL}.
#' @param obsWeights Optional vector for observation weights.
#' @param digits The number of significant digits to use when printing.
#' @param ... additional arguments …
#'
#' @details Summary method for \code{CV.twostageSL}. Calculates the V-fold cross-validated estimate of either the mean squared error or the -2*log(L) depending on the loss function used.
#'
#' @return \code{summary.CV.twostageSL} returns a list with components
#' \item{call}{The function call from \code{CV.twostageSL}.}
#' \item{method}{Describes the loss function used. Currently default is CC.LS.scale.}
#' \item{V}{Number of folds.}
#' \item{Risk.SL}{Risk estimate for the two stage super learner.}
#' \item{Risk.dSL}{Risk estimate for the discrete super learner (the cross-validation selector).}
#' \item{Risk.library}{A matrix with the risk estimates for each algorithm in the library.}
#' \item{Table}{A table with the mean risk estimate and standard deviation across the folds for the two stage super learner and all algorithms in the library.}
#' @export
#' @author Ziyue Wu
#' @seealso \code{\link{CV.twostageSL}}

summary.CV.twostageSL <- function(object, obsWeights = NULL, ...) {

  if ("env" %in% names(object)) {
    env = object$env
  } else {
    env = parent.frame()
  }
  # Default is "method.NNLS".
  method <- if (is.null(as.list(object$call)[["method"]])) {
    method <- "method.CC_LS.scale"
  } else if (is.symbol(as.list(object$call)[["method"]])) {
    method <- get(paste(as.list(object$call)[["method"]]), envir = env)
  } else {
    method <- as.list(object$call)[["method"]]
  }

  library.names <- colnames(coef(object))
  V <- object$V
  n <- length(object$SL.predict)
  if(is.null(obsWeights)) {
    obsWeights <- rep(1, length(object$Y))
  }
  # copy objects locally
  folds <- object$folds
  SL.predict <- object$SL.predict
  discreteSL.predict <- object$discreteSL.predict
  library.predict <- object$library.predict
  Y <- object$Y

  # create placeholders
  Risk.SL <- rep(NA, length = V)
  Risk.dSL <- rep(NA, length = V)
  Risk.library <- matrix(NA, nrow = length(library.names), ncol = V)
  rownames(Risk.library) <- library.names

  # Risk estimate depends on the loss function used
  if (method %in% c("method.NNLS", "method.NNLS2", "method.CC_LS", "method.CC_LS.scale")) {
    for (ii in seq_len(V)) {
      Risk.SL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] - SL.predict[folds[[ii]]])^2)
      Risk.dSL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] - discreteSL.predict[folds[[ii]]])^2)
      Risk.library[, ii] <- apply(library.predict[folds[[ii]], , drop = FALSE], 2, function(x) mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] - x)^2))
    }
    se <- (1 / sqrt(n)) * c(sd(obsWeights * (Y - SL.predict)^2), sd(obsWeights * (Y - discreteSL.predict)^2), apply(library.predict, 2, function(x) sd(obsWeights * (Y - x)^2)))
  } else if (method %in% c("method.NNloglik", "method.CC_nloglik")) {
    for (ii in seq_len(V)) {
      Risk.SL[ii] <- -mean(obsWeights[folds[[ii]]] * ifelse(Y[folds[[ii]]], log(SL.predict[folds[[ii]]]), log(1 - SL.predict[folds[[ii]]])))
      Risk.dSL[ii] <- -mean(obsWeights[folds[[ii]]] * ifelse(Y[folds[[ii]]], log(discreteSL.predict[folds[[ii]]]), log(1 - discreteSL.predict[folds[[ii]]])))
      Risk.library[, ii] <- apply(library.predict[folds[[ii]], , drop = FALSE], 2, function(x) {
        -mean(obsWeights[folds[[ii]]] * ifelse(Y[folds[[ii]]], log(x), log(1 - x)))
      })
    }
    se <- rep.int(NA, (length(library.names) + 2))
  } else if (method %in% c("method.AUC")) {
    requireNamespace("cvAUC") # make sure cvAUC loaded
    # require("ROCR") # cvAUC will load ROCR
    for (ii in seq_len(V)) {
      Risk.SL[ii] <- cvAUC::cvAUC(predictions = SL.predict[folds[[ii]]], labels = Y[folds[[ii]]], folds = NULL)$cvAUC
      Risk.dSL[ii] <- cvAUC::cvAUC(predictions = discreteSL.predict[folds[[ii]]], labels = Y[folds[[ii]]], folds = NULL)$cvAUC
      Risk.library[, ii] <- apply(library.predict[folds[[ii]], , drop = FALSE], 2, function(x) cvAUC::cvAUC(predictions = x, labels = Y[folds[[ii]]], folds = NULL)$cvAUC)
    }
    se <- rep.int(NA, (length(library.names) + 2)) # no se right now?
  } else {
    stop("summary function not available for SuperLearner with loss function/method used")
  }

  Table <- data.frame(Algorithm = c("Two-stage Super Learner", "Discrete SL", library.names), Ave = c(mean(Risk.SL), mean(Risk.dSL), apply(Risk.library, 1, mean)), se = se, Min = c(min(Risk.SL), min(Risk.dSL), apply(Risk.library, 1, min)), Max = c(max(Risk.SL), max(Risk.dSL), apply(Risk.library, 1, max)))
  out <- list(call = object$call, method = method, V = V, Risk.SL = Risk.SL, Risk.dSL = Risk.dSL, Risk.library = Risk.library, Table = Table)
  class(out) <- "summary.CV.twostageSL"
  out
}
#' @rdname summary.CV.twostageSL
#' @export
print.summary.CV.twostageSL <- function(x, digits = max(2, getOption("digits") - 2), ...) {
  cat("\nCall: ", deparse(x$call, width.cutoff = .9*getOption("width")), "\n", fill = getOption("width"))
  cat("Risk is based on: ")
  if(x$method %in% c("method.NNLS", "method.NNLS2", "method.CC_LS","method.CC_LS.scale")) {
    cat("Mean Squared Error")
  } else if (x$method %in% c("method.NNloglik", "method.CC_nloglik")) {
    cat("Negative Log Likelihood (-2*log(L))")
  } else if (x$method %in% c("method.AUC")) {
    cat("Area under ROC curve (AUC)")
  } else {
    stop("summary method not available")
  }
  cat("\n\nAll risk estimates are based on V = ", x$V, "\n\n")
  print(x$Table, digits = digits, row.names = FALSE)
}
