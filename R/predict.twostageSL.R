#' Predict method for twostageSL object
#'
#' Obtains predictions on a new data set from a \code{twostageSL} fit.  May require
#' the original data if one of the library algorithms uses the original data in
#' its predict method.
#'
#' If \code{newdata} is omitted the predicted values from \code{object} are
#' returned.  Each algorithm in the \code{library.2stage} and \code{library.1stage} needs to have a
#' corresponding prediction function with the ``predict.'' prefixed onto the
#' algorithm name (e.g. \code{predict.SL.glm} for \code{SL.glm}).
#'
#' @param object Fitted object from \code{twostageSL}
#' @param newdata New X values for prediction
#' @param X Original data set used to fit \code{object}, if needed by fit object.
#' @param Y Original outcome used to fit \code{object}, if needed by fit object.
#' @param onlySL Logical. If TRUE, only compute predictions for algorithms with
#' non-zero coefficients in the super learner object. Default is FALSE
#' (computes predictions for all algorithms in library).
#' @param \dots Additional arguments passed to the \code{predict.SL.*}
#' functions
#'
#' @return \item{pred}{ Predicted values from twostageSL fit}
#' \item{library.predict}{Predicted values for each algorithm in library}
#' @export
#' @author Ziyue Wu
#'
#' @seealso \code{\link{twostageSL}}
#'
#' @keywords models
predict.twostageSL <- function(object, newdata, X = NULL, Y = NULL,
                                 onlySL = FALSE, ...) {
  if (missing(newdata)) {
    out <- list(pred = object$SL.predict, library.predict = object$library.predict)
    return(out)
  }
  if (!object$control$saveFitLibrary) {
    stop("This SuperLearner fit was created using control$saveFitLibrary = FALSE, so new predictions cannot be made.")
  }
  N <- dim(object$Z)[[1]]
  k.stage1 <- object$library.Num$stage1
  k.stage2 <- object$library.Num$stage2
  k.singlestage <- object$library.Num$stage.single
  if (onlySL) {
    whichLibrary <- which(object$coef > 0)
    # predictions for only algorithm with coef > 0
    predY.stage1 <- matrix(NA, nrow = nrow(newdata), ncol = k.stage1)
    predY.stage2 <- matrix(NA, nrow = nrow(newdata), ncol = k.stage2)
    predY.singlestage <- matrix(NA, nrow = nrow(newdata), ncol = k.singlestage)
    # stage 1
    for (mm in seq(k.stage1)) {
      newdataMM.stage1 <- subset(newdata,
                                 select = object$whichScreen$stage1[object$SL.library$library$twostage[mm*k.stage2, 2], ])
      family.stage1 <- object$family$stage1
      XMM.stage1 <- if (is.null(X)) {
        NULL
      } else {
        subset(X, select = object$whichScreen$stage1[object$SL.library$library$twostage[mm*k.stage2, 2], ])
      }
      predY.stage1[, mm] <- do.call('predict', list(object = object$fitLibrary$stage1[[mm]],
                                                    newdata = newdataMM.stage1,
                                                    family = family.stage1,
                                                    X = XMM.stage1,
                                                    Y = as.numeric(Y==0)))

    }
    # stage 2
    Y.p <- Y[Y>0]
    for (mm in seq(k.stage2)) {
      newdataMM.stage2 <- subset(newdata,
                                 select = object$whichScreen$stage2[object$SL.library$library$twostage[mm, 3], ])
      family.stage2 <- object$family$stage2
      XMM.stage2 <- if (is.null(X)) {
        NULL
      } else {
        dat <- cbind(X,Y)
        X.p <- dat[dat$Y>0,-ncol(dat)]
        subset(X.p, select = object$whichScreen$stage2[object$SL.library$library$twostage[mm, 3], ])
      }
      predY.stage2[, mm] <- do.call('predict', list(object = object$fitLibrary$stage2[[mm]],
                                                    newdata = newdataMM.stage2,
                                                    family = family.stage2,
                                                    X = XMM.stage2,
                                                    Y = Y.p))
    }
    # single stage
    for (mm in seq(k.singlestage)) {
      newdataMM.singestage <- subset(newdata,
                                     select = object$whichScreen$single.stage[object$SL.library$library$singlestage[mm, 2], ])
      family.singlestage <- object$family$stage.single
      XMM.singlestage <- if (is.null(X)) {
        NULL
      } else {
        subset(X, select = object$whichScreen$single.stage[object$SL.library$library$singlestage[mm, 2], ])
      }
      predY.singlestage[, mm] <- do.call('predict', list(object = object$fitLibrary$stage.single[[mm]],
                                                         newdata = newdataMM.singestage,
                                                         family = family.singlestage,
                                                         X = XMM.singlestage,
                                                         Y = Y))
    }
    # get prediction for 2-stage model
    predY <- NULL
    for (i in 1:k.stage1){
      for (j in 1:k.stage2){
        pred <- (1-predY.stage1[,i])*predY.stage2[,j]
        predY <- cbind(predY,pred)
      }
    }
    # combine with prediction from singe-stage model
    predY <- cbind(predY,predY.singlestage)
    colnames(predY) <- object$libraryNames
    getPred <- object$method$computePred(predY = predY, coef = object$coef,
                                         control = object$control)
    predY.onlySL <- predY[,whichLibrary]
    out <- list(pred = getPred, library.predict = predY.onlySL)
    class(out) <- "predict.twostageSL"
    out
  } else {
    # predictions for all algorithms
    predY.stage1 <- matrix(NA, nrow = nrow(newdata), ncol = k.stage1)
    predY.stage2 <- matrix(NA, nrow = nrow(newdata), ncol = k.stage2)
    predY.singlestage <- matrix(NA, nrow = nrow(newdata), ncol = k.singlestage)
    # stage 1
    for (mm in seq(k.stage1)) {
      newdataMM.stage1 <- subset(newdata,
                          select = object$whichScreen$stage1[object$SL.library$library$twostage[mm*k.stage2, 2], ])
      family.stage1 <- object$family$stage1
      XMM.stage1 <- if (is.null(X)) {
        NULL
      } else {
        subset(X, select = object$whichScreen$stage1[object$SL.library$library$twostage[mm*k.stage2, 2], ])
      }
      predY.stage1[, mm] <- do.call('predict', list(object = object$fitLibrary$stage1[[mm]],
                                             newdata = newdataMM.stage1,
                                             family = family.stage1,
                                             X = XMM.stage1,
                                             Y = as.numeric(Y==0)))

    }
    # stage 2
    Y.p <- Y[Y>0]
    for (mm in seq(k.stage2)) {
      newdataMM.stage2 <- subset(newdata,
                                 select = object$whichScreen$stage2[object$SL.library$library$twostage[mm, 3], ])
      family.stage2 <- object$family$stage2
      XMM.stage2 <- if (is.null(X)) {
        NULL
      } else {
        dat <- cbind(X,Y)
        X.p <- dat[dat$Y>0,-ncol(dat)]
        subset(X.p, select = object$whichScreen$stage2[object$SL.library$library$twostage[mm, 3], ])
      }
      predY.stage2[, mm] <- do.call('predict', list(object = object$fitLibrary$stage2[[mm]],
                                                    newdata = newdataMM.stage2,
                                                    family = family.stage2,
                                                    X = XMM.stage2,
                                                    Y = Y.p))
    }
    # single stage
    for (mm in seq(k.singlestage)) {
      newdataMM.singestage <- subset(newdata,
                          select = object$whichScreen$single.stage[object$SL.library$library$singlestage[mm, 2], ])
      family.singlestage <- object$family$stage.single
      XMM.singlestage <- if (is.null(X)) {
        NULL
      } else {
        subset(X, select = object$whichScreen$single.stage[object$SL.library$library$singlestage[mm, 2], ])
      }
      predY.singlestage[, mm] <- do.call('predict', list(object = object$fitLibrary$stage.single[[mm]],
                                             newdata = newdataMM.singestage,
                                             family = family.singlestage,
                                             X = XMM.singlestage,
                                             Y = Y))
    }
    # get prediction for 2-stage model
    predY <- NULL
    for (i in 1:k.stage1){
      for (j in 1:k.stage2){
        pred <- (1-predY.stage1[,i])*predY.stage2[,j]
        predY <- cbind(predY,pred)
      }
    }
    # combine with prediction from singe-stage model
    predY <- cbind(predY,predY.singlestage)
    colnames(predY) <- object$libraryNames
    getPred <- object$method$computePred(predY = predY, coef = object$coef,
                                         control = object$control)
    out <- list(pred = getPred, library.predict = predY)
  }
  class(out) <- "predict.twostageSL"
  cat("Predictions from the two-stage Super Learner:  pred \n\nPrediction for all algorithms in the library:  library.predict\n")
  out
}

