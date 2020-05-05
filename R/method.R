# outline for two-stage Super Learner method
#
# The two-stage SuperLearner method is  the estimation algorithm for the algorithm weights (coefficients) and the model to combine the algorithms
#
# 2 parts need to be included:
#   1) compute coefficients (weights)
#   2) compute predictions

# scaled version of method.CC_LS

method.CC_LS.scale <- function() {
  computeCoef = function(Z, Y, libraryNames, verbose,
                         obsWeights=rep(1, length(Y)),
                         errorsInLibrary = NULL, ...) {
    # compute cvRisk
    cvRisk <- apply(Z, 2, function(x) mean(obsWeights*(x-Y)^2))
    names(cvRisk) <- libraryNames
    # compute coef
    compute <- function(x, y, wt=rep(1, length(y))) {
      wX <- sqrt(wt) * x
      wY <- sqrt(wt) * y
      D <- crossprod(wX)
      d <- crossprod(wX, wY)
      A <- cbind(rep(1, ncol(wX)), diag(ncol(wX)))
      bvec <- c(1, rep(0, ncol(wX)))
      sc <- norm(D,"2")
      # scale D matrix & d vector to aviod inconsistent constraints
      fit <- quadprog::solve.QP(Dmat=D/sc, dvec=d/sc, Amat=A, bvec=bvec, meq=1,
                                factorized = F)
      invisible(fit)
    }
    modZ <- Z
    # check for columns of all zeros. assume these correspond
    # to errors that SuperLearner sets equal to 0. not a robust
    # solution, since in theory an algorithm could predict 0 for
    # all observations (e.g., SL.mean when all Y in training = 0)
    naCols <- which(apply(Z, 2, function(z){ all(z == 0 ) }))
    anyNACols <- length(naCols) > 0
    if(anyNACols){
      # if present, throw warning identifying learners
      warning(paste0(paste0(libraryNames[naCols],collapse = ", "), " have NAs.",
                     "Removing from super learner."))
    }
    # check for duplicated columns
    # set a tolerance level to avoid numerical instability
    tol <- 8
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    if(anyDupCols){
      # if present, throw warning identifying learners
      warning(paste0(paste0(libraryNames[dupCols],collapse = ", "),
                     " are duplicates of previous learners.",
                     " Removing from super learner."))
    }
    # remove from Z if present
    if(anyDupCols | anyNACols){
      rmCols <- unique(c(naCols,dupCols))
      modZ <- Z[,-rmCols]
    }
    # compute coefficients on remaining columns
    fit <- compute(x = modZ, y = Y, wt = obsWeights)
    coef <- fit$solution
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] = 0
    }
    # add in coefficients with 0 weights for algorithms with NAs
    if(anyDupCols | anyNACols){
      ind <- c(seq_along(coef), rmCols - 0.5)
      coef <- c(coef, rep(0, length(rmCols)))
      coef <- coef[order(ind)]
    }
    # Set very small coefficients to 0 and renormalize.
    coef[coef < 1.0e-4] <- 0
    coef <- coef / sum(coef)
    if(!sum(coef) > 0) warning("All algorithms have zero weight", call. = FALSE)
    list(cvRisk = cvRisk, coef = coef, optimizer = fit)
  }

  computePred = function(predY, coef, ...) {
    predY %*% matrix(coef)
  }
  out <- list(require = "quadprog",
              computeCoef = computeCoef,
              computePred = computePred)
  invisible(out)
}
