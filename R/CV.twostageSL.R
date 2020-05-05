#' Function to get V-Fold cross-validated risk estimate for two stage super learner
#'
#' Function to get V-fold cross-validated risk estimate for two stage super learner. This function simply splits the data into V folds and then calls twostageSL. Most of the arguments are passed directly to twostageSL.
#'
#' @param Y The outcome.
#' @param X The covariates.
#' @param V The number of folds for \code{CV.SuperLearner}. This argument will be depreciated and moved into the \code{cvControl}. If Both \code{V} and \code{cvControl} set the number of cross-validation folds, an error message will appear. The recommendation is to use \code{cvControl}. This is not the number of folds for \code{twostageSL}. The number of folds for \code{twostageSL} is controlled with \code{innerCvControl}.
#' @param family.1 Error distribution of the stage 1 outcome for two-stage super learner. Currently only allows \code{binomial} to describe the error distribution. Link function information will be ignored and should be contained in the method argument below.
#' @param family.2 Error distribution of the stage 2 outcome for two-stage super learner. Currently only allows \code{gaussian} to describe the error distribution. Link function information will be ignored and should be contained in the method argument below.
#' @param family.single Error distribution of the outcome for standard super learner. Currently only allows \code{gaussian} to describe the error distribution. Link function information will be ignored and should be contained in the method argument below.
#' @param library.2stage Candidate prediction algorithms in two-stage super learner. A list containing prediction algorithms at stage 1 and stage 2, the prediction algorithms are either a character vector or a list containing character vectors. See details below for examples on the structure. A list of functions included in the SuperLearner package can be found with \code{listWrappers()}.
#' @param library.1stage Candidate prediction algorithms in standard super learner. Either a character vector of prediction algorithms or a list containing character vectors. See details below for examples on the structure. A list of functions included in the SuperLearner package can be found with \code{listWrappers()}.
#' @param twostage logical; TRUE for implementing two-stage super learner; FALSE for implementing standatd super learner
#' @param method Details on estimating the coefficients for the two-stage super learner and the model to combine the individual algorithms in the library. Currently, the built in option is only "method.CC_LS.scale" (default) which is a scaled version of CC_LS. CC_LS.scale uses Goldfarb and Idnani's quadratic programming algorithm to calculate the best convex combination of weights to minimize the squared error loss. In addition, CC_LS.scale divides the quadratic function by a large constant to shrink the huge matrix and vector in quadratic function.
#' @param id Optional cluster identification variable. For the cross-validation splits, \code{id} forces observations in the same cluster to be in the same validation fold. \code{id} is passed to the prediction and screening algorithms in library.2stage and library.1stage, but be sure to check the individual wrappers as many of them ignore the information.
#' @param verbose logical; TRUE for printing progress during the computation (helpful for debugging).
#' @param control A list of parameters to control the estimation process. Parameters include \code{saveFitLibrary} and \code{trimLogit}. See \code{\link{twostageSL.control}} for details.
#' @param cvControl A list of parameters to control the cross-validation process. The outer cross-validation is the sample spliting for evaluating the \code{twostageSL}. Parameters include \code{V}, \code{stratifyCV}, \code{shuffle} and \code{validRows}. See \code{\link{twostageSL.CV.control}} for details.
#' @param innerCvControl A list of lists of parameters to control the inner cross-validation process. It should have \code{V} elements in the list, each a valid \code{cvControl} list. If only a single value, then replicated across all folds. The inner cross-validation are the values passed to each of the \code{V} \code{SuperLearner} calls. Parameters include \code{V}, \code{stratifyCV}, \code{shuffle} and \code{validRows}. See \code{\link{twostageSL.CV.control}} for details.
#' @param obsWeights Optional observation weights variable. As with \code{id} above, \code{obsWeights} is passed to the prediction and screening algorithms, but many of the built in wrappers ignore (or can't use) the information. If you are using observation weights, make sure the library you specify uses the information.
#' @param saveAll Logical; Should the entire \code{twostageSL} object be saved for each fold?
#' @param parallel Options for parallel computation of the V-fold step. Use "seq" (the default) for sequential computation. \code{parallel = 'multicore'} to use \code{mclapply} for the V-fold step (but note that \code{twostageSL()} will still be sequential). The default for \code{mclapply} is to check the \code{mc.cores} option, and if not set to default to 2 cores. Be sure to set \code{options()$mc.cores} to the desired number of cores if you don't want the default. Or \code{parallel} can be the name of a snow cluster and will use \code{parLapply} for the V-fold step. For both multicore and snow, the inner \code{twostageSL} calls will be sequential.
#' @param env Environment containing the learner functions. Defaults to the calling environment.
#'
#' @details The \code{twostageSL} function builds a estimator, but does not contain an estimate on the performance of the estimator. Various methods exist for estimator performance evaluation. If you are familiar with the super learner algorithm, it should be no surprise we recommend using cross-validation to evaluate the honest performance of the two stage super learner estimator. The function \code{CV.twostageSL} computes the usual V-fold cross-validated risk estimate for the two stage super learner (and all algorithms in \code{library.2stage} and \code{library.1stage} for comparison).
#'
#' @return An object of class \code{CV.twostageSL} (a list) with components:
#' \item{call}{The matched call.}
#' \item{AllSL}{If \code{saveAll = TRUE}, a list with output from each call to \code{twostageSL}, otherwise NULL.}
#' \item{SL.predict}{The predicted values from the two stage super learner when each particular row was part of the validation fold.}
#' \item{discreteSL.predict}{The traditional cross-validated selector. Picks the algorithm with the smallest cross-validated risk (in super learner terms, gives that algorithm coefficient 1 and all others 0).}
#' \item{whichDiscreteSL}{A list of length \code{V}. The elements in the list are the algorithm that had the smallest cross-validated risk estimate for that fold.}
#' \item{library.predict}{A matrix with the predicted values from each algorithm in \code{library.2stage} and \code{library.1stage}. The columns are the algorithms in \code{library.2stage} and \code{library.1stage} and the rows represent the predicted values when that particular row was in the validation fold (i.e. not used to fit that estimator).}
#' \item{coef}{A matrix with the coefficients for the two stage super learner on each fold. The columns are the algorithms in \code{library.2stage} and \code{library.1stage} the rows are the folds.}
#' \item{folds}{A list containing the row numbers for each validation fold.}
#' \item{V}{Number of folds for \code{CV.twostageSL}.}
#' \item{number0}{A dataframe indicating the number of zeros in each of the \code{v} fold.}
#' \item{libraryNames}{A character vector with the names of the algorithms in the library. The format is 'predictionAlgorithm_screeningAlgorithm' with '_All' used to denote the prediction algorithm run on all variables in X.}
#' \item{SL.library}{Returns the prediction algorithms and screening algorithms in \code{library.2stage} and \code{library.1stage}.}
#' \item{method}{A list with the method functions.}
#' \item{Y}{The outcome.}
#'
#' @author Ziyue Wu
#' @seealso \code{\link{twostageSL}}.
#'
#' @export
#' @import SuperLearner
#' @import slcost
#'
#' @examples
#' ## simulate data
#' set.seed(12321)
#'
#' ## training set
#' n <- 10000
#' p <- 5
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' colnames(X) <- paste("X", 1:p, sep="")
#' X <- data.frame(X)
#' Y <- rep(NA,n)
#' ## probability of outcome being zero
#' prob <- plogis(1 + X[,1] + X[,2] + X[,1]*X[,2])
#' g <- rbinom(n,1,prob)
#' ## assign zero outcome
#' ind <- g==0
#' Y[ind] <- 0
#' ## assign non-zero outcome
#' ind <- g==1
#' Y[ind] <- 10 + X[ind, 1] + sqrt(abs(X[ind, 2] * X[ind, 3])) + X[ind, 2] - X[ind, 3] + rnorm(sum(ind))
#'
#' ## run the CV.twostageSL
#' cv_sl <- CV.twostageSL(
#'  Y = Y, X = X,
#'  family.1 = binomial,
#'  family.2 = gaussian,
#'  family.single = gaussian,
#'  library.2stage = list(stage1=c("SL.glm","SL.mean","SL.earth"),
#'                        stage2=c("SL.glm","SL.mean","SL.earth")),
#'  library.1stage = c("SL.glm","SL.mean","SL.earth"),
#'  cvControl = list(V = 2),
#'  innerCvControl = list(list(V = 5),
#'                        list(V = 5))
#' )
#' cv_sl$AllSL
#' cv_sl$whichDiscreteSL
#' plot(cv_sl)
#'

CV.twostageSL <- function(Y, X, V = NULL,
                          family.1, family.2, family.single,
                          library.2stage,library.1stage,twostage,
                          method = "method.CC_LS.scale", id = NULL, verbose = FALSE,
                          control = list(saveFitLibrary = FALSE), cvControl = list(),
                          innerCvControl = list(), obsWeights = NULL, saveAll = TRUE,
                          parallel = "seq", env = parent.frame()) {
  call <- match.call()
  # observations & columns for covariates
  N <- dim(X)[1L]
  p <- dim(X)[2L]

  # create CV folds:
  # (1) for comparing different algorithms (include SL)
  if(any(names(cvControl) == "V") & !is.null(V)) {
    stop(paste0("You specified a value for V and a value in the cvControl,
                please only use one, preferably the cvControl"))
  }
  cvControl <- do.call('twostageSL.CV.control', cvControl)
  if(!is.null(V)) {
    # if the user specified V in the function call, override the default in cvControl
    # backward compatibility to not remove the V
    cvControl$V <- V
  }
  # ensure each folds have approximately equal number of obs with y=0
  V <- cvControl$V  # save this because it appears in the output value
  ord <- order(Y)
  cvfold <- rep(c(1:V,V:1),N)[1:N]
  folds <- split(ord, factor(cvfold))
  # check
  tab <- rep(NA,V)
  for (i in 1:V) {
    tab[i] <- sum(Y[folds[[i]]]==0)
  }
  num.0 <- data.frame("fold"=paste0("fold ",c(1:cvControl$V)),"number.of.0"=tab)

  cvControl$validRows = folds

  # (2) for calculating the twostageSL (inner loop)
  if(length(innerCvControl) > 0) {
    if(length(innerCvControl) == 1) {
      warning("Only a single innerCvControl is given, will be replicated across all
              cross-validation split calls to SuperLearner")
      newInnerCvControl <- vector("list", cvControl$V)
      for(ii in seq(cvControl$V)) {
        newInnerCvControl[[ii]] <- unlist(innerCvControl, recursive = FALSE)
      }
      innerCvControl <- newInnerCvControl  # write over previous with replicated list
    }
    if(length(innerCvControl) != cvControl$V) stop("innerCvControl must be a list with
                                                   V cvControl lists")
  } else {
    innerCvControl <- vector("list", cvControl$V)  # if no innerCvControl is given, generate an empty list
    for(ii in seq(cvControl$V)) {
      innerCvControl[[ii]] <- list()
    }
  }

  # put together folds and cvControl (inner loop one) into a list to loop over
  foldsList <- Map(list, folds = folds, cvControl = innerCvControl)

  # check input:
  if(is.null(obsWeights)) {
    obsWeights <- rep(1, N)
  }
  if(!identical(length(obsWeights), N)) {
    stop("obsWeights vector must have the same dimension as Y")
  }

  # check method:
  if(is.character(method)) {
    if(exists(method, mode = 'list')) {
      method <- get(method, mode = 'list')
    } else if(exists(method, mode = 'function')) {
      method <- get(method, mode = 'function')()
    }
  } else if(is.function(method)) {
    method <- method()
  }
  if(!is.list(method)) {
    stop("method is not in the appropriate format. Check out help('method.template')")
  }
  # make some modifications (scale) to the method.CC_LS
  method$computeCoef <- method.CC_LS.scale()$computeCoef

  # family can be either character or function, so these lines put everything together
  #family for stage 1
  if(is.character(family.1))
    family.1 <- get(family.1, mode="function", envir=parent.frame())
  if(is.function(family.1))
    family.1 <- family.1()
  if (is.null(family.1$family)) {
    print(family.1)
    stop("'family' not recognized")
  }
  # family for stage 2
  if(is.character(family.2))
    family.2 <- get(family.2, mode="function", envir=parent.frame())
  if(is.function(family.2))
    family.2 <- family.2()
  if (is.null(family.2$family)) {
    print(family.2)
    stop("'family' not recognized")
  }
  # family for single stage
  if(is.character(family.single))
    family.single <- get(family.single, mode="function", envir=parent.frame())
  if(is.function(family.single))
    family.single <- family.single()
  if (is.null(family.single$family)) {
    print(family.single)
    stop("'family' not recognized")
  }

  # create placeholders:
  # original library of algorithms
  library.stage1 <- library.2stage$stage1
  library.stage2 <- library.2stage$stage2
  library.stage_1 <- SuperLearner:::.createLibrary(library.stage1)
  library.stage_2 <- SuperLearner:::.createLibrary(library.stage2)
  library.stage_single <- SuperLearner:::.createLibrary(library.1stage)

  # number of algorithms and screening algorithms
  k.1 <- nrow(library.stage_1$library)
  k.2 <- nrow(library.stage_2$library)
  k.single <- nrow(library.stage_single$library)
  k.2stage <- k.1*k.2
  k.all <- k.1*k.2+k.single
  kScreen.1 <- length(library.stage_1$screenAlgorithm)
  kScreen.2 <- length(library.stage_2$screenAlgorithm)
  kScreen.single <- length(library.stage_single$screenAlgorithm)
  kscreen.2stage <- kScreen.1*kScreen.2
  kscreen.all <- kScreen.1*kScreen.2+kScreen.single

  # library name
  # stage 1
  libname.stage.1 <- NULL
  lib.stage1 <- library.stage_1$library$predAlgorithm
  lib.stage1.screen <- library.stage_1$screenAlgorithm[library.stage_1$library$rowScreen]
  repname <- function(x) {
    name <- rep(x,k.2)
    return(name)
  }
  libname.stage.1 <- unlist(lapply(lib.stage1,repname), use.names=FALSE)
  libname.stage.1.screen <- unlist(lapply(lib.stage1.screen,repname), use.names=FALSE)
  # stage 2
  libname.stage.2 <- NULL
  lib.stage2 <- library.stage_2$library$predAlgorithm
  lib.stage2.screen <- library.stage_2$screenAlgorithm[library.stage_2$library$rowScreen]
  libname.stage.2 <- rep(lib.stage2,k.1)
  libname.stage.2.screen <- rep(lib.stage2.screen,k.1)
  # single stage
  libname.stage.single <- library.stage_single$library$predAlgorithm
  libname.stage.single.screen <- library.stage_single$screenAlgorithm[library.stage_single$library$rowScreen]

  twostage.library <- paste("S1:",paste(libname.stage.1,libname.stage.1.screen,sep="_"),
                            "+ S2:",paste(libname.stage.2,libname.stage.2.screen,sep="_"))
  wholelibrary <- c(twostage.library,
                    paste("Single:",paste(libname.stage.single,libname.stage.single.screen,sep="_")))

  wholelibrary <- c(wholelibrary,"Standard SuperLearner")
  # library
  # add library for two stages
  lib <- list(twostage=data.frame("predAlgorithm"=paste("S1:",libname.stage.1,
                                                        "+ S2:",libname.stage.2),
                                  "rowScreen.Stage.1"=rep(library.stage_1$library$rowScreen,each=k.2),
                                  "rowScreen.Stage.2"=rep(library.stage_2$library$rowScreen,k.1)),
              singestage=library.stage_single$library)
  library <- list("library"=lib,
                  "screenAlgorithm"=list(stage.1 = library.stage_1$screenAlgorithm,
                                         stage.2 = library.stage_2$screenAlgorithm,
                                         stage.single = library.stage_single$screenAlgorithm))

  # output from each two-part superlearner
  AllSL <- vector('list', V)
  names(AllSL) <- paste("training", 1:V, sep=" ")

  # store the prediction from two-part superlearner
  SL.predict <- rep(NA, N)
  # store the prediction from discrete (single) superlearner
  discreteSL.predict <- rep.int(NA, N)
  # store the name of discrete superlearner
  whichDiscreteSL <- rep.int(NA, V)
  # store the prediction from each of the algorithms from two part superlearner
  library.predict <- matrix(NA, nrow = N, ncol = k.all+1)
  colnames(library.predict) <- wholelibrary

  if(p < 2 & !identical(library$screenAlgorithm, "All")) {
    warning('Screening algorithms specified in combination with single-column X.')
  }

  # run SuperLearner:
  .crossValFun <- function(valid, Y, dataX, family, id, obsWeights, SL.library, method,
                           verbose, control, saveAll) {
    cvLearn <- dataX[-valid[[1]], , drop = FALSE]
    cvOutcome <- Y[-valid[[1]]]
    cvValid <- dataX[valid[[1]], , drop = FALSE]
    cvId <- id[-valid[[1]]]
    cvObsWeights <- obsWeights[-valid[[1]]]

    # fit two-stage SuperLearner
    fit.twoSL <- twostageSL(Y = cvOutcome, X = cvLearn, newX = cvValid,
                            family.1 = family.1,family.2 = family.2,
                            family.single = family.single,
                            library.2stage = library.2stage,
                            library.1stage = library.1stage,
                            twostage = TRUE, method = method, id = cvId,
                            verbose = verbose, control = control, cvControl = valid[[2]],
                            obsWeights = cvObsWeights, env = env)

    # construct one-stage superlearner
    # extract onestage matrix z1
    z1 <- fit.twoSL$Z[,c((k.2stage+1):k.all)]
    onestagename <- colnames(fit.twoSL$library.predict[,c((k.2stage+1):k.all)])
    # get optimum weights for each algorithm in one-stage
    getCoef <- method.CC_LS.scale()$computeCoef(Z=z1,Y=cvOutcome,libraryNames=onestagename,
                                                verbose=FALSE)
    coef.onestage <- getCoef$coef
    # Prediction for each algorithm in one-stage superlearner
    predY.onestage = fit.twoSL$library.predict[,c((k.2stage+1):k.all)]
    # Compute onestage superlearner predictions on newX.
    onestage.pred <- fit.twoSL$method$computePred(predY = predY.onestage, coef = coef.onestage,
                                                     control =fit.twoSL$control)

    # combine the two SL (one-stage & two-stage)
    fit.oneSL <- data.frame("CV.Risk"=getCoef$cvRisk,"Coef"=coef.onestage)
    fit <- list("one stage" = fit.oneSL,
                "two stage" = fit.twoSL)

    lib.predict <- cbind(fit.twoSL$library.predict,onestage.pred)
    fit.coef <- c(fit.twoSL$coef,0)
    names(fit.coef) <- wholelibrary
    # output from two-stage SuperLearner
    out <- list(cvAllSL = if(saveAll) fit, cvSL.predict = fit.twoSL$SL.predict,
                cvdiscreteSL.predict = fit.twoSL$library.predict[, which.min(fit.twoSL$cvRisk)],
                cvwhichDiscreteSL = names(which.min(fit.twoSL$cvRisk)),
                cvlibrary.predict = lib.predict, cvcoef = fit.coef)
    return(out)
  }


  # different parallel option
  if (inherits(parallel, 'cluster')) {
    .SL.require('parallel')
    cvList <- parallel::parLapply(parallel, X = foldsList, fun = .crossValFun, Y = Y,
                                  dataX = X,
                                  method = method, id = id, obsWeights = obsWeights,
                                  verbose = verbose, control = control, saveAll = saveAll)
  } else if (parallel == 'multicore') {
    .SL.require('parallel')
    cvList <- parallel::mclapply(foldsList, FUN = .crossValFun, Y = Y, dataX = X,
                                method = method,
                                 id = id, obsWeights = obsWeights, verbose = verbose,
                                 control = control, saveAll = saveAll, mc.set.seed = FALSE)
  } else if (parallel == "seq") {
    cvList <- lapply(foldsList, FUN = .crossValFun, Y = Y, dataX = X, method = method, id = id,
                     obsWeights = obsWeights, verbose = verbose, control = control,
                     saveAll = saveAll)
  } else {
    stop('parallel option was not recognized, use parallel = "seq" for sequential computation.')
  }

  # check out Biobase::subListExtract to replace the lapply
  AllSL <- lapply(cvList, '[[', 'cvAllSL')

  SL.predict[unlist(folds, use.names = FALSE)] <- unlist(lapply(cvList, '[[', 'cvSL.predict'), use.names = FALSE)

  discreteSL.predict[unlist(folds, use.names = FALSE)] <- unlist(lapply(cvList, '[[', 'cvdiscreteSL.predict'), use.names = FALSE)

  whichDiscreteSL <- lapply(cvList, '[[', 'cvwhichDiscreteSL')

  library.predict[unlist(folds, use.names = FALSE), ] <- do.call('rbind', lapply(cvList, '[[', 'cvlibrary.predict'))

  coef <- do.call('rbind', lapply(cvList, '[[', 'cvcoef'))

  colnames(coef) <- wholelibrary

  # put together output
  out <- list(call = call, AllSL = AllSL, SL.predict = SL.predict,
              discreteSL.predict = discreteSL.predict, whichDiscreteSL = whichDiscreteSL,
              library.predict = library.predict, coef = coef, folds = folds, V = V,number0 = num.0,
              libraryNames = wholelibrary, SL.library = library, method = method, Y = Y)
  class(out) <- 'CV.twostageSL'
  return(out)
}







