#' Two-stage Super Learner Prediction
#'
#' A Prediction Function for the Two-stage Super Learner. The \code{twostageSL} function takes a training set pair (X,Y) and returns the predicted values based on a validation set.
#'
#' @param Y The outcome in the training data set. Must be a numeric vector.
#' @param X The predictor variables in the training data set, usually a data.frame.
#' @param newX The predictor variables in the validation data set. The structure should match X. If missing, uses X for newX.
#' @param library.2stage Candidate prediction algorithms in two-stage super learner. A list containing prediction algorithms at stage 1 and stage 2, the prediction algorithms are either a character vector or a list containing character vectors. See details below for examples on the structure. A list of functions included in the \code{twostageSL} package can be found with \code{twostage_listWrappers}.
#' @param library.1stage Candidate prediction algorithms in standard super learner. Either a character vector of prediction algorithms or a list containing character vectors. See details below for examples on the structure. A list of functions included in the \code{twostageSL} package can be found with \code{twostage_listWrappers}.
#' @param twostage logical; TRUE for implementing two-stage super learner; FALSE for implementing standatd super learner
#' @param family.1 Error distribution of the stage 1 outcome for two-stage super learner. Currently only allows \code{binomial} (default) to describe the error distribution. Link function information will be ignored and should be contained in the method argument below.
#' @param family.2 Error distribution of the stage 2 outcome for two-stage super learner. Currently only allows \code{gaussian} (default) to describe the error distribution. Link function information will be ignored and should be contained in the method argument below.
#' @param family.single Error distribution of the outcome for standard super learner. Currently only allows \code{gaussian} (default) to describe the error distribution. Link function information will be ignored and should be contained in the method argument below.
#' @param method Details on estimating the coefficients for the two-stage super learner and the model to combine the individual algorithms in the library. Currently, the built in option is only "method.CC_LS.scale" (default) which is a scaled version of CC_LS. CC_LS.scale uses Goldfarb and Idnani's quadratic programming algorithm to calculate the best convex combination of weights to minimize the squared error loss. In addition, CC_LS.scale divides the quadratic function by a large constant to shrink the huge matrix and vector in quadratic function.
#' @param id Optional cluster identification variable. For the cross-validation splits, \code{id} forces observations in the same cluster to be in the same validation fold. \code{id} is passed to the prediction and screening algorithms in library.2stage and library.1stage, but be sure to check the individual wrappers as many of them ignore the information.
#' @param verbose logical; TRUE for printing progress during the computation (helpful for debugging).
#' @param control A list of parameters to control the estimation process. Parameters include \code{saveFitLibrary} and \code{trimLogit}. See \code{\link{twostageSL.control}} for details.
#' @param cvControl A list of parameters to control the cross-validation process. Parameters include \code{V}, \code{stratifyCV}, \code{shuffle} and \code{validRows}. See \code{\link{twostageSL.CV.control}} for details.
#' @param obsWeights Optional observation weights variable. As with \code{id} above, \code{obsWeights} is passed to the prediction and screening algorithms, but many of the built in wrappers ignore (or can't use) the information. If you are using observation weights, make sure the library you specify uses the information.
#' @param env Environment containing the learner functions. Defaults to the calling environment.
#' @details \code{twostageSL} fits the two-stage super learner prediction algorithm. The weights for each algorithm in \code{library.2stage} and \code{library.1stage} is estimated, along with the fit of each algorithm.
#'
#' The prescreen algorithms. These algorithms first rank the variables in \code{X} based on either a univariate regression p-value or the \code{randomForest} variable importance. A subset of the variables in \code{X} is selected based on a pre-defined cut-off. With this subset of the X variables, the algorithms in \code{library.2stage} and \code{library.1stage} are then fit.
#'
#' The twostageSL package contains a few prediction and screening algorithm wrappers. The full list of wrappers can be viewed with \code{listWrappers()}. The design of the twostageSL package is such that the user can easily add their own wrappers.
#' @return An object with S3 class \code{twostageSL} containing:
#' \item{call}{The matched call.}
#' \item{libraryNames}{A character vector with the names of the algorithms in the library. The format is 'predictionAlgorithm_screeningAlgorithm' with '_All' used to denote the prediction algorithm run on all variables in X.}
#' \item{library.Num}{Number of prediction algorithms in \code{library.2stage} and \code{library.1stage}.}
#' \item{orig.library}{Returns the prediction algorithms and screening algorithms in each stage of \code{library.2stage} and \code{library.1stage} seperately.}
#' \item{SL.library}{Returns the prediction algorithms and screening algorithms in \code{library.2stage} and \code{library.1stage}.}
#' \item{SL.predict}{The predicted values from the two-stage super learner for the rows in \code{newX}.}
#' \item{coef}{Coefficients for the two-stage super learner.}
#' \item{library.predict}{A matrix with the predicted values from each algorithm in \code{library.2stage} and \code{library.1stage} for the rows in \code{newX}.}
#' \item{Z}{The Z matrix (the cross-validated predicted values for each algorithm in \code{library.2stage} and \code{library.1stage}).}
#' \item{cvRisk}{A numeric vector with the V-fold cross-validated risk estimate for each algorithm in \code{library.2stage} and \code{library.1stage}. Note that this does not contain the CV risk estimate for the two-stage super learner, only the individual algorithms in the library.}
#' \item{family}{Returns the \code{family.1}, \code{family.2} and \code{family.single} value from above}
#' \item{fitLibrary}{A list with the fitted objects for each algorithm in \code{library.2stage} and \code{library.1stage} on the full training data set.}
#' \item{cvfitLibrary}{A list with fitted objects for each algorithm in \code{library.2stage} and \code{library.1stage} on each of \code{v} different training data sets.}
#' \item{varNames}{A character vector with the names of the variables in \code{X}.}
#' \item{validRows}{A list containing the row numbers for the V-fold cross-validation step.}
#' \item{number0}{A dataframe indicating the number of zeros in each of the \code{v} fold.}
#' \item{method}{A list with the method functions.}
#' \item{whichScreen}{A logical matrix indicating which variables passed each screening algorithm.}
#' \item{control}{The \code{control} list.}
#' \item{cvControl}{The \code{cvControl} list.}
#' \item{errorsInCVLibrary}{A logical vector indicating if any algorithms experienced an error within the CV step.}
#' \item{errorsInLibrary}{A logical vector indicating if any algorithms experienced an error on the full data.}
#' \item{data}{The data frame including the predict variables and outcome in the training data set.}
#' \item{env}{Environment passed into function which will be searched to find the learner functions. Defaults to the calling environment.}
#' \item{times}{A list that contains the execution time of the twostageSL, plus separate times for model fitting and prediction.}
#' @author Ziyue Wu
#' @references van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2008) Super Learner, Statistical Applications of Genetics and Molecular Biology, 6, article 25.
#' @seealso \code{\link{SuperLearner}}.
#' @export
#' @import SuperLearner
#' @examples
#' ## simulate data
#' set.seed(123)
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
#' ## test set
#' m <- 1000
#' newX <- matrix(rnorm(m*p), nrow = m, ncol = p)
#' colnames(newX) <- paste("X", 1:p, sep="")
#' newX <- data.frame(newX)
#' newY <- rep(NA,m)
#' ## probability of outcome being zero
#' newprob <- plogis(1 + newX[,1] + newX[,2] + newX[,1]*newX[,2])
#' newg <- rbinom(m,1,newprob)
#' ## assign zero outcome
#' newind <- newg==0
#' newY[newind] <- 0
#' ## assign non-zero outcome
#' newind <- g==1
#' newY[newind] <- 10 + newX[newind, 1] + sqrt(abs(newX[newind, 2] * newX[newind, 3])) + newX[newind, 2] - X[newind, 3] + rnorm(sum(newind))
#'
#' ## generate the Library
#' twostage.library <- list(stage1=c("SL.glm","SL.mean","SL.earth"),
#'                         stage2=c("SL.glm","SL.mean","SL.earth"))
#' onestage.library <- c("SL.glm","SL.mean","SL.earth")
#'
#' ## run the twostage super learner
#' two <- twostageSL(Y=Y,
#'                  X=X,
#'                  newX = newX,
#'                  library.2stage <- twostage.library,
#'                  library.1stage <- onestage.library,
#'                  twostage = TRUE,
#'                  family.1=binomial,
#'                  family.2=gaussian,
#'                  family.single=gaussian,
#'                  cvControl = list(V = 5))
#' two
#' ## run the standard super learner
#' one <- twostageSL(Y=Y,
#'                  X=X,
#'                  newX = newX,
#'                  library.2stage <- twostage.library,
#'                  library.1stage <- onestage.library,
#'                  twostage = FALSE,
#'                  family.1=binomial,
#'                  family.2=gaussian,
#'                  family.single=gaussian,
#'                  cvControl = list(V = 5))
#' one
#'
#' ## library with screening
#' twostage.library <- list(stage1=list(c("SL.glm","screen.glmnet"),
#'                                     c("SL.earth","screen.corP"),
#'                                     c("SL.mean","All")),
#'                         stage2=list(c("SL.glm","screen.glmnet"),
#'                                     c("SL.earth","screen.corP"),
#'                                     c("SL.mean","All")))
#' onestage.library <- list(c("SL.glm","screen.glmnet"),
#'                         c("SL.earth","screen.corP"),
#'                         c("SL.mean","All"))
#'
#' ## run the twostage super learner
#' two <- twostageSL(Y=Y,
#'                  X=X,
#'                  newX = newX,
#'                  library.2stage <- twostage.library,
#'                  library.1stage <- onestage.library,
#'                  twostage = TRUE,
#'                  family.1=binomial,
#'                  family.2=gaussian,
#'                  family.single=gaussian,
#'                  cvControl = list(V = 5))
#' two
#' ## run the standard super learner
#' one <- twostageSL(Y=Y,
#'                  X=X,
#'                  newX = newX,
#'                  library.2stage <- twostage.library,
#'                  library.1stage <- onestage.library,
#'                  twostage = FALSE,
#'                  family.1=binomial,
#'                  family.2=gaussian,
#'                  family.single=gaussian,
#'                  cvControl = list(V = 5))
#' one
#'
#'

twostageSL <- function(Y, X, newX = NULL, library.2stage, library.1stage,twostage,
                       family.1, family.2, family.single, method="method.CC_LS",
                       id=NULL, verbose=FALSE, control = list(), cvControl = list(),
                       obsWeights = NULL, env = parent.frame()){

  # Begin timing how long two-stage SuperLearner takes to execute
  time_start = proc.time()

  # Get details of estimation algorithm for the algorithm weights (coefficients)
  if (is.character(method)) {
    if (exists(method, mode = 'list')) {
      method <- get(method, mode = 'list')
    } else if (exists(method, mode = 'function')) {
      method <- get(method, mode = 'function')()
    }
  } else if (is.function(method)) {
    method <- method()
  }
  # make some modifications (scale) to the superlearner:method.CC_LS
  method$computeCoef <- method.CC_LS.scale()$computeCoef
  if(!is.list(method)) {
    stop("method is not in the appropriate format. Check out help('method.template')")
  }
  if(!is.null(method$require)) {
    sapply(method$require, function(x) require(force(x), character.only = TRUE))
  }

  # get defaults for controls and make sure in correct format
  control <- do.call('SuperLearner.control', control)
  # change the logical for saveCVFitLibrary to TRUE (we are gonna use that)
  control$saveCVFitLibrary <- TRUE

  cvControl <- do.call('SuperLearner.CV.control', cvControl)

  # put together the library
  library.stage1 <- library.2stage$stage1
  library.stage2 <- library.2stage$stage2
  library.stage_1 <- SuperLearner:::.createLibrary(library.stage1)
  library.stage_2 <- SuperLearner:::.createLibrary(library.stage2)
  library.stage_single <- SuperLearner:::.createLibrary(library.1stage)
  SuperLearner:::.check.SL.library(library = c(unique(library.stage_1$library$predAlgorithm),
                                               library.stage_1$screenAlgorithm))
  SuperLearner:::.check.SL.library(library = c(unique(library.stage_2$library$predAlgorithm),
                                               library.stage_2$screenAlgorithm))
  SuperLearner:::.check.SL.library(library = c(unique(library.stage_single$library$predAlgorithm),
                                               library.stage_single$screenAlgorithm))
  call <- match.call(expand.dots = TRUE)

  # should we be checking X and newX for data.frame?
  # data.frame not required, but most of the built-in wrappers assume a data.frame
  if(!inherits(X, 'data.frame')) message('X is not a data frame. Check the algorithms in SL.library to make sure they are compatible with non data.frame inputs')
  varNames <- colnames(X)
  N <- dim(X)[1L]
  p <- dim(X)[2L]
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

  # check if the model use method.AUC
  if (family.1$family != "binomial" & isTRUE("cvAUC" %in% method$require)){
    stop("'method.AUC' is designed for the 'binomial' family only")
  }
  if (family.2$family != "binomial" & isTRUE("cvAUC" %in% method$require)){
    stop("'method.AUC' is designed for the 'binomial' family only")
  }
  if (family.single$family != "binomial" & isTRUE("cvAUC" %in% method$require)){
    stop("'method.AUC' is designed for the 'binomial' family only")
  }

  # chekc whether screen algorithm compatible with number of columns
  if(p < 2 & !identical(library.stage_1$screenAlgorithm, "All")) {
    warning('Screening algorithms specified in combination with single-column X.')
  }
  if(p < 2 & !identical(library.stage_2$screenAlgorithm, "All")) {
    warning('Screening algorithms specified in combination with single-column X.')
  }
  if(p < 2 & !identical(library.stage_single$screenAlgorithm, "All")) {
    warning('Screening algorithms specified in combination with single-column X.')
  }

  # generate library names
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

  # add family for two stages and single stage
  family <- list(stage1 = family.1,stage2 = family.2,
                 stage.single = family.single)

  # add library for two stages
  lib <- list(twostage=data.frame("predAlgorithm"=paste("S1:",libname.stage.1,
                                                        "+ S2:",libname.stage.2),
                                  "rowScreen.Stage.1"=rep(library.stage_1$library$rowScreen,each=k.2),
                                  "rowScreen.Stage.2"=rep(library.stage_2$library$rowScreen,k.1)),
              singlestage=library.stage_single$library)
  library <- list("library"=lib,
                  "screenAlgorithm"=list(stage.1 = library.stage_1$screenAlgorithm,
                                         stage.2 = library.stage_2$screenAlgorithm,
                                         stage.single = library.stage_single$screenAlgorithm))

  # if newX is missing, use X
  if(is.null(newX)) {
    newX <- X
  }

  # Various chekcs for data structure
  if(!identical(colnames(X), colnames(newX))) {
    stop("The variable names and order in newX must be identical to the variable names and order in X")
  }
  if (sum(is.na(X)) > 0 | sum(is.na(newX)) > 0 | sum(is.na(Y)) > 0) {
    stop("missing data is currently not supported. Check Y, X, and newX for missing values")
  }
  if (!is.numeric(Y)) {
    stop("the outcome Y must be a numeric vector")
  }

  # errors records if an algorithm stops either in the CV step and/or in full data
  errorsInCVLibrary <- rep(0, k.all)
  errorsInLibrary <- rep(0, k.all)

  ########################################################################################################
  # Step 0: make valid rows
  # ensure each folds have approximately equal number of obs with y=0
  V <- cvControl$V
  ord <- order(Y)
  cvfold <- rep(c(1:V,V:1),N)[1:N]
  folds <- split(ord, factor(cvfold))
  folds <- lapply(folds,sort,decreasing=FALSE)
  # check
  tab <- rep(NA,V)
  for (i in 1:V) {
    tab[i] <- sum(Y[folds[[i]]]==0)
  }
  num.0 <- data.frame("fold"=paste0("fold ",c(1:cvControl$V)),"number.of.0"=tab)

  cvControl$validRows = folds

  # test id
  if(is.null(id)) {
    id <- seq(N)
  }
  if(!identical(length(id), N)) {
    stop("id vector must have the same dimension as Y")
  }
  # test observation weights
  if(is.null(obsWeights)) {
    obsWeights <- rep(1, N)
  }
  if(!identical(length(obsWeights), N)) {
    stop("obsWeights vector must have the same dimension as Y")
  }

  #########################################################################################################
  # Step 1: fit superlearner for modeling prob of y=0
  time_train_start = proc.time()

  # list all the algorithms considered
  # save cross-validated fits (10) in the control option
  step1.fit <- SuperLearner(Y=as.numeric(Y==0),X=X,family=family.1,
                            SL.library=library.stage1,verbose=verbose,
                            method=method.CC_nloglik,
                            control=list(saveCVFitLibrary=T),
                            cvControl=cvControl)
  # get the cross-validated predicted values for each algorithm in SL.library
  # P(Y=0|X)
  z1 <- step1.fit$Z
  # get the cross-validated fits (10 fits for 10 training set) for each algorithm
  stage1.cvFitLibrary <- step1.fit$cvFitLibrary


  ##########################################################################################################
  # step 2: fit model for E[Y|Y>0,X]
  # create function for the cross-validation step at stage 2:
  .crossValFUN <- function(valid, Y, dataX, predX, id, obsWeights, library, family,
                           kScreen, k, p, libraryNames, saveCVFitLibrary) {
    tempLearn <- dataX[-valid, , drop = FALSE]
    tempOutcome <- Y[-valid]
    tempValid <- predX[valid, , drop = FALSE]
    tempWhichScreen <- matrix(NA, nrow = kScreen, ncol = p)
    tempId <- id[-valid]
    tempObsWeights <- obsWeights[-valid]

    # should this be converted to a lapply also?
    for(s in seq(kScreen)) {
      screen_fn = get(library$screenAlgorithm[s], envir = env)
      testScreen <- try(do.call(screen_fn,
                                list(Y = tempOutcome,
                                     X = tempLearn,
                                     family = family,
                                     id = tempId,
                                     obsWeights = tempObsWeights)))
      if(inherits(testScreen, "try-error")) {
        warning(paste("replacing failed screening algorithm,", library$screenAlgorithm[s], ", with All()", "\n "))
        tempWhichScreen[s, ] <- TRUE
      } else {
        tempWhichScreen[s, ] <- testScreen
      }
      if(verbose) {
        message(paste("Number of covariates in ", library$screenAlgorithm[s], " is: ", sum(tempWhichScreen[s, ]), sep = ""))
      }
    } #end screen

    # should this be converted to a lapply also?
    out <- matrix(NA, nrow = nrow(tempValid), ncol = k)
    if(saveCVFitLibrary){
      model_out <- vector(mode = "list", length = k)
    }else{
      model_out <- NULL
    }

    for(s in seq(k)) {
      pred_fn = get(library$library$predAlgorithm[s], envir = env)
      testAlg <- try(do.call(pred_fn,
                             list(Y = tempOutcome,
                                  X = subset(tempLearn,
                                             select = tempWhichScreen[library$library$rowScreen[s], ],
                                             drop=FALSE),
                                  newX = subset(tempValid,
                                                select = tempWhichScreen[library$library$rowScreen[s], ],
                                                drop=FALSE),
                                  family = family,
                                  id = tempId,
                                  obsWeights = tempObsWeights)))
      if(inherits(testAlg, "try-error")) {
        warning(paste("Error in algorithm", library$library$predAlgorithm[s], "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n" ))
        # errorsInCVLibrary[s] <<- 1
      } else {
        out[, s] <- testAlg$pred
        if(saveCVFitLibrary){
          model_out[[s]] <- testAlg$fit
        }
      }
      if (verbose) message(paste("CV", libraryNames[s]))
    } #end library
    if(saveCVFitLibrary){
      names(model_out) <- libraryNames
    }
    invisible(list(out = out, model_out = model_out))
  }

  # the lapply performs the cross-validation steps to create Z for stage 2
  # additional steps to put things in the correct order
  # rbind unlists the output from lapply
  # need to unlist folds to put the rows back in the correct order

  # create subset with only obs y>0
  pid <- c(1:N)
  dat.p <- cbind(pid,X,Y)
  Y.p <- Y[Y>0]
  X.p <- dat.p[dat.p$Y>0,-c(1,ncol(dat.p))]
  X.test <- dat.p[,-c(1,ncol(dat.p))]
  p.id <- dat.p[dat.p$Y>0,1]
  p.obsWeights <- obsWeights[p.id]

  # find intersect of folds & X.p
  folds.p <- vector(mode = "list", length = V)

  for (i in 1:V){
    folds.p[[i]] <- intersect(folds[[i]],as.numeric(row.names(X.p)))
  }
  names(folds.p) <- seq(1,V)

  crossValFUN_out <- lapply(folds.p, FUN = .crossValFUN,
                            Y = Y.p, dataX = X.p, predX = X.test, id = p.id,
                            obsWeights = p.obsWeights, family = family.2,
                            library = library.stage_2, kScreen = kScreen.2,
                            k = k.2, p = p, libraryNames = library.stage_2$library$predAlgorithm,
                            saveCVFitLibrary = control$saveCVFitLibrary)

  # create matrix to store results
  z2 <- matrix(NA,nrow = N,ncol=k.2)
  z2[unlist(folds.p, use.names = FALSE), ] <- do.call('rbind', lapply(crossValFUN_out, "[[", "out"))
  z2[-unlist(folds.p, use.names = FALSE), ] <- 0

  if(control$saveCVFitLibrary){
    stage2.cvFitLibrary <- lapply(crossValFUN_out, "[[", "model_out")
  }else{
    stage2.cvFitLibrary <- NULL
  }

  # z1 for E[P(Y=0|X)]
  # z2 for E[Y|Y>0,X]
  # multiply (1-z1)*z2 to generate z
  z <- NULL
  for (i in 1:k.1){
    for (j in 1:k.2){
      temp <- rep(0,N)
      for (k in 1:N){
        temp[k] <- (1-z1[k,i])*z2[k,j]
      }
      z <- cbind(z,temp)
    }
  }

  ########################################################################################################
  # step 3: fit the whole model using one stage option (rather than two stages)
  # list all the algorithms considered
  # save cross-validated fits (10) in the control option
  onestage.fit <- SuperLearner(Y=Y,X=X,family=family.single,
                               SL.library=library.1stage,verbose=verbose,
                               method=method.CC_LS.scale,
                               control=list(saveCVFitLibrary=T),
                               cvControl=cvControl)
  # get the cross-validated predicted values for each algorithm in SL.library
  z.single <- onestage.fit$Z
  # get the cross-validated fits (10 fits for 10 training set) for each algorithm
  single.stage.cvFitLibrary <- onestage.fit$cvFitLibrary

  # combine 2 stages output z with 1 stage prediction output z
  z <- cbind(z,z.single)

  # Check for errors. If any algorithms had errors, replace entire column with
  # 0 even if error is only in one fold.
  errorsInCVLibrary <- apply(z, 2, function(x) anyNA(x))
  if (sum(errorsInCVLibrary) > 0) {
    z[, as.logical(errorsInCVLibrary)] <- 0
  }
  if (all(z == 0)) {
    stop("All algorithms dropped from library")
  }

  ########################################################################################################
  # step 4: use cross-validation to calcualte weights for different algorithm at stage 1 & 2
  # using an scaled method.CC_LS in superlearner
  # get optimum weights for each algorithm
  getCoef <- method.CC_LS.scale()$computeCoef(Z=z,Y=Y,libraryNames=wholelibrary,
                                              verbose=verbose)
  coef <- getCoef$coef
  names(coef) <- wholelibrary

  time_train = proc.time() - time_train_start

  # Set a default in case the method does not return the optimizer result.
  if (!("optimizer" %in% names(getCoef))) {
    getCoef["optimizer"] <- NA
  }

  #########################################################################################################
  # step 5: now fit all algorithms in library on entire data set (X) and predict on newX

  .screenFun <- function(fun, list) {
    screen_fn = get(fun, envir = env)
    testScreen <- try(do.call(screen_fn, list))
    if (inherits(testScreen, "try-error")) {
      warning(paste("replacing failed screening algorithm,", fun, ", with All() in full data", "\n "))
      out <- rep(TRUE, ncol(list$X))
    } else {
      out <- testScreen
    }
    return(out)
  }

  time_predict_start = proc.time()

  # stage 1
  whichScreen.stage1 <- sapply(library$screenAlgorithm$stage.1, FUN = .screenFun,
                               list = list(Y = Y, X = X, family = family, id = id, obsWeights = NULL),
                               simplify = FALSE)
  whichScreen.stage1 <- do.call(rbind, whichScreen.stage1)
  # stage 2
  whichScreen.stage2 <- sapply(library$screenAlgorithm$stage.2, FUN = .screenFun,
                               list = list(Y = Y, X = X, family = family, id = id, obsWeights = NULL),
                               simplify = FALSE)
  whichScreen.stage2 <- do.call(rbind, whichScreen.stage2)
  # single stage
  whichScreen.stage.single <- sapply(library$screenAlgorithm$stage.single, FUN = .screenFun,
                                     list = list(Y = Y, X = X, family = family, id = id, obsWeights = NULL),
                                     simplify = FALSE)
  whichScreen.stage.single <- do.call(rbind, whichScreen.stage.single)
  # combine together
  whichScreen <- list(stage1 = whichScreen.stage1,
                      stage2 = whichScreen.stage2,
                      single.stage = whichScreen.stage.single)

  # Prediction for each algorithm
  .predFun <- function(index, lib, Y, dataX, newX, whichScreen, family, id, obsWeights,
                       verbose, control, libraryNames) {
    pred_fn = get(lib$predAlgorithm[index], envir = env)
    testAlg <- try(do.call(pred_fn, list(Y = Y,
                                         X = subset(dataX,
                                                    select = whichScreen[lib$rowScreen[index], ], drop=FALSE),
                                         newX = subset(newX, select = whichScreen[lib$rowScreen[index], ], drop=FALSE),
                                         family = family, id = id, obsWeights = obsWeights)))
    # testAlg <- try(do.call(lib$predAlgorithm[index], list(Y = Y, X = dataX[, whichScreen[lib$rowScreen[index], drop = FALSE]], newX = newX[, whichScreen[lib$rowScreen[index], drop = FALSE]], family = family, id = id, obsWeights = obsWeights)))
    if (inherits(testAlg, "try-error")) {
      warning(paste("Error in algorithm", lib$predAlgorithm[index], " on full data", "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n" ))
      out <- rep.int(NA, times = nrow(newX))
    } else {
      out <- testAlg$pred
      if (control$saveFitLibrary) {
        eval(bquote(fitLibrary[[.(index)]] <- .(testAlg$fit)), envir = fitLibEnv)
      }
    }
    if (verbose) {
      message(paste("full", libraryNames[index]))
    }
    invisible(out)
  }

  # stage 1
  # put fitLibrary at stage 1 in it's own environment to locate later
  fitLibEnv <- new.env()
  assign('fitLibrary', vector('list', length = k.1), envir = fitLibEnv)
  assign('libraryNames', library.stage_1$library$predAlgorithm, envir = fitLibEnv)
  evalq(names(fitLibrary) <- library.stage_1$library$predAlgorithm, envir = fitLibEnv)
  # get prediction for stage 1
  predY.stage1 <- do.call('cbind', lapply(seq(k.1), FUN = .predFun,
                                          lib = library.stage_1$library, Y = as.numeric(Y==0), dataX = X,
                                          newX = newX, whichScreen = whichScreen$stage1,
                                          family = family.1, id = id,
                                          obsWeights = obsWeights, verbose = verbose,
                                          control = control,
                                          libraryNames = library.stage_1$library$predAlgorithm))
  # save fit library for stage 1
  stage1.fitlib <- get("fitLibrary",envir = fitLibEnv)

  # stage 2
  # put fitLibrary at stage 2 in it's own environment to locate later
  fitLibEnv <- new.env()
  assign('fitLibrary', vector('list', length = k.2), envir = fitLibEnv)
  assign('libraryNames', library.stage_2$library$predAlgorithm, envir = fitLibEnv)
  evalq(names(fitLibrary) <- library.stage_2$library$predAlgorithm, envir = fitLibEnv)
  # get prediction for stage 2
  predY.stage2 <- do.call('cbind', lapply(seq(k.2), FUN = .predFun,
                                          lib = library.stage_2$library, Y = Y.p, dataX = X.p,
                                          newX = newX, whichScreen = whichScreen$stage2,
                                          family = family.2, id = p.id,
                                          obsWeights = p.obsWeights, verbose = verbose,
                                          control = control,
                                          libraryNames = library.stage_2$library$predAlgorithm))
  # save fit library for stage 2
  stage2.fitlib <- get("fitLibrary",envir = fitLibEnv)

  # single stage
  # put fitLibrary at single in it's own environment to locate later
  fitLibEnv <- new.env()
  assign('fitLibrary', vector('list', length = k.single), envir = fitLibEnv)
  assign('libraryNames', library.stage_single$library$predAlgorithm, envir = fitLibEnv)
  evalq(names(fitLibrary) <- library.stage_single$library$predAlgorithm, envir = fitLibEnv)
  # save fit library for single stage
  predY.stage.single <- do.call('cbind', lapply(seq(k.single), FUN = .predFun,
                                                lib = library.stage_single$library, Y = Y, dataX = X,
                                                newX = newX, whichScreen = whichScreen$single.stage,
                                                family = family.single, id = id,
                                                obsWeights = obsWeights, verbose = verbose,
                                                control = control,
                                                libraryNames = library.stage_single$library$predAlgorithm))
  # save fit library for single stage
  stage.single.fitlib <- get("fitLibrary",envir = fitLibEnv)

  # get prediction for 2-stage model
  predY <- NULL
  for (i in 1:k.1){
    for (j in 1:k.2){
      pred <- (1-predY.stage1[,i])*predY.stage2[,j]
      predY <- cbind(predY,pred)
    }
  }
  # combine with prediction from singe-stage model
  predY <- cbind(predY,predY.stage.single)

  #generate Fitlibrary
  fitLibrary = list("stage1"=stage1.fitlib,
                    "stage2"=stage2.fitlib,
                    "stage.sinlge"=stage.single.fitlib)

  #generate cross-validation Fitlibrary
  cvfitLibrary <- list("stage1"=stage1.cvFitLibrary,
                       "stage2"=stage2.cvFitLibrary,
                       "stage.single"=single.stage.cvFitLibrary)

  # check for errors
  errorsInLibrary <- apply(predY, 2, function(algorithm) anyNA(algorithm))
  if (sum(errorsInLibrary) > 0) {
    if (sum(coef[as.logical(errorsInLibrary)]) > 0) {
      warning(paste0("Re-running estimation of coefficients removing failed algorithm(s)\n",
                     "Original coefficients are: \n", paste(coef, collapse = ", "), "\n"))
      z[, as.logical(errorsInLibrary)] <- 0
      if (all(z == 0)) {
        stop("All algorithms dropped from library")
      }
      getCoef <- method$computeCoef(Z = z, Y = Y, libraryNames = wholelibrary,
                                    obsWeights = obsWeights, control = control,
                                    verbose = verbose,
                                    errorsInLibrary = errorsInLibrary)
      coef <- getCoef$coef
      names(coef) <- wholelibrary
    } else {
      warning("Coefficients already 0 for all failed algorithm(s)")
    }
  }

  # Compute super learner predictions on newX.
  getPred <- method$computePred(predY = predY, coef = coef, control=control)

  time_predict = proc.time() - time_predict_start

  # Add names of algorithms to the predictions.
  colnames(predY) <- wholelibrary

  # Clean up when errors in library.
  if(sum(errorsInCVLibrary) > 0) {
    getCoef$cvRisk[as.logical(errorsInCVLibrary)] <- NA
  }

  # Finish timing the full SuperLearner execution.
  time_end = proc.time()

  # Compile execution times.
  times = list(everything = time_end - time_start,
               train = time_train,
               predict = time_predict)

  # number of algorithms used in each stage
  library.num <- list(stage1 = k.1,
                      stage2 = k.2,
                      stage.single = k.single)

  #original library for each stage's algorithm
  orig.library <- list(stage1 = library.stage_1,
                       stage2 = library.stage_2,
                       stage.single = library.stage_single)

  # Output whether two-stage resutls or one-stage results
  if (twostage){
    # results
    data.frame("CV.Risk"=getCoef$cvRisk,"Coef"=coef)

    # Put everything together in a list.
    out <- list(
      call = call,
      libraryNames = wholelibrary,
      library.Num = library.num,
      orig.library = orig.library,
      SL.library = library,
      SL.predict = getPred,
      coef = coef,
      library.predict = predY,
      Z = z,
      cvRisk = getCoef$cvRisk,
      family = family,
      fitLibrary = fitLibrary,
      cvfitLibrary = cvfitLibrary,
      varNames = varNames,
      validRows = folds,
      number0 = num.0,
      method = method,
      whichScreen = whichScreen,
      control = control,
      cvControl = cvControl,
      errorsInCVLibrary = errorsInCVLibrary,
      errorsInLibrary = errorsInLibrary,
      metaOptimizer = getCoef$optimizer,
      env = env,
      times = times
    )
    class(out) <- c("SuperLearner")
    out
  } else {
    # results
    onestage.fit
    # Put everything together in a list.
  }
}
