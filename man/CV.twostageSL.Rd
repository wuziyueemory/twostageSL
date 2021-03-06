% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CV.twostageSL.R
\name{CV.twostageSL}
\alias{CV.twostageSL}
\title{Function to get V-Fold cross-validated risk estimate for two stage super learner}
\usage{
CV.twostageSL(
  Y,
  X,
  V = NULL,
  family.1 = binomial,
  family.2 = gaussian,
  family.single = gaussian,
  library.2stage,
  library.1stage,
  twostage,
  method = "method.CC_LS.scale",
  id = NULL,
  verbose = FALSE,
  control = list(saveFitLibrary = FALSE),
  cvControl = list(),
  innerCvControl = list(),
  obsWeights = NULL,
  saveAll = TRUE,
  parallel = "seq",
  env = parent.frame()
)
}
\arguments{
\item{Y}{The outcome.}

\item{X}{The covariates.}

\item{V}{The number of folds for \code{CV.SuperLearner}. This argument will be depreciated and moved into the \code{cvControl}. If Both \code{V} and \code{cvControl} set the number of cross-validation folds, an error message will appear. The recommendation is to use \code{cvControl}. This is not the number of folds for \code{twostageSL}. The number of folds for \code{twostageSL} is controlled with \code{innerCvControl}.}

\item{family.1}{Error distribution of the stage 1 outcome for two-stage super learner. Currently only allows \code{binomial} (default) to describe the error distribution. Link function information will be ignored and should be contained in the method argument below.}

\item{family.2}{Error distribution of the stage 2 outcome for two-stage super learner. Currently only allows \code{gaussian} (default) to describe the error distribution. Link function information will be ignored and should be contained in the method argument below.}

\item{family.single}{Error distribution of the outcome for standard super learner. Currently only allows \code{gaussian} (default) to describe the error distribution. Link function information will be ignored and should be contained in the method argument below.}

\item{library.2stage}{Candidate prediction algorithms in two-stage super learner. A list containing prediction algorithms at stage 1 and stage 2, the prediction algorithms are either a character vector or a list containing character vectors. See details below for examples on the structure. A list of functions included in the SuperLearner package can be found with \code{listWrappers()}.}

\item{library.1stage}{Candidate prediction algorithms in standard super learner. Either a character vector of prediction algorithms or a list containing character vectors. See details below for examples on the structure. A list of functions included in the SuperLearner package can be found with \code{listWrappers()}.}

\item{twostage}{logical; TRUE for implementing two-stage super learner; FALSE for implementing standatd super learner}

\item{method}{Details on estimating the coefficients for the two-stage super learner and the model to combine the individual algorithms in the library. Currently, the built in option is only "method.CC_LS.scale" (default) which is a scaled version of CC_LS. CC_LS.scale uses Goldfarb and Idnani's quadratic programming algorithm to calculate the best convex combination of weights to minimize the squared error loss. In addition, CC_LS.scale divides the quadratic function by a large constant to shrink the huge matrix and vector in quadratic function.}

\item{id}{Optional cluster identification variable. For the cross-validation splits, \code{id} forces observations in the same cluster to be in the same validation fold. \code{id} is passed to the prediction and screening algorithms in library.2stage and library.1stage, but be sure to check the individual wrappers as many of them ignore the information.}

\item{verbose}{logical; TRUE for printing progress during the computation (helpful for debugging).}

\item{control}{A list of parameters to control the estimation process. Parameters include \code{saveFitLibrary} and \code{trimLogit}. See \code{\link{twostageSL.control}} for details.}

\item{cvControl}{A list of parameters to control the cross-validation process. The outer cross-validation is the sample spliting for evaluating the \code{twostageSL}. Parameters include \code{V}, \code{stratifyCV}, \code{shuffle} and \code{validRows}. See \code{\link{twostageSL.CV.control}} for details.}

\item{innerCvControl}{A list of lists of parameters to control the inner cross-validation process. It should have \code{V} elements in the list, each a valid \code{cvControl} list. If only a single value, then replicated across all folds. The inner cross-validation are the values passed to each of the \code{V} \code{SuperLearner} calls. Parameters include \code{V}, \code{stratifyCV}, \code{shuffle} and \code{validRows}. See \code{\link{twostageSL.CV.control}} for details.}

\item{obsWeights}{Optional observation weights variable. As with \code{id} above, \code{obsWeights} is passed to the prediction and screening algorithms, but many of the built in wrappers ignore (or can't use) the information. If you are using observation weights, make sure the library you specify uses the information.}

\item{saveAll}{Logical; Should the entire \code{twostageSL} object be saved for each fold?}

\item{parallel}{Options for parallel computation of the V-fold step. Use "seq" (the default) for sequential computation. \code{parallel = 'multicore'} to use \code{mclapply} for the V-fold step (but note that \code{twostageSL()} will still be sequential). The default for \code{mclapply} is to check the \code{mc.cores} option, and if not set to default to 2 cores. Be sure to set \code{options()$mc.cores} to the desired number of cores if you don't want the default. Or \code{parallel} can be the name of a snow cluster and will use \code{parLapply} for the V-fold step. For both multicore and snow, the inner \code{twostageSL} calls will be sequential.}

\item{env}{Environment containing the learner functions. Defaults to the calling environment.}
}
\value{
An object of class \code{CV.twostageSL} (a list) with components:
\item{call}{The matched call.}
\item{AllSL}{If \code{saveAll = TRUE}, a list with output from each call to \code{twostageSL}, otherwise NULL.}
\item{SL.predict}{The predicted values from the two stage super learner when each particular row was part of the validation fold.}
\item{discreteSL.predict}{The traditional cross-validated selector. Picks the algorithm with the smallest cross-validated risk (in super learner terms, gives that algorithm coefficient 1 and all others 0).}
\item{whichDiscreteSL}{A list of length \code{V}. The elements in the list are the algorithm that had the smallest cross-validated risk estimate for that fold.}
\item{library.predict}{A matrix with the predicted values from each algorithm in \code{library.2stage} and \code{library.1stage}. The columns are the algorithms in \code{library.2stage} and \code{library.1stage} and the rows represent the predicted values when that particular row was in the validation fold (i.e. not used to fit that estimator).}
\item{coef}{A matrix with the coefficients for the two stage super learner on each fold. The columns are the algorithms in \code{library.2stage} and \code{library.1stage} the rows are the folds.}
\item{folds}{A list containing the row numbers for each validation fold.}
\item{V}{Number of folds for \code{CV.twostageSL}.}
\item{number0}{A dataframe indicating the number of zeros in each of the \code{v} fold.}
\item{libraryNames}{A character vector with the names of the algorithms in the library. The format is 'predictionAlgorithm_screeningAlgorithm' with '_All' used to denote the prediction algorithm run on all variables in X.}
\item{SL.library}{Returns the prediction algorithms and screening algorithms in \code{library.2stage} and \code{library.1stage}.}
\item{method}{A list with the method functions.}
\item{Y}{The outcome.}
}
\description{
Function to get V-fold cross-validated risk estimate for two stage super learner. This function simply splits the data into V folds and then calls twostageSL. Most of the arguments are passed directly to twostageSL.
}
\details{
The \code{twostageSL} function builds a estimator, but does not contain an estimate on the performance of the estimator. Various methods exist for estimator performance evaluation. If you are familiar with the super learner algorithm, it should be no surprise we recommend using cross-validation to evaluate the honest performance of the two stage super learner estimator. The function \code{CV.twostageSL} computes the usual V-fold cross-validated risk estimate for the two stage super learner (and all algorithms in \code{library.2stage} and \code{library.1stage} for comparison).
}
\examples{
## simulate data
set.seed(12321)

## training set
n <- 10000
p <- 5
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X) <- paste("X", 1:p, sep="")
X <- data.frame(X)
Y <- rep(NA,n)
## probability of outcome being zero
prob <- plogis(1 + X[,1] + X[,2] + X[,1]*X[,2])
g <- rbinom(n,1,prob)
## assign zero outcome
ind <- g==0
Y[ind] <- 0
## assign non-zero outcome
ind <- g==1
Y[ind] <- 10 + X[ind, 1] + sqrt(abs(X[ind, 2] * X[ind, 3])) + X[ind, 2] - X[ind, 3] + rnorm(sum(ind))

## run the CV.twostageSL
cv_sl <- CV.twostageSL(
 Y = Y, X = X,
 family.1 = binomial,
 family.2 = gaussian,
 family.single = gaussian,
 library.2stage = list(stage1=c("SL.glm","SL.mean","SL.earth"),
                       stage2=c("SL.glm","SL.mean","SL.earth")),
 library.1stage = c("SL.glm","SL.mean","SL.earth"),
 cvControl = list(V = 2),
 innerCvControl = list(list(V = 5),
                       list(V = 5))
)
cv_sl$AllSL
cv_sl$whichDiscreteSL
plot(cv_sl)

}
\seealso{
\code{\link{twostageSL}}.
}
\author{
Ziyue Wu
}
