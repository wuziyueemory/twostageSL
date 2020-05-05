#' Control parameters for the twostageSL
#'
#' Control parameters for the twostageSL
#'
#' @param saveFitLibrary Logical. Should the fit for each algorithm be saved in the output from \code{twostageSL}.
#' @param saveCVFitLibrary Logical. Should cross-validated fits for each algorithm be saved in the output from \code{twostageSL}.
#' @param trimLogit number between 0.0 and 0.5. What level to truncate the logit transformation to maintain a bounded loss function when using the NNloglik method.
#'
#' @return A list containing the control parameters.
#'
twostageSL.control <- function(saveFitLibrary = TRUE,
                                 saveCVFitLibrary = TRUE,
                                 trimLogit = 0.001) {
  if(trimLogit > 0.5) {
    warning('trimLogit must be less than 0.5, will replace with trimLogit = 0.001')
    trimLogit <- 0.001
  }
  list(saveFitLibrary = saveFitLibrary, trimLogit = trimLogit, saveCVFitLibrary = saveCVFitLibrary)
}

