#' @export twostage_listWrappers
twostage_listWrappers <- function(what = "both") {
  everything <- sort(getNamespaceExports("SuperLearner"))
  everything2 <- sort(getNamespaceExports("slcost"))
  if(what == "both") {
    message("Prediction algorithm wrappers in twostageSL from SuperLearner:\n")
    print(everything[grepl(pattern="^[S]L", everything)])
    message("Additional prediction algorithm wrappers in twostageSL from slcost:\n")
    print(everything2[grepl(pattern="^[S]L", everything2)])
    message("\nAll screening algorithm wrappers in SuperLearner:\n")
    print("All")
    print(everything[grepl(pattern="screen", everything)])
  } else if(what == "SL") {
    message("Prediction algorithm wrappers in twostageSL from SuperLearner:\n")
    print(everything[grepl(pattern="^[S]L", everything)])
    message("Additional prediction algorithm wrappers in twostageSL from slcost:\n")
    print(everything2[grepl(pattern="^[S]L", everything2)])
  } else if(what == "screen") {
    message("All screening algorithm wrappers in twostageSL from SuperLearner:\n")
    print("All")
    print(everything[grepl(pattern="screen", everything)])
  } else if(what == 'method') {
    message("Only method in twostageSL package:\n")
    print("method.CC_LS.scale")
  } else {
    stop("Please specify what = 'both', 'SL', or 'screen'")
  }
}
