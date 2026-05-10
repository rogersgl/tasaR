#' @include class-AmpliconSequencing.R
NULL

# --------
# Settings
# --------

#' Retrieve object settings
#'
#' @description
#' Extracts a list of tasaR settings from an S4 object of class tas.object.settings or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.object.settings or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @export
#'
#' @examples
#' getSettings(object)
setGeneric("getSettings", function(object) standardGeneric("getSettings"))

#' @describeIn getSettings Method for class tas.object.settings
#' @export
setMethod("getSettings", signature(object = "tas.object.settings"), function(object) {
  return(as.list(object))
})

#' @describeIn getSettings Method for class AmpliconSequencing
#' @export
setMethod("getSettings", signature(object = "AmpliconSequencing"), function(object) {
  return(as.list(object@Settings))
})
