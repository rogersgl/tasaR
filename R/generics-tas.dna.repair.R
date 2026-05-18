#' @include class-AmpliconSequencing.R
NULL

# ----------------
# DNA Repair Types
# ----------------

#' Retrieve DNA repair pathway inferences
#'
#' @description
#' Extracts a table (data.frame) of DNA repair pathway frequencies from an S4 object of class tas.dna.repair or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.dna.repair or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @export
#'
#' @examples
#' \dontrun{
#'   getMutationTypes(object)
#' }
setGeneric("getMutationTypes", function(object) standardGeneric("getMutationTypes"))

#' @describeIn getMutationTypes Method for class tas.dna.repair
#' @export
setMethod("getMutationTypes", signature(object = "tas.dna.repair"), function(object) {
  if (!isEmpty(object)) {return(as.data.frame(object))} else {return("Empty")}
})

#' @describeIn getMutationTypes Method for class AmpliconSequencing
#' @export
setMethod("getMutationTypes", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@MutationTypes)) {return(as.data.frame(object@MutationTypes))} else {return("Empty")}
})
