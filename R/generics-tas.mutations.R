#' @include class-AmpliconSequencing.R
NULL

# ---------
# Mutations
# ---------

### main definition page ###

#' @title Mutation retrieval methods
#' @description
#' Methods to retrieve mutation information from an object of class tas.mutations or AmpliconSequencing
#'
#' @name getMutations
#' @aliases getMutations,tasaR-method
#' @docType methods
NULL


# -----
#  DNA
# -----

### S4 object tas.mutations ###

#' Retrieve mutation distributions
#'
#'
#' @param object S4 object of class tas.mutations or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @describeIn getMutations Extracts a list containing DNA mutation information
#' @export
#'
#' @aliases getDNAMutations
#'
#' @examples
#' getDNAMutations(object)
setGeneric("getDNAMutations", function(object) standardGeneric("getDNAMutations"))


#' @export
setMethod("getDNAMutations", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(list(AllMutations = object@DNA$AllMutations,
                                     CytosineMutations = object@DNA$CytosineMutations,
                                     NonCytosineMutations = object@DNA$NonCytosineMutations,
                                     MotifSums = object@DNA$MotifSums))
  } else {return("Empty")}
})

#' @export
setMethod("getDNAMutations", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(list(AllMutations = object@Mutations@DNA$AllMutations,
                                               CytosineMutations = object@Mutations@DNA$CytosineMutations,
                                               NonCytosineMutations = object@Mutations@DNA$NonCytosineMutations,
                                               MotifSums = object@Mutations@DNA$MotifSums))
  } else {return("Empty")}
})

### All distribtuion ###

#' @describeIn getMutations Extracts a data.frame quantifying mutation rate at all positions in the sequence
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getMutationDistributionDNA(object)
setGeneric("getMutationDistributionDNA", function(object) standardGeneric("getMutationDistributionDNA"))


#' @export
setMethod("getMutationDistributionDNA", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@DNA$AllMutations)} else {return("Empty")}
})


#' @export
setMethod("getMutationDistributionDNA", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@DNA$AllMutations)} else {return("Empty")}
})


### Cytosine mutations ###

#' @describeIn getMutations Extracts a data.frame quantifying mutation rate at all AID hotspot cytosines in the sequence
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getMutationDistributionCytosine(object)
setGeneric("getMutationDistributionCytosine", function(object) standardGeneric("getMutationDistributionCytosine"))


#' @export
setMethod("getMutationDistributionCytosine", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@DNA$CytosineMutations)} else {return("Empty")}
})


#' @export
setMethod("getMutationDistributionCytosine", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@DNA$CytosineMutations)} else {return("Empty")}
})


### Non-cytosine mutations ###


#' @describeIn getMutations Extracts a data.frame quantifying mutation rate at all bases that are not AID hotspot cytosines
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getMutationDistributionNonCytosine(object)
setGeneric("getMutationDistributionNonCytosine", function(object) standardGeneric("getMutationDistributionNonCytosine"))


#' @export
setMethod("getMutationDistributionNonCytosine", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@DNA$NonCytosineMutations)} else {return("Empty")}
})


#' @export
setMethod("getMutationDistributionNonCytosine", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@DNA$NonCytosineMutations)} else {return("Empty")}
})


### Motif sums ###

#' @describeIn getMutations Extracts a vector summarizing mutations at different motifs and denominators
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getMutationMotifSums(object)
setGeneric("getMutationMotifSums", function(object) standardGeneric("getMutationMotifSums"))


#' @export
setMethod("getMutationMotifSums", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@DNA$MotifSums)} else {return("Empty")}
})


#' @export
setMethod("getMutationMotifSums", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@DNA$MotifSums)} else {return("Empty")}
})

# -----
#  AA
# -----

### All Mutations ###

#' @describeIn getMutations Extracts a data.frame of protein mutation positions and frequencies
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getMutationDistributionAA(object)
setGeneric("getMutationDistributionAA", function(object) standardGeneric("getMutationDistributionAA"))


#' @export
setMethod("getMutationDistributionAA", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@AA$AllMutations)} else {return("Empty")}
})


#' @export
setMethod("getMutationDistributionAA", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@AA$AllMutations)} else {return("Empty")}
})


## AA Mutation Matrix ###

#' @describeIn getMutations Extracts a matrix of protein mutation identities and frequencies
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getMutationMatrixAA(object)
setGeneric("getMutationMatrixAA", function(object) standardGeneric("getMutationMatrixAA"))


#' @export
setMethod("getMutationMatrixAA", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@AA$MutationMatrix)} else {return("Empty")}
})


#' @export
setMethod("getMutationMatrixAA", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@AA$MutationMatrix)} else {return("Empty")}
})


# -----
#  AID
# -----

### Both tables ###

#' Retrieve AID mutation tables
#'
#'
#' @param object S4 object of class tas.mutations or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @describeIn getMutations Extracts a list containing tables for WRCH and WRCY AID hotspot motifs
#' @export
#'
#' @examples
#' getAIDTables(object)
setGeneric("getAIDTables", function(object) standardGeneric("getAIDTables"))


#' @export
setMethod("getAIDTables", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object@AIDTables)) {return(list(WRCH = object@AIDTables$WRCH, WRCY = object@AIDTables$WRCY))} else {return("Empty")}
})


#' @export
setMethod("getAIDTables", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations@AIDTables)) {return(list(WRCH = object@Mutations@AIDTables$WRCH, WRCY = object@Mutations@AIDTables$WRCY))} else {return("Empty")}
})


### WRCH Table ###

#' @describeIn getMutations Extracts a data.frame showing positions and mutations at WRCH AID hotspot motifs
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getWRCHTable(object)
setGeneric("getWRCHTable", function(object) standardGeneric("getWRCHTable"))


#' @export
setMethod("getWRCHTable", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object@AIDTables)) {return(object@AIDTables$WRCH)} else {return("Empty")}
})


#' @export
setMethod("getWRCHTable", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations@AIDTables)) {return(object@Mutations@AIDTables$WRCH)} else {return("Empty")}
})


### WRCY Table ###

#' @describeIn getMutations Extracts a data.frame showing positions and mutations at WRCY AID hotspot motifs
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getWRCYTable(object)
setGeneric("getWRCYTable", function(object) standardGeneric("getWRCYTable"))


#' @export
setMethod("getWRCYTable", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object@AIDTables)) {return(object@AIDTables$WRCY)} else {return("Empty")}
})


#' @export
setMethod("getWRCYTable", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations@AIDTables)) {return(object@Mutations@AIDTables$WRCY)} else {return("Empty")}
})



#' @describeIn getMutations Extracts antibody-specific mutation information
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getAbMutations(object)
setGeneric("getAbMutations", function(object) standardGeneric("getAbMutations"))


#' @export
setMethod("getAbMutations", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@Antibody)} else {return("Empty")}
})


#' @export
setMethod("getAbMutations", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@Antibody)} else {return("Empty")}
})
