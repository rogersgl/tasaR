#' @include class-AmpliconSequencing.R
NULL

# --------------
# Sequence Table
# --------------

### Table ###

#' Retrieve sequence information
#'
#' @description
#' Extracts a table (data.frame) of sequences, including sequence quantification and annotations of DNA and protein mutations detected from an S4 object of class tas.sequences or its parent class AmpliconSequencing.
#'
#'
#'
#' @param object S4 object of class tas.sequences or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @export
#'
#' @examples
#' getSequenceTable(object)
setGeneric("getSequenceTable", function(object) standardGeneric("getSequenceTable"))

#' @describeIn getSequenceTable Method for class tas.sequences
#' @export
setMethod("getSequenceTable", signature(object = "tas.sequences"), function(object) {
  if (!isEmpty(object)) {return(object@Table)} else {return("Empty")}
})

#' @describeIn getSequenceTable Method for class AmpliconSequencing
#' @export
setMethod("getSequenceTable", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences)) {return(object@Sequences@Table)} else {return("Empty")}
})


### Simplify Table ###

#' Retrieve a simplified table of sequences
#'
#' @describeIn getSequenceTableSimplified Extracts a sequence table and simplifies entries for a heterozygous reference sequence.
#' For instance, a control heterozygous sample might have 2 distinct wild-type (WT) alleles. This function would combine those
#' into a single entry and add their proportions, giving the total frequency of each mutation regardless of the underlying allele.
#' Sequence information is removed because it differs between paired alleles.
#'
#' @param object S4 object of class tas.sequences, AmpliconSequencing, or PairedAmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @export
#'
#' @examples
#' getSequenceTableSimplified(object)
setGeneric("getSequenceTableSimplified", function(object) standardGeneric("getSequenceTableSimplified"))

#' @describeIn getSequenceTableSimplified Method for class tas.sequences
#' @export
setMethod("getSequenceTableSimplified", signature (object = "tas.sequences"), function(object){
  dt <- as.data.table(getSequenceTable(object))
  dt.s <- dt[, .(Count = sum(Count), Percent = sum(Percent)), by = .(Indels, BasesChanged)]
  setcolorder(dt.s, c("Count", "Percent", "Indels", "BasesChanged"))
  return(as.data.frame(dt.s))
})

#' @describeIn getSequenceTableSimplified Method for class AmpliconSequencing
#' @export
setMethod("getSequenceTableSimplified", signature (object = "AmpliconSequencing"), function(object){
  dt <- as.data.table(getSequenceTable(object))
  dt.s <- dt[, .(Count = sum(Count), Percent = sum(Percent)), by = .(Indels, BasesChanged)]
  setcolorder(dt.s, c("Count", "Percent", "Indels", "BasesChanged"))
  return(as.data.frame(dt.s))
})

#' @describeIn getSequenceTableSimplified Method for class PairedAmpliconSequencing
#' @export
setMethod("getSequenceTableSimplified", signature (object = "PairedAmpliconSequencing"), function(object){
  out <- list()
  out[slotNames(object)] <- lapply(slotNames(object), function(x) {
    dt <- as.data.table(getSequenceTable(slot(object, x)))
    dt.s <- dt[, .(Count = sum(Count), Percent = sum(Percent)), by = .(Indels, BasesChanged)]
    setcolorder(dt.s, c("Count", "Percent", "Indels", "BasesChanged"))
    return(as.data.frame(dt.s))
  })
  return(out)
})

#' @describeIn getSequenceTable Extracts a character vector of unique DNA sequences
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.sequences or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getSequencesDNA(object)
setGeneric("getSequencesDNA", function(object) standardGeneric("getSequencesDNA"))

#' @describeIn getSequenceTable Method for class tas.sequences
#' @export
setMethod("getSequencesDNA", signature(object = "tas.sequences"), function(object) {
  if (!isEmpty(object)) {return(object@Table$Sequences)} else {return("Empty")}
})

#' @describeIn getSequenceTable Method for class AmpliconSequencing
#' @export
setMethod("getSequencesDNA", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences)) {return(object@Sequences@Table$Sequences)} else {return("Empty")}
})


#' @describeIn getSequenceTable Extracts a character vector of unique protein sequences
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.sequences or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getSequencesAA(object)
setGeneric("getSequencesAA", function(object) standardGeneric("getSequencesAA"))

#' @describeIn getSequenceTable Method for class tas.sequences
#' @export
setMethod("getSequencesAA", signature(object = "tas.sequences"), function(object) {
  if (!isEmpty(object)) {return(object@Table$AA)} else {return("Empty")}
})

#' @describeIn getSequenceTable Method for class AmpliconSequencing
#' @export
setMethod("getSequencesAA", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences)) {return(object@Sequences@Table$AA)} else {return("Empty")}
})

### Supplemental ###


#' @describeIn getSequenceTable Extracts a data.table with additional sequence information and metadata
#' @param object S4 object of class tas.sequences or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @export
#'
#' @examples
#' getSequenceSupplemental(object)
setGeneric("getSequenceSupplemental", function(object) standardGeneric("getSequenceSupplemental"))


#' @describeIn getSequenceTable Method for class tas.sequences
#' @export
setMethod("getSequenceSupplemental", signature(object = "tas.sequences"), function(object) {
  if (!isEmpty(object)) {return(
    cbind(data.table(Sequence = object@Table$Sequences), object@Supplemental[,-"Index"])
  )} else {return("Empty")}
})

#' @describeIn getSequenceTable Method for class AmpliconSequencing
#' @export
setMethod("getSequenceSupplemental", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences)) {return(
    cbind(data.table(Sequence = object@Sequences@Table$Sequences), object@Sequences@Supplemental[,-"Index"])
  )} else {return("Empty")}
})


### List of alignments ###

#' Retrieve pairwise alignments
#'
#' @description
#' Extracts DNA and/or protein (AA) pairwise alignments from an S4 object of class tas.sequences or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.sequences or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @describeIn getAlignments Extracts a list of DNA and AA pairwise alignments
#' @export
#'
#' @examples
#' getAlignments(object)
setGeneric("getAlignments", function(object) standardGeneric("getAlignments"))

#' @describeIn getAlignments Method for class tas.sequences
#' @export
setMethod("getAlignments", signature(object = "tas.sequences"), function(object) {
  if (!isEmpty(object)) {return(list(DNA = object@Alignments$DNA, AA = object@Alignments$AA))} else {return("Empty")}
})

#' @describeIn getAlignments Method for class AmpliconSequencing
#' @export
setMethod("getAlignments", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences@Alignments)) {return(list(DNA = object@Sequences@Alignments$DNA, AA = object@Sequences@Alignments$AA))} else {return("Empty")}
})




### DNA ###

#' @describeIn getAlignments Extracts a pairwise alignment of DNA sequences
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.sequences or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getDNAalign(object)
setGeneric("getDNAalign", function(object) standardGeneric("getDNAalign"))

#' @describeIn getAlignments Method for class tas.sequences
#' @export
setMethod("getDNAalign", signature(object = "tas.sequences"), function(object) {
  if (!isEmpty(object@Alignments$DNA)) {return(object@Alignments$DNA)} else {return("Empty")}
})

#' @describeIn getAlignments Method for class AmpliconSequencing
#' @export
setMethod("getDNAalign", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences@Alignments$DNA)) {return(object@Sequences@Alignments$DNA)} else {return("Empty")}
})




### Protein ###

#' @describeIn getAlignments Extracts a pairwise alignment of protein sequences
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.sequences or AmpliconSequencing
#' @returns NULL
#' @examples
#' getAAalign(object)
setGeneric("getAAalign", function(object) standardGeneric("getAAalign"))

#' @describeIn getAlignments Method for class tas.sequences
#' @export
setMethod("getAAalign", signature(object = "tas.sequences"), function(object) {
  if (!isEmpty(object@Alignments$AA)) {return(object@Alignments$AA)} else {return("Empty")}
})

#' @describeIn getAlignments Method for class AmpliconSequencing
#' @export
setMethod("getAAalign", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences@Alignments$AA)) {return(object@Sequences@Alignments$AA)} else {return("Empty")}
})




### Read Counts ###

setGeneric("getFilterCounts", function(object) standardGeneric("getFilterCounts"))

setMethod("getFilterCounts", signature(object = "tas.sequences"), function(object) {
  if (!all(is.na(object@ReadCounts))) {return(object@ReadCounts)} else {return("Empty")}
})

setMethod("getFilterCounts", signature(object = "AmpliconSequencing"), function(object) {
  if (!all(is.na(object@Sequences@ReadCounts))) {return(object@Sequences@ReadCounts)} else {return("Empty")}
})
