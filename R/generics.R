###############################################################################
# Getter methods
# -------------


# ----------
# Alignments
# ----------

### List of alignments ###

#' Retrieve pairwise alignments
#'
#' @description
#' Extracts DNA and/or protein (AA) pairwise alignments from an S4 object of class tas.alignment or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.alignment or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @describeIn getAlignments Extracts a list of DNA and AA pairwise alignments
#' @export
#'
#' @examples
#' getAlignments(object)
setGeneric("getAlignments", function(object) standardGeneric("getAlignments"))

#' @describeIn getAlignments Method for class tas.alignment
#' @export
setMethod("getAlignments", signature(object = "tas.alignment"), function(object) {
  if (!isEmpty(object)) {return(list(DNA = object@DNA, AA = object@AA))}
  else {return("Empty")}
})

#' @describeIn getAlignments Method for class AmpliconSequencing
#' @export
setMethod("getAlignments", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Alignment)) {return(list(DNA = object@Alignment@DNA, AA = object@Alignment@AA))}
  else {return("Empty")}
})


### DNA ###

#' @describeIn getAlignments Extracts a pairwise alignment of DNA sequences
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.alignment or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getDNAalign(object)
setGeneric("getDNAalign", function(object) standardGeneric("getDNAalign"))

#' @describeIn getAlignments Method for class tas.alignment
#' @export
setMethod("getDNAalign", signature(object = "tas.alignment"), function(object) {
  if (!isEmpty(object@DNA)) {return(object@DNA)}
  else {return("Empty")}
})

#' @describeIn getAlignments Method for class AmpliconSequencing
#' @export
setMethod("getDNAalign", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Alignment@DNA)) {return(object@Alignment@DNA)}
  else {return("Empty")}
})


### Protein ###

#' @describeIn getAlignments Extracts a pairwise alignment of protein sequences
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.alignment or AmpliconSequencing
#' @returns NULL
#' @examples
#' getAAalign(object)
setGeneric("getAAalign", function(object) standardGeneric("getAAalign"))

#' @describeIn getAlignments Method for class tas.alignment
#' @export
setMethod("getAAalign", signature(object = "tas.alignment"), function(object) {
  if (!isEmpty(object@AA)) {return(object@AA)}
  else {return("Empty")}
})

#' @describeIn getAlignments Method for class AmpliconSequencing
#' @export
setMethod("getAAalign", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Alignment@AA)) {return(object@Alignment@AA)}
  else {return("Empty")}
})



# --------------
# Sequence Table
# --------------

### Table ###

#' Retrieve table of sequences
#'
#' @description
#' Extracts a table (data.frame) of sequences, including sequence quantification and annotations of DNA and protein mutations detected from an S4 object of class tas.sequences or its parent class AmpliconSequencing.
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
  if (!isEmpty(object)) {return(object@Table)}
  else {return("Empty")}
})

#' @describeIn getSequenceTable Method for class AmpliconSequencing
#' @export
setMethod("getSequenceTable", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences)) {return(object@Sequences@Table)}
  else {return("Empty")}
})


### Supplemental ###


#' @describeIn getSequenceTable Extracts a data.table linking sequences, Illumina read IDs, and UMIs (if applicable)
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
  if (!isEmpty(object)) {return(object@Supplemental)}
  else {return("Empty")}
})

#' @describeIn getSequenceTable Method for class AmpliconSequencing
#' @export
setMethod("getSequenceSupplemental", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences)) {return(object@Sequences@Supplemental)}
  else {return("Empty")}
})

# ---------
# Mutations
# ---------

### S4 object tas.mutations ###

#' Retrieve mutation distributions
#'
#' @description
#' Extracts quantifications of mutations from an S4 object of class tas.mutations or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.mutations or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @describeIn getMutations Extracts a list containing mutation information
#' @export
#'
#' @examples
#' getMutations(object)
setGeneric("getMutations", function(object) standardGeneric("getMutations"))

#' @describeIn getMutations Method for class tas.mutations
#' @export
setMethod("getMutations", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(list(AllMutations = object@AllMutations,
                                     CytosineMutations = object@CytosineMutations,
                                     NonCytosineMutations = object@NonCytosineMutations,
                                     MotifSums = object@MotifSums))}
  else {return("Empty")}
})

#' @describeIn getMutations Method for class AmpliconSequencing
#' @export
setMethod("getMutations", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(list(AllMutations = object@Mutations@AllMutations,
                                               CytosineMutations = object@Mutations@CytosineMutations,
                                               NonCytosineMutations = object@Mutations@NonCytosineMutations,
                                               MotifSums = object@Mutations@MotifSums))}
  else {return("Empty")}
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
#' getMutationDistribution(object)
setGeneric("getMutationDistribution", function(object) standardGeneric("getMutationDistribution"))

#' @describeIn getMutations Method for class tas.mutations
#' @export
setMethod("getMutationDistribution", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@AllMutations)}
  else {return("Empty")}
})

#' @describeIn getMutations Method for class AmpliconSequencing
#' @export
setMethod("getMutationDistribution", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@AllMutations)}
  else {return("Empty")}
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

#' @describeIn getMutations Method for class tas.mutations
#' @export
setMethod("getMutationDistributionCytosine", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@CytosineMutations)}
  else {return("Empty")}
})

#' @describeIn getMutations Method for class AmpliconSequencing
#' @export
setMethod("getMutationDistributionCytosine", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@CytosineMutations)}
  else {return("Empty")}
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

#' @describeIn getMutations Method for class tas.mutations
#' @export
setMethod("getMutationDistributionNonCytosine", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@NonCytosineMutations)}
  else {return("Empty")}
})

#' @describeIn getMutations Method for class AmpliconSequencing
#' @export
setMethod("getMutationDistributionNonCytosine", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@NonCytosineMutations)}
  else {return("Empty")}
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

#' @describeIn getMutations Method for class tas.mutations
#' @export
setMethod("getMutationMotifSums", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@MotifSums)}
  else {return("Empty")}
})

#' @describeIn getMutations Method for class AmpliconSequencing
#' @export
setMethod("getMutationMotifSums", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@MotifSums)}
  else {return("Empty")}
})


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
#' getMutationTypes(object)
setGeneric("getMutationTypes", function(object) standardGeneric("getMutationTypes"))

#' @describeIn getMutationTypes Method for class tas.dna.repair
#' @export
setMethod("getMutationTypes", signature(object = "tas.dna.repair"), function(object) {
  if (!isEmpty(object)) {return(as.data.frame(object))}
  else {return("Empty")}
})

#' @describeIn getMutationTypes Method for class AmpliconSequencing
#' @export
setMethod("getMutationTypes", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@MutationTypes)) {return(as.data.frame(object@MutationTypes))}
  else {return("Empty")}
})




# ----------
# AID Tables
# ----------

### Both tables ###

#' Retrieve AID mutation tables
#'
#' @description
#' Extracts tables summarizing positions and mutations at AID hotspots from an S4 object of class tas.aid.tables or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.aid.tables or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @describeIn getAIDTables Extracts a list containing tables for WRCH and WRCY AID hotspot motifs
#' @export
#'
#' @examples
#' getAIDTables(object)
setGeneric("getAIDTables", function(object) standardGeneric("getAIDTables"))

#' @describeIn getAIDTables Method for class tas.aid.tables
#' @export
setMethod("getAIDTables", signature(object = "tas.aid.tables"), function(object) {
  if (!isEmpty(object)) {return(list(WRCH = object@WRCH, WRCY = object@WRCY))}
  else {return("Empty")}
})

#' @describeIn getAIDTables Method for class AmpliconSequencing
#' @export
setMethod("getAIDTables", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@AIDTables)) {return(list(WRCH = object@AIDTables@WRCH, WRCY = object@AIDTables@WRCY))}
  else {return("Empty")}
})


### WRCH Table ###

#' @describeIn getAIDTables Extracts a data.frame showing positions and mutations at WRCH AID hotspot motifs
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.aid.tables or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getWRCHTable(object)
setGeneric("getWRCHTable", function(object) standardGeneric("getWRCHTable"))

#' @describeIn getAIDTables Method for class tas.aid.tables
#' @export
setMethod("getWRCHTable", signature(object = "tas.aid.tables"), function(object) {
  if (!isEmpty(object)) {return(object@WRCH)}
  else {return("Empty")}
})

#' @describeIn getAIDTables Method for class AmpliconSequencing
#' @export
setMethod("getWRCHTable", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@AIDTables)) {return(object@AIDTables@WRCH)}
  else {return("Empty")}
})


### WRCY Table ###

#' @describeIn getAIDTables Extracts a data.frame showing positions and mutations at WRCY AID hotspot motifs
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.aid.tables or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getWRCYTable(object)
setGeneric("getWRCYTable", function(object) standardGeneric("getWRCYTable"))

#' @describeIn getAIDTables Method for class tas.aid.tables
#' @export
setMethod("getWRCYTable", signature(object = "tas.aid.tables"), function(object) {
  if (!isEmpty(object)) {return(object@WRCY)}
  else {return("Empty")}
})

#' @describeIn getAIDTables Method for class AmpliconSequencing
#' @export
setMethod("getWRCYTable", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@AIDTables)) {return(object@AIDTables@WRCY)}
  else {return("Empty")}
})



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

