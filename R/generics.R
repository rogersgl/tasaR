###############################################################################
# Getter methods
# -------------


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
  if (!isEmpty(object)) {return(object@Table)} else {return("Empty")}
})

#' @describeIn getSequenceTable Method for class AmpliconSequencing
#' @export
setMethod("getSequenceTable", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences)) {return(object@Sequences@Table)} else {return("Empty")}
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
  if (!is.na(object@Sequences@ReadCounts)) {return(object@Sequences@ReadCounts)} else {return("Empty")}
})




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
  if (!isEmpty(object@AA)) {return(object@AA$AllMutations)} else {return("Empty")}
})


#' @export
setMethod("getMutationDistributionAA", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations@AA)) {return(object@Mutations@AA$AllMutations)} else {return("Empty")}
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
  if (!isEmpty(object@AA)) {return(object@AA$MutationMatrix)} else {return("Empty")}
})


#' @export
setMethod("getMutationMatrixAA", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations@AA)) {return(object@Mutations@AA$MutationMatrix)} else {return("Empty")}
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
  if (!isEmpty(object)) {return(as.data.frame(object))} else {return("Empty")}
})

#' @describeIn getMutationTypes Method for class AmpliconSequencing
#' @export
setMethod("getMutationTypes", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@MutationTypes)) {return(as.data.frame(object@MutationTypes))} else {return("Empty")}
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






###############################################################################
# Graphing functions
# ------------------

#' Graph results of tasaR analysis
#'
#' @description
#' Graphs a variety of data outputs from tasaR analysis. Graph type is set based on
#' the 'output' parameter.
#'
#' @param results S4 object of class AmpliconSequencing, tas.mutation, tas.sequences, or tas.dna.repair
#' @param output Character vector specifying the type of graph to make
#' @param labelCDRs (Optional) Logical vector specifying whether the CDR regions (if an antibody) should be labeled on the resulting graphs
#' @param mutationTypes (Optional) Character vector specifying the type of DNA mutation summary graph to make
#'
#'
#' @section Output Types:
#' \describe{
#'  \item{\code{dna.positions}}{A bar chart showing the frequency of mutations at each nucleotide position.}
#'  \item{\code{dna.positions.labeled}}{A bar chart showing the frequency of mutations at each nucleotide position. AID cytosines are highlighted in red.}
#'  \item{\code{aa.positions}}{A bar chart showing the frequency of mutations at each amino acid position.}
#'  \item{\code{aa.positions.labeled}}{A bar chart showing the frequency of mutations at each amino acid position. Amino acids whose codon includes an AID cytosine are highlighted in red.}
#'  \item{\code{aid.boxplot}}{A box plot comparing the mutation frequency at AID cytosines to all other bases.}
#'  \item{\code{aa.mutations}}{A bar chart showing the frequency and identity of all amino acid mutations at each amino acid position.}
#'  \item{\code{mutation.types.stacked}}{Charts showing the frequency of inferred DNA repair types. Using the mutationTypes argument, can generate a "bar" (stacked bar graph), "pie" chart, or "donut" chart.}
#'  \item{\code{histogram}}{Histogram showing the count of sequences with differing numbers of mutations.}
#'  \item{\code{dna.align}}{Sequence alignment of top DNA sequences against the reference sequence. Default 10 sequences, change with optional settings seq.number = x.}
#'  \item{\code{aa.align}}{Sequence alignment of top protein sequences against the reference sequence. Default 10 sequences, change with optional settings seq.number = x.}
#' }
#'
#'
#' @section Input data types for results:
#' \describe{
#'  \item{\code{AmpliconSequencing}}{A complete AmpliconSequencing object contains the data to generate each type of graph.}
#'  \item{\code{tas.mutations}}{A complete tas.mutations object contains the data to generate all graphs except for stacked mutation types and histogram. If labeling CDRs or AID cytosines for position graphs, a settings object of tas.object.settings or an appropriately formatted list or data.frame must also be provided.}
#'  \item{\code{tas.sequences}}{Contains the data to generate a histogram.}
#'  \item{\code{tas.dna.repair}}{Contains the data to generate mutation type graphs.}
#' }
#'
#'
#'
#' @usage NULL
#' @returns A ggplot2 graph as specified in the arguments
#' @export
setGeneric("graphResults", function(results, output, ...) standardGeneric("graphResults"))

#' @export
setMethod("graphResults", signature = "AmpliconSequencing", function(results, output = NULL, labelCDRs = TRUE, mutationTypes = "bar", seq.number = 10, ...) {
  if (is.null(output)) {stop("Graph output type must be provided.")}

  match.arg(output, c("dna.positions", "dna.positions.labeled",
                      "aa.positions", "aa.positions.labeled",
                      "aid.boxplot", "aa.mutations",
                      "mutation.types", "histogram",
                      "dna.align", "aa.align"
                      ))

  match.arg(mutationTypes, c("bar", "pie", "donut"))

  o <- switch(output, "dna.positions" = gg.dna.pos(results),
         "dna.positions.labeled" = gg.dna.pos.cyt(results),
         "aa.positions" = gg.aa.pos(results),
         "aa.positions.labeled" = gg.aa.pos.cyt(results),
         "aid.boxplot" = gg.aid.box(results),
         "aa.mutations" = gg.aa.muts.stacked(results),
         "mutation.types" = gg.dna.mut.types(results, mutationTypes),
         "histogram" = gg.dna.mut.count.hist(results),
         "dna.align" = gg.dna.align(results, seq.number),
         "aa.align" = gg.aa.align(results, seq.number),
         NULL)
  return(o)
})

#' @export
setMethod("graphResults", signature = "tas.mutations", function(results, output = NULL, settings = NULL, labelCDRs = TRUE, mutationTypes = "bar", ...) {
  if (is.null(output)) {stop("Graph output type must be provided.")}

  match.arg(output, c("dna.positions", "dna.positions.labeled",
                      "aa.positions", "aa.positions.labeled",
                      "aid.boxplot", "aa.mutations"
                      ))

  o <- switch(output, "aid.boxplot" = gg.aid.box(results),
              "aa.mutations" = gg.aa.muts.stacked(results),
              NULL)
  if (!is.null(o)) {return(o)}

  if (labelCDRs == FALSE) {
      o <- switch(output, "dna.positions" = gg.dna.pos(results),
                          "aa.positions" = gg.aa.pos(results),
                           NULL)
  }
  if (!is.null(o)) {return(o)}


  if (!class(settings) %in% c("tas.object.settings", "list", "data.frame")) {
    stop("For output ", output, ", settings must be a tas.object.settings object, or an appropriately formatted list or data.frame")
  }

  if (class(settings) != "tas.object.settings") {
    ns <- readSettings(settings)
  } else {
    ns <- settings
  }

  as <- new("AmpliconSequencing", Mutations = results, Settings = ns)

  o <- switch(output, "dna.positions" = gg.dna.pos(as),
                      "dna.positions.labeled" = gg.dna.pos.cyt(as),
                      "aa.positions" = gg.aa.pos(as),
                      "aa.positions.labeled" = gg.aa.pos.cyt(as),
                      NULL)

  return(o)
})

#' @export
setMethod("graphResults", signature = "tas.sequences", function(results, output = NULL, settings = NULL, labelCDRs = TRUE, mutationTypes = "bar", ...) {
  if (is.null(output)) {stop("Graph output type must be provided.")}
  match.arg(output, c("histogram"))

  o <- gg.dna.mut.count.hist(results)
  return(o)
})

#' @export
setMethod("graphResults", signature = "tas.dna.repair", function(results, output = NULL, settings = NULL, labelCDRs = TRUE, mutationTypes = "bar", ...) {
  if (is.null(output)) {stop("Graph output type must be provided.")}
  match.arg(output, c("mutation.types"))
  match.arg(mutationTypes, c("bar", "pie", "donut"))
  o <- gg.dna.mut.types(results, mutationTypes)
  return(o)
})

