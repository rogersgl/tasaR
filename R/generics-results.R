#' @include class-AmpliconSequencing.R
#' @include utils-plots.R
#' @include utils-plots-batch.R
#' @include utils-alignments.R
NULL

###############################################################################
# Results functions
# -----------------


# --------
# Graphing
# --------


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
#'  \item{\code{all}}{Return a list of all of the below graph types.}
#'  \item{\code{dna.positions}}{A bar chart showing the frequency of mutations at each nucleotide position.}
#'  \item{\code{dna.positions.labeled}}{A bar chart showing the frequency of mutations at each nucleotide position. AID cytosines are highlighted in red.}
#'  \item{\code{aa.positions}}{A bar chart showing the frequency of mutations at each amino acid position.}
#'  \item{\code{aa.positions.labeled}}{A bar chart showing the frequency of mutations at each amino acid position. Amino acids whose codon includes an AID cytosine are highlighted in red.}
#'  \item{\code{aid.boxplot}}{A box plot comparing the mutation frequency at AID cytosines to all other bases.}
#'  \item{\code{aa.mutations}}{A bar chart showing the frequency and identity of all amino acid mutations at each amino acid position.}
#'  \item{\code{mutation.types.stacked}}{Charts showing the frequency of inferred DNA repair types. Using the mutationTypes argument, can generate a "bar" (stacked bar graph), "pie" chart, or "donut" chart.}
#'  \item{\code{histogram}}{Histogram showing the count of sequences with differing numbers of mutations.}
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
setMethod("graphResults", signature = "AmpliconSequencing", function(results, output = "all", labelCDRs = TRUE, mutationTypes = "bar", ...) {
  if (is.null(output)) {stop("Graph output type must be provided.")}

  match.arg(output, c("all", "dna.positions", "dna.positions.labeled",
                      "aa.positions", "aa.positions.labeled",
                      "aid.boxplot", "aa.mutations",
                      "mutation.types", "histogram"
                      ))

  match.arg(mutationTypes, c("bar", "pie", "donut"))

  if (output == "all"){
    o <- list(gg.dna.pos(results),
              gg.dna.pos.cyt(results),
              gg.aa.pos(results),
              gg.aa.pos.cyt(results),
              gg.aid.box(results),
              gg.aa.muts.stacked(results),
              gg.dna.mut.types(results, mutationTypes),
              gg.dna.mut.count.hist(results))
  } else {
  o <- list(switch(output, "dna.positions" = gg.dna.pos(results),
         "dna.positions.labeled" = gg.dna.pos.cyt(results),
         "aa.positions" = gg.aa.pos(results),
         "aa.positions.labeled" = gg.aa.pos.cyt(results),
         "aid.boxplot" = gg.aid.box(results),
         "aa.mutations" = gg.aa.muts.stacked(results),
         "mutation.types" = gg.dna.mut.types(results, mutationTypes),
         "histogram" = gg.dna.mut.count.hist(results),
         NULL))
  }
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





# -------------
# Export tables
# -------------




#' Export tables of tasaR analysis
#'
#' @description
#' Exports tables of the results of tasaR analyses to a series of .csv files, and copies them to a path specified by the user.
#'
#' @param results S4 object of class AmpliconSequencing, tas.mutation, tas.sequences, or tas.dna.repair
#' @param path Character vector specifying the file path to write the .csv files to.
#'
#'
#' @section Output Tables:
#' \describe{
#'  \item{\code{Sequence Table.csv}}{Overview table showing seqences, counts, mutations, and translations. From tas.mutations.}
#'  \item{\code{Read Counts.csv}}{Table showing the read counts at different stages of filtering.}
#'  \item{\code{All DNA Mutations.csv}}{Table showing the % DNA mutation at each nt position.}
#'  \item{\code{Cytosine DNA Mutations.csv}}{Table showing the % DNA mutation at each AID cytosine position.}
#'  \item{\code{Non-Cytosine DNA Mutations.csv}}{Table showing the % DNA mutation at each nt position that is not an AID cytosine.}
#'  \item{\code{Motif Sums.csv}}{Table summarizing the frequency of mutation at different motifs and denominators.}
#'  \item{\code{All Protein Mutations.csv}}{Table showing the % protein mutation at each amino acid position.}
#'  \item{\code{Protein Mutation Matrix.csv}}{Table showing the frequency of each amino acid occurring at each position in the protein sequence.}
#'  \item{\code{WRCH Table.csv}}{Table summarizing the position of AID hotspots (WRCH motif) and mutation frequency at the cytosine.}
#'  \item{\code{WRCY Table.csv}}{Table summarizing the position of AID hotspots (WRCY motif) and mutation frequency at the cytosine.}
#'  \item{\code{DNA Repair Types.csv}}{Table showing the inferred frequency of different DNA repair pathways..}
#' }
#'
#' @usage NULL
#' @returns NULL
#' @export
setGeneric("exportTables", function(results, path, ...) standardGeneric("exportTables"))


#' @export
setMethod("exportTables", signature = c("AmpliconSequencing", "character"), function(results, path, ...) {

  if (!dir.exists(file.path(tempdir(), "export"))) {
    dir.create(file.path(tempdir(), "export"))
  }

  write.csv(getSequenceTable(results), file = file.path(tempdir(), "export", "Sequence Table.csv"))
  write.csv(getFilterCounts(results), file = file.path(tempdir(), "export", "Read Counts.csv"))
  write.csv(getMutationDistributionDNA(results), file = file.path(tempdir(), "export", "All DNA Mutations.csv"))
  write.csv(getMutationDistributionCytosine(results), file = file.path(tempdir(), "export", "Cytosine DNA Mutations.csv"))
  write.csv(getMutationDistributionNonCytosine(results), file = file.path(tempdir(), "export", "Non-Cytosine DNA Mutations.csv"))
  write.csv(getMutationMotifSums(results), file = file.path(tempdir(), "export", "Motif Sums.csv"))
  write.csv(getMutationDistributionAA(results), file = file.path(tempdir(), "export", "All Protein Mutations.csv"))
  write.csv(getMutationMatrixAA(results), file = file.path(tempdir(), "export", "Protein Mutation Matrix.csv"))
  write.csv(getWRCHTable(results), file = file.path(tempdir(), "export", "WRCH Table.csv"))
  write.csv(getWRCYTable(results), file = file.path(tempdir(), "export", "WRCY Table.csv"))
  write.csv(getMutationTypes(results), file = file.path(tempdir(), "export", "DNA Repair Types.csv"))

  if (!dir.exists(path)) {
    dir.create(path)
  }

  files <- list.files(file.path(tempdir(), "export"), full.names = TRUE)
  file.copy(files, path, recursive = TRUE, overwrite = TRUE)
  file.remove(files)

})

#' @export
setMethod("exportTables", signature = c("tas.sequences", "character"), function(results, path, ...) {
  if (!dir.exists(file.path(tempdir(), "export"))) {dir.create(file.path(tempdir(), "export"))}
  write.csv(getSequenceTable(results), file = file.path(tempdir(), "export", "Sequence Table.csv"))
  write.csv(getFilterCounts(results), file = file.path(tempdir(), "export", "Read Counts.csv"))
  if (!dir.exists(path)) {dir.create(path)}
  files <- list.files(file.path(tempdir(), "export"), full.names = TRUE)
  file.copy(files, path, recursive = TRUE, overwrite = TRUE)
  file.remove(files)
})

#' @export
setMethod("exportTables", signature = c("tas.mutations", "character"), function(results, path, ...) {
  if (!dir.exists(file.path(tempdir(), "export"))) {dir.create(file.path(tempdir(), "export"))}
  write.csv(getMutationDistributionDNA(results), file = file.path(tempdir(), "export", "All DNA Mutations.csv"))
  write.csv(getMutationDistributionCytosine(results), file = file.path(tempdir(), "export", "Cytosine DNA Mutations.csv"))
  write.csv(getMutationDistributionNonCytosine(results), file = file.path(tempdir(), "export", "Non-Cytosine DNA Mutations.csv"))
  write.csv(getMutationMotifSums(results), file = file.path(tempdir(), "export", "Motif Sums.csv"))
  write.csv(getMutationDistributionAA(results), file = file.path(tempdir(), "export", "All Protein Mutations.csv"))
  write.csv(getMutationMatrixAA(results), file = file.path(tempdir(), "export", "Protein Mutation Matrix.csv"))
  write.csv(getWRCHTable(results), file = file.path(tempdir(), "export", "WRCH Table.csv"))
  write.csv(getWRCYTable(results), file = file.path(tempdir(), "export", "WRCY Table.csv"))
  if (!dir.exists(path)) {dir.create(path)}
  files <- list.files(file.path(tempdir(), "export"), full.names = TRUE)
  file.copy(files, path, recursive = TRUE, overwrite = TRUE)
  file.remove(files)
})

#' @export
setMethod("exportTables", signature = c("tas.dna.repair", "character"), function(results, path, ...) {
  if (!dir.exists(file.path(tempdir(), "export"))) {dir.create(file.path(tempdir(), "export"))}
  write.csv(getMutationTypes(results), file = file.path(tempdir(), "export", "DNA Repair Types.csv"))
  if (!dir.exists(path)) {dir.create(path)}
  files <- list.files(file.path(tempdir(), "export"), full.names = TRUE)
  file.copy(files, path, recursive = TRUE, overwrite = TRUE)
  file.remove(files)
})





# --------------
# Summary Report
# --------------

# TODO: might need some kind of LaTeX parser or something to write this?






