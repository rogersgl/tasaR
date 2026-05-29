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
#' @param settings (Optional) tasaR settings as a tas.object.settings object, list, or data.frame. Only used for input object signatures that are not AmpliconSequencing.
#' @param ... Additional parameters
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
#' @usage NULL
#' @returns A list of ggplot2 graph(s) as specified in the arguments
#' @export
setGeneric("graphResults", function(results, output, ...) standardGeneric("graphResults"))

#' @rdname graphResults
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
    o <- list("DNA Mutation Distribution" = gg.dna.pos(results),
              "DNA Mutation Distribution labeled cytosines" = gg.dna.pos.cyt(results),
              "AA Mutation Distribution" = gg.aa.pos(results),
              "AA Mutation Distribution labeled cytosines" = gg.aa.pos.cyt(results),
              "Mutation at cytosines boxplot" = gg.aid.box(results),
              "All AA Mutations" = gg.aa.muts.stacked(results),
              "DNA Repair Types" = gg.dna.mut.types(results, mutationTypes),
              "Mutations per read histogram" = gg.dna.mut.count.hist(results))
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

#' @rdname graphResults
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

  if (!is(settings, "tas.object.settings")) {
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

#' @rdname graphResults
#' @export
setMethod("graphResults", signature = "tas.sequences", function(results, output = NULL, settings = NULL, labelCDRs = TRUE, mutationTypes = "bar", ...) {
  if (is.null(output)) {stop("Graph output type must be provided.")}
  match.arg(output, c("histogram"))

  o <- gg.dna.mut.count.hist(results)
  return(o)
})

#' @rdname graphResults
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
#' @param ... Additional parameters
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

#' @rdname exportTables
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

#' @rdname exportTables
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

#' @rdname exportTables
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

#' @rdname exportTables
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







# ---------------------
# PairedAnalyzeAmplicon
# ---------------------

# TODO: refactor these functions after moving the nuclease gap correction into buildSequenceTable

setGeneric("exportNucleaseAnalysis", function(paired.seq.results, ...) standardGeneric("exportNucleaseAnalysis"))




setMethod("exportNucleaseAnalysis", signature("PairedAmpliconSequencing"), function(paired.seq.results,
                                                                                    export.dir,
                                                                                    gRNA.seq = NULL,
                                                                                    manual.cut.site = NULL,
                                                                                    window.size = 60,
                                                                                    seq.count = 10){
  if (is.null(gRNA.seq) && is.null(manual.cut.site)){stop("Cut site(s) must be specified either by providing the gRNA sequence or manually as a numeric vector. See help at ?exportNucleaseAnalysis.")}
  if (!is.null(gRNA.seq) && !is.null(manual.cut.site)){warning("Manual cut sites and gRNA sequecne detected. Manual settings will override gRNA sequence matching.")}

  ref.seq <- lapply(getSettings(paired.seq.results@Control)$ReferenceSequence, Biostrings::DNAStringSet)

  # determine expected cut site depending on input type
  if (!is.null(gRNA.seq) && is.null(manual.cut.site)){
    manual.cut.site <- .findCutSitesFromSequence(gRNA.seq, ref.seq)
  }

  # get coordinates for window around cut site
  if (is.null(manual.cut.site) || !is.na(manual.cut.site)){
    gRNA.window <- (min(manual.cut.site) - window.size/2):(max(manual.cut.site) + window.size/2)
    if(min(gRNA.window) < 1) {
      gRNA.window <- gRNA.window[gRNA.window >= 1]
    }
    if (max(gRNA.window) > Biostrings::nchar(ref.seq[[1]])) {
      gRNA.window <- gRNA.window[gRNA.window <= Biostrings::nchar(ref.seq[[1]])]
    }
  } else {
    gRNA.window <- 1:nchar(ref.seq[[1]])
  }

  if (!dir.exists(file.path(tempdir(), "paired/"))) {
    dir.create(file.path(tempdir(), "paired/"))
  }

  # build and export control alignment
  residues_per_line <- ifelse(length(gRNA.window) < 100, length(gRNA.window), 100)

  ctrl.seq <- .alignment.seqs.ampseq(paired.seq.results@Control, gRNA.window, seq.count)
  ctrl.aln <- .nucleaseAlignmentCorrection(msa::msaClustalW(ctrl.seq, order = "input"), cut.site = which(gRNA.window == manual.cut.site), reference = reference)
  ctrl.files <- list.files(path = file.path(tempdir(), "paired/"),
                           pattern = "^ctrl-aln\\..*$",
                           full.names = TRUE)
  invisible(file.remove(ctrl.files[file.exists(ctrl.files)]))
  ctrl.tex <- make_texshade_from_dna(ctrl.aln,
                                     file.path(tempdir(),"paired/ctrl-aln.tex"),
                                     reference = 1L,
                                     residues_per_line = residues_per_line,
                                     cut.sites = which(gRNA.window == manual.cut.site))
  invisible(file.remove("ctrl-aln.fasta"))

  # build and export experimental alignment
  expt.seq <- .alignment.seqs.ampseq(paired.seq.results@Experimental, gRNA.window, seq.count)
  expt.aln <- .nucleaseAlignmentCorrection(msa::msaClustalW(expt.seq, order = "input"), cut.site = which(gRNA.window == manual.cut.site), reference = reference)
  expt.files <- list.files(path = file.path(tempdir(), "paired/"),
                           pattern = "^expt-aln\\..*$",
                           full.names = TRUE)
  invisible(file.remove(expt.files[file.exists(expt.files)]))
  expt.tex <- make_texshade_from_dna(expt.aln,
                                     file.path(tempdir(),"paired/expt-aln.tex"),
                                     reference = 1L,
                                     residues_per_line = residues_per_line,
                                     cut.sites = which(gRNA.window == manual.cut.site))
  invisible(file.remove("expt-aln.fasta"))

  # overlay graph of mutation frequency at each base in the window

  p.pair <- gg.dna.pos.paired(paired.seq.results, gRNA.window)
  ggplot2::ggsave("Paired Mutations at Cut Site.pdf", p.pair, path = file.path(tempdir(), "paired"))

})



setMethod("exportNucleaseAnalysis", signature("list"), function(paired.seq.results,
                                                                export.dir,
                                                                gRNA.seq = NULL,
                                                                manual.cut.site = NULL,
                                                                window.size = 60,
                                                                seq.count = 10){
  if (!any(sapply(paired.seq.results, is, "AmpliconSequencing"))) {stop("All elements in the input list should be of class AmpliconSequencing.")}
  if (!"Control" %in% names(paired.seq.results)){stop("Control AmpliconSequencing object should be in a list element named 'Control'.")}
  if (is.null(gRNA.seq) && is.null(manual.cut.site)){stop("Cut site(s) must be specified either by providing the gRNA sequence or manually as a numeric vector. See help at ?exportNucleaseAnalysis.")}

  sample.number <- length(paired.seq.results) - 1
  if (!is.null(gRNA.seq) && sample.number != length(gRNA.seq)){stop("Must provide a gRNA sequence for each experimental element as a character vector.")}
  if (!is.null(manual.cut.site) && sample.number != length(manual.cut.site)){stop("Must provide a vector with a cut site coordinate for each experimental element.")}

  # get refence sequence list from shared control element
  ref.seq <- lapply(getSettings(paired.seq.results$Control)$ReferenceSequence, Biostrings::DNAStringSet)

  if (!dir.exists(file.path(tempdir(), "paired/"))) {
    dir.create(file.path(tempdir(), "paired/"))
  }

  for (n in names(paired.seq.results)[names(paired.seq.results) != "Control"]) {

    # determine expected cut site depending on input type
    if (!is.null(gRNA.seq) && is.null(manual.cut.site)){
      manual.cut.site <- .findCutSitesFromSequence(gRNA.seq[which(names(paired.seq.results) == n)], ref.seq)
    }

    # get coordinates for window around cut site
    if (is.null(manual.cut.site) || !is.na(manual.cut.site)){
      gRNA.window <- (min(manual.cut.site) - window.size/2):(max(manual.cut.site) + window.size/2)
      if(min(gRNA.window) < 1) {
        gRNA.window <- gRNA.window[gRNA.window >= 1]
      }
      if (max(gRNA.window) > Biostrings::nchar(ref.seq[[1]])) {
        gRNA.window <- gRNA.window[gRNA.window <= Biostrings::nchar(ref.seq[[1]])]
      }
    } else {
      gRNA.window <- 1:nchar(ref.seq[[1]])
    }

    # build and export control alignment
    residues_per_line <- ifelse(length(gRNA.window) < 100, length(gRNA.window), 100)

    ctrl.seq <- .alignment.seqs.ampseq(paired.seq.results$Control, gRNA.window, seq.count)
    ctrl.aln <- .nucleaseAlignmentCorrection(msa::msaClustalW(ctrl.seq, order = "input"), cut.site = which(gRNA.window == manual.cut.site), reference = reference)
    ctrl.files <- list.files(path = file.path(tempdir(), "paired/"),
                             pattern = stringr::str_c("^ctrl-aln-", n, "\\..*$"),
                             full.names = TRUE)
    invisible(file.remove(ctrl.files[file.exists(ctrl.files)]))
    ctrl.tex <- make_texshade_from_dna(ctrl.aln,
                                       file.path(tempdir(), stringr::str_c("paired/ctrl-aln-", n, ".tex")),
                                       reference = 1L,
                                       residues_per_line = residues_per_line,
                                       cut.sites = which(gRNA.window == manual.cut.site))
    invisible(file.remove(stringr::str_c("ctrl-aln-", n, ".fasta")))


    # build and export experimental alignment
    expt.seq <- .alignment.seqs.ampseq(paired.seq.results[[n]], gRNA.window, seq.count)
    expt.aln <- .nucleaseAlignmentCorrection(msa::msaClustalW(expt.seq, order = "input"), cut.site = which(gRNA.window == manual.cut.site), reference = reference)
    expt.files <- list.files(path = file.path(tempdir(), "paired/"),
                             pattern = stringr::str_c("^expt-aln-", n, "\\..*$"),
                             full.names = TRUE)
    invisible(file.remove(expt.files[file.exists(expt.files)]))
    expt.tex <- make_texshade_from_dna(expt.aln,
                                       file.path(tempdir(), stringr::str_c("paired/expt-aln-", n, ".tex")),
                                       reference = 1L,
                                       residues_per_line = residues_per_line,
                                       cut.sites = which(gRNA.window == manual.cut.site))
    invisible(file.remove(stringr::str_c("expt-aln-", n, ".fasta")))

    p.pair <- gg.dna.pos.paired(methods::new("PairedAmpliconSequencing",
                                             Control = paired.seq.results$Control,
                                             Experimental = paired.seq.results[[n]]),
                                gRNA.window)
    ggplot2::ggsave(stringr::str_c("Paired Mutations at Cut Site - ", n, ".pdf"), p.pair, path = file.path(tempdir(), "paired"))

  }

})


