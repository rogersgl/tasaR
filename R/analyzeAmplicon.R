#' @include class-AmpliconSequencing.R
NULL

#' Analyze NGS amplicons
#'
#' @description
#' Takes a set of tasaR settings and performs a complete analysis pipeline. _input.settings_ can be a tas.object.settings object, a compliant R list or data.frame, or a file.path to a .csv file.
#' See also \code{\link{buildSettings}}, \code{\link{readSettings}}, and \code{\link{makeSettingsCSV}}.
#'
#' @param input.settings A set of tasaR settings. Can be either a directory path to a .csv file (see makeSettingsCSV)
#' @param ... additional arguments passed through to sub-functions
#' @param progress_callback Optional function receiving structured progress events during analysis
#'
#' @returns An S4 object of class AmpliconSequencing
#'
#' @section Additional arguments:
#' \describe{
#'   \item{`row = NULL`}{set the row number if reading settings from a .csv file}
#'  }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   analyzeAmplicon('/.../folder/settings.csv', row = 2)
#'   analyzeAmplicon(data.frame)
#'   analyzeAmplicon(list)
#' }
analyzeAmplicon <- function(input.settings, ..., progress_callback = NULL) {
  if (!is(input.settings, "tas.object.settings")) {
    settings <- readSettings(input.settings, ...)
  } else {
    settings <- input.settings
  }

  sequence.table <- buildSequenceTable(settings, progress_callback = progress_callback)
  .emit_analysis_progress(progress_callback, settings@Name, 5L, 5L, "Measuring mutations")
  mut.pos <- measureMutations(sequence.table, settings)
  mut.dna <- measureDNArepair(sequence.table, settings)

  new("AmpliconSequencing", Sequences = sequence.table,
                            Mutations = mut.pos,
                            MutationTypes = mut.dna,
                            Settings = settings)
}



#' Summarize amplicon results
#'
#' @description
#' Convenience wrapper to graph all results from a single AmpliconSequencing
#' object and save them along with extracted tables to a separate list. Can
#' also export all results to a specified directory.
#'
#'
#' @param results S4 object of class AmpliconSequencing
#' @param export Logical specifying whether to export results or not
#' @param path (Optional) Character, File path to directory to export results to
#' @param suppressConsoleOutput Logical, specify whether to suppress output to console if appropriate
#'
#' @returns A list containing Graphs and Tables for the sequencing object
#' @export
#'
#' @examples
#' \dontrun{
#' Summarize(AmpliconSequencing, export = TRUE, path = file.path(tempdir(), export))
#' }
Summarize <- function(results, export = FALSE, path = NULL, suppressConsoleOutput = FALSE) {
  if (!dir.exists(path)) dir.create(path)
  this.path <- file.path(path, results@Settings@Name)
  if (!dir.exists(this.path)) dir.create(this.path)
  graphs <- graphResults(results, output = "all")

  if (export) {
    exportTables(results, this.path)
    gn <- names(graphs)
    suppressMessages(invisible(lapply(seq_len(length(graphs)), function(x) {
      ggplot2::ggsave(filename = stringr::str_c(gn[[x]], ".pdf"),
                      device = "pdf",
                      plot = graphs[[x]],
                      path = this.path)
    })))
  }

  out <- list(Graphs = graphs,
       Tables = list("Sequence Table" = getSequenceTable(results),
                     "Read Counts" = getFilterCounts(results),
                     "All DNA Mutations" = getMutationDistributionDNA(results),
                     "Cytosine DNA Mutations" = getMutationDistributionCytosine(results),
                     "Non-Cytosine DNA Mutations" = getMutationDistributionNonCytosine(results),
                     "Motif Sums" = getMutationMotifSums(results),
                     "All Protein Mutations" = getMutationDistributionAA(results),
                     "Protein Mutation Matrix" = getMutationMatrixAA(results),
                     "WRCH Table" = getWRCHTable(results),
                     "WRCY Table" = getWRCYTable(results),
                     "DNA Repair Types" = getMutationTypes(results)
                     )
       )

  if (suppressConsoleOutput) {
    return(invisible(out))
  } else {
    return(out)
  }
}
