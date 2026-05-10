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
#' analyzeAmplicon('/.../folder/settings.csv', row = 2)
#' analyzeAmplicon(data.frame)
#' analyzeAmplicon(list)
analyzeAmplicon <- function(input.settings, ...) {
  if (class(input.settings) != "tas.object.settings") {
    settings <- readSettings(input.settings, ...)
  } else {
    settings <- input.settings
  }

  sequence.table <- buildSequenceTable(settings)
  mut.pos <- measureMutations(sequence.table, settings)
  mut.dna <- measureDNArepair(sequence.table, settings)

  new("AmpliconSequencing", Sequences = sequence.table,
                            Mutations = mut.pos,
                            MutationTypes = mut.dna,
                            Settings = settings)
}
