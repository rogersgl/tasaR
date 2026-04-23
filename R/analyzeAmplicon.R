

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



#' @describeIn analyzeAmplicon Analyze multiple amplicons in sequence
#' @param settings.list A list of settings
#' @param rows Optional vector of row numbers to analyze in a .csv input
#' @section settings.list formats:
#' \describe{
#'    \item{`A list of tas.object.settings objects`}{}
#'    \item{`A file path to a .csv file and a vector of row numbers`}{}
#'    \item{`A list of settings as R lists`}{}
#'    \item{`A data.frame of settings, each row is a sample`}{}
#'  }
#' @returns A list of AmpliconSequencing analysis objects
#' @export
#' @examples
#' vanalyzeAmplicon(list(tas.object.settings))
#' vanalyzeAmplicon('/.../folder/settings.csv', rows = 2:11)
#' vanalyzeAmplicon(list(list))
#' #' vanalyzeAmplicon(data.frame)
vanalyzeAmplicon <- function(settings.list, rows = NULL){

  # list of tas.object.settings
  if (class(settings.list) == "list" && all(sapply(settings.list, is, "tas.object.settings"))) {
    return(lapply(settings.list, analyzeAmplicon))
  }

  # a file.path and vector of row numbers
  if (class(settings.list) == "character" && !is.null(rows) && file.exists(settings.list)) {
    return(lapply(rows, function(x) {
            analyzeAmplicon(settings.list, x)
          }))
  }

  # a list of settings lists
  if (class(settings.list) == "list" && all(sapply(settings.list, is, "list"))) {
    return(lapply(settings.list, analyzeAmplicon))
  }

  #a data.frame of settings
  if (class(settings.list) == "data.frame") {
    return(lapply(1:nrow(settings.list), function(x) {
      analyzeAmplicon(settings.list[x,])
    }))
  }

}









