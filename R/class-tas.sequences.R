#' @include helper-fns.R
NULL

# -------------
# tas.sequences
# -------------


#' @title S4 class tas.sequences
#' @description S4 object of analyzed sequence table
#' @slot Table A data.frame with a variable number of rows and 7 named columns:
#' \describe{
#'   \item{\code{Sequences}}{character vector of DNA sequences}
#'   \item{\code{Count}}{numeric vector of sequence counts}
#'   \item{\code{Percent}}{numeric vector of sequence frequencies}
#'   \item{\code{Indels}}{character vector annotating insertions and deletions in the sequence}
#'   \item{\code{BasesChanged}}{integer vector counting the number of mutated DNA bases in the sequence}
#'   \item{\code{AA}}{character vector of protein sequences}
#'   \item{\code{ProteinMutation}}{character vector annotating protein mutations in the sequence}
#' }
#'
#' @slot Alignments A list of DNA and protein (AA) pairwise alignments
#' @slot Supplemental A data.table of DNA sequences paired with UMIs and Illumina sequence IDs
#' @slot ReadCounts Integer vector showing number of reads after different filtering steps
#'
#' @export
setClass("tas.sequences", slots = list(Table = "data.frame",
                                       Supplemental = "data.table",
                                       Alignments = "list",
                                       ReadCounts = "integer"),
         prototype = list(Table = data.frame(Sequences = "",
                                             Count = NA_real_,
                                             Percent = NA_real_,
                                             Indels = "",
                                             BasesChanged = NA_real_,
                                             AA = "",
                                             ProteinMutation = ""),
                          Supplemental = data.table::data.table(Index = NA_integer_,
                                                                UMIs = list(),
                                                                IDs = list(),
                                                                IndelStart = NA_real_,
                                                                IndelType = NA_character_,
                                                                RefIdx = NA_real_),
                          Alignments = list(DNA = list(.empty.pass()),
                                            AA = list(.empty.pass())),
                          ReadCounts = NA_integer_))


setValidity("tas.sequences", function(object) {
  if (ncol(object@Table) != 7) {return("tas.sequences must contain exactly 7 columns.")}
  if (any(colnames(object@Table) != c("Sequences", "Count", "Percent", "Indels", "BasesChanged", "AA", "ProteinMutation"))) {return("Column names in @Table are incorrect.")}
  if (any(S4Vectors::grepl("[^ACGTMRWSYKVHDBN-]", object@Table$Sequences, ignore.case = TRUE))) {return("Sequences must be valid DNA sequences.")}
  if (any(S4Vectors::grepl("[^ARNDCQEGHILKMFPSTWYVUOBJZX*+-]", object@Table$AA, ignore.case = TRUE))) {return("AA must be valid protein sequences.")}
  if (any(colnames(object@Supplemental) != c("Index", "UMIs", "IDs", "IndelStart", "IndelType", "RefIdx"))) {return("Column names in @Supplemental are incorrect.")}

  if (all(names(object@Alignments) != c("DNA", "AA"))) {return("Invalid list names in @Alignments.")}
  if (unique(sapply(object@Alignments$DNA, class)) != "PairwiseAlignmentsSingleSubject" || unique(sapply(object@Alignments$AA, class)) != "PairwiseAlignmentsSingleSubject") {return("Class of elements in @Alignments should be 'PairwiseAlignmentsSingleSubject'")}
  if (any(sapply(object@Alignments$DNA, length) != sapply(object@Alignments$AA, length))) {return("DNA and AA alignments are of unequal lengths.")}
  return(TRUE)
})
