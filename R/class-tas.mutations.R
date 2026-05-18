#' @include helper-fns.R
NULL

# -------------
# tas.mutations
# -------------

#' @title S4 class tas.mutations
#' @description S4 object of sequence mutations
#'
#' @slot DNA a list of DNA mutation information
#' @slot AA a list of protein mutation information
#' @slot AIDTables a list of summary information of mutations at AID hotspots
#' @slot Antibody a list of antibody-specific mutation information
#'
#' @section DNA:
#' \describe{
#'  \item{\code{AllMutations}}{data.frame showing position and frequency of mutations at all positions}
#'  \item{\code{CytosineMutations}}{data.frame showing position and frequency of mutations at AID cytosines}
#'  \item{\code{NonCytosineMutations}}{data.frame showing position and frequency of mutations at all positions that are not AID cytosines}
#'  \item{\code{MotifSums}}{numeric vector summarizing the frequency of mutations among different motifs and denominators}
#' }
#'
#'
#' @section AA:
#' \describe{
#'  \item{\code{AllMutations}}{data.frame showing position and frequency of mutations at all positions}
#'  \item{\code{MutationMatrix}}{matrix showing frequency and identity of all protein mutations (WT omitted)}
#' }
#'
#'
#' @section AIDTables Columns:
#' \describe{
#'   \item{\code{Motif}}{character vector DNA motif (top strand = WRCH/WRCY, bottom strand = DGYW/RGYW)}
#'   \item{\code{Start}}{numeric vector of starting coordinates of AID motifs}
#'   \item{\code{End}}{numeric vector of ending coordinates of AID motifs}
#'   \item{\code{Cytosine}}{numeric vector of coordinates of AID cytosines}
#'   \item{\code{CytosineMutationFrequency}}{numeric vector of the mutation frequency at the AID cytosine}
#' }
#'
#' @section Antibody:
#' \describe{
#'   \item{\code{RegionMutations}}{numeric vector of percent deviation from expected mutation rate for each antibody region, assuming an equal distribution across length}
#' }
#'
#'
#' @export
setClass("tas.mutations", slots = list(DNA = "list",
                                       AA = "list",
                                       AIDTables = "list",
                                       Antibody = "list"),
         prototype = list(DNA = list(AllMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                                     CytosineMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                                     NonCytosineMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                                     MotifSums = c(AverageAllMutations = NA_real_,
                                                   AverageCytosineMutations = NA_real_,
                                                   AverageNonCytosineMutations = NA_real_,
                                                   FrequencyOfAllMutationsAtCytosines = NA_real_,
                                                   FrequencyOfAllMutationsAtNonCytosines = NA_real_)),
                          AA = list(AllMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                                    MutationMatrix = matrix()),
                          AIDTables = list(WRCH = data.frame(Motif = character(), Start = integer(), End = integer(), Cytosine = integer(), CytosineMutationFrequency = numeric()),
                                           WRCY = data.frame(Motif = character(), Start = integer(), End = integer(), Cytosine = integer(), CytosineMutationFrequency = numeric())),
                          Antibody = list(RegionMutations = c(FR1 = NA_real_,
                                                              CDR1 = NA_real_,
                                                              FR2 = NA_real_,
                                                              CDR2 = NA_real_,
                                                              FR3 = NA_real_,
                                                              CDR3 = NA_real_,
                                                              FR4 = NA_real_))
         ))


setValidity("tas.mutations", function(object) {
  if (any(slotNames(object) != c("DNA", "AA", "AIDTables", "Antibody"))) {return("Invalid slot names.")}
  if (!is(object@DNA, "list") ||
      any(names(object@DNA) != c("AllMutations", "CytosineMutations", "NonCytosineMutations", "MotifSums")) ||
      any(!sapply(object@DNA, inherits, c("data.frame", "data.frame", "data.frame", "numeric")))) {return("Invalid @DNA.")}
  if (!is(object@AA, "list") ||
      any(names(object@AA) != c("AllMutations", "MutationMatrix")) ||
      any(!BiocGenerics::unlist(sapply(object@AA, inherits, c("data.frame", "matrix", "array"))))) {return("Invalid @AA.")}
  if (!is(object@AIDTables, "list") ||
      any(names(object@AIDTables) != c("WRCH", "WRCY")) ||
      any(!sapply(object@AIDTables, inherits, c("data.frame", "data.frame")))) {return("Invalid @AIDTables.")}
  if (!is(object@Antibody, "list") ||
      any(sapply(object@Antibody, names) != c("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")) ||
      any(!sapply(object@Antibody, inherits, "numeric"))) {return("Invalid @Antibody")}
  dfs <- list(object@DNA$AllMutations, object@DNA$CytosineMutations, object@DNA$NonCytosineMutations, object@AA$AllMutations)
  for (i in dfs) {
    if (ncol(i) != 2) {return(stringr::str_c("@",i," must have exactly 2 columns."))}
    if (any(!sapply(i, is.numeric))) {return(stringr::str_c("Columns in @",i," must be numeric vectors."))}
    if (any(colnames(i) != c("Position", "MutationFrequency"))) {return(stringr::str_c("Columns in @",i," must be named 'Position' and 'MutationFrequency'."))}
    if (length(i$Position) != length(i$MutationFrequency)) {return(stringr::str_c("Columns in @",i," must be of equal length."))}
  }
  if (!all(sapply(object@DNA$MotifSums, inherits, "numeric"))) {return("MotifSums must a numeric vector.")}
  if (!all(is.na(object@DNA$MotifSums)) && any(object@DNA$MotifSums < 0)) {return("MotifSums values must be >= 0.")}

  for (i in object@AIDTables) {
    if (any(colnames(i) != c("Motif", "Start", "End", "Cytosine", "CytosineMutationFrequency"))) {return("Invalid column names in @AIDTables")}
    if (!is(i$Motif, "character") || any(!sapply(i[,-1], is.numeric))) {return("Invalid column types in @AIDTables")}
  }
  return(TRUE)
})
