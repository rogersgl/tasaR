#' @include helper-fns.R
#' @include class-tas.object.settings.R
#' @include class-tas.sequences.R
#' @include class-tas.mutations.R
#' @include class-tas.dna.repair.R
NULL

###############################################################################
# Map of AmpliconSequencing, the primery S4 class in tasaR. Contains 4 different
# S4 classes that enumerate different aspects of the analysis.
#
# {tas.sequences} contains sequence information as shown
# {tas.mutations} contains information about distribution of mutations
# {tas.dna.repair} contains information about types of DNA repair mutations
# {tas.object.settings} contains settings used for analysis
#
# -----------------------------------------------------------------------------
# AmpliconSequencing                                                          |
#   |- Sequences (S4) [tas.sequences]                                         |
#       |- Table (data.frame)                                                 |
#         |- Sequences (character)                                            |
#         |- Count (numeric)                                                  |
#         |- Percent (numeric)                                                |
#         |- Indels (character)                                               |
#         |- BasesChanged (character)                                         |
#         |- AA (character)                                                   |
#         |- ProteinMutation (character)                                      |
#       |- Supplemental (data.table)                                          |
#         |- Index (integer)                                                  |
#         |- UMIs (list)                                                      |
#         |- IDs (list)                                                       |
#         |- IndelStart (integer)                                             |
#         |- IndelType (character)                                            |
#         |- RefIdx (numeric)                                                 |
#       |- Alignment (list)                                                   |
#          |- DNA (list)                                                      |
#             |- (PairwiseAlignmentsSingleSubject)                            |
#          |- AA (list)                                                       |
#             |- (PairwiseAlignmentsSingleSubject)                            |
#       |- ReadCounts (integer)                                               |
#   |- Mutations (S4) [tas.mutations]                                         |
#       |- DNA (list)                                                         |
#          |- AllMutations (data.frame)                                       |
#          |- CytosineMutations (data.frame)                                  |
#          |- NonCytosineMutations (data.frame)                               |
#          |- MotifSums (numeric)                                             |
#             |- AverageAllMutations (numeric)                                |
#             |- AverageCytosineMutations (numeric)                           |
#             |- AverageNonCytosineMutations (numeric)                        |
#             |- FrequencyOfAllMutationsAtCytosines (numeric)                 |
#             |- FrequencyOfAllMutationsAtNonCytosines (numeric)              |
#       |- AA (list)                                                          |
#          |- AllMutations (data.frame)                                       |
#          |- MutationMatrix (matrix)                                         |
#       |- AIDTables (list)                                                   |
#          |- WRCH (data.frame)                                               |
#             |- Motif (character)                                            |
#             |- Start (integer)                                              |
#             |- End (integer)                                                |
#             |- Cytosine (numeric)                                           |
#             |- CytosineMutationFrequency (numeric)                          |
#         |- WRCY (data.frame)                                                |
#             |- Motif (character)                                            |
#             |- Start (integer)                                              |
#             |- End (integer)                                                |
#             |- Cytosine (numeric)                                           |
#             |- CytosineMutationFrequency (numeric)                          |
#       |- Anitbody (list)                                                    |
#         |- RegionMutations (numeric)                                        |
#   |- MutationTypes (S4) [tas.dna.repair]                                    |
#       |- WT (numeric)                                                       |
#       |- NHEJ (numeric)                                                     |
#       |- MMEJ (numeric)                                                     |
#       |- BaseChange (numeric)                                               |
#       |- IndelBaseChange (numeric)                                          |
#       |- Other (numeric)                                                    |
#   |- Settings (S4) [tas.object.settings]                                    |
#       |- Name (character)                                                   |
#       |- IsAntibody (logical)                                               |
#       |- MergedFASTQPath (character)                                        |
#       |- ReferenceSequence (list)                                           |
#       |- ForwardExtensionType (character)                                   |
#       |- ForwardExtension (character)                                       |
#       |- ForwardPrimer (character)                                          |
#       |- ReverseExtensionType (character)                                   |
#       |- ReverseExtension (character)                                       |
#       |- ReversePrimer (character)                                          |
#       |- AmpliconLength (character)                                         |
#       |- InsertStart (character)                                            |
#       |- InsertEnd (character)                                              |
#       |- AntibodyRegions (integer)                                          |
#           |- FR1Start (integer)                                             |
#           |- CDR1Start (integer)                                            |
#           |- FR2Start (integer)                                             |
#           |- CDR2Start (integer)                                            |
#           |- FR3Start (integer)                                             |
#           |- CDR3Start (integer)                                            |
#           |- FR4Start (integer)                                             |
#           |- FR4End (integer)                                               |
# -----------------------------------------------------------------------------




# ------------------
# AmpliconSequencing
# ------------------

#' @title S4 class AmpliconSequencing
#'
#' @description S4 class that holds analysis results from tasaR
#'
#' @slot Sequences S4 object of class \code{\linkS4class{tas.sequences}}
#' @slot Mutations S4 object of class \code{\linkS4class{tas.mutations}}
#' @slot MutationTypes S4 object of class \code{\linkS4class{tas.dna.repair}}
#' @slot Settings S4 object of class \code{\linkS4class{tas.object.settings}}
#'
#'
#' @section List of accessor functions:
#'
#' \describe{
#'   \item{`getAlignments()`}{returns list of DNA and AA alignments}
#'   \item{`getDNAalign()`}{returns S4 object of DNA pairwise alignment}
#'   \item{`getAAalign()`}{returns S4 object of protein pairwise alignment}
#'   \item{`getSequenceTable()`}{returns data.frame table of sequences}
#'   \item{`getMutations()`}{returns list comprising the Mutations slot}
#'   \item{`getMutationDistribution()`}{returns data.frame table of positions and frequencies all mutations}
#'   \item{`getMutationDistributionCytosine()`}{returns data.frame table of positions and frequences of all mutations at AID cytosines}
#'   \item{`getMutationDistributionNonCytosine()`}{returns data.frame table of positions and frequences of all mutations at all bases that are not AID cytosines}
#'   \item{`getMutationMotifSums()`}{returns numeric vector summarizing mutations across different patterns and denominators}
#'   \item{`getMutationTypes()`}{returns numeric vector predicting DNA repair pathway usage}
#'   \item{`getAIDTables()`}{returns list of AID cytosine mutation summary tables}
#'   \item{`getWRCHTable()`}{returns data.frame table of positions and mutations of AID cytosines and their motifs for pattern WRCH}
#'   \item{`getWRCYTable()`}{returns data.frame table of positions and mutations of AID cytosines and their motifs for pattern WRCY}
#'   \item{`getSettings()`}{returns list of object-specific tasaR settings used for analysis}
#' }
#'
#'
#' @export
setClass("AmpliconSequencing", slots = list(Sequences = "tas.sequences",
                                            Mutations = "tas.mutations",
                                            MutationTypes = "tas.dna.repair",
                                            Settings = "tas.object.settings"),
         prototype = list(Sequences = new("tas.sequences"),
                          Mutations = new("tas.mutations"),
                          MutationTypes = new("tas.dna.repair"),
                          Settings = new("tas.object.settings")))


setValidity("AmpliconSequencing", function(object) {
  if (length(slotNames(object)) != 4) {
    return("Incorrect number of slots.")
  }
  for (i in slotNames(object)) {
    if (!validObject(slot(object, i))) {
      return(stringr::str_c("Invalid slot ",i," detected."))
    }
  }
  return(TRUE)
})
