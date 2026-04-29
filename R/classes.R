#' @include helper-fns.R
NULL

###############################################################################
# S4 class definitions and hierarchy
# -------
# AmpliconSequencing
#   |- Sequences (S4) [tas.sequences]
#       |- Table (data.frame)
#         |- Sequences (character)
#         |- Count (numeric)
#         |- Percent (numeric)
#         |- Indels (character)
#         |- BasesChanged (character)
#         |- AA (character)
#         |- ProteinMutation (character)
#       |- Supplemental (data.table)
#         |- Sequence (character)
#         |- UMIs (list)
#         |- IDs (list)
#       |- Alignment (list)
#          |- DNA (PairwiseAlignmentsSingleSubject)
#          |- AA (PairwiseAlignmentsSingleSubject)
#   |- Mutations (S4) [tas.mutations]
#       |- DNA (list)
#          |- AllMutations (data.frame)
#          |- CytosineMutations (data.frame)
#          |- NonCytosineMutations (data.frame)
#          |- MotifSums (numeric)
#             |- AverageAllMutations (numeric)
#             |- AverageCytosineMutations (numeric)
#             |- AverageNonCytosineMutations (numeric)
#             |- FrequencyOfAllMutationsAtCytosines (numeric)
#             |- FrequencyOfAllMutationsAtNonCytosines (numeric)
#       |- AA (list)
#          |- AllMutations (data.frame)
#          |- MutationMatrix (matrix)
#       |- AIDTables (list)
#          |- WRCH (data.frame)
#             |- Motif (character)
#             |- Start (integer)
#             |- End (integer)
#             |- Cytosine (numeric)
#             |- CytosineMutationFrequency (numeric)
#         |- WRCY (data.frame)
#             |- Motif (character)
#             |- Start (integer)
#             |- End (integer)
#             |- Cytosine (numeric)
#             |- CytosineMutationFrequency (numeric)
#   |- MutationTypes (S4) [tas.dna.repair]
#       |- WT (numeric)
#       |- NHEJ (numeric)
#       |- MMEJ (numeric)
#       |- BaseChange (numeric)
#       |- IndelBaseChange (numeric)
#       |- Other (numeric)
#   |- Settings (S4) [tas.object.settings]
#       |- Name (character)
#       |- ReferenceSequence (character)
#       |- IsAntibody (logical)
#       |- MeasureSHM (logical)
#       |- MeasureDNARepairTypes (logical)
#       |- MergedFASTQPath (character)
#       |- ForwardExtensionType (character)
#       |- ForwardExtension (character)
#       |- ForwardPrimer (character)
#       |- ReverseExtensionType (character)
#       |- ReverseExtension (character)
#       |- ReversePrimer (character)
#       |- AmpliconLength (character)
#       |- InsertStart (character)
#       |- InsertEnd (character)
#       |- AntibodyRegions (integer)
#           |- FR1Start (integer)
#           |- CDR1Start (integer)
#           |- FR2Start (integer)
#           |- CDR2Start (integer)
#           |- FR3Start (integer)
#           |- CDR3Start (integer)
#           |- FR4Start (integer)
#           |- FR4End (integer)
# ------



# -------------------
# tas.object.settings
# -------------------


#' @title S4 class tas.object.settings
#' @description S4 object containing analysis settings for tasaR
#' @slot Name character vector of sequence name
#' @slot IsAntibody logical vector if sequence is an antibody
#' @slot MergedFASTQPath character vector of file path to the merged paired-end FASTQ file to be analyzed
#' @slot ReferenceSequence character vector of the reference (expected) amplicon sequence to be analyzed (do not include primers)
#' @slot ForwardExtensionType character vector ("barcode", "umi", or "") of forward primer extension type (if any)
#' @slot ForwardExtension character vector of DNA sequence of the forward extension
#' @slot ForwardPrimer character vector of DNA sequence of the forward primer (only sequence that binds to the target DNA)
#' @slot ReverseExtensionType character vector ("barcode", "umi", or "") of reverse primer extension type (if any)
#' @slot ReverseExtension character vector of DNA sequence of the reverse extension
#' @slot ReversePrimer character vector of DNA sequence of the forward primer (only sequence that binds to the target DNA)
#' @slot AmpliconLength integer vector of the total amplicon length, including primers and extensions (exclude any Illumina sequencing adapters)
#' @slot InsertStart integer vector pointing to the coordinate of the first nt of the ReferenceSequence within the total amplicon
#' @slot InsertEnd NEGATIVE integer vector pointing to the coordinate of the last nt of the ReferenceSequence within the total amplicon (count from the end)
#' @slot AntibodyRegions integer vector containing coordinates of the first nt of each region of an antibody (or all NA if IsAntibody = FALSE)
#'
#' @export
setClass("tas.object.settings", slots = list(Name = "character",
                                             IsAntibody = "logical",
                                             MergedFASTQPath = "character",
                                             ReferenceSequence = "character",
                                             ForwardExtensionType = "character",
                                             ForwardExtension = "character",
                                             ForwardPrimer = "character",
                                             ReverseExtensionType = "character",
                                             ReverseExtension = "character",
                                             ReversePrimer = "character",
                                             AmpliconLength = "integer",
                                             InsertStart = "integer",
                                             InsertEnd = "integer",
                                             AntibodyRegions = "integer"),
         prototype = list(Name = "unknown",
                          IsAntibody = FALSE,
                          MergedFASTQPath = "",
                          ReferenceSequence = "",
                          ForwardExtensionType = "",
                          ForwardExtension = "",
                          ForwardPrimer = "",
                          ReverseExtensionType = "",
                          ReverseExtension = "",
                          ReversePrimer = "",
                          AmpliconLength = NA_integer_,
                          InsertStart = NA_integer_,
                          InsertEnd = NA_integer_,
                          AntibodyRegions = c(FR1Start = NA_integer_,
                                              CDR1Start = NA_integer_,
                                              FR2Start = NA_integer_,
                                              CDR2Start = NA_integer_,
                                              FR3Start = NA_integer_,
                                              CDR3Start = NA_integer_,
                                              FR4Start = NA_integer_,
                                              FR4End = NA_integer_)))






setValidity("tas.object.settings", function(object) {
  # check vector lengths
  for (i in slotNames(object)[1:(length(slotNames(object)) - 1)]) {
    if (length(slot(object, i)) != 1) {return("@",i," must be a vector of length = 1.")}
  }

  if (length(object@AntibodyRegions) != 8) {return("@AntibodyRegions must be an integer vector of length = 8.")}

  # sanity checks for input values
  if (object@ReferenceSequence != "" && nchar(object@ReferenceSequence) > object@AmpliconLength) {return("@ReferenceSequence cannot be longer than the AmpliconLength.")}
  if (any(S4Vectors::grepl("[^ACGTMRWSYKVHDBN-]", object@ReferenceSequence, ignore.case = TRUE))) {return("@ReferenceSequence must be a DNA sequence. Invalid characters detected.")}

  if (any(S4Vectors::grepl("[^ACGTMRWSYKVHDBN-]", object@ForwardExtension, ignore.case = TRUE))) {return("@ForwardExtension must be a DNA sequence. Invalid characters detected.")}
  if (any(S4Vectors::grepl("[^ACGTMRWSYKVHDBN-]", object@ForwardPrimer, ignore.case = TRUE))) {return("@ForwardPrimer must be a DNA sequence. Invalid characters detected.")}
  if (any(S4Vectors::grepl("[^ACGTMRWSYKVHDBN-]", object@ReverseExtension, ignore.case = TRUE))) {return("@ReverseExtension must be a DNA sequence. Invalid characters detected.")}
  if (any(S4Vectors::grepl("[^ACGTMRWSYKVHDBN-]", object@ReversePrimer, ignore.case = TRUE))) {return("@ReversePrimer must be a DNA sequence. Invalid characters detected.")}
  if (!any(c("barcode", "umi", "") %in% tolower(object@ForwardExtensionType))) {return("@ForwardExtensionType must be either 'barcode' or 'umi'.")}
  if (!any(c("barcode", "umi", "") %in% tolower(object@ReverseExtensionType))) {return("@ReverseExtensionType must be either 'barcode' or 'umi'.")}
  if (!is.na(object@InsertStart) && !object@InsertStart >= 1) {return("@InsertStart must be a positive number >= 1.")}
  if (!is.na(object@InsertEnd) &&!object@InsertEnd < 0) {return("@InsertEnd must be a negative number. Count backwards from the 3' end of the amplicon.")}

  if ((object@MergedFASTQPath != "") && !stringr::str_detect(object@MergedFASTQPath, ".fastq")) {
    return("@MergedFASTQPath should be a file path to a .fastq file.")
  }

  region_starts <- c(FR1 = object@AntibodyRegions["FR1Start"],
                     CDR1 = object@AntibodyRegions["CDR1Start"],
                     FR2  = object@AntibodyRegions["FR2Start"],
                     CDR2 = object@AntibodyRegions["CDR2Start"],
                     FR3  = object@AntibodyRegions["FR3Start"],
                     CDR3 = object@AntibodyRegions["CDR3Start"],
                     FR4  = object@AntibodyRegions["FR4Start"],
                     END = object@AntibodyRegions["FR4End"])
  if (!all(is.na(region_starts)) && any(diff(region_starts) <= 0)) {
    stop("@AntibodyRegion start positions must be strictly increasing.")
  }

  return(TRUE)
})





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
                                                                IndelType = NA_character_),
                          Alignments = list(DNA = empty.pass(),
                                            AA = empty.pass()),
                          ReadCounts = NA_integer_))


setValidity("tas.sequences", function(object) {
  if (ncol(object@Table) != 7) {return("tas.sequences must contain exactly 7 columns.")}
  if (any(colnames(object@Table) != c("Sequences", "Count", "Percent", "Indels", "BasesChanged", "AA", "ProteinMutation"))) {return("Column names in @Table are incorrect.")}
  if (any(S4Vectors::grepl("[^ACGTMRWSYKVHDBN-]", object@Table$Sequences, ignore.case = TRUE))) {return("Sequences must be valid DNA sequences.")}
  if (any(S4Vectors::grepl("[^ARNDCQEGHILKMFPSTWYVUOBJZX*+-]", object@Table$AA, ignore.case = TRUE))) {return("AA must be valid protein sequences.")}
  if (any(colnames(object@Supplemental) != c("Index", "UMIs", "IDs", "IndelStart", "IndelType"))) {return("Column names in @Supplemental are incorrect.")}

  if (all(names(object@Alignments) != c("DNA", "AA"))) {return("Invalid list names in @Alignments.")}
  if (class(object@Alignments$DNA) != "PairwiseAlignmentsSingleSubject" || class(object@Alignments$AA) != "PairwiseAlignmentsSingleSubject") {return("Class of elements in @Alignments should be 'PairwiseAlignmentsSingleSubject'")}
  if (length(object@Alignments$DNA) != length(object@Alignments$AA)) {return("DNA and AA alignments are of unequal lengths.")}
  return(TRUE)
})




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
  if (class(object@DNA) != "list" ||
      any(names(object@DNA) != c("AllMutations", "CytosineMutations", "NonCytosineMutations", "MotifSums")) ||
      any(sapply(object@DNA, class) != c("data.frame", "data.frame", "data.frame", "numeric"))) {return("Invalid @DNA.")}
  if (class(object@AA) != "list" ||
      any(names(object@AA) != c("AllMutations", "MutationMatrix")) ||
      any(BiocGenerics::unlist(sapply(object@AA, class)) != c("data.frame", "matrix", "array"))) {return("Invalid @AA.")}
  if (class(object@AIDTables) != "list" ||
      any(names(object@AIDTables) != c("WRCH", "WRCY")) ||
      any(sapply(object@AIDTables, class) != c("data.frame", "data.frame"))) {return("Invalid @AIDTables.")}
  if (class(object@Antibody) != "list" ||
      any(sapply(object@Antibody, names) != c("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")) ||
      any(sapply(object@Antibody, class) != c("numeric"))) {return("Invalid @Antibody")}
  dfs <- list(object@DNA$AllMutations, object@DNA$CytosineMutations, object@DNA$NonCytosineMutations, object@AA$AllMutations)
  for (i in dfs) {
    if (ncol(i) != 2) {return(stringr::str_c("@",i," must have exactly 2 columns."))}
    if (any(!sapply(i, is.numeric))) {return(stringr::str_c("Columns in @",i," must be numeric vectors."))}
    if (any(colnames(i) != c("Position", "MutationFrequency"))) {return(stringr::str_c("Columns in @",i," must be named 'Position' and 'MutationFrequency'."))}
    if (length(i$Position) != length(i$MutationFrequency)) {return(stringr::str_c("Columns in @",i," must be of equal length."))}
  }
  if (!all(sapply(object@DNA$MotifSums, class)=="numeric")) {return("MotifSums must a numeric vector.")}
  if (!all(is.na(object@DNA$MotifSums)) && any(object@DNA$MotifSums < 0)) {return("MotifSums values must be >= 0.")}

  for (i in object@AIDTables) {
    if (any(colnames(i) != c("Motif", "Start", "End", "Cytosine", "CytosineMutationFrequency"))) {return("Invalid column names in @AIDTables")}
    if (class(i$Motif) != "character" || any(!sapply(i[,-1], is.numeric))) {return("Invalid column types in @AIDTables")}
  }


  return(TRUE)
})


# --------------
# tas.dna.repair
# --------------

#' @title S4 class tas.dna.repair
#' @description S4 object showing the inferred frequency of different types of DNA repair
#'
#' @slot WT numeric vector of percent of unmutated wild-type (WT) sequences
#' @slot NHEJ numeric vector of percent of non-homologous end joining (NHEJ) repair (insertion, -1, or -2)
#' @slot MMEJ numeric vector of percent of microhomology-mediated end joining (MMEJ) repair (< -2)
#' @slot BaseChange numeric vector of percent of sequences with mutated bases but no indels
#' @slot IndelBaseChange numeric vector of percent of sequences with mutated bases and indels
#' @slot Other numeric vector of percent of sequences not fitting into any other category above
#'
#'
#' @export
setClass("tas.dna.repair", slots = list(WT = "numeric",
                                        NHEJ = "numeric",
                                        MMEJ = "numeric",
                                        BaseChange = "numeric",
                                        IndelBaseChange = "numeric",
                                        Other = "numeric"),
         prototype = list(WT = NA_real_,
                          NHEJ = NA_real_,
                          MMEJ = NA_real_,
                          BaseChange = NA_real_,
                          IndelBaseChange = NA_real_,
                          Other = NA_real_))

setValidity("tas.dna.repair", function(object) {
  for (i in slotNames(object)) {
    if (length(slot(object, i)) != 1) {return("@",i," must be a numeric vector of length = 1.")}
  }
  if (all(slotNames(object) != c("WT", "NHEJ", "MMEJ", "BaseChange", "IndelBaseChange", "Other"))) {
    return("Column names are incorrect.")
  }
  if (!isEmpty(object)) {
    if (abs(sum(object@WT, object@NHEJ, object@MMEJ, object@BaseChange, object@IndelBaseChange, object@Other) - 100) > 1e7) {
      return("Frequencies do not add up to 100%.")
    }
  }
  return(TRUE)
})


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
#' @aliases alias
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
