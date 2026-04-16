###############################################################################
# S4 class definitions and hierarchy
# -------
# AmpliconSequencing
#   |- Alignment (S4)
#       |- DNA (PairwiseAlignmentsSingleSubject)
#       |- AA (PairwiseAlignmentsSingleSubject)
#   |- Sequences (S4)
#       |- Table
#         |- Sequences (character)
#         |- Count (numeric)
#         |- Percent (numeric)
#         |- Indels (character)
#         |- BasesChanged (character)
#         |- AA (character)
#         |- ProteinMutation (character)
#       |- Supplemental
#         |- Sequence (character)
#         |- UMIs (list)
#         |- IDs (list)
#   |- Mutations (S4)
#       |- AllMutations (data.frame)
#       |- CytosineMutations (data.frame)
#       |- NonCytosineMutations (data.frame)
#       |- MotifSums (numeric)
#           |- TotalPercentMutated (numeric)
#           |- TotalPercentMutatedInMotif (numeric)
#           |- TotalPercentMutatedInCytosine (numeric)
#           |- MotifMutationAverage (numeric)
#           |- MotifMutationFrequencyofTotal (numeric)
#           |- MotifMutationNonCAverage (numeric)
#           |- MotifMutationNonCFrequencyOfTotal (numeric)
#           |- CytosineMutationAverage (numeric)
#           |- NonCytosineMutationAverage (numeric)
#           |- CytosineMutationFrequencyOfTotal (numeric)
#           |- CytosineMutationFrequencyOfMotif (numeric)
#   |- MutationTypes (S4)
#       |- WT (numeric)
#       |- NHEJ (numeric)
#       |- MMEJ (numeric)
#       |- BaseChange (numeric)
#       |- IndelBaseChange (numeric)
#       |- Other (numeric)
#   |- AIDTables (S4)
#       |- WRCH (data.frame)
#           |- Motif (character)
#           |- Start (integer)
#           |- End (integer)
#           |- Cytosine (integer)
#           |- MotifMutagenesis (numeric)
#           |- CytosineMutagenesis (numeric)
#       |- WRCY (data.frame)
#           |- Motif (character)
#           |- Start (integer)
#           |- End (integer)
#           |- Cytosine (integer)
#           |- MotifMutagenesis (numeric)
#           |- CytosineMutagenesis (numeric)
#   |- Settings (S4)
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

setRefClass("tas.global.settings", fields = list(max.deletion = "integer",
                                                 max.insertion = "integer",
                                                 read.frequency.limit = "numeric"))
tasGlobalSettings <- new("tas.global.settings",
                         max.deletion = 25L,
                         max.insertion = 25L,
                         read.frequency.limit = 0.01)



# -------------------
# tas.object.settings
# -------------------


#' @title S4 class tas.object.settings
#' @description S4 object containing analysis settings for tasaR
#' @slot Name character vector of sequence name
#' @slot IsAntibody logical vector if sequence is an antibody
#' @slot MeasureSHM logical vector defining whether to measure somatic hypermutation (SHM)
#' @slot MeasureDNARepairTypes logical vector defining whether to measure DNA repair types
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
                                             MeasureSHM = "logical",
                                             MeasureDNARepairTypes = "logical",
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
                          MeasureSHM = FALSE,
                          MeasureDNARepairTypes = FALSE,
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
  if (any(grepl("[^ACGTMRWSYKVHDBN-]", object@ReferenceSequence, ignore.case = TRUE))) {return("@ReferenceSequence must be a DNA sequence. Invalid characters detected.")}

  if (any(grepl("[^ACGTMRWSYKVHDBN-]", object@ForwardExtension, ignore.case = TRUE))) {return("@ForwardExtension must be a DNA sequence. Invalid characters detected.")}
  if (any(grepl("[^ACGTMRWSYKVHDBN-]", object@ForwardPrimer, ignore.case = TRUE))) {return("@ForwardPrimer must be a DNA sequence. Invalid characters detected.")}
  if (any(grepl("[^ACGTMRWSYKVHDBN-]", object@ReverseExtension, ignore.case = TRUE))) {return("@ReverseExtension must be a DNA sequence. Invalid characters detected.")}
  if (any(grepl("[^ACGTMRWSYKVHDBN-]", object@ReversePrimer, ignore.case = TRUE))) {return("@ReversePrimer must be a DNA sequence. Invalid characters detected.")}
  if (!any(c("barcode", "umi", "") %in% tolower(object@ForwardExtensionType))) {return("@ForwardExtensionType must be either 'barcode' or 'umi'.")}
  if (!any(c("barcode", "umi", "") %in% tolower(object@ReverseExtensionType))) {return("@ReverseExtensionType must be either 'barcode' or 'umi'.")}
  if (!is.na(object@InsertStart) && !object@InsertStart >= 1) {return("@InsertStart must be a positive number >= 1.")}
  if (!is.na(object@InsertEnd) &&!object@InsertEnd < 0) {return("@InsertEnd must be a negative number. Count backwards from the 3' end of the amplicon.")}

  if ((object@MergedFASTQPath != "") && !str_detect(object@MergedFASTQPath, ".fastq")) {
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
# tas.mutations
# -------------

#' @title S4 class tas.mutations
#' @description S4 object of sequence mutations by position
#'
#' @slot AllMutations data.frame showing position and frequency of mutations at all positions
#' @slot CytosineMutations data.frame showing position and frequency of mutations at AID cytosines
#' @slot NonCytosineMutations data.frame showing position and frequency of mutations at all positions that are not AID cytosines
#' @slot MotifSums numeric vector summarizing the frequency of mutations among different motifs and denominators
#'
#' @importFrom S4Vectors isEmpty
#'
#' @export
setClass("tas.mutations", slots = list(AllMutations = "data.frame",
                                       CytosineMutations = "data.frame",
                                       NonCytosineMutations = "data.frame",
                                       MotifSums = "numeric"),
         prototype = list(AllMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                          CytosineMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                          NonCytosineMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                          MotifSums = c(TotalPercentMutated = NA_real_,
                                        TotalPercentMutatedInMotif = NA_real_,
                                        TotalPercentMutatedInCytosine = NA_real_,
                                        MotifMutationAverage = NA_real_,
                                        MotifMutationFrequencyofTotal = NA_real_,
                                        MotifMutationNonCAverage = NA_real_,
                                        MotifMutationNonCFrequencyOfTotal = NA_real_,
                                        CytosineMutationAverage = NA_real_,
                                        NonCytosineMutationAverage = NA_real_,
                                        CytosineMutationFrequencyOfTotal = NA_real_,
                                        CytosineMutationFrequencyOfMotif = NA_real_)))


setValidity("tas.mutations", function(object) {
  for (i in slotNames(object)[1:3]) {
    tdf <- slot(object,i)
    if (ncol(tdf) != 2) {return(str_c("@",i," must have exactly 2 columns."))}
    if (any(sapply(tdf, class) != "numeric")) {return(str_c("Columns in @",i," must be numeric vectors."))}
    if (any(colnames(tdf) != c("Position", "MutationFrequency"))) {return(str_c("Columns in @",i," must be named 'Position' and 'MutationFrequency'."))}
    if (length(tdf$Position) != length(tdf$MutationFrequency)) {return(str_c("Columns in @",i," must be of equal length."))}
  }
  if (!all(sapply(object@MotifSums, class)=="numeric") || !all(sapply(object@MotifSums, length)==1)) {return("MotifSums must a numeric vector and each entry must be of length = 1.")}
  if (!all(is.na(object@MotifSums[i])) && any(object@MotifSums < 0)) {return("MotifSums values must be >= 0.")}

  return(TRUE)
})



# -------------
# tas.sequences
# -------------


#' @title S4 class tas.sequences
#' @description S4 object of analyzed sequence table
#' @slot Table A data.frame with a variable number of rows and 7 named columns:
#' \describe{
#'   \item{Sequences}{character vector of DNA sequences}
#'   \item{Count}{numeric vector of sequence counts}
#'   \item{Percent}{numeric vector of sequence frequencies}
#'   \item{Indels}{character vector annotating insertions and deletions in the sequence}
#'   \item{BasesChanged}{integer vector counting the number of mutated DNA bases in the sequence}
#'   \item{AA}{character vector of protein sequences}
#'   \item{ProteinMutation}{character vector annotating protein mutations in the sequence}
#' }
#'
#' @importFrom S4Vectors isEmpty
#' @import data.table
#'
#' @export
setClass("tas.sequences", slots = list(Table = "data.frame",
                                       Supplemental = "data.table"),
         prototype = list(Table = data.frame(Sequences = "",
                                             Count = NA_real_,
                                             Percent = NA_real_,
                                             Indels = "",
                                             BasesChanged = NA_real_,
                                             AA = "",
                                             ProteinMutation = ""),
                          Supplemental = data.table::data.table(Sequence = "",
                                                                UMIs = list(),
                                                                IDs = list())))

setValidity("tas.sequences", function(object) {
  if (ncol(object@Table) != 7) {return("tas.sequences must contain exactly 7 columns.")}
  if (any(colnames(object@Table) != c("Sequences", "Count", "Percent", "Indels", "BasesChanged", "AA", "ProteinMutation"))) {return("Column names in @Table are incorrect.")}
  if (any(grepl("[^ACGTMRWSYKVHDBN-]", object@Table$Sequences, ignore.case = TRUE))) {return("Sequences must be valid DNA sequences.")}
  if (any(grepl("[^ARNDCQEGHILKMFPSTWYVUOBJZX*-]", object@Table$AA, ignore.case = TRUE))) {return("AA must be valid protein sequences.")}
  if (any(colnames(object@Supplemental) != c("Sequence", "UMIs", "IDs"))) {return("Column names in @Supplemental are incorrect.")}
  if (any(grepl("[^ACGTMRWSYKVHDBN-]", object@Supplemental$Sequence, ignore.case = TRUE))) {return("Sequence entries in @Supplemental must be valid DNA sequences.")}

  return(TRUE)
})



# --------------
# tas.aid.tables
# --------------

#' @title S4 class tas.aid.tables
#' @description S4 object of tables summarizing mutations at AID cytosine motifs. Separate tables for the broader WRCH motif and the more specific WRCY motif.
#'
#' @slot WRCH data.frame. Coordinates and mutations at WRCH AID cytosines.
#' @slot WRCY data.frame. Coordinates and mutations at WRCY AID cytosines.
#'
#' @section Table Columns:
#' \describe{
#'   \item{Motif}{character vector DNA motif (top strand = WRCH/WRCY, bottom strand = DGYW/RGYW)}
#'   \item{Start}{numeric vector of starting coordinates of AID motifs}
#'   \item{End}{numeric vector of ending coordinates of AID motifs}
#'   \item{Cytosine}{numeric vector of coordinates of AID cytosines}
#'   \item{MotifMutagenesis}{numeric vector of the total frequency of mutagenesis within the motif (max 400\%)}
#'   \item{CytosineMutagenesis}{numeric vector of the frequency of mutagenesis at the AID cytosine}
#' }
#' }
#'
#' @importFrom S4Vectors isEmpty
#'
#' @export
setClass("tas.aid.tables", slots = list(WRCH = "data.frame",
                                        WRCY = "data.frame"),
         prototype = list(WRCH = data.frame(Motif = "",
                                            Start = NA_integer_,
                                            End = NA_integer_,
                                            Cytosine = NA_integer_,
                                            MotifMutagenesis = NA_real_,
                                            CytosineMutagenesis = NA_real_),
                          WRCY = data.frame(Motif = "",
                                            Start = NA_integer_,
                                            End = NA_integer_,
                                            Cytosine = NA_integer_,
                                            MotifMutagenesis = NA_real_,
                                            CytosineMutagenesis = NA_real_)))


setValidity("tas.aid.tables", function(object) {
  for (i in slotNames(object)) {
    if (all(class(colnames(slot(object,i))) != c("character", "integer", "integer", "integer", "numeric", "numeric"))){
      stop(str_c("Column data types in ",i," are incorrect."))
    }
    if (!all(colnames(slot(object, i)) == c("Motif", "Start", "End", "Cytosine", "MotifMutagenesis", "CytosineMutagenesis"))) {
      stop(str_c("Column names in ",i," are incorrect."))
    }
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
#' @importFrom S4Vectors isEmpty
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



# -------------
# tas.alignment
# -------------

#' @title S4 class tas.alignment
#' @description S4 class containing pairwise alignments of sequences against the reference. Contains both DNA and protein (AA) alignments.
#'
#' @slot DNA S4 object of class PairwiseAlignmentsSingleSubject for DNA sequences
#' @slot AA S4 object of class PairwiseAlignmentsSingleSubject for protein sequences
#'
#' @importFrom S4Vectors isEmpty
#'
#' @export
setClass("tas.alignment", slots = list(DNA = "PairwiseAlignmentsSingleSubject",
                                       AA = "PairwiseAlignmentsSingleSubject"),
         prototype = list(DNA = empty.pass(),
                          AA = empty.pass()))


setValidity("tas.alignment", function(object) {
  n <- slotNames(object)
  if (length(slot(object, n[1])) != length(slot(object, n[2]))) {
    return("Length of DNA and AA pairwise alignments are not equal.")
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
#' @slot Alignment S4 object of class \linkS4class{tas.alignment}
#' @slot Sequences S4 object of class \linkS4class{tas.sequences}
#' @slot Mutations S4 object of class \linkS4class{tas.mutations}
#' @slot MutationTypes S4 object of class \linkS4class{tas.dna.repair}
#' @slot AIDTables S4 object of class \linkS4class{tas.aid.tables}
#' @slot Settings S4 object of class \linkS4class{tas.object.settings}
#'
#'
#' @section List of accessor functions:
#'
#' \describe{
#'   \item{`getAlignments()`}{returns S4 object of class tas.alignment}
#'   \item{`getDNAalign()`}{returns S4 object of DNA pairwise alignment}
#'   \item{`getAAalign()`}{returns S4 object of protein pairwise alignment}
#'   \item{`getSequenceTable()`}{returns data.frame table of sequences}
#'   \item{`getMutations()`}{returns S4 object of class tas.mutations}
#'   \item{`getMutationDistribution()`}{returns data.frame table of positions and frequencies all mutations}
#'   \item{`getMutationDistributionCytosine()`}{returns data.frame table of positions and frequences of all mutations at AID cytosines}
#'   \item{`getMutationDistributionNonCytosine()`}{returns data.frame table of positions and frequences of all mutations at all bases that are not AID cytosines}
#'   \item{`getMutationMotifSums()`}{returns numeric vector summarizing mutations across different patterns and denominators}
#'   \item{`getMutationTypes()`}{returns numeric vector predicting DNA repair pathway usage}
#'   \item{`getAIDTables()`}{returns S4 object of class tas.aid.tables}
#'   \item{`getWRCHTable()`}{returns data.frame table of positions and mutations of AID cytosines and their motifs for pattern WRCH}
#'   \item{`getWRCYTable()`}{returns data.frame table of positions and mutations of AID cytosines and their motifs for pattern WRCY}
#'   \item{`getSettings()`}{returns S4 object of class tas.object.settings}
#' }
#'
#' @importFrom S4Vectors isEmpty
#'
#' @export
setClass("AmpliconSequencing", slots = list(Alignment = "tas.alignment",
                                            Sequences = "tas.sequences",
                                            Mutations = "tas.mutations",
                                            MutationTypes = "tas.dna.repair",
                                            AIDTables = "tas.aid.tables",
                                            Settings = "tas.object.settings"),
         prototype = list(Alignment = new("tas.alignment"),
                          Sequences = new("tas.sequences"),
                          Mutations = new("tas.mutations"),
                          MutationTypes = new("tas.dna.repair"),
                          AIDTables = new("tas.aid.tables"),
                          Settings = new("tas.object.settings")))


setValidity("AmpliconSequencing", function(object) {
  if (length(slotNames(object)) != 6) {
    return("Incorrect number of slots.")
  }
  for (i in slotNames(object)) {
    if (!validObject(slot(object, i))) {
      return(str_c("Invalid slot ",i," detected."))
    }
  }
  return(TRUE)
})
