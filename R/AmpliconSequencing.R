# So far, have been building a script, not a package
# Everything needs to be re-factored to better reflect the package framework
#   - functions designed for 1 sample
#   - how to manage settings?
#   - split up export functions into individual graph types

# Option to import from a .csv file




###############################################################################
# Helper functions
# ------

# ------------------
# Exported functions
# ------------------

#' Write a .CSV file to input settings
#'
#' @param filedest Destination for the .csv file to be written to. May be a folder or file path ("./" or "./X.csv"). If folder, will be automatically named 'settings.csv'.
#' @export
#'
#' @examples
#' makeSettingsCSV('/.../folder')
#' makeSettingsCSV('/.../folder/file.csv')
makeSettingsCSV <- function(filedest) {
  titles <- slotNames("tas.object.settings")
  abr <- new("tas.object.settings")@AntibodyRegions
  df <- data.frame(matrix("", 1, (length(titles) - 1 + length(abr))))
  colnames(df) <- c(titles[1:(length(titles)-1)], names(abr))
  if (str_detect(filedest, ".csv")){
    write.csv(df, file = filedest, quote = FALSE, row.names = FALSE)
    message("Created file ", filedest, "\n")
  } else {
    write.csv(df, file = file.path(filedest, "settings.csv"), quote = FALSE, row.names = FALSE)
    message("Created file ", file.path(filedest, "settings.csv"),"\n", sep = "")
  }
}



# ---------------
# Hidden funcions
# ---------------


#' @title Function empty.pass
#' @description Helper function to create an empty PairwiseAlignmentsSingleSubject object. Used to initialize S4 class "tas.alignment".
#'
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom pwalign pairwiseAlignment
#'
#' @returns An empty PairwiseAlignmentsSingleSubject object.
empty.pass <- function() {
  empty_pattern <- DNAStringSet(character(0))
  subject <- DNAString("ACGT")
  pairwiseAlignment(pattern = empty_pattern, subject = subject)
}


#' @title Function charDisplayTrim
#' @description Helper function to trim longer character vectors to 12 characters and (...). Used in "show" methods for more compact display with long sequences.
#' @param string character vector to be trimmed
#' @returns a truncated string
charDisplayTrim <- function(string) {
  str_trunc(string, width = 23, side = "right")
}

#' @title Define an S4 object of class tas.object.settings
#' @description Helper function for readSettings().
#'
#' @param df A data.frame
#' @param snames A character vector of slot names for tas.object.settings defined by parent function readSettings().
#'
#' @returns An S4 object of class tas.object.settings
defineSettingsObject <- function(df, snames) {
  if (class(df) != "data.frame") {stop("Invalid input class. Must be a data.frame.")}
  if (nrow(df) != 1) {stop("Invalid number of rows, df must have 1 row only.")}
  if (any(snames != colnames(df))) {stop("Invalid input column names.")}

  tos <- new("tas.object.settings")
  for (i in 1:(length(slotNames("tas.object.settings"))-1)) {
    cl <- class(slot(new("tas.object.settings"), snames[i]))
    slot(tos, snames[i]) <- as(df[,i], cl)
  }
  for (i in (length(snames) - 7):(length(snames))) {
    tos@AntibodyRegions[snames[i]] <- as.integer(df[,i])
  }
  return(tos)
}


#function to determine optimal sequence for each umi and filter non-passing umis as "Rejected"

#' @title Merge data.table by UMIs
#'
#' @description
#' Helper function for tas_sequence_table. Determines the optimal sequence for each UMI and filters non-passing UMIs as "Rejected".
#'
#' @param x Input vector of sequences to compare.
#' @import data.table
#' @returns Either an optimal consensus sequence for the UMI, or "Rejected".
umi_pileup <- function(x){
  t <- table(x)
  d <- data.table::data.table(seq = names(t),count=as.numeric(t))
  setorder(d,-count)
  if((length(d$count)==1) || (d$count[1]>=(d$count[2]+2) && (d$count[2]<d$count[1]*0.7))){
    return(d$seq[1])
  }else{
    return("Rejected")
  }
}




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
#' @description S4 object containing analysis settings for tasa
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


#' Coerce tas.object.settings to list
#'
#' @param x A tas.object.settings object
#' @param ... Passed to data.frame
#' @export
#' @method as.data.frame tas.dna.repair
as.list.tas.object.settings <- function(x) {
  list(Name = x@Name,
       IsAntibody = x@IsAntibody,
       MeasureSHM = x@MeasureSHM,
       MeasureDNARepairTypes = x@MeasureDNARepairTypes,
       MergedFASTQPath = x@MergedFASTQPath,
       ReferenceSequence = x@ReferenceSequence,
       ForwardExtensionType = x@ForwardExtensionType,
       ForwardExtension = x@ForwardExtension,
       ForwardPrimer = x@ForwardPrimer,
       ReverseExtensionType = x@ReverseExtensionType,
       ReverseExtension = x@ReverseExtension,
       ReversePrimer = x@ReversePrimer,
       AmpliconLength = x@AmpliconLength,
       InsertStart = x@InsertStart,
       InsertEnd = x@InsertEnd,
       AntibodyRegions = c(x@AntibodyRegions["FR1Start"],
                           x@AntibodyRegions["CDR1Start"],
                           x@AntibodyRegions["FR2Start"],
                           x@AntibodyRegions["CDR2Start"],
                           x@AntibodyRegions["FR3Start"],
                           x@AntibodyRegions["CDR3Start"],
                           x@AntibodyRegions["FR4Start"],
                           x@AntibodyRegions["FR4End"]))
}




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


#' Coerce tas.dna.repair to data.frame
#'
#' @param x A tas.dna.repair object
#' @param ... Passed to data.frame
#' @export
#' @method as.data.frame tas.dna.repair
as.data.frame.tas.dna.repair <- function(x, ..., row.names = NULL) {
  data.frame(
    WT = x@WT,
    NHEJ = x@NHEJ,
    MMEJ = x@MMEJ,
    BaseChange = x@BaseChange,
    IndelBaseChange = x@IndelBaseChange,
    Other = x@Other,
    row.names = row.names
  )
}


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
#' @description S4 class that holds analysis results from tasa
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




###############################################################################
# Getter methods
# -------------



# ----------
# Alignments
# ----------

### List of alignments ###

#' Retrieve pairwise alignments
#'
#' @description
#' Extracts DNA and/or protein (AA) pairwise alignments from an S4 object of class tas.alignment or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.alignment or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @describeIn getAlignments Extracts a list of DNA and AA pairwise alignments
#' @export
#'
#' @examples
#' getAlignments(object)
setGeneric("getAlignments", function(object) standardGeneric("getAlignments"))

#' @describeIn getAlignments Method for class tas.alignment
#' @export
setMethod("getAlignments", signature(object = "tas.alignment"), function(object) {
  if (!isEmpty(object)) {return(list(DNA = object@DNA, AA = object@AA))}
  else {return("Empty")}
})

#' @describeIn getAlignments Method for class AmpliconSequencing
#' @export
setMethod("getAlignments", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Alignment)) {return(list(DNA = object@Alignment@DNA, AA = object@Alignment@AA))}
  else {return("Empty")}
})


### DNA ###

#' @describeIn getAlignments Extracts a pairwise alignment of DNA sequences
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.alignment or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getDNAalign(object)
setGeneric("getDNAalign", function(object) standardGeneric("getDNAalign"))

#' @describeIn getAlignments Method for class tas.alignment
#' @export
setMethod("getDNAalign", signature(object = "tas.alignment"), function(object) {
  if (!isEmpty(object@DNA)) {return(object@DNA)}
  else {return("Empty")}
})

#' @describeIn getAlignments Method for class AmpliconSequencing
#' @export
setMethod("getDNAalign", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Alignment@DNA)) {return(object@Alignment@DNA)}
  else {return("Empty")}
})


### Protein ###

#' @describeIn getAlignments Extracts a pairwise alignment of protein sequences
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.alignment or AmpliconSequencing
#' @returns NULL
#' @examples
#' getAAalign(object)
setGeneric("getAAalign", function(object) standardGeneric("getAAalign"))

#' @describeIn getAlignments Method for class tas.alignment
#' @export
setMethod("getAAalign", signature(object = "tas.alignment"), function(object) {
  if (!isEmpty(object@AA)) {return(object@AA)}
  else {return("Empty")}
})

#' @describeIn getAlignments Method for class AmpliconSequencing
#' @export
setMethod("getAAalign", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Alignment@AA)) {return(object@Alignment@AA)}
  else {return("Empty")}
})



# --------------
# Sequence Table
# --------------

### Table ###

#' Retrieve table of sequences
#'
#' @description
#' Extracts a table (data.frame) of sequences, including sequence quantification and annotations of DNA and protein mutations detected from an S4 object of class tas.sequences or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.sequences or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @export
#'
#' @examples
#' getSequenceTable(object)
setGeneric("getSequenceTable", function(object) standardGeneric("getSequenceTable"))

#' @describeIn getSequenceTable Method for class tas.sequences
#' @export
setMethod("getSequenceTable", signature(object = "tas.sequences"), function(object) {
  if (!isEmpty(object)) {return(object@Table)}
  else {return("Empty")}
})

#' @describeIn getSequenceTable Method for class AmpliconSequencing
#' @export
setMethod("getSequenceTable", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences)) {return(object@Sequences@Table)}
  else {return("Empty")}
})


### Supplemental ###


#' @describeIn getSequenceTable Extracts a data.table linking sequences, Illumina read IDs, and UMIs (if applicable)
#' @param object S4 object of class tas.sequences or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @export
#'
#' @examples
#' getSequenceSupplemental(object)
setGeneric("getSequenceSupplemental", function(object) standardGeneric("getSequenceSupplemental"))


#' @describeIn getSequenceTable Method for class tas.sequences
#' @export
setMethod("getSequenceSupplemental", signature(object = "tas.sequences"), function(object) {
  if (!isEmpty(object)) {return(object@Supplemental)}
  else {return("Empty")}
})

#' @describeIn getSequenceTable Method for class AmpliconSequencing
#' @export
setMethod("getSequenceSupplemental", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Sequences)) {return(object@Sequences@Supplemental)}
  else {return("Empty")}
})

# ---------
# Mutations
# ---------

### S4 object tas.mutations ###

#' Retrieve mutation distributions
#'
#' @description
#' Extracts quantifications of mutations from an S4 object of class tas.mutations or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.mutations or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @describeIn getMutations Extracts a list containing mutation information
#' @export
#'
#' @examples
#' getMutations(object)
setGeneric("getMutations", function(object) standardGeneric("getMutations"))

#' @describeIn getMutations Method for class tas.mutations
#' @export
setMethod("getMutations", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(list(AllMutations = object@AllMutations,
                                     CytosineMutations = object@CytosineMutations,
                                     NonCytosineMutations = object@NonCytosineMutations,
                                     MotifSums = object@MotifSums))}
  else {return("Empty")}
})

#' @describeIn getMutations Method for class AmpliconSequencing
#' @export
setMethod("getMutations", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(list(AllMutations = object@Mutations@AllMutations,
                                               CytosineMutations = object@Mutations@CytosineMutations,
                                               NonCytosineMutations = object@Mutations@NonCytosineMutations,
                                               MotifSums = object@Mutations@MotifSums))}
  else {return("Empty")}
})

### All distribtuion ###

#' @describeIn getMutations Extracts a data.frame quantifying mutation rate at all positions in the sequence
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getMutationDistribution(object)
setGeneric("getMutationDistribution", function(object) standardGeneric("getMutationDistribution"))

#' @describeIn getMutations Method for class tas.mutations
#' @export
setMethod("getMutationDistribution", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@AllMutations)}
  else {return("Empty")}
})

#' @describeIn getMutations Method for class AmpliconSequencing
#' @export
setMethod("getMutationDistribution", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@AllMutations)}
  else {return("Empty")}
})


### Cytosine mutations ###

#' @describeIn getMutations Extracts a data.frame quantifying mutation rate at all AID hotspot cytosines in the sequence
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getMutationDistributionCytosine(object)
setGeneric("getMutationDistributionCytosine", function(object) standardGeneric("getMutationDistributionCytosine"))

#' @describeIn getMutations Method for class tas.mutations
#' @export
setMethod("getMutationDistributionCytosine", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@CytosineMutations)}
  else {return("Empty")}
})

#' @describeIn getMutations Method for class AmpliconSequencing
#' @export
setMethod("getMutationDistributionCytosine", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@CytosineMutations)}
  else {return("Empty")}
})


### Non-cytosine mutations ###


#' @describeIn getMutations Extracts a data.frame quantifying mutation rate at all bases that are not AID hotspot cytosines
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getMutationDistributionNonCytosine(object)
setGeneric("getMutationDistributionNonCytosine", function(object) standardGeneric("getMutationDistributionNonCytosine"))

#' @describeIn getMutations Method for class tas.mutations
#' @export
setMethod("getMutationDistributionNonCytosine", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@NonCytosineMutations)}
  else {return("Empty")}
})

#' @describeIn getMutations Method for class AmpliconSequencing
#' @export
setMethod("getMutationDistributionNonCytosine", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@NonCytosineMutations)}
  else {return("Empty")}
})


### Motif sums ###

#' @describeIn getMutations Extracts a vector summarizing mutations at different motifs and denominators
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.mutations or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getMutationMotifSums(object)
setGeneric("getMutationMotifSums", function(object) standardGeneric("getMutationMotifSums"))

#' @describeIn getMutations Method for class tas.mutations
#' @export
setMethod("getMutationMotifSums", signature(object = "tas.mutations"), function(object) {
  if (!isEmpty(object)) {return(object@MotifSums)}
  else {return("Empty")}
})

#' @describeIn getMutations Method for class AmpliconSequencing
#' @export
setMethod("getMutationMotifSums", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@Mutations)) {return(object@Mutations@MotifSums)}
  else {return("Empty")}
})


# ----------------
# DNA Repair Types
# ----------------


#' Retrieve DNA repair pathway inferences
#'
#' @description
#' Extracts a table (data.frame) of DNA repair pathway frequencies from an S4 object of class tas.dna.repair or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.dna.repair or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @export
#'
#' @examples
#' getMutationTypes(object)
setGeneric("getMutationTypes", function(object) standardGeneric("getMutationTypes"))

#' @describeIn getMutationTypes Method for class tas.dna.repair
#' @export
setMethod("getMutationTypes", signature(object = "tas.dna.repair"), function(object) {
  if (!isEmpty(object)) {return(as.data.frame(object))}
  else {return("Empty")}
})

#' @describeIn getMutationTypes Method for class AmpliconSequencing
#' @export
setMethod("getMutationTypes", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@MutationTypes)) {return(as.data.frame(object@MutationTypes))}
  else {return("Empty")}
})




# ----------
# AID Tables
# ----------

### Both tables ###

#' Retrieve AID mutation tables
#'
#' @description
#' Extracts tables summarizing positions and mutations at AID hotspots from an S4 object of class tas.aid.tables or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.aid.tables or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @describeIn getAIDTables Extracts a list containing tables for WRCH and WRCY AID hotspot motifs
#' @export
#'
#' @examples
#' getAIDTables(object)
setGeneric("getAIDTables", function(object) standardGeneric("getAIDTables"))

#' @describeIn getAIDTables Method for class tas.aid.tables
#' @export
setMethod("getAIDTables", signature(object = "tas.aid.tables"), function(object) {
  if (!isEmpty(object)) {return(list(WRCH = object@WRCH, WRCY = object@WRCY))}
  else {return("Empty")}
})

#' @describeIn getAIDTables Method for class AmpliconSequencing
#' @export
setMethod("getAIDTables", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@AIDTables)) {return(list(WRCH = object@AIDTables@WRCH, WRCY = object@AIDTables@WRCY))}
  else {return("Empty")}
})


### WRCH Table ###

#' @describeIn getAIDTables Extracts a data.frame showing positions and mutations at WRCH AID hotspot motifs
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.aid.tables or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getWRCHTable(object)
setGeneric("getWRCHTable", function(object) standardGeneric("getWRCHTable"))

#' @describeIn getAIDTables Method for class tas.aid.tables
#' @export
setMethod("getWRCHTable", signature(object = "tas.aid.tables"), function(object) {
  if (!isEmpty(object)) {return(object@WRCH)}
  else {return("Empty")}
})

#' @describeIn getAIDTables Method for class AmpliconSequencing
#' @export
setMethod("getWRCHTable", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@AIDTables)) {return(object@AIDTables@WRCH)}
  else {return("Empty")}
})


### WRCY Table ###

#' @describeIn getAIDTables Extracts a data.frame showing positions and mutations at WRCY AID hotspot motifs
#'
#' @usage NULL
#'
#' @param object An S4 object of class tas.aid.tables or AmpliconSequencing
#' @returns NULL
#' @export
#'
#' @examples
#' getWRCYTable(object)
setGeneric("getWRCYTable", function(object) standardGeneric("getWRCYTable"))

#' @describeIn getAIDTables Method for class tas.aid.tables
#' @export
setMethod("getWRCYTable", signature(object = "tas.aid.tables"), function(object) {
  if (!isEmpty(object)) {return(object@WRCY)}
  else {return("Empty")}
})

#' @describeIn getAIDTables Method for class AmpliconSequencing
#' @export
setMethod("getWRCYTable", signature(object = "AmpliconSequencing"), function(object) {
  if (!isEmpty(object@AIDTables)) {return(object@AIDTables@WRCY)}
  else {return("Empty")}
})



# --------
# Settings
# --------

#' Retrieve object settings
#'
#' @description
#' Extracts a list of tasa settings from an S4 object of class tas.object.settings or its parent class AmpliconSequencing.
#'
#' @param object S4 object of class tas.object.settings or AmpliconSequencing
#' @usage NULL
#' @returns NULL
#' @export
#'
#' @examples
#' getSettings(object)
setGeneric("getSettings", function(object) standardGeneric("getSettings"))

#' @describeIn getSettings Method for class tas.object.settings
#' @export
setMethod("getSettings", signature(object = "tas.object.settings"), function(object) {
  return(as.list(object))
})

#' @describeIn getSettings Method for class AmpliconSequencing
#' @export
setMethod("getSettings", signature(object = "AmpliconSequencing"), function(object) {
  return(as.list(object@Settings))
})



###############################################################################
# --------------
# Setter methods
# --------------

#' Read settings for tasa
#'
#' @description
#' The tasa package has a relatively complex set of settings that must be defined for each sample. This function reads those settings from either a .csv file on disk or an R object (list/data.frame) and generates an appropriate S4 object of class tas.object.settings.
#'
#' @param x The object to be analyzed. May be a file path ("./X.csv"), or a data.frame or list with the appropriate names and classes for the object. Function should automatically detect which type input 'x' is.
#' @param row If using a .csv input, the row number in the sheet that contains the information for this sample. Use the **displayed number**, the function will account for and remove the header row during processing.
#'
#' @returns An S4 object of class tas.object.settings
#' @export
#'
#' @examples
#' readSettings('/.../folder/settings.csv', row = 2)
#' readSettings(data.frame)
#' readSettings(list)
readSettings <- function(x, row = NULL){
  snames <- c(slotNames("tas.object.settings")[1:(length(slotNames("tas.object.settings")) - 1)], names(new("tas.object.settings")@AntibodyRegions))
  if (length(x) == 1 && str_detect(x, ".csv")) {
    message("Processing 'x' as a .csv file.\n")
    if (!file.exists(x)) {stop(str_c("File not found at ",x,"\n"))}
    if (is.null(row)) {stop("Row number NULL is invalid.")}
    t <- suppressWarnings(read.csv(x))
    if (any(colnames(t) != snames)) {stop("Invalid column names.")}
    out <- t[(row-1),]
  } else if (length(x) > 1 || (class(x) == "data.frame" && ncol(x) > 1)) {
    message("Processing 'x' as an R object.\n")
    if (class(x) == "list" && length(x) != length(snames)) {stop("Invalid length for list 'x'.")}
    if (class(x) == "data.frame" && ncol(x) != length(snames)) {stop("Invalid column number for data.frame 'x'.")}
    out <- as.data.frame(x)
  } else {
    stop("Invalid value for x.")
  }
  output <- defineSettingsObject(out, snames)
  return(output)
}



#' Make a labeled table of sequences
#'
#' @description
#' Uses the settings specified by a tas.settings object, including a file path to a merged .fastq file, to generate a labeled table of sequences, with counts, frequencies, DNA mutations, and AA mutations annotated.
#'
#' Will use barcodes for sequence filtering/demultiplexing if provided. Will bin sequences and count by UMIs if provided (only supports 1 UMI currently). Otherwise, will bin sequences and filter based on the global setting tasGlobalSettings$read.frequency.limit (the minimum % to accept a sequence).
#'
#' @param settings An object of S4 class tas.settings
#'
#' @returns An S4 object of class tas.sequences
#' @export
#'
#' @examples
#' buildSequenceTable(settings)
buildSequenceTable <- function(settings) {
  filter.counts <- numeric()
  if (!file.exists(settings@MergedFASTQPath)) {stop("FASTQ file not found.")}
  cat("Reading .fastq file...\n")
  reads <- ShortRead::readFastq(settings@MergedFASTQPath)
  filter.counts <- c(filter.counts, Merged = length(reads))
  cat("Filtering reads...\n")
  reads.filtered <- tas_filter2(reads, settings)
  filter.counts <- c(filter.counts, Filtered = length(reads.filtered))
  reads <- NULL
  cat("Building sequence table...\n")
  tas_sequence_table2(reads.filtered, settings)
}






###############################################################################
# Extended functions
# ------

# ---------
# isEmpty()
# ---------

### main definition page ###

#' @title isEmpty methods for package tasa
#' @description Extends the S4 function isEmpty to detect whether S4 classes in tasa are empty, holding only initialized values.
#' @param x An S4 object from tasa
#' @importMethodsFrom S4Vectors isEmpty
#' @name isEmpty,tasa-method
NULL

### class tas.object.settings ###

#' @aliases NULL
#' @describeIn isEmpty,tasa-method Method for class tas.object.settings
#' @importMethodsFrom S4Vectors isEmpty
#' @export
setMethod("isEmpty", "tas.object.settings", function(x) {
  n <- slotNames(x)
  n <- setdiff(n, c("Name", "AntibodyRegions"))
  empty <- logical()
  for (i in n) {
    if (class(slot(x,i)) == "numeric" || class(slot(x,i)) == "integer"){
      empty <- c(empty, is.na(slot(x,i)))
    } else if (class(slot(x,i)) == "character") {
      empty <- c(empty, (slot(x,i) == ""))
      }
  }
  for (i in names(x@AntibodyRegions)) {
    empty <- c(empty, is.na(x@AntibodyRegions[i]))
  }
  if (all(empty)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})


### class tas.mutations ###

#' @aliases NULL
#' @describeIn isEmpty,tasa-method Methods for class tas.mutations
#' @importMethodsFrom S4Vectors isEmpty
#' @export
setMethod("isEmpty", "tas.mutations", function(x) {
  n <- slotNames(x)
  df_empty <- logical()
  for (i in n[1:3]){
    df_empty <- c(df_empty, isEmpty(slot(x, i)))
  }
  ms <- all(is.na((x@MotifSums)))
  if (all(df_empty) && ms){
    return(TRUE)
  } else {
    return(FALSE)
  }
})


### class tas.sequences ###

#' @aliases NULL
#' @describeIn isEmpty,tasa-method Methods for class tas.sequences
#' @importMethodsFrom S4Vectors isEmpty
#' @export
setMethod("isEmpty", "tas.sequences", function(x) {
  df <- x@Table
  if (nrow(df) < 2 && all(is.na(df[1,c(2,3,5)])) && all(df[1,c(1,4,6,7)] == "")){
    return(TRUE)
  } else {
    return(FALSE)
  }
})


### class tas.aid.tables ###

#' @aliases NULL
#' @describeIn isEmpty,tasa-method Methods for class tas.aid.tables
#' @importMethodsFrom S4Vectors isEmpty
#' @export
setMethod("isEmpty", "tas.aid.tables", function(x) {
  empty.all <- logical()
  for (i in slotNames(x)){
    if (nrow(slot(x, i)) > 1){
      empty.all <- c(empty.all, FALSE)
    } else {
      empty <- logical()
      empty <- c(empty, slot(x, i)[1,1] == "")
      for (j in 2:ncol(slot(x, i))) {
        empty <- c(empty, is.na(slot(x, i)[1,j]))
      }
      if (all(empty)) {
        empty.all <- c(empty.all, TRUE)
      } else {
        empty.all <- c(empty.all, FALSE)
      }
    }
  }
  if (all(empty.all)){
    return(TRUE)
  } else {
    return(FALSE)
  }
})


### class tas.dna.repair ###

#' @aliases NULL
#' @describeIn isEmpty,tasa-method Methods for class tas.dna.repair
#' @importMethodsFrom S4Vectors isEmpty
#' @export
setMethod("isEmpty", "tas.dna.repair", function(x) {
  empty <- logical()
  for (i in slotNames(x)) {
    empty <- c(empty, is.na(slot(x, i)))
  }
  if (all(empty)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})


### class tas.alignment ###

#' @aliases NULL
#' @describeIn isEmpty,tasa-method Methods for class tas.alignment
#' @importMethodsFrom S4Vectors isEmpty
#' @export
setMethod("isEmpty", "tas.alignment", function(x) {
  n <- slotNames(x)
  len <- numeric(length(n))
  for (i in n){
    len <- c(len, length(slot(x, i)) )
  }
  if (sum(len) == 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
})


### class AmpliconSequencing ###

#' @aliases NULL
#' @describeIn isEmpty,tasa-method Methods for class AmpliconSequencing
#' @importMethodsFrom S4Vectors isEmpty
#' @export
setMethod("isEmpty", "AmpliconSequencing", function(x) {
  empty <- logical()
  for (i in slotNames(x)) {
    empty <- c(empty, isEmpty(slot(x, i)))
  }
  if (all(empty)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})





# ------
# show()
# ------


### class tas.object.settings ###

#' @export
setMethod("show", "tas.object.settings", function(object) {
  cat("Settings for ", object@Name, ":\n",
      "IsAntibody : ", object@IsAntibody, "\n",
      "MeasureSHM : ", object@MeasureSHM, "\n",
      "MeasureDNARepairTypes : ", object@MeasureDNARepairTypes, "\n", sep = "")
  for (i in slotNames(object)[6:12]) {
    s <- slot(object, i)
    if (s != ""){
      s <- charDisplayTrim(s)
      cat("",i,": ", s, "\n", sep = "")}
  }
  for (i in slotNames(object)[13:15]) {
    if (!is.na(slot(object, i))){cat("",i,": ", slot(object, i), "\n", sep = "")}
  }
  if (!all(is.na(object@AntibodyRegions))){
    cat("Antibody region coordinates:\n",
        " FR1: ", object@AntibodyRegions["FR1Start"], "\n",
        "CDR1: ", object@AntibodyRegions["CDR1Start"], "\n",
        " FR2: ", object@AntibodyRegions["FR2Start"], "\n",
        "CDR2: ", object@AntibodyRegions["CDR2Start"], "\n",
        " FR3: ", object@AntibodyRegions["FR3Start"], "\n",
        "CDR3: ", object@AntibodyRegions["CDR3Start"], "\n",
        " FR4: ", object@AntibodyRegions["FR4Start"], "\n",
        " End: ", object@AntibodyRegions["FR4End"], "\n", sep = "")
  }
})


### class tas.mutations ###

#' @export
setMethod("show", "tas.mutations", function(object) {
  cat("DNA mutation tables:\n")
  for (i in slotNames(object)[1:3]) {
    df <- slot(object, i)
    if (!isEmpty(df)) {
      cat("",i,":\n", sep = "")
      if (nrow(df) <= 6){
        show(df)
      } else {
        show(df[1:3,])
        cat("\n          ..........          \n\n")
        show(df[(nrow(df)-2):nrow(df),])
      }
    }
    cat("\n")
  }
  cat("Percent mutation for different motifs:\n")
  for (i in names(object@MotifSums)) {
    if (!is.na(object@MotifSums[i])) {
      cat("",i,": ", object@MotifSums[i],"%\n", sep = "")
    }
  }
})


### class tas.sequences ###

#' @export
setMethod("show", "tas.sequences", function(object) {
  df <- object@Table
  if (nrow(df) <= 6){
    for (i in colnames(df)) {
      if (class(df[,i]) == "character"){
        df[,i] <- sapply(df[,i], charDisplayTrim, USE.NAMES = FALSE)
      }
    }
    show(df)
  } else {
    df <- rbind(df[1:3,], df[(nrow(df)-2):nrow(df),])
    for (i in colnames(df)) {
      if (class(df[,i]) == "character"){
        df[,i] <- sapply(df[,i], charDisplayTrim, USE.NAMES = FALSE)
      }
    }
    show(df[1:3,])
    cat("----------\n")
    show(df[(nrow(df)-2):nrow(df),])
  }
})


### class tas.aid.tables ###

#' @export
setMethod("show", "tas.aid.tables", function(object) {
  for (i in slotNames(object)) {
    cat("Slot",i,":\n", sep = "")

    df <- slot(object, i)
    if (nrow(df) <= 6){
      show(df)
    } else {
      df <- rbind(df[1:3,], df[(nrow(df)-2):nrow(df),])
      show(df[1:3,])
      cat("\n               .........................               \n\n")
      show(df[(nrow(df)-2):nrow(df),])
    }

    cat("\n\n")
  }
})


### class tas.dna.repair ###

#' @export
setMethod("show", "tas.dna.repair", function(object) {
  cat("Inferred DNA repair pathways:\n")
  for (i in slotNames(object)) {
    if (!is.na(slot(object, i))){
      cat(i,": ", slot(object, i), "%\n", sep = "")
    }
  }
})


### class AmpliconSequencing ###

#' @export
setMethod("show", "AmpliconSequencing", function(object) {
  cat("An S4 object of class AmpliconSequencing:\n\n")
  for (i in slotNames(object)[1:5]) {
    cat(i,":\n", sep = "")
    o <- slot(object, i)
    if (isEmpty(o)) {
      cat("Empty\n\n")
    } else {
      cat("Access using @",i,"\n\n")
    }
  }
  show(object@Settings)
})



