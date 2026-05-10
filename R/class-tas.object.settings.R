#' @include helper-fns.R
NULL

# -------------------
# tas.object.settings
# -------------------

#' @title S4 class tas.object.settings
#'
#' @description S4 object containing analysis settings for tasaR
#'
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
                                             ReferenceSequence = "list",
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
                          ReferenceSequence = list(""),
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
  snames <- slotNames(object)[1:(length(slotNames(object)) - 1)]
  snames <- snames[-4]
  for (i in snames) {
    if (length(slot(object, i)) != 1) {return(stringr::str_c("@",i," must be a vector of length = 1."))}
  }

  if (length(object@AntibodyRegions) != 8) {return("@AntibodyRegions must be an integer vector of length = 8.")}

  # sanity checks for input values
  if (all(object@ReferenceSequence != "") && any(nchar(object@ReferenceSequence) > object@AmpliconLength)) {return("@ReferenceSequence cannot be longer than the AmpliconLength.")}
  if (any(S4Vectors::grepl("[^ACGTMRWSYKVHDBN-]", unlist(object@ReferenceSequence), ignore.case = TRUE))) {return("@ReferenceSequence must be a DNA sequence. Invalid characters detected.")}

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
