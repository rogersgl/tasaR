###############################################################################
# readSettings()
# --------------


#' Read settings for tasaR
#'
#' @description
#' The tasaR package has a relatively complex set of settings that must be defined for each sample. This function reads those settings from either a .csv file on disk or an R object (list/data.frame) and generates an appropriate S4 object of class tas.object.settings.
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



#' Make a settings object for tasaR
#'
#' @param Name Character vector of sample name
#' @param IsAntibody Logical vector of if sequence is an antibody
#' @param MergedFASTQPath Character vector of file path to the .fastq file to analyze
#' @param ReferenceSequence Character vector of the expected amplicon sequence
#' @param ForwardPrimer Character vector of the forward primer DNA sequence that _binds to the amplicon_
#' @param ReversePrimer Character vector of the reverse primer DNA sequence that _binds to the amplicon_
#' @param AmpliconLength Total length of the amplicon (exluding sequencing adapters)
#' @param InsertStart Integer of the first nt of the reference sequence
#' @param InsertEnd Integer of the last nt of the reference sequence. Should be a negative number, count _backwards_ from the 3' end.
#' @param ForwardExtensionType Character vector of the type of extension on the forward primer. Can be ("", "Barcode", or "UMI").
#' @param ForwardExtension Character vector of the DNA sequence of the forward extension. Accepts standard ambiguity codes.
#' @param ReverseExtensionType Character vector of the type of extension on the forward primer. Can be ("", "Barcode", or "UMI").
#' @param ReverseExtension Character vector of the DNA sequence of the reverse extension. Accepts standard ambiguity codes.
#' @param FR1Start Integer of first nt of the FR1 region if sequence is an antibody. Defaults to NA.
#' @param CDR1Start Integer of first nt of the CDR1 region if sequence is an antibody. Defaults to NA.
#' @param FR2Start Integer of first nt of the FR2 region if sequence is an antibody. Defaults to NA.
#' @param CDR2Start Integer of first nt of the CDR2 region if sequence is an antibody. Defaults to NA.
#' @param FR3Start Integer of first nt of the FR3 region if sequence is an antibody. Defaults to NA.
#' @param CDR3Start Integer of first nt of the CDR3 region if sequence is an antibody. Defaults to NA.
#' @param FR4Start Integer of first nt of the FR4 region if sequence is an antibody. Defaults to NA.
#' @param FR4End Integer of last nt of the FR4 region if sequence is an antibody. Defaults to NA.
#'
#' @returns An S4 object of class tas.object.settings
#' @export
#'
#' @examples
#' buildSettings(Name = "test", IsAntibody = TRUE, MergedFASTQPath = file.path(tempdir(), "test-merged.fastq"), ReferenceSequence = "GTTCAACTGGTGGAAAGCGG...",
#'               ForwardPrimer = "GTAAAACGACGGCCAGT", ReversePrimer = "CAGGAAACAGCTATGAC", AmpliconLength = 395L, InsertStart = 25L, InsertEnd = -30L,
#'               ForwardExtensionType = "Barcode", ForwardExtension = "GCTAGCC",  ReverseExtensionType = "UMI", ReverseExtension = "NNNYRNNNYRNN",
#'               FR1Start = 1L, CDR1Start = 73L, FR2Start = 97L, CDR2Start = 148L, FR3Start = 172L, CDR3Start = 286L, FR4Start = 310L, FR4End = 342L)
#'

buildSettings <- function(Name = 'unknown',
                          IsAntibody = FALSE,
                          MergedFASTQPath,
                          ReferenceSequence,
                          ForwardPrimer,
                          ReversePrimer,
                          AmpliconLength,
                          InsertStart,
                          InsertEnd,
                          ForwardExtensionType = "",
                          ForwardExtension = "",
                          ReverseExtensionType = "",
                          ReverseExtension = "",
                          FR1Start = NA_integer_,
                          CDR1Start = NA_integer_,
                          FR2Start = NA_integer_,
                          CDR2Start = NA_integer_,
                          FR3Start = NA_integer_,
                          CDR3Start = NA_integer_,
                          FR4Start = NA_integer_,
                          FR4End = NA_integer_) {

  new("tas.object.settings", Name = Name,
      IsAntibody = IsAntibody,
      MergedFASTQPath = MergedFASTQPath,
      ReferenceSequence = ReferenceSequence,
      ForwardPrimer = ForwardPrimer,
      ReversePrimer = ReversePrimer,
      AmpliconLength = AmpliconLength,
      InsertStart = InsertStart,
      InsertEnd = InsertEnd,
      ForwardExtensionType = ForwardExtensionType,
      ForwardExtension = ForwardExtension,
      ReverseExtensionType = ReverseExtensionType,
      ReverseExtension = ReverseExtension,
      AntibodyRegions = c(FR1Start = FR1Start,
                          CDR1Start = CDR1Start,
                          FR2Start = FR2Start,
                          CDR2Start = CDR2Start,
                          FR3Start = FR3Start,
                          CDR3Start = CDR3Start,
                          FR4Start = FR4Start,
                          FR4End = FR4End
      ))
}




###############################################################################
# Helper functions
# ------

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








