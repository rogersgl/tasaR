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





