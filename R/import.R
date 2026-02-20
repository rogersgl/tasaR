#' @name tas_load_dependencies
#' @title Checks installation and loads dependencies.
#' @description
#' Checks installation and silently loads all package dependencies.
#' @returns No returns within R.
#' @export
tas_load_dependencies <- function(){

  # if(!Operating.System %in% c("Linux","MacOS","Windows")){
  #   stop("Operating System must be one of: 'Linux', 'MacOS', or 'Windows'.")
  # }

  OS <- Sys.info()["sysname"]
  #NGS_Config <- c(as.list(NGS_Config),
  #list(Functions=data.frame(Modules=c("Pandaseq","Barcodes","UMIs","AID","MSA"),Status=c(1,1,1,1,10),row.names = c("Pandaseq","Barcodes","UMIs","AID","MSA"))))
  #Modules are planned to allow use of or omit features of analysis, but are not currently implemented 1/14/26
  if (OS=="Darwin"){
    OS <- "MacOS"
  }

  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", quietly = TRUE)
    library("BiocManager", quietly = TRUE)}

  if (!require("R.utils", quietly = TRUE)){
    install.packages("R.utils", quietly = TRUE)
    library("R.utils", quietly = TRUE)}

  if (!require("zip", quietly = TRUE)){
    install("zip", quietly = TRUE)
    library("zip", quietly = TRUE)}

  if (!require("data.table", quietly = TRUE)){
    install("data.table", quietly = TRUE)
    library("data.table", quietly = TRUE)}

  if (!require("openxlsx", quietly = TRUE)){
    install.packages("openxlsx", quietly = TRUE)
    library("openxlsx", quietly = TRUE)}

  if (!require("stringr", quietly = TRUE)){
    install.packages("stringr", quietly = TRUE)
    library("stringr", quietly = TRUE)}

  if (!require("ggplot2", quietly = TRUE)){
    install.packages("ggplot2", quietly = TRUE)
    library("ggplot2", quietly = TRUE)}

  if (!require("ggseqlogo", quietly = TRUE)){
    install.packages("ggseqlogo", quietly = TRUE)
    library("ggseqlogo", quietly = TRUE)}

  if (!require("parallel", quietly = TRUE)){
    install.packages("parallel", quietly = TRUE)
    library("parallel", quietly = TRUE)}

  if (!require("parallelly", quietly = TRUE)){
    install.packages("parallelly", quietly = TRUE)
    library("parallelly", quietly = TRUE)}

  if (!require("seqinr", quietly = TRUE)){
    install.packages("seqinr", quietly = TRUE)
    library("seqinr", quietly = TRUE)}

  if (!require("ape", quietly = TRUE)){
    install.packages("ape", quietly = TRUE)
    library("ape", quietly = TRUE)}

  if (!require("Biostrings", quietly = TRUE)){
    BiocManager::install("Biostrings", quietly = TRUE)
    library("Biostrings", quietly = TRUE)}

  if (!require("ShortRead", quietly = TRUE)){
    BiocManager::install("ShortRead", quietly = TRUE)
    library("ShortRead", quietly = TRUE)}

  if (!require("pwalign", quietly = TRUE)){
    BiocManager::install("pwalign", quietly = TRUE)
    library("pwalign", quietly = TRUE)}

  if (!require("msa", quietly = TRUE)){
    BiocManager::install("msa", quietly = TRUE)
    library("msa", quietly = TRUE)}

  if (!require("ggmsa", quietly = TRUE)){
    BiocManager::install("ggmsa", quietly = TRUE)
    library("ggmsa", quietly = TRUE)}
}




#' @name tas_import
#' @title Read input .xlsx file.
#' @description
#' Reads the input .csv file and imports all parameters that will be used for downstream analysis.
#' @param Input.File Complete file path of the import .csv file.
#' @param config A character vector provided by parent funtion tas_analyze. Can be 'sheet', 'manual', or 'shiny'.
#' @returns A data frame containing all input parameters and configurations required for analysis.
#' @export
tas_import <- function(Input.File){

  #import processing key file "Input.csv"
  if (file.exists(Input.File)==TRUE){
    NGS_Input <- suppressWarnings(read.csv(Input.File))
    NGS_Input[is.na(NGS_Input)] <- ""
  }else{
    if (shiny.env){
      spsComps::shinyCatch(stop("Input spreadsheet not found."),position = "top-center")
    }else {stop("Input spreadsheet not found. Check filepath and format (.csv).")}
  }

    OS <- Sys.info()["sysname"]
    if (OS=="Darwin"){
      OS <- "MacOS"
    }

    NGS_Config <- list(WorkingDirectory=tempdir(),
                       OperatingSystem=OS,
                       nCores = availableCores())

  #check operating system syntax
  if (!NGS_Config$OperatingSystem %in% c("Windows","MacOS","Linux")){
    stop("The host operating system could not be correctly identified.")
  }
  #Testing has not been performed to validate OS identification outside of MacOS
  NGS <- list()
  NGS <- list(Input = NGS_Input,Config = NGS_Config)
  return(NGS)
}



#' @name tas_check
#' @title Check for presence and formatting of all input files.
#' @description
#' Checks whether InputDataFrame and Config.List contain all expected columns/entries. Also checks whether all unpaired read files specified have been provided.
#'
#' @param InputFilePath The complete file path of the input .csv file. Make sure all R1 and R2 files are in the same directory.
#' @param Input.DataFrame A data frame imported from the .csv file specifying details of each sample.
#' @param Config.List A list containing configuration parameters for tasAnalyzer.
#' @param shiny.env A logical (TRUE/FALSE) identifying whether processing is being perfomed by the Shiny app. Assigned by parent functions.
#' @param shiny.fileTable A data.table with information about uploaded files in the Shiny app. Assigned automatically by parent functions.
#' @param WD A character vector indicating the file path to the working directory containing input files.
#'
#' @returns No returns within R.
#' @export

tas_check <- function(InputFilePath,Input.DataFrame, Config.List, shiny.env, shiny.fileTable, WD){
  Input.Columns <- c("SampleName",
                     "ForwardFASTQFileName",
                     "ReverseFASTQFileName",
                     "MergedFASTQFileName",
                     "ForwardExtensionType",
                     "ForwardExtension",
                     "ForwardPrimer",
                     "ReverseExtensionType",
                     "ReverseExtension",
                     "ReversePrimer",
                     "ReferenceSequence",
                     "AmpliconLength",
                     "MaxDeletion",
                     "MaxInsertion",
                     "InsertStart",
                     "InsertEnd",
                     "Antibody",
                     "FR1Start",
                     "CDR1Start",
                     "FR2Start",
                     "CDR2Start",
                     "FR3Start",
                     "CDR3Start",
                     "FR4Start")

  Config.Entries <- c("WorkingDirectory",
                      "OperatingSystem",
                      "nCores",
                      "merge.reads",
                      "measure.shm",
                      "sequence.alignment.count",
                      "read.frequency.limit",
                      "protein.mutations",
                      "PhyloTree",
                      "dna.repair.pathways",
                      "multicore")

  if(FALSE %in% (Input.Columns %in% colnames(Input.DataFrame))){
    if (shiny.env){spsComps::shinyCatch(stop(c("The following columns were not found in the input file: ",
                                               str_flatten(Input.Columns[!Input.Columns %in% colnames(Input.DataFrame)],collapse = ", "))),
                                        position = "top-center")
    }else{
      stop(c("The following columns were not found in the input file: ",
             str_flatten(Input.Columns[!Input.Columns %in% colnames(Input.DataFrame)],collapse = ", ")))
    }
  }

  if(FALSE %in% (Config.Entries %in% names(Config.List))){
    if (shiny.env){spsComps::shinyCatch(stop(c("The following columns were not found in the input file: ",
                                               str_flatten(Config.Entries[!Config.Entries %in% names(Config.List)],collapse = ", "))),
                                        position = "top-center")
    }else{
      stop(c("The following columns were not found in the input file: ",
             str_flatten(Config.Entries[!Config.Entries %in% names(Config.List)],collapse = ", ")))
    }
  }

  if (!isSingleNumber(Config.List$sequence.alignment.count)){
    stop("sequence.alignment.count must be a single number")
  }
  if (Config.List$sequence.alignment.count<0){
    stop("sequence.alignment.count must be > 0")
  }
  if (ceiling(Config.List$sequence.alignment.count)!=Config.List$sequence.alignment.count){
    stop("sequence.alignment.count must be a whole number")
  }
  if (!isSingleNumber(Config.List$read.frequency.limit)){
    stop("read.frequency.limit must be a single number")
  }
  if (Config.List$read.frequency.limit<0){
    stop("read.frequency.limit must be > 0")
  }

  filenames <- c(Input.DataFrame$ForwardFASTQFileName,Input.DataFrame$ReverseFASTQFileName)

  if (shiny.env){
    shiny.filenames <- shiny.fileTable$name
    if (all(filenames %in% shiny.filenames)==FALSE){
      spsComps::shinyCatch(stop(print(c("These files were not found: ",
                                      filenames[!filenames %in% shiny.filenames]))),
                           position = "top-center")
      }
    } else {
      input.dir <- substr(InputFilePath,1,str_locate(InputFilePath,"Input.csv")[,"start"]-1)

      file.copy(str_c(input.dir,filenames)[file.exists(str_c(input.dir,filenames))],
                str_c(WD,filenames)[file.exists(str_c(input.dir,filenames))])
      dirfiles <- dir(WD)
    if (all(filenames %in% dirfiles)==FALSE){
      stop(print(c("These files were not found: ",
                 filenames[!filenames %in% dirfiles])))
    }
  }

}
