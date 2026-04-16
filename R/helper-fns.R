# So far, have been building a script, not a package
# Everything needs to be re-factored to better reflect the package framework
#   - functions designed for 1 sample
#   - how to manage settings?
#   - split up export functions into individual graph types

# Option to import from a .csv file




###############################################################################
# Helper functions
# ------



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
# --------------
# Setter methods
# --------------








