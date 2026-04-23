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


# build empty pairwiseAlignments for object initialization

#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom pwalign pairwiseAlignment
empty.pass <- function() {
  empty_pattern <- DNAStringSet(character(0))
  subject <- DNAString("ACGT")
  pairwiseAlignment(pattern = empty_pattern, subject = subject)
}


# trim long strings to 20 characters, used for show() functions

charDisplayTrim <- function(string) {
  str_trunc(string, width = 23, side = "right")
}

#function to determine optimal sequence for each umi and filter non-passing umis as "Rejected"

#' @import data.table
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








