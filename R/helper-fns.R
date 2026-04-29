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


empty.pass <- function() {
  empty_pattern <- Biostrings::DNAStringSet(character(0))
  subject <- Biostrings::DNAString("ACGT")
  pwalign::pairwiseAlignment(pattern = empty_pattern, subject = subject)
}


# trim long strings to 20 characters, used for show() functions

charDisplayTrim <- function(string) {
  stringr::str_trunc(string, width = 23, side = "right")
}

#function to determine optimal sequence for each umi and filter non-passing umis as "Rejected"

umi_pileup <- function(x){
  t <- table(x)
  d <- data.table(seq = names(t),count=as.numeric(t))
  setorder(d,-count)
  if((length(d$count)==1) || (d$count[1]>=(d$count[2]+2) && (d$count[2]<d$count[1]*0.7))){
    return(d$seq[1])
  }else{
    return("Rejected")
  }
}















###############################################################################
# --------------
# Constants
# --------------

.pt <- 72.27 / 25.4  # points per mm






