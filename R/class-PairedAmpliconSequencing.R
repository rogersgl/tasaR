#' @include class-AmpliconSequencing.R
NULL

# ------------------------
# PairedAmpliconSequencing
# ------------------------

#' S4 class PairedAmpliconSequencing
#'
#' @description
#' S4 class that holds results of paired tasaR analysis. Incorporates a control and
#' an experimental condition for comparison.
#'
#'
#' @slot Control S4 object of class AmpliconSequencing
#' @slot Experimental S4 object of class AmpliconSequencing
#'
#' @export
setClass("PairedAmpliconSequencing", slots = list(Control = "AmpliconSequencing",
                                                  Experimental = "AmpliconSequencing"),
         prototype = list(Control = methods::new("AmpliconSequencing"),
                          Experimental = methods::new("AmpliconSequencing")))

setValidity("PairedAmpliconSequencing", function(object) {
  if (length(slotNames(object)) != 2) {
    return("Incorrect number of slots.")
  }
  for (i in slotNames(object)) {
    if (!validObject(slot(object, i))) {
      return(stringr::str_c("Invalid slot ",i," detected."))
    }
  }
  return(TRUE)
})

