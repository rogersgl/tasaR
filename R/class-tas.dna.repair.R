#' @include helper-fns.R
NULL

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
