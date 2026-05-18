#' @include class-AmpliconSequencing.R
NULL


###############################################################################
# Coercion methods
# ------

# ------------
# as.numeric()
# ------------

#' Coerce tas.dna.repair to numeric vector
#'
#' @param x A tas.dna.repair object
#' @aliases as.numeric,tasaR-method
#' @export
setMethod("as.numeric", signature = "tas.dna.repair", function(x) {
  stats::setNames(c(x@WT, x@NHEJ, x@MMEJ, x@BaseChange, x@IndelBaseChange, x@Other),
           c("WT", "NHEJ", "MMEJ", "BaseChange", "IndelBaseChange", "Other"))
})

# ---------
# as.list()
# ---------

### class tas.object.settings ###

#' Coerce tas.object.settings to list
#'
#' @param x A tas.object.settings object
#' @param ... Additional parameters
#' @aliases as.list,tasaR-method
#' @export
#' @method as.list tas.object.settings
as.list.tas.object.settings <- function(x, ...) {
  list(Name = x@Name,
       IsAntibody = x@IsAntibody,
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



# ---------------
# as.data.frame()
# ---------------

### class tas.dna.repair ###

#' Coerce tas.dna.repair to data.frame
#'
#' @param x A tas.dna.repair object
#' @param ... Passed to data.frame
#' @param row.names (Optional) character vector to give row names ot data.frame. Default is NULL.
#' @aliases as.data.frame,tasaR-method
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


