#' @include class-AmpliconSequencing.R
NULL

# ---------
# isEmpty()
# ---------

### main definition page ###

#' @title isEmpty methods for package tasaR
#' @description
#' Extends the isEmpty function to S4 objects in package tasaR to detect whether S4 classes in tasa are empty, holding only initialized values.
#'
#' @name isEmpty-tasaR
#' @aliases isEmpty,tasaR-method
#' @param x object to check for uninitialized values
#' @docType methods
NULL

### class tas.object.settings ###


#' @rdname isEmpty-tasaR
#' @export
setMethod("isEmpty", "tas.object.settings", function(x) {
  n <- methods::slotNames(x)
  n <- setdiff(n, c("Name", "AntibodyRegions"))
  empty <- logical()
  for (i in n) {
    if (class(methods::slot(x,i)) == "numeric" || class(methods::slot(x,i)) == "integer"){
      empty <- c(empty, is.na(methods::slot(x,i)))
    } else if (class(methods::slot(x,i)) == "character") {
      empty <- c(empty, (methods::slot(x,i) == ""))
    } else if (class(methods::slot(x, i)) == "list") {
      empty <- c(empty, unlist((methods::slot(x,i)) == ""))
    }
  }
  for (i in names(x@AntibodyRegions)) {
    empty <- c(empty, is.na(x@AntibodyRegions[i]))
  }
  if (all(empty)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})


### class tas.mutations ###

#' @rdname isEmpty-tasaR
#' @export
setMethod("isEmpty", "tas.mutations", function(x) {
  dfs <- c(x@DNA$AllMutations, x@DNA$CytosineMutations, x@DNA$NonCytosineMutations, x@AA$AllMutations, x@AIDTables$WRCH, x@AIDTables$WRCY)
  df_empty <- logical()
  for (i in dfs){
    df_empty <- c(df_empty, isEmpty(i))
  }
  ms <- all(is.na((x@DNA$MotifSums)))
  mm <- is.na(x@AA$MutationMatrix)
  if (all(df_empty) && ms && mm){
    return(TRUE)
  } else {
    return(FALSE)
  }
})


### class tas.sequences ###

#' @rdname isEmpty-tasaR
#' @export
setMethod("isEmpty", "tas.sequences", function(x) {
  df <- x@Table
  if (nrow(df) < 2 && all(is.na(df[1,c(2,3,5)])) && all(df[1,c(1,4,6,7)] == "")){
    return(TRUE)
  } else {
    return(FALSE)
  }
})


### class tas.dna.repair ###

#' @rdname isEmpty-tasaR
#' @export
setMethod("isEmpty", "tas.dna.repair", function(x) {
  empty <- logical()
  for (i in methods::slotNames(x)) {
    empty <- c(empty, is.na(methods::slot(x, i)))
  }
  if (all(empty)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})



### class AmpliconSequencing ###

#' @rdname isEmpty-tasaR
#' @export
setMethod("isEmpty", "AmpliconSequencing", function(x) {
  empty <- logical()
  for (i in methods::slotNames(x)) {
    empty <- c(empty, isEmpty(methods::slot(x, i)))
  }
  if (all(empty)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})

### class PairedAmpliconSequencing ###

#' @rdname isEmpty-tasaR
#' @export
setMethod("isEmpty", "PairedAmpliconSequencing", function(x) {
  empty <- logical()
  for (i in methods::slotNames(x)) {
    for (j in methods::slotNames(slot(x, i)))
      empty <- c(empty, isEmpty(methods::slot(methods::slot(x, i), j)))
  }
  if (all(empty)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})


### class PairwiseAlignmentsSingleSubject ###

#' @rdname isEmpty-tasaR
#' @export
setMethod("isEmpty", "PairwiseAlignmentsSingleSubject", function(x) {
  o <- utils::capture.output(x)
  if (length(o) == 1 && stringr::str_detect(o, "Empty")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})


### class DNAMultipleAlignment ###
#' @rdname isEmpty-tasaR
#' @export
setMethod("isEmpty", "DNAMultipleAlignment", function(x) {
  if (nrow(x) == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})
