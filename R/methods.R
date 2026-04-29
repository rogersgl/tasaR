###############################################################################
# Access methods
# ------

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
#' @docType methods
#' @slot x object to check for uninitialized values
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


### class PairwiseAlignmentsSingleSubject ###

setMethod("isEmpty", "PairwiseAlignmentsSingleSubject", function(x) {
  o <- utils::capture.output(x)
  if (length(o) == 1 && stringr::str_detect(o, "Empty")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})

# ------
# show()
# ------


### class tas.object.settings ###

#' @export
setMethod("show", "tas.object.settings", function(object) {
  cat("Settings for ", object@Name, ":\n",
      "IsAntibody : ", object@IsAntibody, "\n")
  for (i in methods::slotNames(object)[3:10]) {
    s <- methods::slot(object, i)
    if (s != ""){
      s <- charDisplayTrim(s)
      cat("",i,": ", s, "\n", sep = "")}
  }
  for (i in methods::slotNames(object)[11:13]) {
    if (!is.na(methods::slot(object, i))){cat("",i,": ", methods::slot(object, i), "\n", sep = "")}
  }
  if (!all(is.na(object@AntibodyRegions))){
    cat("Antibody region coordinates:\n",
        " FR1: ", object@AntibodyRegions["FR1Start"], "\n",
        "CDR1: ", object@AntibodyRegions["CDR1Start"], "\n",
        " FR2: ", object@AntibodyRegions["FR2Start"], "\n",
        "CDR2: ", object@AntibodyRegions["CDR2Start"], "\n",
        " FR3: ", object@AntibodyRegions["FR3Start"], "\n",
        "CDR3: ", object@AntibodyRegions["CDR3Start"], "\n",
        " FR4: ", object@AntibodyRegions["FR4Start"], "\n",
        " End: ", object@AntibodyRegions["FR4End"], "\n", sep = "")
  }
})


### class tas.mutations ###

#' @export
setMethod("show", "tas.mutations", function(object) {
  cat("DNA mutation tables:\n\n")
  for (i in names(object@DNA)[1:3]) {
    df <- object@DNA[[i]]
    if (!isEmpty(df)) {
      cat("",i,":\n", sep = "")
      if (nrow(df) <= 6){
        show(df)
        cat("\n-------------------------\n")
      } else {
        show(df[1:3,])
        cat("\n          ..........          \n\n")
        show(df[(nrow(df)-2):nrow(df),])
        cat("\n-------------------------\n")
      }
    }
    cat("\n")
  }
  cat("Percent mutation for different motifs:\n")
  for (i in names(object@DNA$MotifSums)) {
    if (!all(is.na(object@DNA$MotifSums[i]))) {
      cat("",i,": ", object@DNA$MotifSums[i],"%\n", sep = "")
    }
  }
  cat("\nProtein mutation tables:\n")
  if (!isEmpty(object@AA$AllMutations)) {
    df <- object@AA$AllMutations
    if (nrow(df) <= 6){
      show(df)
      cat("\n-------------------------\n")
    } else {
      show(df[1:3,])
      cat("\n          ..........          \n\n")
      show(df[(nrow(df)-2):nrow(df),])
      cat("\n-------------------------\n")
    }
  }
  cat("\n")

  if (!all(is.na(object@AA$MutationMatrix))) {
    m <- object@AA$MutationMatrix
    cat("MutationMatrix: maxrix of protein mutations for ", ncol(m), " amino acid positions.\n", sep = "")
  }
  cat("\nAID summary tables:\n")
  for (df in object@AIDTables) {
    if (!isEmpty(df)) {
      if (nrow(df) <= 6){
        show(df)
        cat("\n-------------------------\n")
      } else {
        show(df[1:3,])
        cat("\n          ..........          \n\n")
        show(df[(nrow(df)-2):nrow(df),])
        cat("\n-------------------------\n")
      }
    }
  }

})


### class tas.sequences ###

#' @export
setMethod("show", "tas.sequences", function(object) {
  df <- object@Table
  if (nrow(df) <= 6){
    for (i in colnames(df)) {
      if (class(df[,i]) == "character"){
        df[,i] <- sapply(df[,i], charDisplayTrim, USE.NAMES = FALSE)
      }
    }
    show(df)
  } else {
    df <- rbind(df[1:3,], df[(nrow(df)-2):nrow(df),])
    for (i in colnames(df)) {
      if (class(df[,i]) == "character"){
        df[,i] <- sapply(df[,i], charDisplayTrim, USE.NAMES = FALSE)
      }
    }
    show(df[1:3,])
    cat("----------\n")
    show(df[(nrow(df)-2):nrow(df),])
  }
})


### class tas.dna.repair ###

#' @export
setMethod("show", "tas.dna.repair", function(object) {
  cat("Inferred DNA repair pathways:\n")
  for (i in methods::slotNames(object)) {
    if (!is.na(methods::slot(object, i))){
      cat(i,": ", methods::slot(object, i), "%\n", sep = "")
    }
  }
})


### class AmpliconSequencing ###

#' @export
setMethod("show", "AmpliconSequencing", function(object) {
  for (i in methods::slotNames(object)) {
    cat(rep("-", nchar(i)+1),
        "\n", i,":\n",
        rep("-", nchar(i)+1),"\n\n", sep = "")
    o <- methods::slot(object, i)
    if (isEmpty(o)) {
      cat("Empty\n\n")
    } else {
      show(o)
    }
    cat("\n\n\n")
  }
})



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
#' @aliases as.numeric,tasaR-method
#' @export
#' @method as.list tas.object.settings
as.list.tas.object.settings <- function(x) {
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
#' @aliases as.numeric,tasaR-method
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


