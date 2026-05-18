#' @include class-AmpliconSequencing.R
#' @include class-PairedAmpliconSequencing.R
NULL

# ------
# show()
# ------

#' @title show methods for package tasaR
#' @description
#' Extends the show function to print to console contents of S4 objects in the tasaR package.
#'
#' @name show-tasaR
#' @aliases show,tasaR-method
#' @param object S4 object to print to console
#' @docType methods
NULL


### class tas.object.settings ###

#' @rdname show-tasaR
#' @export
setMethod("show", "tas.object.settings", function(object) {
  cat("Settings for ", object@Name, ":\n",
      "IsAntibody : ", object@IsAntibody, "\n")
  for (i in methods::slotNames(object)[3:10]) {
    s <- methods::slot(object, i)
    if (any(s != "") && !is(s, "list")){
      s <- .charDisplayTrim(s)
      cat("",i,": ", s, "\n", sep = "")
    } else if (is(s, "list")) {
      for (j in 1:length(s)) {
        cat("",i," [", j, "]: ", s[[j]], "\n", sep = "")
      }
    }
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

#' @rdname show-tasaR
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

#' @rdname show-tasaR
#' @export
setMethod("show", "tas.sequences", function(object) {
  df <- object@Table
  if (nrow(df) <= 6){
    for (i in colnames(df)) {
      if (is(df[,i], "character")){
        df[,i] <- sapply(df[,i], .charDisplayTrim, USE.NAMES = FALSE)
      }
    }
    show(df)
  } else {
    df <- rbind(df[1:3,], df[(nrow(df)-2):nrow(df),])
    for (i in colnames(df)) {
      if (is(df[,i], "character")){
        df[,i] <- sapply(df[,i], .charDisplayTrim, USE.NAMES = FALSE)
      }
    }
    show(df[1:3,])
    cat("----------\n")
    show(df[(nrow(df)-2):nrow(df),])
  }
})


### class tas.dna.repair ###

#' @rdname show-tasaR
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

#' @rdname show-tasaR
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

### class PairedAmpliconSequencing ###

#' @rdname show-tasaR
#' @export
setMethod("show", "PairedAmpliconSequencing", function(object) {
  for (i in methods::slotNames(object)) {
    show(slot(object, i))
  }
})
