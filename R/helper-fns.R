###############################################################################
# Helper functions
# ------

# ---------------
# Hidden funcions
# ---------------


# build empty pairwiseAlignments for object initialization

.empty.pass <- function() {
  empty_pattern <- Biostrings::DNAStringSet(character(0))
  subject <- Biostrings::DNAString("ACGT")
  pwalign::pairwiseAlignment(pattern = empty_pattern, subject = subject)
}


# trim long strings to 20 characters, used for show() functions

.charDisplayTrim <- function(string) {
  stringr::str_trunc(string, width = 23, side = "right")
}


#function to determine optimal sequence for each umi and filter non-passing umis as "Rejected"

.umi_pileup <- function(x){
  t <- table(x)
  d <- data.table(seq = names(t),count=as.numeric(t))
  setorder(d,-count)
  if((length(d$count)==1) || (d$count[1]>=(d$count[2]+2) && (d$count[2]<d$count[1]*0.7))){
    return(d$seq[1])
  }else{
    return("Rejected")
  }
}


# adjust mutation frequencies of separate reference alignments based on their proportions and sum corrected %'s

.mutationFrequencyMerge <- function(mf, id, ref.idx) {
  sum(mf*length(which(id == ref.idx))/length(ref.idx))
}


# find cut sites from a gRNA sequence

.findCutSitesFromSequence <- function(gRNA.seq, ref.seq) {

  gRNA.seq <- DNAString(gRNA.seq)
  gRNA.coords.c <- suppressWarnings(as.data.frame(sapply(ref.seq, function(x){
    Biostrings::vmatchPattern(gRNA.seq, x)
  })))

  gRNA.rc <- Biostrings::reverseComplement(gRNA.seq)
  gRNA.coords.rc <- suppressWarnings(as.data.frame(sapply(ref.seq, function(x){
    Biostrings::vmatchPattern(gRNA.rc, x)
  })))

  w <- which(!sapply(list(gRNA.coords.c, gRNA.coords.rc), isEmpty))
  if (w == 1) {
    gRNA.coords <- gRNA.coords.c
  } else if (w == 2) {
    gRNA.coords <- gRNA.coords.rc
    gRNA.coords$start <- gRNA.coords.rc$end
    gRNA.coords$end <- gRNA.coords.rc$start
  } else if (length(w) == 2){
    stop("gRNA matches found in both strands.")
  }

  if (isEmpty(gRNA.coords)){stop("Specified gRNA sequence not found in reference sequence.")}
  if (nrow(gRNA.coords) > length(ref.seq)){stop("Multiple gRNA matches found in reference sequence.")}
  if (nrow(gRNA.coords) < length(ref.seq)){stop("Did not find gRNA match in all reference sequences.")}

  if (nuclease.type == "Cas9") {
    return(gRNA.coords$End - 3.5)
  }

  if (nuclease.type == "Cas12a") {
    return(c(gRNA.coords$start + 18.5, gRNA.coords$start + 22.5))
  }
}


### check if latex is installed ###

.isTinytexOk <- function(required_pkgs = character()) {

  # Best case: actual TinyTeX
  if (isTRUE(tinytex::is_tinytex())){return(TRUE)}

  # Need both a TeX engine and tlmgr for TinyTeX-style recovery
  pdflatex <- Sys.which("pdflatex")
  tlmgr    <- Sys.which("tlmgr")
  if (!nzchar(pdflatex) || !nzchar(tlmgr)){
    stop("No system-level LaTeX instalation detected. Please manually run the function tinytex::install_tinytex() to enable writing alignments to pdf.")
  }

  # Confirm this is TeX Live, not just any pdflatex wrapper
  tlver <- tinytex::tlmgr_version()
  if (!any(sapply(tlver, function(x){grepl("TeX Live", x)}))){
    stop("TeX Live installation not detected, unable to ensure compatability for writing alignments to pdf.")
  }

  # Check if the tlmgr directory has write permission and if not that the required packages are already installed
  if (file.access("/Library/TeX/texbin/pdflatex", mode = 2) != 0 && length(required_pkgs) > 0) {
    ok <- tryCatch(
      all(tinytex::check_installed(required_pkgs)),
      error = function(e) FALSE
    )
    if (!ok){
      stop(stringr::str_c("LaTeX installation detected in a directory without write permission,
                          and the follwing required packages are missing: ",
                          str_flatten(required_pkgs, collapse = ", ", "\n")))
      }
  }
  TRUE
}



# check in R whether paths provided to pandaseq are viable before calling the C function
.validate_pandaseq_paths <- function(forward_fastq,
                                     reverse_fastq,
                                     output_fastq,
                                     log_file = NULL,
                                     create_output_dir = TRUE,
                                     create_log_dir = TRUE) {
  stopifnot(length(forward_fastq) == 1L)
  stopifnot(length(reverse_fastq) == 1L)
  stopifnot(length(output_fastq) == 1L)

  forward_fastq <- normalizePath(forward_fastq, mustWork = FALSE)
  reverse_fastq  <- normalizePath(reverse_fastq, mustWork = FALSE)
  output_fastq   <- normalizePath(output_fastq, mustWork = FALSE)

  if (!file.exists(forward_fastq)) {
    stop("Forward FASTQ not found: ", forward_fastq, call. = FALSE)
  }
  if (!file.exists(reverse_fastq)) {
    stop("Reverse FASTQ not found: ", reverse_fastq, call. = FALSE)
  }

  output_dir <- dirname(output_fastq)
  if (!dir.exists(output_dir)) {
    if (isTRUE(create_output_dir)) {
      ok <- dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      if (!ok && !dir.exists(output_dir)) {
        stop("Output directory could not be created: ", output_dir, call. = FALSE)
      }
    } else {
      stop("Output directory does not exist: ", output_dir, call. = FALSE)
    }
  }

  if (file.access(output_dir, 2) != 0) {
    stop("Output directory is not writable: ", output_dir, call. = FALSE)
  }

  if (!is.null(log_file)) {
    stopifnot(length(log_file) == 1L)
    log_file <- normalizePath(log_file, mustWork = FALSE)
    log_dir <- dirname(log_file)

    if (!dir.exists(log_dir)) {
      if (isTRUE(create_log_dir)) {
        ok <- dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
        if (!ok && !dir.exists(log_dir)) {
          stop("Log directory could not be created: ", log_dir, call. = FALSE)
        }
      } else {
        stop("Log directory does not exist: ", log_dir, call. = FALSE)
      }
    }

    if (file.access(log_dir, 2) != 0) {
      stop("Log directory is not writable: ", log_dir, call. = FALSE)
    }
  }

  invisible(list(
    forward_fastq = forward_fastq,
    reverse_fastq = reverse_fastq,
    output_fastq = output_fastq,
    log_file = log_file
  ))
}


# resolve a set of DNA bases into an IUPAC ambiguity code


# make_ambiguity_code <- function(x) {
#   if (length(x) == 1) {return(x)}
#   b <- stringr::str_flatten(sort(unlist(x, use.names = FALSE)))
#   switch(b,
#          "AC" = "M",
#          "AG" = "R",
#          "AT" = "W",
#          "CT" = "Y",
#          "CG" = "S",
#          "GT" = "K",
#          "ACG" = "V",
#          "ACT" = "H",
#          "AGT" = "D",
#          "CGT" = "B",
#          "ACGT" = "N"
#          )
# }
# M = A/C
# R = A/G
# W = A/T
# Y = C/T
# S = C/G
# K = G/T

# V = A/C/G
# H = A/C/T
# D = A/G/T
# B = C/G/T









###############################################################################
# --------------
# Constants
# --------------

.pt <- 72.27 / 25.4  # points per mm






