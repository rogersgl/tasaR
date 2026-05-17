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

  gRNA.seq <- Biostrings::DNAString(gRNA.seq)
  gRNA.coords.c <- suppressWarnings(BiocGenerics::as.data.frame(sapply(ref.seq, function(x){
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

  colnames(gRNA.coords) <- c("group", "group_name", "start", "end", "width")

  return(gRNA.coords$end - 3)

  # if (nuclease.type == "Cas9") {
  #   return(gRNA.coords$end - 3)
  # }

  # if (nuclease.type == "Cas12a") {
  #   return(c(gRNA.coords$start + 17, gRNA.coords$start + 21))
  # }
}


### extract and name sequences for alignment
.alignment.seqs <- function(seqs, refs, seq.indels, seq.basesChanged, seq.percent, gRNA.window, seq.count) {
  seq.count <- min(seq.count, length(seqs))
  ref.seq <- character(length(refs))
  for (i in seq_along(refs)) {
    ref.seq[i] <- stringr::str_sub(refs[[i]], min(gRNA.window), max(gRNA.window))
    ref.seq <- setNames(ref.seq[i], stringr::str_c("Reference [", i, "]"))
  }
  seqs <- sapply(seqs, function(x){
    stringr::str_sub(x, min(gRNA.window), max(gRNA.window))
  }, USE.NAMES = FALSE)
  # seq_tbl <- getSequenceTable(amplicon.seq)[1:seq.count,]
  n <- stringr::str_c(ifelse(nzchar(seq.indels), seq.indels, ""),
                      ifelse(seq.basesChanged != 0,
                             ifelse(nzchar(seq.indels), ", ", "") |>
                               stringr::str_c(seq.basesChanged, " SNV"), ""),
                      stringr::str_c(" – ", round(seq.percent, 2), "%")
  )
  seqs <- setNames(seqs, n)
  return(Biostrings::DNAStringSet(c(ref.seq, seqs)))
}


.alignment.seqs.ampseq <- function(amplicon.seq, gRNA.window, seq.count) {
  seq.count <- min(seq.count, nrow(getSequenceTable(amplicon.seq)))
  refs <- getSettings(amplicon.seq)$ReferenceSequence
  ref.seq <- character(length(refs))
  for (i in seq_along(refs)) {
    ref.seq[i] <- stringr::str_sub(refs[[i]], min(gRNA.window), max(gRNA.window))
    ref.seq <- setNames(ref.seq[i], stringr::str_c("Reference [", i, "]"))
  }
  seqs <- getSequencesDNA(amplicon.seq)[1:seq.count]
  seqs <- sapply(seqs, function(x){
    stringr::str_sub(x, min(gRNA.window), max(gRNA.window))
  }, USE.NAMES = FALSE)
  seq_tbl <- getSequenceTable(amplicon.seq)[1:seq.count,]
  n <- stringr::str_c(ifelse(nzchar(seq_tbl$Indels), seq_tbl$Indels, ""),
                              ifelse(seq_tbl$BasesChanged != 0,
                                 ifelse(nzchar(seq_tbl$Indels), ", ", "") |>
                                     stringr::str_c(seq_tbl$BasesChanged, " SNV"), ""),
                              stringr::str_c(" – ", round(seq_tbl$Percent, 2), "%")
                             )
  seqs <- setNames(seqs, n)
  return(Biostrings::DNAStringSet(c(ref.seq, seqs)))
}

### correct gap placement relative to Cas9 cut site in alignments ###

.nucleaseAlignmentCorrection_find_gap_runs_not_touching_cut <- function(seq, cut_col) {
  # cut_col = aligned column immediately left of the cut
  loc <- stringr::str_locate_all(seq, "-+")[[1]]

  if (nrow(loc) == 0) {
    return(data.frame(start = integer(), end = integer(), length = integer()))
  }

  out <- data.frame(
    start  = loc[, 1],
    end    = loc[, 2],
    length = loc[, 2] - loc[, 1] + 1L
  )

  out$touches_cut <- out$start <= (cut_col + 1L) & out$end >= cut_col

  out[!out$touches_cut, c("start", "end", "length"), drop = FALSE]
}

.nucleaseAlignmentCorrection_ref_pos_to_aln_cut <- function(ref_aln, ref_pos) {
  ref_chars <- strsplit(ref_aln, "", fixed = TRUE)[[1]]
  nongap_cols <- which(ref_chars != "-")

  if (ref_pos < 1L || ref_pos > length(nongap_cols)) {
    stop("ref_pos out of range")
  }

  nongap_cols[ref_pos]
}


# gap_runs: data.frame with start/end/length from find_gap_runs_not_touching_cut()
# cut.site:  ungapped reference position immediately left of the cut
# cut_col:   aligned-column position of that same cut in the reference alignment
#    gap_runs <- .nucleaseAlignmentCorrection_can_cut_align_runs(seqs[i], ref_aln, cut.site, cut_col, runs)
.nucleaseAlignmentCorrection_can_cut_align_runs <- function(seq_aln, ref_aln, cut.site, cut_col, gap_runs) {
  if (nrow(gap_runs) == 0) {
    return(gap_runs)
  }

  ins <- stringr::str_locate_all(ref_aln, "-")[[1]]
  seq_u <- seq_aln
  if (!isEmpty(ins)){
    for (i in 1:nrow(ins)){
      seq_u <- stringr::`str_sub<-`(seq_u, start = ins[i,1], end = ins[i,2], value = "-")
    }
  }
  ref_u <- stringr::str_remove_all(ref_aln, "-")
  # seq_u <- stringr::str_remove_all(seq_u, "\\*")
  n <- nchar(ref_u)

  gap_runs$side <- ifelse(gap_runs$end < cut_col, "left", "right")

  gap_runs$candidate_start <- ifelse(
    gap_runs$side == "right",
    cut.site + 1L,                      # delete bases immediately right of cut
    cut.site - gap_runs$length + 1L     # delete bases immediately left of cut
  )

  gap_runs$can_align <- apply(gap_runs, 1, function(gap_rows) {
    k <- as.numeric(gap_rows["length"])
    pos <- as.numeric(gap_rows["candidate_start"])
    start <- as.numeric(gap_rows["start"])
    end <- as.numeric(gap_rows["end"])
    len <- as.numeric(gap_rows["length"])

    if (pos < 1L || (pos + k - 1L) > n) return(FALSE)

    candidate <- paste0(
      substr(ref_u, 1L, pos - 1L),
      substr(ref_u, pos + k, n)
    )

    # for each gap_run, generate a reference sequence that ignores the other gaps
    # (fills in with reference sequence) to check each one individually
    if (nrow(gap_runs) > 1){
      seq_u_tmp <- stringr::`str_sub<-`(seq_u,
                                        start = gap_runs$start[gap_runs$start != start],
                                        end = gap_runs$end[gap_runs$end != end],
                                        value = stringr::str_sub(ref_aln,
                                                                 start = gap_runs$start[gap_runs$start != start],
                                                                 end = gap_runs$end[gap_runs$end != end]))
    } else {
      seq_u_tmp <- unname(seq_u)
    }
    seq_u_tmp <- stringr::str_remove_all(seq_u_tmp, "-")
    identical(seq_u_tmp, candidate)
  })

  gap_runs
}



  # gap_runs$can_align <- mapply(function(k, pos) {
  #   if (pos < 1L || (pos + k - 1L) > n) return(FALSE)
  #
  #   candidate <- paste0(
  #     substr(ref_u, 1L, pos - 1L),
  #     substr(ref_u, pos + k, n)
  #   )
  #
  #   # goal: for each row in gap_runs, remove any other gaps from other rows and replace with the sequence from ref_aln
  #   # hypothesis: will fix calling of can_align that is messed up by another gap and will always fail with original implementation
  #   seq_u_tmp <- stringr::`str_sub<-`(seq_u,
  #                                     start = gap_runs$start[gap_runs$start != pos],
  #                                     end = gap_runs$end[gap_runs$end != (pos + k)],
  #                                     value = stringr::str_sub(ref_u,
  #                                                              start = gap_runs$start[gap_runs$start != pos],
  #                                                              end = gap_runs$end[gap_runs$end != (pos + k)]))
  #
  #   identical(unname(seq_u), candidate)



.nucleaseAlignmentCorrection_ref_pos_to_aln_cols <- function(ref_aln, ref_pos, k) {
  ref_chars <- strsplit(ref_aln, "", fixed = TRUE)[[1]]
  ref_cols <- which(ref_chars != "-")

  if (ref_pos < 1L || (ref_pos + k - 1L) > length(ref_cols)) {
    stop("candidate_start out of range")
  }

  ref_cols[ref_pos:(ref_pos + k - 1L)]
}

.nucleaseAlignmentCorrection_apply_cut_align_runs <- function(seq_aln, ref_aln, gap_runs) {
  s <- strsplit(seq_aln, "", fixed = TRUE)[[1]]

  if (S4Vectors::isEmpty(gap_runs)) return(seq_aln)
  keep <- which(gap_runs$can_align)
  if (length(keep) == 0L) return(seq_aln)

  # process longer runs first, then nearer-to-cut runs if you want
  keep <- keep[order(-gap_runs$length[keep])]

  for (i in keep) {
    old_cols <- gap_runs$start[i]:gap_runs$end[i]
    k <- gap_runs$length[i]
    new_cols <- .nucleaseAlignmentCorrection_ref_pos_to_aln_cols(ref_aln, gap_runs$candidate_start[i], k)

    # move the block:
    # old gap positions get the sequence that was at the target,
    # target positions become gaps
    tmp <- s[old_cols]
    s[old_cols] <- s[new_cols]
    s[new_cols] <- tmp
  }

  paste0(s, collapse = "")
}

.nucleaseAlignmentCorrection_fillInsertEnd <- function(seqs, idx.ins, end.coords, ref_aln, reference) {
  sapply(1:nrow(end.coords), function(x) {
    stringr::`str_sub<-`(seqs[idx.ins[x]],
                         start = end.coords[x,"start"],
                         end = end.coords[x,"end"],
                         value = stringr::str_sub(seqs[reference],
                                                  start = end.coords[x,"start"],
                                                  end = end.coords[x,"end"]
                         )
    )
  })
}

.nucleaseAlignmentCorrection <- function(aln, cut.site, reference) {
  seqs <- as.character(aln)
  ref_aln <- seqs[reference]
  cut_col <- .nucleaseAlignmentCorrection_ref_pos_to_aln_cols(ref_aln, cut.site, 1L)[1]

  # fill in insertion ends
  idx.ins <- which(stringr::str_detect(BiocGenerics::rownames(aln), "\\+"))
  end.coords <- stringr::str_locate(seqs[idx.ins], "-+$")
  if (!isEmpty(idx.ins)){
    seqs[idx.ins] <- .nucleaseAlignmentCorrection_fillInsertEnd(seqs, idx.ins, end.coords, ref_aln, reference)
  }



  idx.del <- which(stringr::str_detect(BiocGenerics::rownames(aln), "-"))

  for (i in idx.del) {
    runs <- .nucleaseAlignmentCorrection_find_gap_runs_not_touching_cut(seqs[i], cut_col)
    gap_runs <- .nucleaseAlignmentCorrection_can_cut_align_runs(seqs[i], ref_aln, cut.site, cut_col, runs)
    seqs[i] <- .nucleaseAlignmentCorrection_apply_cut_align_runs(seqs[i], ref_aln, gap_runs)
  }

  Biostrings::DNAMultipleAlignment(seqs, use.names = TRUE)
}


# .nucleaseAlignmentCorrection(aln, cut.sites, reference)
# this is working for the test case
# TODO: build a robust test suite to thoroughly evaluate this workflow
# parameters to test:
#   - fix an alignment issue on the left side of the sequence?
#   - adjustment of > 1 nt (microhomology 2 nt, 3 nt, 4 nt?)
#   - can it handle multiple gaps correctly?
#   - can it handle multiple cut sites e.g. Cas12a?
#   - sequence with both ins and del?


### check if LaTeX is installed ###

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






