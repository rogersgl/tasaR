
#' Merge paired end reads
#'
#' @description
#' A function to merge paired end reads using an included installation of the software PANDAseq.
#' Writes a merged fastq file to disk with the given arguments.
#'
#'
#' @param forward_fastq Character vector, file path to the forward (R1) read
#' @param reverse_fastq Character vector, file path to the reverse (R2) read
#' @param output_fastq Character vector, file path to write the merged read to. Sub-directories will be automatically created.
#' @param log_file Character vector, file path for the log file to be written to. Sub-directories will be automatically created.
#' @param min_length Numeric vector, minimum sequence length to allow
#' @param max_length Numeric vector, maximum sequence length to allow
#' @param extra_args Character vector, additional arguments are set with reasonable defaults. See PANDAseq manual for more details.
#'
#' @returns NULL
#' @export
pandaseq_merge_files <- function(forward_fastq,
                                 reverse_fastq,
                                 output_fastq,
                                 log_file = NULL,
                                 min_length = NULL,
                                 max_length = NULL,
                                 extra_args = c("-F", "-d", "bFSrk"),
                                 verbose = FALSE) {
  paths <- .validate_pandaseq_paths(
    forward_fastq = forward_fastq,
    reverse_fastq = reverse_fastq,
    output_fastq = output_fastq,
    log_file = log_file
  )

  args <- c(
    "pandaseq",
    "-f", paths$forward_fastq,
    "-r", paths$reverse_fastq,
    "-w", paths$output_fastq,
    if (!is.null(paths$log_file)) c("-g", paths$log_file) else NULL,
    if (!is.null(min_length)) c("-l", as.character(min_length)) else NULL,
    if (!is.null(max_length)) c("-L", as.character(max_length)) else NULL,
    extra_args
  )

  res <- .Call("pandaseq_merge_tokens", args, PACKAGE = "tasaR")
  res$args <- args
  class(res) <- c("pandaseq_result", "list")

  if (isTRUE(verbose) && nzchar(res$stderr)) {
    message(res$stderr)
  }

  invisible(res)
}

#' @export
print.pandaseq_result <- function(x, ...) {
  if (isTRUE(x$ok)) {
    cat("<pandaseq_result> OK:", x$output_fastq, "\n")
  } else {
    cat("<pandaseq_result> FAIL phase:", x$phase, "\n")
    if (nzchar(x$stderr)) cat(x$stderr, "\n")
  }
  invisible(x)
}


