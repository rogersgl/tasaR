test_that("Nuclease alignment correction handles cases as expected", {

  cut.site <- 26

  # GGACGACGGCAACTACAAGACCCGCG[|]CCGAGGTGAAGTTCGAGGGCGACAC #

  ### right mismatch ###
  input.seqs <- c("GGACGACGGCAACTACAAGACCCGCG-CCGAGGTGAAGTTCGAGGGCGACAC",
                  "GGACGACGGCAACTACAAGACCCGCG-C-GAGGTGAAGTTCGAGGGCGACAC")
  input.seqs <- setNames(input.seqs, c("Reference", "-2 – 100%"))
  expect_output(aln <- msa::msaClustalW(Biostrings::DNAStringSet(input.seqs, use.names = TRUE), order = "input"))
  out <- .nucleaseAlignmentCorrection(aln, cut.site, 1L)
  expect_equal(unname(as.character(out)[2]),
               "GGACGACGGCAACTACAAGACCCGCG--CGAGGTGAAGTTCGAGGGCGACAC")

  ### left mismatch ###
  input.seqs <- c("GGACGACGGCAACTACAAGACCCGCGACGAGGTGAAGTTCGAGGGCGACAC",
                  "GGACGACGGCAACTACAAGACC--CGACGAGGTGAAGTTCGAGGGCGACAC")
  input.seqs <- setNames(input.seqs, c("Reference", "-2 – 100%"))
  expect_output(aln <- msa::msaClustalW(Biostrings::DNAStringSet(input.seqs, use.names = TRUE), order = "input"))
  out <- .nucleaseAlignmentCorrection(aln, cut.site, 1L)
  expect_equal(unname(as.character(out)[2]),
               "GGACGACGGCAACTACAAGACCCG--ACGAGGTGAAGTTCGAGGGCGACAC")

  ### 3 nt microhomology ###
  input.seqs <- c("GGACGACGGCAACTACAAGACCCGCGAGGAGGTGAAGTTCGAGGGCGACAC",
                  "GGACGACGGCAACTACAAGACCCGCGAGG---TGAAGTTCGAGGGCGACAC")
  input.seqs <- setNames(input.seqs, c("Reference", "-3 – 100%"))
  expect_output(aln <- msa::msaClustalW(Biostrings::DNAStringSet(input.seqs, use.names = TRUE), order = "input"))
  out <- .nucleaseAlignmentCorrection(aln, cut.site, 1L)
  expect_equal(unname(as.character(out)[2]),
               "GGACGACGGCAACTACAAGACCCGCG---AGGTGAAGTTCGAGGGCGACAC")

  ### sequence with both an insertion and a deletion to its right ###
  input.seqs <- c("GGACGACGGCAACTACAAGACCCGCG--CCGAGGTGAAGTTCGAGGGCGACAC",
                  "GGACGACGGCAACTACAAGACCCGCGATC-GAGGTGAAGTTCGAGGGCGAC--")
  input.seqs <- setNames(input.seqs, c("Reference", "+2, -1 – 100%"))
  expect_output(aln <- msa::msaClustalW(Biostrings::DNAStringSet(input.seqs, use.names = TRUE), order = "input"))
  out <- .nucleaseAlignmentCorrection(aln, cut.site, reference = 1L)
  expect_equal(unname(as.character(out)[2]),
               "GGACGACGGCAACTACAAGACCCGCGAT-CGAGGTGAAGTTCGAGGGCGACAC")

  ### sequence with both an insertion and a deletion to its left ###
  input.seqs <- c("GGACGACGGCAACTACAAGACCCGCG--CCGAGGTGAAGTTCGAGGGCGACAC",
                  "GGACGACGGCAACTACAAGACC--CGATCCGAGGTGAAGTTCGAGGGCGAC--")
  input.seqs <- setNames(input.seqs, c("Reference", "+2, -1 – 100%"))
  expect_output(aln <- msa::msaClustalW(Biostrings::DNAStringSet(input.seqs, use.names = TRUE), order = "input"))
  out <- .nucleaseAlignmentCorrection(aln, cut.site, reference = 1L)
  expect_equal(unname(as.character(out)[2]),
               "GGACGACGGCAACTACAAGACCCG--ATCCGAGGTGAAGTTCGAGGGCGACAC")

  ### sequence with multiple gaps ###
  input.seqs <- c("GGACGACGGCAACTACAAGACCCGCG-CCGAGGTGAAGTTCGAGGGCGACAC",
                  "GGACGACGGC-ACTACAAGACCCGCG-C-GAGGTGAAGTTCGAGGGCGACAC")
  input.seqs <- setNames(input.seqs, c("Reference", "-1, -1 – 100%"))
  expect_output(aln <- msa::msaClustalW(Biostrings::DNAStringSet(input.seqs, use.names = TRUE), order = "input"))
  out <- .nucleaseAlignmentCorrection(aln, cut.site, reference = 1L)
  expect_equal(unname(as.character(out)[2]),
               "GGACGACGGC-ACTACAAGACCCGCG--CGAGGTGAAGTTCGAGGGCGACAC")

  ### sequence with mutiple gaps on the same side of the cut site ###
  input.seqs <- c("GGACGACGGCAACTACAAGACCCGCG-CCGAGGTGAAGTTCGAGGGCGACAC",
                  "GGACGACGGCAACTACAAGACCCGCG-C-GAGGTGAAGTTC-AGGGCGACAC")
  input.seqs <- setNames(input.seqs, c("Reference", "-1, -1 – 100%"))
  expect_output(aln <- msa::msaClustalW(Biostrings::DNAStringSet(input.seqs, use.names = TRUE), order = "input"))
  out <- .nucleaseAlignmentCorrection(aln, cut.site, reference = 1L)
  expect_equal(unname(as.character(out)[2]),
               "GGACGACGGCAACTACAAGACCCGCG--CGAGGTGAAGTTC-AGGGCGACAC")

  ### an ambiguously placed alignment with multiple solutions after correction ###
  # moves next to gap but does not give result spanning the gap
  # either could be correct, so don't think this needs to be changed
  # just note that software favors 1-sided deletions in output alignments
  input.seqs <- c("GGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACAC",
                  "GGACGACGGCAACTACAAGACC--CGCCGAGGTGAAGTTCGAGGGCGACAC"
                  )
  input.seqs <- setNames(input.seqs, c("Reference", "-2 – 100%"))
  expect_output(aln <- msa::msaClustalW(Biostrings::DNAStringSet(input.seqs, use.names = TRUE), order = "input"))
  out <- .nucleaseAlignmentCorrection(aln, cut.site, reference = 1L)
  expect_equal(unname(as.character(out)[2]),
               "GGACGACGGCAACTACAAGACCCG--CCGAGGTGAAGTTCGAGGGCGACAC")


  ### adjustment with multiple cut sites ###
  # TODO: maybe add support for multiple cut sites (eg Cas12a or ZFN/TALENs)
  # currently shelving because looks like it will be complex for an edge case


})

test_that("analysis progress callback errors are non-fatal", {
  callback <- function(event) {
    stop("callback failed")
  }

  expect_no_error(.emit_analysis_progress(callback, "sample_1", 1L, 5L, "Reading FASTQ"))
})

test_that("analysis progress cancellation conditions are rethrown", {
  callback <- function(event) {
    stop(structure(
      list(message = "Mutation analysis canceled.", call = NULL),
      class = c("tasaR_analysis_canceled", "error", "condition")
    ))
  }

  expect_error(
    .emit_analysis_progress(callback, "sample_1", 1L, 5L, "Reading FASTQ"),
    class = "tasaR_analysis_canceled"
  )
})
