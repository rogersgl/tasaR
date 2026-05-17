test_that(".validate_pandaseq_paths fails clearly on missing input files", {
  forward <- tempfile("missing-forward-", fileext = ".fastq.gz")
  reverse <- tempfile("missing-reverse-", fileext = ".fastq.gz")
  output  <- tempfile("merged-", fileext = ".fastq")

  expect_error(
    .validate_pandaseq_paths(
      forward_fastq = forward,
      reverse_fastq = reverse,
      output_fastq = output
    ),
    regexp = "Forward FASTQ not found",
    fixed = FALSE
  )
})

test_that(".validate_pandaseq_paths creates output directory when requested", {
  td <- tempfile("pandaseq-outdir-")
  forward <- tempfile("forward-", fileext = ".fastq.gz")
  reverse <- tempfile("reverse-", fileext = ".fastq.gz")
  output  <- file.path(td, "nested", "merged.fastq")

  dir.create(dirname(forward), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(reverse), recursive = TRUE, showWarnings = FALSE)
  writeLines(character(), forward)
  writeLines(character(), reverse)

  out <- .validate_pandaseq_paths(
    forward_fastq = forward,
    reverse_fastq = reverse,
    output_fastq = output,
    create_output_dir = TRUE
  )

  expect_true(dir.exists(dirname(out$output_fastq)))
})


test_that("pandaseq_merge_files runs and merges simple test files", {

  p <- pandaseq_merge_files(forward_fastq = file.path(tempdir(), "R1-test.fastq.gz"),
                         reverse_fastq = file.path(tempdir(), "R2-test.fastq.gz"),
                         output_fastq = file.path(tempdir(), "merged-test.fastq.gz"),
                         log_file = file.path(tempdir(), "log/test-log.txt"),
                         min_length = 370,
                         max_length = 420,
                         extra_args = c("-F", "-d", "bFSrk"),
                         verbose = FALSE)
  expect_true(p$ok)

  expect_true(file.exists(file.path(tempdir(), "merged-test.fastq.gz")))

  expect_equal(as.character(ShortRead::sread(ShortRead::readFastq(file.path(tempdir(), "merged-test.fastq.gz")))),
               rep("GCTAGCCGTAAAACGACGGCCAGTGTTCAACTGGTGGAAAGCGGCGGTGCTCTGGTACAACCGGGCGGTAGTCTGCGCCTGAGCTGTGCCGCAAGCGGTTTCCCAGTCAACCGCTACTCTATGCGTTGGTATCGCCAGGCGCCTGGTAAAGAACGTGAATGGGTTGCCGGCATGAGCAGTGCGGGCGATCGTTCTAGTTACGAGGACTCTGTTAAAGGTCGTTTTACAATTAGCCGTGATGATGCGCGCAATACCGTGTATCTGCAAATGAACAGTCTGAAGCCGGAGGACACCGCAGTATATTATTGCAATGTCAACGTGGGGTTTGAATATTGGGGCCAGGGGACTCAGGTGACGGTGAGCTCTGTCATAGCTGTTTCCTGAGCGAATTAGAG", 3
                   )
               )

  # TODO: add failure tests

})
