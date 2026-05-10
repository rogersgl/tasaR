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

# need a builder for unpaired reads to make this test work
# Need to implement, important test for the package
# test_that("pandaseq_merge_files runs without crashing", {
#   skip("Add small test FASTQ files to inst/extdata first")
#
#   forward <- system.file("extdata", "test_R1.fastq.gz", package = "tasaR")
#   reverse <- system.file("extdata", "test_R2.fastq.gz", package = "tasaR")
#   output  <- tempfile(fileext = ".fastq")
#
#   expect_true(
#     pandaseq_merge_files(forward, reverse, output)
#   )
#
#   expect_true(file.exists(output))
# })
