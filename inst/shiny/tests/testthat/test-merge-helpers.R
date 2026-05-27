test_that("merge manifest validation requires expected columns", {
  valid <- data.frame(
    "Sample Name" = "sample_1",
    "Left (R1) File" = "sample_1_R1.fastq.gz",
    "Right (R2) File" = "sample_1_R2.fastq.gz",
    check.names = FALSE
  )
  expect_true(.validate_batch_manifest(valid))

  invalid <- valid[, c("Sample Name", "Left (R1) File"), drop = FALSE]
  expect_error(.validate_batch_manifest(invalid), "Missing required columns")
})

test_that("merge status detects ready samples and duplicates", {
  df <- data.frame(
    "Sample Name" = c("sample_1", "sample_2"),
    "Left (R1) File" = c("sample_1_R1.fastq.gz", "sample_2_R1.fastq.gz"),
    "Right (R2) File" = c("sample_1_R2.fastq.gz", "sample_2_R2.fastq.gz"),
    check.names = FALSE
  )
  expect_equal(.sample_table_with_status(df)$Status, c("Ready", "Ready"))

  duplicated <- df
  duplicated[2, "Sample Name"] <- "sample_1"
  expect_match(.sample_table_with_status(duplicated)$Status[1], "Duplicate name")
})

test_that("single merge input preparation returns ordered pipeline inputs", {
  queued <- data.frame(
    "Queue ID" = c("a", "b"),
    "Sample Name" = c("sample_a", "sample_b"),
    "R1 Path" = c("/tmp/a_R1.fastq.gz", "/tmp/b_R1.fastq.gz"),
    "R2 Path" = c("/tmp/a_R2.fastq.gz", "/tmp/b_R2.fastq.gz"),
    "Status" = c("Ready", "Ready"),
    "Source" = c("single", "single"),
    check.names = FALSE
  )

  out <- .prepare_merge_inputs(mode = "single", queued_df = .order_sample_df(queued, c("b", "a"), "Queue ID"))
  expect_equal(out$Sample_Name, c("sample_b", "sample_a"))
  expect_equal(out$R1_Path, c("/tmp/b_R1.fastq.gz", "/tmp/a_R1.fastq.gz"))
})

test_that("merge job rows assemble queue ids with prepared inputs", {
  ordered_queue <- data.frame(
    "Queue ID" = c("a", "b"),
    "Sample Name" = c("sample_a", "sample_b"),
    "R1 Path" = c("/tmp/a_R1.fastq.gz", "/tmp/b_R1.fastq.gz"),
    "R2 Path" = c("/tmp/a_R2.fastq.gz", "/tmp/b_R2.fastq.gz"),
    "Status" = c("Ready", "Ready"),
    "Source" = c("single", "single"),
    check.names = FALSE
  )
  merge_inputs <- .prepare_merge_inputs(mode = "single", queued_df = ordered_queue)
  settings <- list(merge_prefix = "merge", merge_filetype = "FASTQ")
  output_dir <- file.path(tempdir(), "merge", "session_1")

  expect_no_error({
    job_rows <- .merge_job_rows(
      ordered_queue = ordered_queue,
      merge_inputs = merge_inputs,
      settings = settings,
      output_dir = output_dir
    )
  })

  expect_equal(job_rows[["Queue ID"]], c("a", "b"))
  expect_equal(job_rows[["Sample_Name"]], c("sample_a", "sample_b"))
  expect_equal(
    job_rows[["Output_Path"]],
    file.path(output_dir, c("merge-sample_a.fastq", "merge-sample_b.fastq"))
  )
  expect_equal(
    job_rows[["Log_Path"]],
    file.path(output_dir, c("sample_a-log.txt", "sample_b-log.txt"))
  )
})

test_that("merge output paths use configured prefix and file type", {
  settings <- list(merge_prefix = "merge", merge_filetype = "FASTQ")
  paths <- .merge_output_paths("sample_a", settings, output_dir = "/tmp/merge")

  expect_equal(paths$output_path, "/tmp/merge/merge-sample_a.fastq")
  expect_equal(paths$log_path, "/tmp/merge/sample_a-log.txt")

  settings$merge_prefix <- ""
  settings$merge_filetype <- "FASTA"
  paths <- .merge_output_paths("sample_a", settings, output_dir = "/tmp/merge")
  expect_equal(paths$output_path, "/tmp/merge/sample_a.fasta")
})

test_that("pandaseq extra args use primer values", {
  settings <- list(
    min_overlap = 20,
    merge_similarity = 0.9,
    forward_primer = "ACGT",
    reverse_primer = "TGCA",
    merge_seq_build_log = FALSE,
    merge_file_errors = FALSE,
    merge_statistics = FALSE,
    merge_seq_reconstruction_log = FALSE,
    merge_kmer = FALSE,
    merge_discard_ambig = TRUE
  )

  args <- .pandaseq_extra_args(settings)
  expect_true(all(c("-p", "ACGT", "-q", "TGCA", "-N") %in% args))
  expect_false("settings$forward_primer" %in% args)
  expect_false("settings$reverse_primer" %in% args)
})

test_that("single merge sample handles pandaseq result status", {
  row <- data.frame(
    Sample_Name = "sample_a",
    R1_Path = "/tmp/r1.fastq.gz",
    R2_Path = "/tmp/r2.fastq.gz",
    Source = "single",
    check.names = FALSE
  )
  settings <- list(
    merge_prefix = "merge",
    merge_filetype = "FASTQ",
    min_length = 100,
    max_length = 200,
    min_overlap = 20,
    merge_similarity = 0.9,
    forward_primer = "",
    reverse_primer = "",
    merge_seq_build_log = FALSE,
    merge_file_errors = FALSE,
    merge_statistics = FALSE,
    merge_seq_reconstruction_log = FALSE,
    merge_kmer = FALSE,
    merge_discard_ambig = FALSE
  )

  ok_fun <- function(...) {
    structure(list(ok = TRUE, output_fastq = "out.fastq"), class = c("pandaseq_result", "list"))
  }
  fail_fun <- function(...) {
    structure(
      list(ok = FALSE, phase = "run", status = 2L, stderr = "bad overlap"),
      class = c("pandaseq_result", "list")
    )
  }

  expect_true(.run_single_merge_sample(row, settings, merge_fun = ok_fun)$ok)
  expect_error(
    .run_single_merge_sample(row, settings, merge_fun = fail_fun),
    "PANDAseq merge failed for sample_a"
  )
  expect_error(
    .run_single_merge_sample(row, settings, merge_fun = fail_fun),
    "bad overlap"
  )
})

test_that("single merge sample passes explicit output and log paths unchanged", {
  row <- data.frame(
    Sample_Name = "sample_a",
    R1_Path = "/tmp/r1.fastq.gz",
    R2_Path = "/tmp/r2.fastq.gz",
    Source = "single",
    Output_Path = "/tmp/merge/session_1/merge-sample_a.fastq",
    Log_Path = "/tmp/merge/session_1/sample_a-log.txt",
    check.names = FALSE
  )
  settings <- list(
    merge_prefix = "ignored",
    merge_filetype = "FASTA",
    min_length = 100,
    max_length = 200,
    min_overlap = 20,
    merge_similarity = 0.9,
    forward_primer = "",
    reverse_primer = "",
    merge_seq_build_log = FALSE,
    merge_file_errors = FALSE,
    merge_statistics = FALSE,
    merge_seq_reconstruction_log = FALSE,
    merge_kmer = FALSE,
    merge_discard_ambig = FALSE
  )

  seen <- NULL
  capture_fun <- function(...) {
    seen <<- list(...)
    structure(list(ok = TRUE, output_fastq = seen$output_fastq), class = c("pandaseq_result", "list"))
  }

  .run_single_merge_sample(row, settings, merge_fun = capture_fun)

  expect_equal(seen$output_fastq, row[["Output_Path"]])
  expect_equal(seen$log_file, row[["Log_Path"]])
})
