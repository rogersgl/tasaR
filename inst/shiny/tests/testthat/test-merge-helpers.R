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
