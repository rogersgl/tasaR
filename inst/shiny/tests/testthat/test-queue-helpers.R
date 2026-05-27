test_that("sample order reconciliation preserves existing order and appends new ids", {
  expect_equal(
    .reconcile_sample_order(c("a", "b", "c"), c("b", "z", "a", "b")),
    c("b", "a", "c")
  )
})

test_that("sample data frames can be ordered and moved", {
  df <- data.frame(
    id = c("a", "b", "c"),
    value = 1:3,
    check.names = FALSE
  )

  expect_equal(.order_sample_df(df, c("c", "a"), "id")$id, c("c", "a", "b"))
  expect_equal(.move_order_item(c("a", "b", "c"), "b", "up"), c("b", "a", "c"))
  expect_equal(.move_order_item(c("a", "b", "c"), "b", "down"), c("a", "c", "b"))
  expect_equal(.move_order_item(c("a", "b", "c"), "a", "up"), c("a", "b", "c"))
})

test_that("rows can be removed by id", {
  df <- data.frame(
    "Queue ID" = c("a", "b"),
    "Sample Name" = c("sample_a", "sample_b"),
    check.names = FALSE
  )

  res <- .remove_batch_row_by_id(df, "b", id_col = "Queue ID")
  expect_true(res$success)
  expect_equal(res$removed[["Sample Name"]], "sample_b")
  expect_equal(res$df[["Queue ID"]], "a")

  missing <- .remove_batch_row_by_id(df, "z", id_col = "Queue ID")
  expect_false(missing$success)
  expect_equal(missing$df, df)
})

test_that("queue action HTML keeps expected classes and labels", {
  html <- .queue_actions_html("row_1", "remove_id", "move_up_id", "move_down_id")

  expect_match(html, "queue-move-btn", fixed = TRUE)
  expect_match(html, "queue-remove-btn", fixed = TRUE)
  expect_match(html, "aria-label='Move up'", fixed = TRUE)
  expect_match(html, "aria-label='Move down'", fixed = TRUE)
  expect_match(html, "aria-label='Remove'", fixed = TRUE)
})

test_that("queue progress HTML renders processing and complete states", {
  html <- .sample_progress_html(c("processing", "done", ""))

  expect_match(html[[1]], "queue-progress-spinner", fixed = TRUE)
  expect_match(html[[2]], "queue-progress-done", fixed = TRUE)
  expect_match(html[[3]], "queue-progress-empty", fixed = TRUE)
})

test_that("queue progress HTML renders row download links for completed files", {
  html <- .sample_progress_html(
    "done",
    download_href = "merge-results-session/merged.fastq",
    download_filename = "merged.fastq"
  )

  expect_match(html[[1]], "queue-progress-complete", fixed = TRUE)
  expect_match(html[[1]], "queue-progress-download", fixed = TRUE)
  expect_match(html[[1]], "fa-download", fixed = TRUE)
  expect_match(html[[1]], "href='merge-results-session/merged.fastq'", fixed = TRUE)
  expect_match(html[[1]], "download='merged.fastq'", fixed = TRUE)
})

test_that("queue small text column defs target the expected visible columns", {
  merge_df <- .make_sample_table_data(
    data.frame(
      "Sample Name" = "sample_1",
      "Left (R1) File" = "r1.fastq.gz",
      "Right (R2) File" = "r2.fastq.gz",
      "Source" = "single",
      "Status" = "",
      check.names = FALSE
    ),
    id_col = "Sample Name",
    remove_input_id = "remove_id"
  )
  merge_defs <- .queue_small_text_column_defs(
    merge_df,
    c("Left (R1) File", "Right (R2) File")
  )
  expect_equal(merge_defs[[1]]$targets, c(3L, 4L))
  expect_false("Source" %in% names(merge_df))

  mutation_df <- .make_mutation_table_data(
    within(
      data.frame(as.list(.mutation_default_settings()), stringsAsFactors = FALSE, check.names = FALSE),
      {
        `Mutation Queue ID` <- "mut_1"
        `Sample Name` <- "sample_1"
        `Merged File` <- "merged.fastq.gz"
        `Merged Path` <- "/tmp/merged.fastq.gz"
        Source <- "single"
        Status <- ""
        Progress <- ""
      }
    ),
    id_col = "Mutation Queue ID",
    remove_input_id = "remove_id"
  )
  mutation_defs <- .queue_small_text_column_defs(
    mutation_df,
    c("Merged File")
  )
  expect_equal(mutation_defs[[1]]$targets, 3L)
  expect_false("Source" %in% names(mutation_df))

  nuclease_df <- .make_nuclease_table_data(
    within(
      data.frame(as.list(.nuclease_default_settings()), stringsAsFactors = FALSE, check.names = FALSE),
      {
        `Nuclease Queue ID` <- "nuc_1"
        `Sample Name` <- "sample_1"
        `Control File` <- "control.fastq.gz"
        `Control Path` <- "/tmp/control.fastq.gz"
        `Experimental File` <- "experimental.fastq.gz"
        `Experimental Path` <- "/tmp/experimental.fastq.gz"
        `gRNA Sequence` <- "ACGT"
        `Cut Position` <- "10"
        Source <- "single"
        Status <- ""
        Progress <- ""
      }
    ),
    id_col = "Nuclease Queue ID",
    remove_input_id = "remove_id"
  )
  nuclease_defs <- .queue_small_text_column_defs(
    nuclease_df,
    c("Control File", "Experimental File")
  )
  expect_equal(nuclease_defs[[1]]$targets, c(3L, 4L))
  expect_false("Source" %in% names(nuclease_df))
})

test_that("queue table can target progress column for centered styling", {
  table_df <- .make_sample_table_data(
    data.frame(
      "Queue ID" = "queue_1",
      "Sample Name" = "sample_1",
      "Left (R1) File" = "r1.fastq.gz",
      "Right (R2) File" = "r2.fastq.gz",
      "Source" = "single",
      "Status" = "Ready",
      "Progress" = "processing",
      check.names = FALSE
    ),
    id_col = "Queue ID",
    remove_input_id = "remove_id"
  )

  progress_col_index <- match("Progress", names(table_df)) - 1L
  progress_col_def <- if (!is.na(progress_col_index)) {
    list(list(className = "queue-progress-cell", targets = progress_col_index))
  } else {
    list()
  }

  expect_equal(progress_col_def[[1]]$className, "queue-progress-cell")
  expect_equal(progress_col_def[[1]]$targets, 7L)
})
