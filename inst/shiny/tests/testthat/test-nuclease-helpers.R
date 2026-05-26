valid_nuclease_row <- function() {
  row <- as.data.frame(as.list(.nuclease_default_settings()), stringsAsFactors = FALSE, check.names = FALSE)
  row[["Sample Name"]] <- "sample_1"
  row[["Control File"]] <- "control.fastq.gz"
  row[["Control Path"]] <- "/tmp/control.fastq.gz"
  row[["Experimental File"]] <- "experimental.fastq.gz"
  row[["Experimental Path"]] <- "/tmp/experimental.fastq.gz"
  row[["gRNA Sequence"]] <- "ACGTACGTACGTACGTACGT"
  row[["Cut Position"]] <- "10"
  row[["Source"]] <- "single"
  row[["Forward Primer"]] <- "ACGT"
  row[["Reverse Primer"]] <- "TGCA"
  row[["Status"]] <- ""
  row[["Progress"]] <- ""
  row
}

test_that("nuclease defaults and normalization fill missing values", {
  defaults <- .nuclease_default_settings()
  expect_equal(defaults[["Forward Extension Type"]], "None")
  expect_equal(defaults[["Reverse Extension Type"]], "None")

  row <- .normalize_nuclease_row(data.frame("Cut Position" = NA, check.names = FALSE))
  expect_equal(row[["Cut Position"]], "")
  expect_equal(row[["Forward Primer"]], "")
})

test_that("nuclease status validates required settings", {
  row <- valid_nuclease_row()
  expect_equal(.nuclease_status_for_row(row, all_samples = "sample_1"), "Ready")

  row[["Cut Position"]] <- "not_numeric"
  expect_match(.nuclease_status_for_row(row, all_samples = "sample_1"), "Invalid cut position")

  row <- valid_nuclease_row()
  row[["Forward Extension Type"]] <- "UMI"
  expect_match(.nuclease_status_for_row(row, all_samples = "sample_1"), "Missing forward extension sequence")
})

test_that("nuclease queue table data includes stable control columns", {
  row <- valid_nuclease_row()
  row[["Nuclease Queue ID"]] <- "nuc_1"
  row <- row[, names(.nuclease_empty_queue()), drop = FALSE]
  row <- .nuclease_queue_with_status(row)

  table <- .make_nuclease_table_data(row, "Nuclease Queue ID", "remove_nuc")
  expect_true(all(c("::", "Display Order", "Settings", "Actions", "Sample ID") %in% names(table)))
  expect_equal(table[["Status"]], "<div class='queue-status-text queue-status-ready'>Ready</div>")
})
