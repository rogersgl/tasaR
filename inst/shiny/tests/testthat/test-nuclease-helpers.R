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

  row <- .normalize_nuclease_row(data.frame(
    "gRNA Sequence" = "acguacguacguacguacgu",
    "Cut Position" = NA,
    check.names = FALSE
  ))
  expect_equal(row[["Cut Position"]], "")
  expect_equal(row[["Forward Primer"]], "")
  expect_equal(row[["gRNA Sequence"]], "ACGTACGTACGTACGTACGT")
})

test_that("nuclease gRNA parser and validator enforce sequence rules", {
  expect_equal(.parse_nuclease_grna(" acguacguacguacguacgu "), "ACGTACGTACGTACGTACGT")
  expect_true(.valid_nuclease_grna("acguacguacguacguacgu"))
  expect_false(.valid_nuclease_grna("ACGTACGTACGTACGTACGN"))
  expect_false(.valid_nuclease_grna("ACGTACGTACGTACGTACG"))
  expect_false(.valid_nuclease_grna("ACGTACGTACGTACGTACGTA"))
})

test_that("nuclease status validates required settings", {
  row <- valid_nuclease_row()
  expect_equal(.nuclease_status_for_row(row, all_samples = "sample_1"), "Ready")

  row <- valid_nuclease_row()
  row[["Cut Position"]] <- ""
  expect_equal(.nuclease_status_for_row(row, all_samples = "sample_1"), "Ready")

  row <- valid_nuclease_row()
  row[["gRNA Sequence"]] <- "acguacguacguacguacgu"
  row[["Cut Position"]] <- ""
  expect_equal(.nuclease_status_for_row(row, all_samples = "sample_1"), "Ready")

  row <- valid_nuclease_row()
  row[["gRNA Sequence"]] <- "ACGTACGTACGTACGTACGN"
  row[["Cut Position"]] <- ""
  expect_match(.nuclease_status_for_row(row, all_samples = "sample_1"), "Invalid gRNA sequence")

  row <- valid_nuclease_row()
  row[["gRNA Sequence"]] <- "ACGTACGTACGTACGTACG"
  row[["Cut Position"]] <- ""
  expect_match(.nuclease_status_for_row(row, all_samples = "sample_1"), "Invalid gRNA sequence")

  row <- valid_nuclease_row()
  row[["gRNA Sequence"]] <- "ACGTACGTACGTACGTACGN"
  expect_equal(.nuclease_status_for_row(row, all_samples = "sample_1"), "Ready")

  row <- valid_nuclease_row()
  row[["gRNA Sequence"]] <- ""
  expect_equal(.nuclease_status_for_row(row, all_samples = "sample_1"), "Ready")

  row <- valid_nuclease_row()
  row[["gRNA Sequence"]] <- ""
  row[["Cut Position"]] <- ""
  expect_match(.nuclease_status_for_row(row, all_samples = "sample_1"), "Missing gRNA sequence or cut position")

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
  expect_true(all(c("::", "Display Order", "Target", "Actions", "Sample ID") %in% names(table)))
  expect_false(any(c("Settings", "gRNA Sequence", "Cut Position") %in% names(table)))
  expect_equal(table[["Status"]], "<div class='queue-status-text queue-status-ready'>Ready</div>")
})

test_that("nuclease target display prefers manual cut position", {
  row <- valid_nuclease_row()
  expect_equal(.nuclease_target_display(row), "10")

  row[["Cut Position"]] <- ""
  expect_equal(.nuclease_target_display(row), "ACGTACGTACGTACGTACGT")

  row[["gRNA Sequence"]] <- "acguacguacguacguacgu"
  expect_equal(.nuclease_target_display(row), "ACGTACGTACGTACGTACGT")

  row[["gRNA Sequence"]] <- ""
  expect_equal(.nuclease_target_display(row), "auto")

  row[["gRNA Sequence"]] <- "ACGTACGTACGTACGTACGT"
  row[["Cut Position"]] <- "not_numeric"
  expect_equal(.nuclease_target_display(row), "not_numeric")
})

test_that("nuclease target cell includes edit action", {
  row <- valid_nuclease_row()
  cell <- .nuclease_target_cell(row, "nuc_1", "nuclease_open_settings")

  expect_match(cell, "10", fixed = TRUE)
  expect_match(cell, "queue-settings-edit-btn", fixed = TRUE)
  expect_match(cell, "nuclease_open_settings", fixed = TRUE)
  expect_match(cell, "nuc_1", fixed = TRUE)
})
