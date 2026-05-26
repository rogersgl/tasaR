valid_mutation_row <- function() {
  row <- as.data.frame(as.list(.mutation_default_settings()), stringsAsFactors = FALSE, check.names = FALSE)
  row[["Sample Name"]] <- "sample_1"
  row[["Merged File"]] <- "sample_1.fastq.gz"
  row[["Merged Path"]] <- "/tmp/sample_1.fastq.gz"
  row[["Source"]] <- "single"
  row[["Reference Sequence"]] <- "ACGTACGT"
  row[["Forward Primer"]] <- "ACGT"
  row[["Reverse Primer"]] <- "TGCA"
  row[["Amplicon Length"]] <- "100"
  row[["Insert Start Coordinate"]] <- "1"
  row[["Insert End Coordinate"]] <- "-1"
  row[["Status"]] <- ""
  row[["Progress"]] <- ""
  row
}

test_that("mutation defaults and normalization fill missing values", {
  defaults <- .mutation_default_settings()
  expect_false(defaults[["Antibody"]])
  expect_equal(defaults[["Forward Extension Type"]], "None")

  row <- .normalize_mutation_row(data.frame("Antibody" = NA, check.names = FALSE))
  expect_false(row[["Antibody"]])
  expect_equal(row[["Reference Sequence"]], "")
})

test_that("mutation status validates required settings", {
  row <- valid_mutation_row()
  expect_equal(.mutation_status_for_row(row, all_samples = "sample_1"), "Ready")

  row[["Forward Extension Type"]] <- "Barcode"
  expect_match(.mutation_status_for_row(row, all_samples = "sample_1"), "Missing forward extension sequence")

  row <- valid_mutation_row()
  row[["Antibody"]] <- TRUE
  expect_match(.mutation_status_for_row(row, all_samples = "sample_1"), "Missing fr1 start")
})

test_that("mutation queue table data includes stable control columns", {
  row <- valid_mutation_row()
  row[["Mutation Queue ID"]] <- "mut_1"
  row <- row[, names(.mutation_empty_queue()), drop = FALSE]
  row <- .mutation_queue_with_status(row)

  table <- .make_mutation_table_data(row, "Mutation Queue ID", "remove_mut")
  expect_true(all(c("::", "Display Order", "Settings", "Actions", "Sample ID") %in% names(table)))
  expect_equal(table[["Status"]], "<div class='queue-status-text queue-status-ready'>Ready</div>")
})
