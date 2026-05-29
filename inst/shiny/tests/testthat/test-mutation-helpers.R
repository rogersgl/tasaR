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
  row[["Result Path"]] <- ""
  row[["Output Path"]] <- ""
  row[["Download Href"]] <- ""
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

test_that("mutation table renders completed downloads", {
  row <- valid_mutation_row()
  row[["Mutation Queue ID"]] <- "mut_1"
  row[["Progress"]] <- "done"
  row[["Output Path"]] <- "/tmp/sample_1.zip"
  row[["Download Href"]] <- "mutation-results-token/sample_1.zip"
  row <- row[, names(.mutation_empty_queue()), drop = FALSE]

  table <- .make_mutation_table_data(row, "Mutation Queue ID", "remove_mut")
  expect_match(table[["Progress"]], "queue-progress-download", fixed = TRUE)
  expect_match(table[["Progress"]], "mutation-results-token/sample_1.zip", fixed = TRUE)
})

test_that("mutation analyzed count tracks completed progress states", {
  row <- valid_mutation_row()
  rows <- row[rep(1, 5), , drop = FALSE]
  rows[["Progress"]] <- c("", "processing", "done", "completed", "success")

  expect_equal(.mutation_analyzed_count(rows), 3L)
  expect_equal(.mutation_analyzed_count(.mutation_empty_queue()), 0L)
})

test_that("prep mutation settings emits tasaR-compatible names", {
  row <- valid_mutation_row()
  row[["Mutation Queue ID"]] <- "mut_1"
  row[["Antibody"]] <- TRUE
  row[["CDR3 Start"]] <- "30"
  row <- row[, names(.mutation_empty_queue()), drop = FALSE]

  settings <- .prep_tas_settings(row)
  expect_true("IsAntibody" %in% names(settings))
  expect_false("IsAntiboy" %in% names(settings))
  expect_equal(settings[["IsAntibody"]], TRUE)
  expect_equal(settings[["CDR3Start"]], "30")
})

test_that("mutation settings df preserves ordered rows for batch analysis", {
  row <- valid_mutation_row()
  rows <- row[rep(1, 2), , drop = FALSE]
  rows[["Sample Name"]] <- c("sample_1", "sample_2")
  rows[["Merged Path"]] <- c("/tmp/sample_1.fastq.gz", "/tmp/sample_2.fastq.gz")

  settings <- .mutation_settings_df(rows)

  expect_equal(settings[["Name"]], c("sample_1", "sample_2"))
  expect_equal(settings[["MergedFASTQPath"]], c("/tmp/sample_1.fastq.gz", "/tmp/sample_2.fastq.gz"))
  expect_true("IsAntibody" %in% names(settings))
})

test_that("zip directory preserves relative nested paths", {
  root <- tempfile("mutation-results-")
  dir.create(file.path(root, "sample_1", "nested", "deeper"), recursive = TRUE)
  writeLines("x", file.path(root, "sample_1", "nested", "result.txt"))
  writeLines("y", file.path(root, "sample_1", "nested", "deeper", "result.txt"))
  zipfile <- tempfile(fileext = ".zip")

  output <- utils::capture.output(.zip_directory(zipfile, file.path(root, "sample_1"), compression_level = 1L))
  entries <- utils::unzip(zipfile, list = TRUE)$Name

  expect_true(file.exists(zipfile))
  expect_false(any(grepl("adding:", output, fixed = TRUE)))
  expect_true("nested/result.txt" %in% entries)
  expect_true("nested/deeper/result.txt" %in% entries)
  expect_false("result.txt" %in% entries)
})
