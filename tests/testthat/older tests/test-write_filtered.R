################################################################################
############################ TAS_WRITE_FILTERED ################################
################################################################################


skip_if_not(exists("tas_write_filtered", mode = "function"), "tas_write_filtered() not found in this session")
skip_if_not_installed("ShortRead")

# Helper: copy function and make mclapply deterministic (use lapply)
make_fcopy_mclapply_lapply <- function(fn) {
  fcopy <- fn
  env_mock <- new.env(parent = environment(fn))
  env_mock$mclapply <- function(X, FUN, mc.cores = 1, ...) lapply(X, FUN, ...)
  environment(fcopy) <- env_mock
  fcopy
}

# Helper: create a ShortReadQ from character sequences
make_srq <- function(seq_chars, id_prefix = "r") {
  seqs <- DNAStringSet(seq_chars)
  qs <- BStringSet(vapply(width(seqs), function(len) paste(rep("I", len), collapse = ""), character(1)))
  ids <- BStringSet(paste0(id_prefix, seq_along(seqs)))
  ShortReadQ(sread = seqs, quality = qs, id = ids)
}

test_that("tas_write_filtered writes filtered FASTQ.gz files into WD/filtered", {
  fcopy <- make_fcopy_mclapply_lapply(tas_write_filtered)

  # Setup WD and sample data
  tmpwd <- file.path(tempdir(), "tas_write_filtered_test1")
  if (dir.exists(tmpwd)) unlink(tmpwd, recursive = TRUE)
  dir.create(tmpwd, recursive = TRUE)
  WD <- paste0(normalizePath(tmpwd, mustWork = TRUE), "/")

  # ensure logs directory exists because function may write logs when replacing files in other tests
  dir.create(file.path(WD, "logs"), showWarnings = FALSE)

  Sample.Names <- c("S1", "S2")
  # Create simple ShortReadQ per sample
  srq1 <- make_srq(c("ACGTACGT"))
  srq2 <- make_srq(c("TTTTGGGG"))
  Reads.Filtered.List <- list(srq1, srq2)

  cfg <- list(nCores = 1)

  # Call the function
  expect_silent(
    fcopy(Sample.Names = Sample.Names,
          Reads.Filtered.List = Reads.Filtered.List,
          Config.List = cfg,
          WD = WD)
  )

  # Validate filtered directory created and files present
  filtered_dir <- file.path(WD, "filtered")
  expect_true(dir.exists(filtered_dir))
  expect_true(file.exists(file.path(filtered_dir, "S1-filtered.fastq.gz")))
  expect_true(file.exists(file.path(filtered_dir, "S2-filtered.fastq.gz")))

  # Clean up
  try(unlink(tmpwd, recursive = TRUE), silent = TRUE)
})

test_that("tas_write_filtered replaces existing filtered file and logs replacement", {
  fcopy <- make_fcopy_mclapply_lapply(tas_write_filtered)

  tmpwd <- file.path(tempdir(), "tas_write_filtered_test2")
  if (dir.exists(tmpwd)) unlink(tmpwd, recursive = TRUE)
  dir.create(tmpwd, recursive = TRUE)
  WD <- paste0(normalizePath(tmpwd, mustWork = TRUE), "/")

  # create logs dir (function logs to WD/logs/tasAnalyzer logs.txt)
  dir.create(file.path(WD, "logs"), recursive = TRUE, showWarnings = FALSE)

  Sample.Names <- c("S1")
  srq1 <- make_srq(c("ACGTACGT"))
  Reads.Filtered.List <- list(srq1)

  cfg <- list(nCores = 1)

  filtered_dir <- file.path(WD, "filtered")
  dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
  existing_file <- file.path(filtered_dir, "S1-filtered.fastq.gz")

  # Create a dummy existing file to trigger replacement branch
  cat("OLDCONTENT", file = existing_file)
  expect_true(file.exists(existing_file))

  # Call function; it should detect existing file, append log message, remove file, and write new filtered file
  expect_silent(
    fcopy(Sample.Names = Sample.Names,
          Reads.Filtered.List = Reads.Filtered.List,
          Config.List = cfg,
          WD = WD)
  )

  # After function, the new file should exist (overwriting the old one)
  expect_true(file.exists(existing_file))

  # Check that the log file contains the replacement message
  log_file <- file.path(WD, "logs", "tasAnalyzer logs.txt")
  expect_true(file.exists(log_file))
  log_contents <- readLines(log_file, warn = FALSE)
  expect_true(any(grepl("The following file was replaced", log_contents)))

  # Clean up
  try(unlink(tmpwd, recursive = TRUE), silent = TRUE)
})
