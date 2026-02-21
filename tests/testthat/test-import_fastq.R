
################################################################################
############################# TAS_IMPORT_FASTQ #################################
################################################################################


skip_if_not(exists("tas_import_fastq", mode = "function"),
            "tas_import_fastq() not found in this session")

# Helper to construct Fastq.File.Path exactly as in your package
make_fastq_paths <- function(Sample.Names, WD) {
  sapply(Sample.Names, function(x) {
    str_c(WD, "merged/", x, "-merged.fastq.gz")
  }, USE.NAMES = TRUE)
}

test_that("tas_import_fastq (mocked readFastq) returns named list with names from Sample.Names and paths constructed from WD", {
  skip_if_not(exists("tas_import_fastq", mode = "function"))

  # copy function so we can attach a mock environment
  fcopy <- tas_import_fastq

  # create mock environment with readFastq stub
  env_mock <- new.env(parent = environment(tas_import_fastq))
  env_mock$readFastq_calls <- character(0)
  env_mock$readFastq <- function(path) {
    env_mock$readFastq_calls <<- c(env_mock$readFastq_calls, path)
    paste0("MOCK_FASTQ_OBJECT_FOR_", basename(path))
  }
  environment(fcopy) <- env_mock

  # Build sample names and WD and paths using the exact sapply recipe
  Sample.Names <- c("S1", "S2")
  WD <- paste0(normalizePath(tempdir(), mustWork = TRUE), "/")
  dir.create(file.path(WD, "merged"), recursive = TRUE, showWarnings = FALSE)

  Fastq.File.Path <- make_fastq_paths(Sample.Names, WD)

  # Ensure the vector is named by sample names (sapply does this)
  expect_named(Fastq.File.Path, Sample.Names)

  result <- fcopy(Sample.Names = Sample.Names, Fastq.File.Path = Fastq.File.Path)

  # Structure checks
  expect_type(result, "list")
  expect_named(result, Sample.Names)

  # Content checks (mock return values)
  expect_equal(result$S1, paste0("MOCK_FASTQ_OBJECT_FOR_", basename(Fastq.File.Path["S1"])))
  expect_equal(result$S2, paste0("MOCK_FASTQ_OBJECT_FOR_", basename(Fastq.File.Path["S2"])))

  # Ensure readFastq called once per file with the exact constructed paths
  expect_equal(length(env_mock$readFastq_calls), 2)
  expect_equal(env_mock$readFastq_calls, as.vector(Fastq.File.Path))
})

test_that("tas_import_fastq propagates readFastq errors (mocked)", {
  skip_if_not(exists("tas_import_fastq", mode = "function"))

  fcopy <- tas_import_fastq
  env_mock <- new.env(parent = environment(tas_import_fastq))
  env_mock$readFastq <- function(path) stop("readFastq failed for path: ", path)
  environment(fcopy) <- env_mock

  Sample.Names <- c("S1")
  WD <- paste0(normalizePath(tempdir(), mustWork = TRUE), "/")
  dir.create(file.path(WD, "merged"), recursive = TRUE, showWarnings = FALSE)
  Fastq.File.Path <- make_fastq_paths(Sample.Names, WD)

  expect_error(
    fcopy(Sample.Names = Sample.Names, Fastq.File.Path = Fastq.File.Path),
    regexp = "readFastq failed for path",
    fixed = FALSE
  )
})

test_that("tas_import_fastq integration test reads real gzipped FASTQ files (ShortRead)", {
  skip_if_not_installed("ShortRead")
  skip_if_not(exists("tas_import_fastq", mode = "function"))

  # Setup WD and merged directory
  Sample.Names <- c("S1", "S2")
  WD <- paste0(normalizePath(tempdir(), mustWork = TRUE), "/")
  merged_dir <- file.path(WD, "merged")
  if (!dir.exists(merged_dir)) dir.create(merged_dir, recursive = TRUE)

  # Create two valid FASTQ records and write them to gzipped files at the expected locations
  fastq_content <- "@SEQ_ID\nACGT\n+\n!!!!\n"  # exactly 4 non-empty lines per record

  paths <- make_fastq_paths(Sample.Names, WD)
  # write valid gzipped FASTQ using gzfile()
  for (p in paths) {
    con <- gzfile(p, "w")
    cat(fastq_content, file = con)
    close(con)
    expect_true(file.exists(p))
  }

  # Call real function (uses ShortRead::readFastq internally)
  result <- tas_import_fastq(Sample.Names = Sample.Names, Fastq.File.Path = paths)

  expect_named(result, Sample.Names)
  # Verify returned objects are ShortReadQ
  expect_s4_class(result$S1, "ShortReadQ")
  expect_s4_class(result$S2, "ShortReadQ")

  # Clean up test files (best-effort)
  try(unlink(merged_dir, recursive = TRUE), silent = TRUE)
})

test_that("Fastq.File.Path naming: sapply creates named vector matching Sample.Names", {
  Sample.Names <- c("S1", "S2", "sample-3")
  WD <- paste0(normalizePath(tempdir(), mustWork = TRUE), "/")
  Fastq.File.Path <- make_fastq_paths(Sample.Names, WD)
  expect_named(Fastq.File.Path, Sample.Names)
  # Check that constructed paths end with the expected suffix
  expect_true(all(grepl("-merged.fastq.gz$", as.vector(Fastq.File.Path))))
})
