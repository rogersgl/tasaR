
################################################################################
########################## TAS_PANDASEQ_TERMINAL ###############################
################################################################################

skip_if_not(exists("tas_pandaseq_terminal", mode = "function"), "tas_pandaseq_terminal() not found in this session")

# Helper: minimal valid Input.DataFrame and Config.List for tas_pandaseq_terminal
make_valid_input_df_for_pandaseq <- function(forward, reverse, amp_len = 300,
                                             max_del = 25, max_ins = 25, os = "Linux") {
  data.frame(
    SampleName = "S1",
    ForwardFASTQFileName = forward,
    ReverseFASTQFileName = reverse,
    AmpliconLength = amp_len,
    MaxDeletion = max_del,
    MaxInsertion = max_ins,
    OperatingSystem = os, # note: used by pandaseq check error path
    stringsAsFactors = FALSE
  )
}

make_valid_config <- function() {
  list(
    WorkingDirectory = tempdir(),
    OperatingSystem = "Linux",
    nCores = 1,
    merge.reads = TRUE
  )
}

test_that("input column validation for tas_pandaseq_terminal", {
  skip_if_not(exists("tas_pandaseq_terminal", mode = "function"))
  # Missing required column (ReverseFASTQFileName)
  df <- make_valid_input_df_for_pandaseq("a.fastq", "b.fastq")
  df$ReverseFASTQFileName <- NULL
  cfg <- make_valid_config()
  expect_error(
    tas_pandaseq_terminal(Sample.Names = "S1",
                          Input.DataFrame = df,
                          Config.List = cfg,
                          WD = paste0(tempdir(), "/"),
                          shiny.env = FALSE,
                          shiny.fileTable = NULL),
    regexp = "The following columns were not found in the input file"
  )
})

test_that("config entries validation for tas_pandaseq_terminal", {
  skip_if_not(exists("tas_pandaseq_terminal", mode = "function"))
  df <- make_valid_input_df_for_pandaseq("a.fastq", "b.fastq")
  cfg <- make_valid_config()
  cfg$merge.reads <- NULL
  expect_error(
    tas_pandaseq_terminal(Sample.Names = "S1",
                          Input.DataFrame = df,
                          Config.List = cfg,
                          WD = paste0(tempdir(), "/"),
                          shiny.env = FALSE,
                          shiny.fileTable = NULL),
    regexp = "The following columns were not found in the input file"
  )
})

test_that("pandaseq missing triggers OS-specific install message (MacOS)", {
  skip_if_not(exists("tas_pandaseq_terminal", mode = "function"))
  # Make a copy of the function and mock system2 to return something other than the expected check string
  fcopy <- tas_pandaseq_terminal
  body_text <- paste(deparse(body(fcopy)), collapse = "\n")
  # Replace R.utils::gzip with gzip so we can mock it in the local environment of the copy
  body_text <- gsub("R.utils::gzip", "gzip", body_text, fixed = TRUE)
  body(fcopy) <- parse(text = body_text)[[1]]

  # Prepare env with a system2 mock that returns something else (to simulate pandaseq not present)
  env_mock <- new.env(parent = environment(tas_pandaseq_terminal))
  env_mock$system2_calls <- list()
  env_mock$system2 <- function(command, args = NULL, stdout = FALSE, stderr = FALSE, ...) {
    # Record the call
    env_mock$system2_calls[[length(env_mock$system2_calls) + 1]] <<- list(command = command, args = args)
    # Return text that is *not* the expected "You must supply both forward and reverse reads."
    return("pandaseq not found")
  }
  # Provide trivial mclapply and gzip mocks to avoid runtime errors if execution reaches them (it shouldn't here)
  env_mock$mclapply <- function(X, FUN, mc.cores = 1, ...) lapply(X, FUN, ...)
  env_mock$gzip <- function(...) invisible(NULL)

  # Place fcopy into that environment
  environment(fcopy) <- env_mock

  df <- make_valid_input_df_for_pandaseq("fwd.fastq", "rev.fastq", os = "MacOS")
  cfg <- make_valid_config()
  # Provide WD; function will run the check and should stop with MacOS specific message
  expect_error(
    fcopy(Sample.Names = "S1",
          Input.DataFrame = df,
          Config.List = cfg,
          WD = paste0(tempdir(), "/"),
          shiny.env = FALSE,
          shiny.fileTable = NULL),
    regexp = "Please install the PANDAseq.pkg",
    ignore.case = TRUE
  )
  # Also assert that our system2 mock was called for the check
  expect_true(length(env_mock$system2_calls) >= 1)
})

test_that("pandaseq missing triggers OS-specific install message (Linux)", {
  skip_if_not(exists("tas_pandaseq_terminal", mode = "function"))
  fcopy <- tas_pandaseq_terminal
  body_text <- paste(deparse(body(fcopy)), collapse = "\n")
  body_text <- gsub("R.utils::gzip", "gzip", body_text, fixed = TRUE)
  body(fcopy) <- parse(text = body_text)[[1]]

  env_mock <- new.env(parent = environment(tas_pandaseq_terminal))
  env_mock$system2_calls <- list()
  env_mock$system2 <- function(command, args = NULL, stdout = FALSE, stderr = FALSE, ...) {
    env_mock$system2_calls[[length(env_mock$system2_calls) + 1]] <<- list(command = command, args = args)
    return("pandaseq not found")
  }
  env_mock$mclapply <- function(X, FUN, mc.cores = 1, ...) lapply(X, FUN, ...)
  env_mock$gzip <- function(...) invisible(NULL)
  environment(fcopy) <- env_mock

  df <- make_valid_input_df_for_pandaseq("fwd.fastq", "rev.fastq", os = "Linux")
  cfg <- make_valid_config()
  expect_error(
    fcopy(Sample.Names = "S1",
          Input.DataFrame = df,
          Config.List = cfg,
          WD = paste0(tempdir(), "/"),
          shiny.env = FALSE,
          shiny.fileTable = NULL),
    regexp = "Please install as appropriate for your distribution",
    ignore.case = TRUE
  )
  expect_true(length(env_mock$system2_calls) >= 1)
})

test_that("successful run path: directories created, files moved, pandaseq and gzip invoked (mocked)", {
  skip_if_not(exists("tas_pandaseq_terminal", mode = "function"))
  # Copy function and prepare it to use mocks (replace R.utils::gzip text)
  fcopy <- tas_pandaseq_terminal
  body_text <- paste(deparse(body(fcopy)), collapse = "\n")
  body_text <- gsub("R.utils::gzip", "gzip", body_text, fixed = TRUE)
  body(fcopy) <- parse(text = body_text)[[1]]

  # Environment to mock side effects and capture calls
  env_mock <- new.env(parent = environment(tas_pandaseq_terminal))
  env_mock$system2_calls <- list()
  env_mock$system2 <- function(command, args = NULL, stdout = FALSE, stderr = FALSE, ...) {
    # If args is NULL or missing => initial check; return the exact expected check string
    env_mock$system2_calls[[length(env_mock$system2_calls) + 1]] <<- list(command = command, args = args)
    if (is.null(args) || identical(args, character(0))) {
      return("You must supply both forward and reverse reads.")
    } else {
      # simulate successful pandaseq run: create the merged file so later gzip call can see something if needed
      # but we won't rely on real gzip: the gzip mock will be recorded instead.
      return("") # normal command execution returns ""
    }
  }

  # Capture gzip invocations
  env_mock$gzip_calls <- list()
  env_mock$gzip <- function(src, destname = NULL, remove = TRUE, overwrite = TRUE, ...) {
    env_mock$gzip_calls[[length(env_mock$gzip_calls) + 1]] <<- list(src = src, destname = destname, remove = remove, overwrite = overwrite)
    # Optionally, create the destname file to simulate compression outcome (not strictly necessary)
    if (!is.null(destname)) {
      # create an empty gz file as marker
      cat("", file = destname)
    }
    invisible(TRUE)
  }

  # Mock mclapply to call the function sequentially and return results
  env_mock$mclapply <- function(X, FUN, mc.cores = 1, ...) {
    lapply(X, FUN, ...)
  }

  # Attach other function dependencies via parent env so stringr::str_* etc are available
  environment(fcopy) <- env_mock

  # Prepare a temp WD and two sample FASTQ files in WD to simulate the non-shiny path
  tmpwd <- file.path(tempdir(), "tas_pandaseq_test_wd")
  dir.create(tmpwd, showWarnings = FALSE)
  WD <- paste0(normalizePath(tmpwd, mustWork = TRUE), "/")

  fwd <- "sample1_R1.fastq"
  rev <- "sample1_R2.fastq"
  # Create the fastq files in WD (the function will rename them into WD/unpaired/)
  file.create(file.path(WD, fwd))
  file.create(file.path(WD, rev))

  # Build Input.DataFrame and Config.List
  df <- make_valid_input_df_for_pandaseq(fwd, rev, amp_len = 250, max_del = 5, max_ins = 10, os = "Linux")
  cfg <- make_valid_config()
  cfg$nCores <- 1

  # Run the mocked function
  expect_silent(
    fcopy(Sample.Names = "S1",
          Input.DataFrame = df,
          Config.List = cfg,
          WD = WD,
          shiny.env = FALSE,
          shiny.fileTable = NULL)
  )

  # Check directories created
  expect_true(dir.exists(file.path(WD, "unpaired")))
  expect_true(dir.exists(file.path(WD, "merged")))
  expect_true(dir.exists(file.path(WD, "logs")))

  # Check that the original files were moved into unpaired
  expect_false(file.exists(file.path(WD, fwd)))
  expect_false(file.exists(file.path(WD, rev)))
  expect_true(file.exists(file.path(WD, "unpaired", fwd)))
  expect_true(file.exists(file.path(WD, "unpaired", rev)))

  # Check that system2 was called at least once for the initial check and then for pandaseq invocation
  expect_true(length(env_mock$system2_calls) >= 1)
  # Check that gzip got called for the merged file (our mock records calls)
  expect_true(length(env_mock$gzip_calls) >= 1)
  # Validate the destname shape is what function constructs (ends with -merged.fastq.gz)
  dests <- vapply(env_mock$gzip_calls, function(x) basename(x$destname), character(1))
  expect_true(any(grepl("-merged.fastq.gz$", dests)))
})



test_that("shiny.env = TRUE branch moves uploaded files from shiny.fileTable and runs pandaseq/gzip (mocked)", {
  skip_if_not(exists("tas_pandaseq_terminal", mode = "function"))

  # Make a copy of the function and replace R.utils::gzip with gzip for mocking
  fcopy <- tas_pandaseq_terminal
  body_text <- paste(deparse(body(fcopy)), collapse = "\n")
  body_text <- gsub("R.utils::gzip", "gzip", body_text, fixed = TRUE)
  body(fcopy) <- parse(text = body_text)[[1]]

  # Prepare environment to mock system2, mclapply and gzip and to capture calls
  env_mock <- new.env(parent = environment(tas_pandaseq_terminal))
  env_mock$system2_calls <- list()
  env_mock$system2 <- function(command, args = NULL, stdout = FALSE, stderr = FALSE, ...) {
    env_mock$system2_calls[[length(env_mock$system2_calls) + 1]] <<- list(command = command, args = args)
    # Return the expected pandaseq check string for the initial (no-args) check
    if (is.null(args) || identical(args, character(0))) {
      return("You must supply both forward and reverse reads.")
    } else {
      return("") # simulate successful pandaseq invocation for subsequent calls
    }
  }

  env_mock$gzip_calls <- list()
  env_mock$gzip <- function(src, destname = NULL, remove = TRUE, overwrite = TRUE, ...) {
    env_mock$gzip_calls[[length(env_mock$gzip_calls) + 1]] <<- list(src = src, destname = destname, remove = remove, overwrite = overwrite)
    # simulate creating gz file marker
    if (!is.null(destname)) cat("", file = destname)
    invisible(TRUE)
  }

  env_mock$mclapply <- function(X, FUN, mc.cores = 1, ...) lapply(X, FUN, ...)

  # Put the function into the mocked environment
  environment(fcopy) <- env_mock

  # Prepare a temp WD and two temporary upload files to simulate shiny uploads
  tmpwd <- file.path(tempdir(), "tas_pandaseq_shiny_wd")
  if (!dir.exists(tmpwd)) dir.create(tmpwd, recursive = TRUE, showWarnings = FALSE)
  WD <- paste0(normalizePath(tmpwd, mustWork = TRUE), "/")

  # Create two temp files that act as uploaded FASTQ files (datapath)
  upload1 <- tempfile(pattern = "upload1_", tmpdir = tempdir(), fileext = ".fastq")
  upload2 <- tempfile(pattern = "upload2_", tmpdir = tempdir(), fileext = ".fastq")
  cat("fake", file = upload1)
  cat("fake", file = upload2)

  # Build shiny.fileTable with names and datapath
  shiny.fileTable <- data.frame(
    name = c("shiny_forward.fastq", "shiny_reverse.fastq"),
    datapath = c(upload1, upload2),
    stringsAsFactors = FALSE
  )

  # Input.DataFrame with matching filenames in the 'name' column (function uses Input.DataFrame$ForwardFASTQFileName etc.)
  df <- data.frame(
    SampleName = "S1",
    ForwardFASTQFileName = "shiny_forward.fastq",
    ReverseFASTQFileName = "shiny_reverse.fastq",
    AmpliconLength = 250,
    MaxDeletion = 5,
    MaxInsertion = 10,
    OperatingSystem = "Linux",
    stringsAsFactors = FALSE
  )

  cfg <- list(
    WorkingDirectory = WD,
    OperatingSystem = "Linux",
    nCores = 1,
    merge.reads = TRUE
  )

  # Call the function (should be silent because we mock external commands)
  expect_silent(
    fcopy(Sample.Names = "S1",
          Input.DataFrame = df,
          Config.List = cfg,
          WD = WD,
          shiny.env = TRUE,
          shiny.fileTable = shiny.fileTable)
  )

  # After call, uploaded files should have been moved into WD/unpaired with the names from shiny.fileTable$name
  expect_true(file.exists(file.path(WD, "unpaired", "shiny_forward.fastq")))
  expect_true(file.exists(file.path(WD, "unpaired", "shiny_reverse.fastq")))

  # pandaseq initial check and invocation(s) should have been attempted (system2_calls recorded)
  expect_true(length(env_mock$system2_calls) >= 1)

  # gzip should have been invoked at least once (mock recorded)
  expect_true(length(env_mock$gzip_calls) >= 1)

})


################################################################################
################################ TAS_IMPORT ####################################
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
