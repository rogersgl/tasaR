################################################################################
########################## TAS_LOAD_DEPENDENCIES ###############################
################################################################################

test_that("Required packages are loaded", {
  packages <- c("BiocManager",
                "R.utils",
                "zip",
                "data.table",
                "openxlsx",
                "stringr",
                "ggplot2",
                "ggseqlogo",
                "parallel",
                "parallelly",
                "seqinr",
                "ape",
                "Biostrings",
                "ShortRead",
                "pwalign",
                "msa",
                "ggmsa")

  #check all
  expect_all_true(packages %in% .packages())

  #check 1 by 1
  for (i in packages){
    expect_true(i %in% .packages())
    }
})


################################################################################
################################ TAS_IMPORT ####################################
################################################################################


test_that("tas_import returns a data.frame and replaces NA with empty string when file exists", {
  # preserve previous shiny.env and restore on exit
  old_shiny <- if (exists("shiny.env", envir = .GlobalEnv)) get("shiny.env", envir = .GlobalEnv) else NULL
  assign("shiny.env", FALSE, envir = .GlobalEnv)
  on.exit({
    if (is.null(old_shiny)) {
      if (exists("shiny.env", envir = .GlobalEnv)) rm("shiny.env", envir = .GlobalEnv)
    } else {
      assign("shiny.env", old_shiny, envir = .GlobalEnv)
    }
  }, add = TRUE)

  # create a temporary CSV with some NA values
  tmp <- tempfile(fileext = ".csv")
  dat <- data.frame(
    id = 1:3,
    name = c("Alice", NA, "Charlie"),
    score = c(NA, 85, 92),
    stringsAsFactors = FALSE
  )
  write.csv(dat, tmp, row.names = FALSE, na = "") # write with empty strings for NA for easy round-trip
  # overwrite one value to be explicit NA in file to simulate NA in read.csv:
  # write.csv wrote empty strings for NA; to test replacement, write one explicit NA text row
  # simpler: just write using utils::write.table with NA
  # write.table(dat, tmp, sep = ",", row.names = FALSE, col.names = TRUE, na = "", quote = TRUE)

  # call function
  result <- tas_import(tmp)

  # expectations
  expect_s3_class(result$Input, "data.frame")
  # ensure no NA remain (converted to "")
  expect_false(any(is.na(result)))
  # check that previously NA entries became empty string
  expect_equal(as.character(result$Input$name[2]), "")
  expect_equal(as.character(result$Input$score[1]), "")
})


#################################################################################


test_that("tas_import throws informative error when file is missing and shiny.env is FALSE", {
  # ensure shiny.env is FALSE
  old_shiny <- if (exists("shiny.env", envir = .GlobalEnv)){get("shiny.env", envir = .GlobalEnv)} else NULL
  assign("shiny.env", FALSE, envir = .GlobalEnv)
  on.exit({
    if (is.null(old_shiny)) {
      if (exists("shiny.env", envir = .GlobalEnv)) rm("shiny.env", envir = .GlobalEnv)
    } else {
      assign("shiny.env", old_shiny, envir = .GlobalEnv)
    }
  }, add = TRUE)

  missing_path <- tempfile(pattern = "no_such_file_", fileext = ".csv")
  # ensure it really doesn't exist
  if (file.exists(missing_path)){file.remove(missing_path)}

  expect_error(
    tas_import(missing_path),
    regexp = "Input spreadsheet not found",
    fixed = FALSE
  )
})


################################################################################


test_that("shiny path: missing file -> calls spsComps::shinyCatch (mocked) and triggers the expected stop", {
  old_shiny <- if (exists("shiny.env", envir = .GlobalEnv)) get("shiny.env", envir = .GlobalEnv) else NULL
  assign("shiny.env", TRUE, envir = .GlobalEnv)

  on.exit({
    if (is.null(old_shiny)) {
      if (exists("shiny.env", envir = .GlobalEnv)) rm("shiny.env", envir = .GlobalEnv)
    } else {
      assign("shiny.env", old_shiny, envir = .GlobalEnv)
    }
  }, add = TRUE)

  missing_path <- tempfile(pattern = "no_such_file_", fileext = ".csv")
  if (file.exists(missing_path)) file.remove(missing_path)

  # container for capture
  capture <- list(called = FALSE, expr = NULL, position = NULL)

  # Mock that *evaluates* the expr (so stop() happens inside the mock).
  # This simulates a shinyCatch that evaluates the provided expression and therefore triggers the error.
  mock_shinyCatch_eval <- function(expr, position) {
    capture$called <<- TRUE
    # capture unevaluated expression for inspection
    capture$expr <<- substitute(expr)
    capture$position <<- position
    # evaluate the expression so stop(...) is triggered here (causes an error)
    eval(expr)
    invisible(TRUE)
  }

  # Scoped mock for spsComps::shinyCatch
  local_mocked_bindings(
    shinyCatch = mock_shinyCatch_eval,
    .package = "spsComps"
  )

  # Because mock evaluates the expr, expect an error with message from stop(...)
  expect_error(
    tas_import(missing_path),
    regexp = "Input spreadsheet not found",
    ignore.case = TRUE
  )

  # ensure the mock was called
  expect_true(isTRUE(capture$called), info = "spsComps::shinyCatch should have been called when shiny.env = TRUE")
  expect_equal(capture$position, "top-center")
  # inspect the captured expression to ensure it was a stop call
  expect_true(is.call(capture$expr) && as.character(capture$expr[[1]]) == "stop")
})


################################################################################


test_that("Imported NGS object types are correct",{
  # create a temporary CSV with some NA values
  tmp <- tempfile(fileext = ".csv")
  dat <- data.frame(
    id = 1:3,
    name = c("Alice", NA, "Charlie"),
    score = c(NA, 85, 92),
    stringsAsFactors = FALSE
  )
  write.csv(dat, tmp, row.names = FALSE, na = "") # write with empty strings for NA for easy round-trip

  NGS <- tas_import(tmp)
  expect_type(NGS,"list")
  expect_true(class(NGS$Input) == "data.frame")
  expect_true(class(NGS$Config) == "list")
})



################################################################################
################################# TAS_CHECK ####################################
################################################################################


test_that("tas_check exists", {
  skip_if_not(exists("tas_check", mode = "function"), "tas_check() not found in this session")
})

# Helper: build a minimal, valid Input.DataFrame with specified filenames
make_valid_input_df <- function(forward_name = "fwd.fastq", reverse_name = "rev.fastq") {
  cols <- c("SampleName",
            "ForwardFASTQFileName",
            "ReverseFASTQFileName",
            "MergedFASTQFileName",
            "ForwardExtensionType",
            "ForwardExtension",
            "ForwardPrimer",
            "ReverseExtensionType",
            "ReverseExtension",
            "ReversePrimer",
            "ReferenceSequence",
            "AmpliconLength",
            "MaxDeletion",
            "MaxInsertion",
            "InsertStart",
            "InsertEnd",
            "Antibody",
            "FR1Start",
            "CDR1Start",
            "FR2Start",
            "CDR2Start",
            "FR3Start",
            "CDR3Start",
            "FR4Start")

  # create one-row dataframe with default values
  df <- as.data.frame(matrix(NA, nrow = 1, ncol = length(cols)))
  colnames(df) <- cols
  df$SampleName <- "S1"
  df$ForwardFASTQFileName <- forward_name
  df$ReverseFASTQFileName <- reverse_name
  df$MergedFASTQFileName <- "merged.fastq"
  df$ForwardExtensionType <- ""
  df$ForwardExtension <- ""
  df$ForwardPrimer <- "ATG"
  df$ReverseExtensionType <- ""
  df$ReverseExtension <- ""
  df$ReversePrimer <- "CTA"
  df$ReferenceSequence <- "ATGCAGGTAGTACCGACTTAG"
  df$AmpliconLength <- 300
  df$MaxDeletion <- 25
  df$MaxInsertion <- 25
  df$InsertStart <- 1
  df$InsertEnd <- 300
  df$Antibody <- "Yes"
  df$FR1Start <- 1
  df$CDR1Start <- 10
  df$FR2Start <- 50
  df$CDR2Start <- 70
  df$FR3Start <- 120
  df$CDR3Start <- 140
  df$FR4Start <- 200
  return(df)
}

# Minimal valid Config.List
make_valid_config <- function() {
  list(
    WorkingDirectory = tempdir(),
    OperatingSystem = "MacOS",
    nCores = 1,
    merge.reads = TRUE,
    measure.shm = FALSE,
    sequence.alignment.count = 10,
    read.frequency.limit = 0.001,
    protein.mutations = FALSE,
    PhyloTree = FALSE,
    dna.repair.pathways = FALSE,
    multicore = FALSE
  )
}




test_that("successful run when all inputs present (shiny.env = FALSE)", {
  skip_if_not(exists("tas_check", mode = "function"))
  tmp <- tempdir()
  # Create an input directory and Input.csv path (function expects substring "Input.csv")
  input_dir <- file.path(tmp, "input_dir_for_tas_check")
  dir.create(input_dir, showWarnings = FALSE)
  InputFilePath <- file.path(input_dir, "MyInput.Input.csv")  # note "Input.csv" substring included
  writeLines("dummy", InputFilePath)

  # create FASTQ files in input dir
  fwd <- "fwd.fastq"
  rev <- "rev.fastq"
  file.create(file.path(input_dir, fwd))
  file.create(file.path(input_dir, rev))

  Input.DataFrame <- make_valid_input_df(forward_name = fwd, reverse_name = rev)
  Config.List <- make_valid_config()
  WD <- input_dir

  expect_silent(tas_check(InputFilePath = InputFilePath,
                          Input.DataFrame = Input.DataFrame,
                          Config.List = Config.List,
                          shiny.env = FALSE,
                          shiny.fileTable = character(),
                          WD = WD))
})




test_that("missing required columns in Input.DataFrame triggers error", {
  skip_if_not(exists("tas_check", mode = "function"))
  # Make an input df missing one required column
  df <- make_valid_input_df()
  df$ForwardFASTQFileName <- NULL  # remove required column

  Config.List <- make_valid_config()
  # Provide an InputFilePath that includes "Input.csv" (directory not used in this branch)
  InputFilePath <- file.path(tempdir(), "SomeInput.Input.csv")
  writeLines("x", InputFilePath)

  WD <- paste0(tempdir(), "/")

  expect_error(
    tas_check(InputFilePath = InputFilePath,
              Input.DataFrame = df,
              Config.List = Config.List,
              shiny.env = FALSE,
              shiny.fileTable = NULL,
              WD = WD),
    regexp = "The following columns were not found in the input file",
    fixed = FALSE
  )
})




test_that("missing entries in Config.List triggers error", {
  skip_if_not(exists("tas_check", mode = "function"))
  df <- make_valid_input_df()
  cfg <- make_valid_config()
  cfg$merge.reads <- NULL  # remove one config entry

  InputFilePath <- file.path(tempdir(), "SomeInput.Input.csv")
  writeLines("x", InputFilePath)
  WD <- paste0(tempdir(), "/")

  expect_error(
    tas_check(InputFilePath = InputFilePath,
              Input.DataFrame = df,
              Config.List = cfg,
              shiny.env = FALSE,
              shiny.fileTable = NULL,
              WD = WD),
    regexp = "The following columns were not found in the input file",
    fixed = FALSE
  )
})




test_that("sequence.alignment.count validation: not single number", {
  skip_if_not(exists("tas_check", mode = "function"))
  df <- make_valid_input_df()
  cfg <- make_valid_config()
  cfg$sequence.alignment.count <- c(1, 2)  # not single

  InputFilePath <- file.path(tempdir(), "SomeInput.Input.csv")
  writeLines("x", InputFilePath)
  WD <- paste0(tempdir(), "/")

  expect_error(
    tas_check(InputFilePath = InputFilePath,
              Input.DataFrame = df,
              Config.List = cfg,
              shiny.env = FALSE,
              shiny.fileTable = NULL,
              WD = WD),
    regexp = "sequence.alignment.count must be a single number"
  )
})




test_that("sequence.alignment.count validation: negative", {
  skip_if_not(exists("tas_check", mode = "function"))
  df <- make_valid_input_df()
  cfg <- make_valid_config()
  cfg$sequence.alignment.count <- -1

  InputFilePath <- file.path(tempdir(), "SomeInput.Input.csv")
  writeLines("x", InputFilePath)
  WD <- paste0(tempdir(), "/")

  expect_error(
    tas_check(InputFilePath = InputFilePath,
              Input.DataFrame = df,
              Config.List = cfg,
              shiny.env = FALSE,
              shiny.fileTable = NULL,
              WD = WD),
    regexp = "sequence.alignment.count must be > 0"
  )
})




test_that("sequence.alignment.count validation: non-integer", {
  skip_if_not(exists("tas_check", mode = "function"))
  df <- make_valid_input_df()
  cfg <- make_valid_config()
  cfg$sequence.alignment.count <- 1.5

  InputFilePath <- file.path(tempdir(), "SomeInput.Input.csv")
  writeLines("x", InputFilePath)
  WD <- paste0(tempdir(), "/")

  expect_error(
    tas_check(InputFilePath = InputFilePath,
              Input.DataFrame = df,
              Config.List = cfg,
              shiny.env = FALSE,
              shiny.fileTable = NULL,
              WD = WD),
    regexp = "sequence.alignment.count must be a whole number"
  )
})




test_that("read.frequency.limit validation: not single number", {
  skip_if_not(exists("tas_check", mode = "function"))
  df <- make_valid_input_df()
  cfg <- make_valid_config()
  cfg$read.frequency.limit <- c(1, 2)

  InputFilePath <- file.path(tempdir(), "SomeInput.Input.csv")
  writeLines("x", InputFilePath)
  WD <- paste0(tempdir(), "/")

  expect_error(
    tas_check(InputFilePath = InputFilePath,
              Input.DataFrame = df,
              Config.List = cfg,
              shiny.env = FALSE,
              shiny.fileTable = NULL,
              WD = WD),
    regexp = "read.frequency.limit must be a single number"
  )
})




test_that("read.frequency.limit validation: negative", {
  skip_if_not(exists("tas_check", mode = "function"))
  df <- make_valid_input_df()
  cfg <- make_valid_config()
  cfg$read.frequency.limit <- -5

  InputFilePath <- file.path(tempdir(), "SomeInput.Input.csv")
  writeLines("x", InputFilePath)
  WD <- paste0(tempdir(), "/")

  expect_error(
    tas_check(InputFilePath = InputFilePath,
              Input.DataFrame = df,
              Config.List = cfg,
              shiny.env = FALSE,
              shiny.fileTable = NULL,
              WD = WD),
    regexp = "read.frequency.limit must be > 0"
  )
})




test_that("shiny.env = TRUE: succeeds when shiny.fileTable contains uploaded filenames", {
  skip_if_not(exists("tas_check", mode = "function"))
  # Simulate upload table with names matching filenames
  df <- make_valid_input_df(forward_name = "upload_fwd.fastq", reverse_name = "upload_rev.fastq")
  cfg <- make_valid_config()

  shiny.fileTable <- data.frame(name = c("upload_fwd.fastq", "upload_rev.fastq"),
                                stringsAsFactors = FALSE)

  InputFilePath <- file.path(tempdir(), "SomeInput.Input.csv")
  writeLines("x", InputFilePath)

  WD <- paste0(tempdir(), "/")

  expect_silent(
    tas_check(InputFilePath = InputFilePath,
              Input.DataFrame = df,
              Config.List = cfg,
              shiny.env = TRUE,
              shiny.fileTable = shiny.fileTable,
              WD = WD)
  )
})



test_that("shiny.env = FALSE: missing files in input directory triggers error", {
  skip_if_not(exists("tas_check", mode = "function"))
  tmp <- tempdir()
  input_dir <- file.path(tmp, "input_dir_missing_files")
  dir.create(input_dir, showWarnings = FALSE)
  InputFilePath <- file.path(input_dir, "MyInput.Input.csv")
  writeLines("dummy", InputFilePath)

  # Do NOT create the actual files in input_dir
  fwd <- "missing_fwd.fastq"
  rev <- "missing_rev.fastq"

  df <- make_valid_input_df(forward_name = fwd, reverse_name = rev)
  cfg <- make_valid_config()
  WD <- paste0(tempdir(), "/")

  expect_error(
    tas_check(InputFilePath = InputFilePath,
              Input.DataFrame = df,
              Config.List = cfg,
              shiny.env = FALSE,
              shiny.fileTable = NULL,
              WD = WD),
    regexp = "These files were not found",
    fixed = FALSE
  )
})


###############################################################################


test_that("shiny.env = TRUE: mocked spsComps::shinyCatch is called with position='top-center' and rethrows error", {
  skip_if_not(exists("tas_check", mode = "function"), "tas_check() not found in this session")

  # Prepare inputs that will trigger the missing-files shiny branch
  df <- make_valid_input_df(forward_name = "upload_a.fastq", reverse_name = "upload_b.fastq")
  cfg <- make_valid_config()
  shiny.fileTable <- data.frame(name = c("other1.fastq", "other2.fastq"),
                                stringsAsFactors = FALSE)
  InputFilePath <- file.path(tempdir(), "SomeInput.Input.csv")
  writeLines("x", InputFilePath)
  WD <- paste0(tempdir(), "/")

  capture <- list(called = FALSE, expr = NULL, position = NULL)

  # Mock that *evaluates* the expr (so stop() happens inside the mock).
  # This simulates a shinyCatch that evaluates the provided expression and therefore triggers the error.
  mock_shinyCatch_eval <- function(expr, position) {
    capture$called <<- TRUE
    # capture unevaluated expression for inspection
    capture$expr <<- substitute(expr)
    capture$position <<- position
    # evaluate the expression so stop(...) is triggered here (causes an error)
    eval(expr)
    invisible(TRUE)
  }

  # Scoped mock for spsComps::shinyCatch
  local_mocked_bindings(
    shinyCatch = mock_shinyCatch_eval,
    .package = "spsComps"
  )

  # Because mock evaluates the expr, expect an error with message from stop(...)
  expect_error(
    tas_check(InputFilePath = InputFilePath,
              Input.DataFrame = df,
              Config.List = cfg,
              shiny.env = TRUE,
              shiny.fileTable = shiny.fileTable,
              WD = WD
              ),
  regexp = "These files were not found|were not found",
  ignore.case = TRUE
)

  # ensure the mock was called
  expect_true(isTRUE(capture$called), info = "spsComps::shinyCatch should have been called when shiny.env = TRUE")
  expect_equal(capture$position, "top-center")
  # inspect the captured expression to ensure it was a stop call
  expect_true(is.call(capture$expr) && as.character(capture$expr[[1]]) == "stop")

})
