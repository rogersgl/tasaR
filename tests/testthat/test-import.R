
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
