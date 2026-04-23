# ---------------
# makeSettingsCSV
# ---------------

test_that("Template .csv file can be created.", {

  # default name
  expect_message(makeSettingsCSV(tempdir()), "Created file")
  suppressMessages(makeSettingsCSV(tempdir()))
  expect_true(file.exists(file.path(tempdir(),"settings.csv")))

  #specified name
  expect_message(makeSettingsCSV(file.path(tempdir(), "test.csv")), "Created file")
  suppressMessages(makeSettingsCSV(file.path(tempdir(), "test.csv")))
  expect_true(file.exists(file.path(tempdir(),"test.csv")))

  # check column names match what the S4 object expects
  cnames <- colnames(read.csv(file.path(tempdir(),"settings.csv")))
  snames <- slotNames("tas.object.settings")
  expect_all_true(cnames[1:(length(snames)-1)] == snames[-length(snames)])
  expect_all_true(cnames[length(snames):length(cnames)] == names(new("tas.object.settings")@AntibodyRegions))
})


# ------------
# readSettings
# ------------

test_that("Settings can be read from a .csv file.", {
  dfset <- make.set.df()
  write.csv(dfset, file = file.path(tempdir(), "settings.csv"), row.names = FALSE)
  expect_message(o <- readSettings(file.path(tempdir(), "settings.csv"), row = 2), regexp = ".csv file")
  # o <- suppressMessages(readSettings(file.path(tempdir(), "settings.csv"), row = 2))
  expect_true(class(o) == "tas.object.settings")
  expect_true(validObject(o))

  onames <- slotNames(o)
  snames <- slotNames("tas.object.settings")
  expect_all_true(onames == snames)
  expect_all_true(names(o@AntibodyRegions) == names(new("tas.object.settings")@AntibodyRegions))

  expect_error(suppressMessages((readSettings(file.path(tempdir(), "nofile.csv", row = 2)))), regexp = "File not found")
  expect_error(suppressMessages((readSettings(file.path(tempdir(), "settings.csv")))), regexp = "Row number NULL is invalid")
})

test_that("Settings can be read from an R object.", {
  dfset <- make.set.df()
  expect_message(readSettings(dfset), regexp = "R object")
  o <- suppressMessages(readSettings(dfset))
  expect_true(class(o) == "tas.object.settings")
  expect_true(validObject(o))

  onames <- slotNames(o)
  snames <- slotNames("tas.object.settings")
  expect_all_true(onames == snames)
  expect_all_true(names(o@AntibodyRegions) == names(new("tas.object.settings")@AntibodyRegions))
})

test_that("Settings can be built with buildSettings", {
  o <- buildSettings(Name = "test",
                     IsAntibody = TRUE,
                     MergedFASTQPath = file.path(tempdir(), "test-merged.fastq.gz"),
                     ReferenceSequence = "GTTCAACTGGTGGAAAGCGGCGGTGCTCTGGTACAACCGGGCGGTAGTCTGCGCCTGAGCTGTGCCGCAAGCGGTTTCCCAGTCAACCGCTACTCTATGCGTTGGTATCGCCAGGCGCCTGGTAAAGAACGTGAATGGGTTGCCGGCATGAGCAGTGCGGGCGATCGTTCTAGTTACGAGGACTCTGTTAAAGGTCGTTTTACAATTAGCCGTGATGATGCGCGCAATACCGTGTATCTGCAAATGAACAGTCTGAAGCCGGAGGACACCGCAGTATATTATTGCAATGTCAACGTGGGGTTTGAATATTGGGGCCAGGGGACTCAGGTGACGGTGAGCTCT",
                     ForwardExtensionType = "Barcode",
                     ForwardExtension = "GCTAGCC",
                     ForwardPrimer = "GTAAAACGACGGCCAGT",
                     ReverseExtensionType = "UMI",
                     ReverseExtension = "NNNYRNNNYRNN",
                     ReversePrimer = "CAGGAAACAGCTATGAC",
                     AmpliconLength = 395L,
                     InsertStart = 25L,
                     InsertEnd = -30L,
                     FR1Start = 1L,
                     CDR1Start = 73L,
                     FR2Start = 97L,
                     CDR2Start = 148L,
                     FR3Start = 172L,
                     CDR3Start = 286L,
                     FR4Start = 310L,
                     FR4End = 342L)
  expect_true(validObject(o))
  onames <- slotNames(o)
  snames <- slotNames("tas.object.settings")
  expect_all_true(onames == snames)
  expect_all_true(names(o@AntibodyRegions) == names(new("tas.object.settings")@AntibodyRegions))

})
