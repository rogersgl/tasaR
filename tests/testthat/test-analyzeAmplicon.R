


test_that("analyzeAmplicon parses input data types as expected.", {
  expect_output(expect_no_error(amp.s4 <- analyzeAmplicon(test.settings)))
  expect_output(expect_message(expect_no_error(amp.r <- analyzeAmplicon(make.set.df()))))
  expect_message(expect_error(analyzeAmplicon(getSettings(test.settings))))
  expect_message(expect_no_error(fset <- readSettings(file.path(tempdir(), "settings.csv"), row = 2)), regexp = ".csv file")
  expect_output(amp.f <- analyzeAmplicon(fset))

  amp.a <- expect_output(new("AmpliconSequencing", Sequences = test.seq,
                                                   Mutations = measureMutations(test.seq, test.settings),
                                                   MutationTypes = measureDNArepair(test.seq, test.settings),
                                                   Settings = test.settings))

  amp.a@Sequences@Supplemental$IDs <- rep(NA, length(amp.a@Sequences@Supplemental$IDs))
  amp.s4@Sequences@Supplemental$IDs <- rep(NA, length(amp.s4@Sequences@Supplemental$IDs))
  amp.r@Sequences@Supplemental$IDs <- rep(NA, length(amp.r@Sequences@Supplemental$IDs))
  amp.f@Sequences@Supplemental$IDs <- rep(NA, length(amp.f@Sequences@Supplemental$IDs))

  expect_equal(amp.s4, amp.a)
  expect_equal(amp.r, amp.a)
  expect_equal(amp.f, amp.a)
})
