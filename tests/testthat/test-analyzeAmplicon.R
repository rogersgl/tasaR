


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



test_that("analyzeAmplicon does not fail with low-diversity inputs", {
  sett <- test.settings
  sett@MergedFASTQPath <- file.path(tempdir(), "test-merged-lowdiv.fastq.gz")
  expect_output(expect_no_error(amp <- analyzeAmplicon(sett)))


  if (file.exists(file.path(tempdir(), "test-merged-fem.fastq.gz"))) {file.remove(file.path(tempdir(), "test-merged-fem.fastq.gz"))}
  writeFastq(makeFakeSRQ(DNAStringSet(rep("GCTgaCCGTAAAACGACGGCCAGTGTTCAACTGGTGGAAAGCGGCGGTGCTCTGGTACAACCGGGCGGTAGTCTGCGCCTGAGCTGTGCCGCAAGCGGTTTCCCAGTCAACCGCTACTCTATGCGTTGGTATCGCCAGGCGCCTGGTAAAGAACGTGAATGGGTTGCCGGCATGAGCAGTGCGGGCGATCGTTCTAGTTACGAGGACTCTGTTAAAGGTCGTTTTACAATTAGCCGTGATGATGCGCGCAATACCGTGTATCTGCAAATGAACAGTCTGAAGCCGGAGGACACCGCAGTATATTATTGCAATGTCAACGTGGGGTTTGAATATTGGGGCCAGGGGACTCAGGTGACGGTGAGCTCTGTCATAGCTGTTTCCTGAGCGAATTAGAG", 33))),
             file.path(tempdir(), "test-merged-fem.fastq.gz"))
  sett@MergedFASTQPath <- file.path(tempdir(), "test-merged-fem.fastq.gz")
  expect_output(expect_error(buildSequenceTable(sett), "No reads"))

  if (file.exists(file.path(tempdir(), "test-merged-fpm.fastq.gz"))) {file.remove(file.path(tempdir(), "test-merged-fpm.fastq.gz"))}
  writeFastq(makeFakeSRQ(DNAStringSet(rep("GCTAGCCGTAAAACGAgcGCCAGTGTTCAACTGGTGGAAAGCGGCGGTGCTCTGGTACAACCGGGCGGTAGTCTGCGCCTGAGCTGTGCCGCAAGCGGTTTCCCAGTCAACCGCTACTCTATGCGTTGGTATCGCCAGGCGCCTGGTAAAGAACGTGAATGGGTTGCCGGCATGAGCAGTGCGGGCGATCGTTCTAGTTACGAGGACTCTGTTAAAGGTCGTTTTACAATTAGCCGTGATGATGCGCGCAATACCGTGTATCTGCAAATGAACAGTCTGAAGCCGGAGGACACCGCAGTATATTATTGCAATGTCAACGTGGGGTTTGAATATTGGGGCCAGGGGACTCAGGTGACGGTGAGCTCTGTCATAGCTGTTTCCTGAGCGAATTAGAG", 33))),
             file.path(tempdir(), "test-merged-fpm.fastq.gz"))
  sett@MergedFASTQPath <- file.path(tempdir(), "test-merged-fpm.fastq.gz")
  expect_output(expect_error(buildSequenceTable(sett), "No reads"))

  if (file.exists(file.path(tempdir(), "test-merged-rem.fastq.gz"))) {file.remove(file.path(tempdir(), "test-merged-rem.fastq.gz"))}
  writeFastq(makeFakeSRQ(DNAStringSet(rep("GCTAGCCGTAAAACGACGGCCAGTGTTCAACTGGTGGAAAGCGGCGGTGCTCTGGTACAACCGGGCGGTAGTCTGCGCCTGAGCTGTGCCGCAAGCGGTTTCCCAGTCAACCGCTACTCTATGCGTTGGTATCGCCAGGCGCCTGGTAAAGAACGTGAATGGGTTGCCGGCATGAGCAGTGCGGGCGATCGTTCTAGTTACGAGGACTCTGTTAAAGGTCGTTTTACAATTAGCCGTGATGATGCGCGCAATACCGTGTATCTGCAAATGAACAGTCTGAAGCCGGAGGACACCGCAGTATATTATTGCAATGTCAACGTGGGGTTTGAATATTGGGGCCAGGGGACTCAGGTGACGGTGAGCTCTGTCATAGCTGTTTCCTGAGCGAATatGAG", 33))),
             file.path(tempdir(), "test-merged-rem.fastq.gz"))
  sett@MergedFASTQPath <- file.path(tempdir(), "test-merged-rem.fastq.gz")
  expect_output(expect_error(buildSequenceTable(sett), "No reads"))

  if (file.exists(file.path(tempdir(), "test-merged-rpm.fastq.gz"))) {file.remove(file.path(tempdir(), "test-merged-rpm.fastq.gz"))}
  writeFastq(makeFakeSRQ(DNAStringSet(rep("GCTAGCCGTAAAACGACGGCCAGTGTTCAACTGGTGGAAAGCGGCGGTGCTCTGGTACAACCGGGCGGTAGTCTGCGCCTGAGCTGTGCCGCAAGCGGTTTCCCAGTCAACCGCTACTCTATGCGTTGGTATCGCCAGGCGCCTGGTAAAGAACGTGAATGGGTTGCCGGCATGAGCAGTGCGGGCGATCGTTCTAGTTACGAGGACTCTGTTAAAGGTCGTTTTACAATTAGCCGTGATGATGCGCGCAATACCGTGTATCTGCAAATGAACAGTCTGAAGCCGGAGGACACCGCAGTATATTATTGCAATGTCAACGTGGGGTTTGAATATTGGGGCCAGGGGACTCAGGTGACGGTGAGCTCTGTCATAGCTtgTTCCTGAGCGAATTAGAG", 33))),
             file.path(tempdir(), "test-merged-rpm.fastq.gz"))
  sett@MergedFASTQPath <- file.path(tempdir(), "test-merged-rpm.fastq.gz")
  expect_output(expect_error(buildSequenceTable(sett), "No reads"))

})
