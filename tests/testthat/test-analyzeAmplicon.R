


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

test_that("analyzeAmplicon emits structured progress events", {
  events <- list()
  callback <- function(event) {
    events[[length(events) + 1L]] <<- event
  }

  expect_output(expect_no_error(analyzeAmplicon(test.settings, progress_callback = callback)))

  labels <- vapply(events, `[[`, character(1), "label")
  expect_true(all(c(
    "Reading FASTQ",
    "Filtering reads",
    "Labeling mutations",
    "Measuring mutations"
  ) %in% labels))
  expect_true(any(labels %in% c("Binning sequences", "Binning UMIs")))
  expect_true(all(vapply(events, `[[`, integer(1), "total") == 5L))
  expect_true(all(vapply(events, `[[`, character(1), "sample_name") == test.settings@Name))
})



test_that("analyzeAmplicon does not fail with low-diversity inputs", {
  sett <- test.settings
  sett@MergedFASTQPath <- file.path(tempdir(), "test-merged-lowdiv.fastq.gz")
  expect_output(expect_no_error(amp <- analyzeAmplicon(sett)))


  if (file.exists(file.path(tempdir(), "test-merged-fem.fastq.gz"))) {file.remove(file.path(tempdir(), "test-merged-fem.fastq.gz"))}
  writeFastq(makeFakeSRQ(Biostrings::DNAStringSet(rep("GCTgaCCGTAAAACGACGGCCAGTGTTCAACTGGTGGAAAGCGGCGGTGCTCTGGTACAACCGGGCGGTAGTCTGCGCCTGAGCTGTGCCGCAAGCGGTTTCCCAGTCAACCGCTACTCTATGCGTTGGTATCGCCAGGCGCCTGGTAAAGAACGTGAATGGGTTGCCGGCATGAGCAGTGCGGGCGATCGTTCTAGTTACGAGGACTCTGTTAAAGGTCGTTTTACAATTAGCCGTGATGATGCGCGCAATACCGTGTATCTGCAAATGAACAGTCTGAAGCCGGAGGACACCGCAGTATATTATTGCAATGTCAACGTGGGGTTTGAATATTGGGGCCAGGGGACTCAGGTGACGGTGAGCTCTGTCATAGCTGTTTCCTGAGCGAATTAGAG", 33))),
             file.path(tempdir(), "test-merged-fem.fastq.gz"))
  sett@MergedFASTQPath <- file.path(tempdir(), "test-merged-fem.fastq.gz")
  expect_output(expect_error(buildSequenceTable(sett), "No reads"))

  if (file.exists(file.path(tempdir(), "test-merged-fpm.fastq.gz"))) {file.remove(file.path(tempdir(), "test-merged-fpm.fastq.gz"))}
  writeFastq(makeFakeSRQ(Biostrings::DNAStringSet(rep("GCTAGCCGTAAAACGAgcGCCAGTGTTCAACTGGTGGAAAGCGGCGGTGCTCTGGTACAACCGGGCGGTAGTCTGCGCCTGAGCTGTGCCGCAAGCGGTTTCCCAGTCAACCGCTACTCTATGCGTTGGTATCGCCAGGCGCCTGGTAAAGAACGTGAATGGGTTGCCGGCATGAGCAGTGCGGGCGATCGTTCTAGTTACGAGGACTCTGTTAAAGGTCGTTTTACAATTAGCCGTGATGATGCGCGCAATACCGTGTATCTGCAAATGAACAGTCTGAAGCCGGAGGACACCGCAGTATATTATTGCAATGTCAACGTGGGGTTTGAATATTGGGGCCAGGGGACTCAGGTGACGGTGAGCTCTGTCATAGCTGTTTCCTGAGCGAATTAGAG", 33))),
             file.path(tempdir(), "test-merged-fpm.fastq.gz"))
  sett@MergedFASTQPath <- file.path(tempdir(), "test-merged-fpm.fastq.gz")
  expect_output(expect_error(buildSequenceTable(sett), "No reads"))

  if (file.exists(file.path(tempdir(), "test-merged-rem.fastq.gz"))) {file.remove(file.path(tempdir(), "test-merged-rem.fastq.gz"))}
  writeFastq(makeFakeSRQ(Biostrings::DNAStringSet(rep("GCTAGCCGTAAAACGACGGCCAGTGTTCAACTGGTGGAAAGCGGCGGTGCTCTGGTACAACCGGGCGGTAGTCTGCGCCTGAGCTGTGCCGCAAGCGGTTTCCCAGTCAACCGCTACTCTATGCGTTGGTATCGCCAGGCGCCTGGTAAAGAACGTGAATGGGTTGCCGGCATGAGCAGTGCGGGCGATCGTTCTAGTTACGAGGACTCTGTTAAAGGTCGTTTTACAATTAGCCGTGATGATGCGCGCAATACCGTGTATCTGCAAATGAACAGTCTGAAGCCGGAGGACACCGCAGTATATTATTGCAATGTCAACGTGGGGTTTGAATATTGGGGCCAGGGGACTCAGGTGACGGTGAGCTCTGTCATAGCTGTTTCCTGAGCGAATatGAG", 33))),
             file.path(tempdir(), "test-merged-rem.fastq.gz"))
  sett@MergedFASTQPath <- file.path(tempdir(), "test-merged-rem.fastq.gz")
  expect_output(expect_error(buildSequenceTable(sett), "No reads"))

  if (file.exists(file.path(tempdir(), "test-merged-rpm.fastq.gz"))) {file.remove(file.path(tempdir(), "test-merged-rpm.fastq.gz"))}
  writeFastq(makeFakeSRQ(Biostrings::DNAStringSet(rep("GCTAGCCGTAAAACGACGGCCAGTGTTCAACTGGTGGAAAGCGGCGGTGCTCTGGTACAACCGGGCGGTAGTCTGCGCCTGAGCTGTGCCGCAAGCGGTTTCCCAGTCAACCGCTACTCTATGCGTTGGTATCGCCAGGCGCCTGGTAAAGAACGTGAATGGGTTGCCGGCATGAGCAGTGCGGGCGATCGTTCTAGTTACGAGGACTCTGTTAAAGGTCGTTTTACAATTAGCCGTGATGATGCGCGCAATACCGTGTATCTGCAAATGAACAGTCTGAAGCCGGAGGACACCGCAGTATATTATTGCAATGTCAACGTGGGGTTTGAATATTGGGGCCAGGGGACTCAGGTGACGGTGAGCTCTGTCATAGCTtgTTCCTGAGCGAATTAGAG", 33))),
             file.path(tempdir(), "test-merged-rpm.fastq.gz"))
  sett@MergedFASTQPath <- file.path(tempdir(), "test-merged-rpm.fastq.gz")
  expect_output(expect_error(buildSequenceTable(sett), "No reads"))

})

test_that("Summarize works as expected", {
  if (!dir.exists(file.path(tempdir(), "single_export_test"))) dir.create(file.path(tempdir(), "single_export_test"))
  expect_output(expect_no_error(amp <- analyzeAmplicon(test.settings)))
  summary_list <- Summarize(results = amp, export = TRUE, path = file.path(tempdir(), "single_export_test"))
  expect_true(is(summary_list, "list"))
  expect_equal(length(summary_list), 2)
  expect_equal(names(summary_list), c("Graphs", "Tables"))
  expect_equal(names(summary_list$Graphs), c("DNA Mutation Distribution",
                                             "DNA Mutation Distribution labeled cytosines",
                                             "AA Mutation Distribution",
                                             "AA Mutation Distribution labeled cytosines",
                                             "Mutation at cytosines boxplot",
                                             "All AA Mutations",
                                             "DNA Repair Types",
                                             "Mutations per read histogram"))
  expect_equal(names(summary_list$Tables), c("Sequence Table",
                                             "Read Counts",
                                             "All DNA Mutations",
                                             "Cytosine DNA Mutations",
                                             "Non-Cytosine DNA Mutations",
                                             "Motif Sums",
                                             "All Protein Mutations",
                                             "Protein Mutation Matrix",
                                             "WRCH Table",
                                             "WRCY Table",
                                             "DNA Repair Types"))
  expect_all_true(file.exists(file.path(tempdir(), "single_export_test/test/", stringr::str_c(c("DNA Mutation Distribution",
                                                                                        "DNA Mutation Distribution labeled cytosines",
                                                                                        "AA Mutation Distribution",
                                                                                        "AA Mutation Distribution labeled cytosines",
                                                                                        "Mutation at cytosines boxplot",
                                                                                        "All AA Mutations",
                                                                                        "DNA Repair Types",
                                                                                        "Mutations per read histogram"),
                                                                                      ".pdf"))))
  expect_all_true(file.exists(file.path(tempdir(), "single_export_test/test/", stringr::str_c(c("Sequence Table",
                                                                                        "Read Counts",
                                                                                        "All DNA Mutations",
                                                                                        "Cytosine DNA Mutations",
                                                                                        "Non-Cytosine DNA Mutations",
                                                                                        "Motif Sums",
                                                                                        "All Protein Mutations",
                                                                                        "Protein Mutation Matrix",
                                                                                        "WRCH Table",
                                                                                        "WRCY Table",
                                                                                        "DNA Repair Types"),
                                                                                      ".csv"))))
})
