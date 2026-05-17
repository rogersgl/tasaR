





test_that("pairedAnalyzeAmplicon works with a homozygous reference and various input types.", {
  path.ctrl <- file.path(tempdir(), "pair-homo-ctrl-merged.fastq.gz")
  path.expt <- file.path(tempdir(), "pair-homo-expt-merged.fastq.gz")
  expect_output(pas <- pairedAnalyzeAmplicon(file.control = path.ctrl,
                                             file.experimental = path.expt,
                                             settings = full.settings.ctrl,
                                             with.nuclease = TRUE,
                                             gRNA.seq = "TACAAGACCCGCGCCGAGGT",
                                             manual.cut.site = NULL))

  expect_equal(length(slotNames(pas)), 2)
  expect_equal(unlist(getSettings(pas@Control)$ReferenceSequence), unlist(getSettings(pas@Experimental)$ReferenceSequence))
  expect_equal(nrow(getSequenceTable(pas@Control)), 1)
  expect_equal(nrow(getSequenceTable(pas@Experimental)), 4)
  expect_equal(getSequenceTable(pas@Control)$Count, 2)
  expect_equal(getSequenceTable(pas@Experimental)$Count, c(1, 1, 1, 1))
  expect_equal(getSequenceTable(pas@Control)$Indels, "WT")
  expect_equal(getSequenceTable(pas@Experimental)$Indels, c("WT", "-1", "+1", "-5"))
  expect_equal(getMutationTypes(pas@Control), data.frame(WT = 100, NHEJ = 0, MMEJ = 0, BaseChange = 0, IndelBaseChange = 0, Other = 0))
  expect_equal(getMutationTypes(pas@Experimental), data.frame(WT = 25, NHEJ = 50, MMEJ = 25, BaseChange = 0, IndelBaseChange = 0, Other = 0))
  expect_equal(getMutationDistributionDNA(pas@Control), data.frame(Position = 1:360, MutationFrequency = rep(0, 360)))
  expect_equal(getMutationDistributionDNA(pas@Experimental), data.frame(Position = c(1:331, 331.5, 332:360), MutationFrequency = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 25, 25, 25, 50, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)))

  # with a data.frame input
  suppressMessages(expect_output(pas.l <- pairedAnalyzeAmplicon(path.ctrl,
                                                                path.expt,
                                                                make.pair.set.df("ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTG"),
                                                                with.nuclease = TRUE,
                                                                gRNA.seq = "TACAAGACCCGCGCCGAGGT",
                                                                manual.cut.site = NULL
                                                                )))
  expect_equal(pas.l, pas)

  # with a manual input
  expect_output(pas.m <- pairedAnalyzeAmplicon(path.ctrl, path.expt,
                                             fwd.primer = "GTAAAACGACGGCCAGT",
                                             fwd.exten.type = "Barcode",
                                             fwd.exten.seq = "GCTAGCC",
                                             rev.primer = "CAGGAAACAGCTATGAC",
                                             rev.exten.type = "UMI",
                                             rev.exten.seq = "NNNYRNNNYRNN",
                                             with.nuclease = TRUE,
                                             gRNA.seq = "TACAAGACCCGCGCCGAGGT",
                                             manual.cut.site = NULL))
  expect_equal(pas.m, pas)
})



test_that("pairedAnalyzeAmplicon works with a heterozygous reference sequence and various input types.", {
  path.ctrl <- file.path(tempdir(), "pair-het-ctrl-merged.fastq.gz")
  path.expt <- file.path(tempdir(), "pair-het-expt-merged.fastq.gz")
  config.ctrl <- full.settings.ctrl
  expect_output(pas <- pairedAnalyzeAmplicon(path.ctrl,
                                             path.expt,
                                             config.ctrl,
                                             with.nuclease = TRUE,
                                             gRNA.seq = "TACAAGACCCGCGCCGAGGT",
                                             manual.cut.site = NULL))

  expect_equal(length(slotNames(pas)), 2)
  expect_equal(unlist(getSettings(pas@Control)$ReferenceSequence), unlist(getSettings(pas@Experimental)$ReferenceSequence))
  expect_equal(nrow(getSequenceTable(pas@Control)), 2)
  expect_equal(nrow(getSequenceTable(pas@Experimental)), 8)
  expect_equal(getSequenceTable(pas@Control)$Count, c(1, 1))
  expect_equal(getSequenceTable(pas@Experimental)$Count, c(1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(getSequenceTable(pas@Control)$Indels, c("WT", "WT"))
  expect_equal(getSequenceTable(pas@Experimental)$Indels, c("WT", "-1", "+1", "-5", "WT", "-1", "+1", "-5"))
  expect_equal(getMutationTypes(pas@Control), data.frame(WT = 100, NHEJ = 0, MMEJ = 0, BaseChange = 0, IndelBaseChange = 0, Other = 0))
  expect_equal(getMutationTypes(pas@Experimental), data.frame(WT = 25, NHEJ = 50, MMEJ = 25, BaseChange = 0, IndelBaseChange = 0, Other = 0))
  expect_equal(getMutationDistributionDNA(pas@Control), data.frame(Position = 1:360, MutationFrequency = rep(0, 360)))
  expect_equal(getMutationDistributionDNA(pas@Experimental), data.frame(Position = c(1:331, 331.5, 332:360), MutationFrequency = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 25, 25, 25, 50, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)))

})
