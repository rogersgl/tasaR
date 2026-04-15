# -------
# isEmpty
# -------

test_that("Validate tasAnalyzer extensions for function isEmpty.", {

  #tas.object.settings
  obj.set <- new("tas.object.settings", Name = "unknown",
                 IsAntibody = FALSE,
                 MeasureSHM = FALSE,
                 MeasureDNARepairTypes = FALSE,
                 MergedFASTQPath = "",
                 ReferenceSequence = "",
                 ForwardExtensionType = "",
                 ForwardExtension = "",
                 ForwardPrimer = "",
                 ReverseExtensionType = "",
                 ReverseExtension = "",
                 ReversePrimer = "",
                 AmpliconLength = NA_integer_,
                 InsertStart = NA_integer_,
                 InsertEnd = NA_integer_,
                 AntibodyRegions = c(FR1Start = NA_integer_,
                                     CDR1Start = NA_integer_,
                                     FR2Start = NA_integer_,
                                     CDR2Start = NA_integer_,
                                     FR3Start = NA_integer_,
                                     CDR3Start = NA_integer_,
                                     FR4Start = NA_integer_,
                                     FR4End = NA_integer_))

  expect_true(validObject(obj.set))
  expect_true(isEmpty(obj.set))
  obj.set@ForwardPrimer <- "ATGCAG"
  expect_true(validObject(obj.set))
  expect_false(isEmpty(obj.set))

  # tas.mutations
  obj.mut <- new("tas.mutations", AllMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                                  CytosineMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                                  NonCytosineMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                                  MotifSums = c(TotalPercentMutated = NA_real_,
                                                TotalPercentMutatedInMotif = NA_real_,
                                                TotalPercentMutatedInCytosine = NA_real_,
                                                MotifMutationAverage = NA_real_,
                                                MotifMutationFrequencyofTotal = NA_real_,
                                                MotifMutationNonCAverage = NA_real_,
                                                MotifMutationNonCFrequencyOfTotal = NA_real_,
                                                CytosineMutationAverage = NA_real_,
                                                NonCytosineMutationAverage = NA_real_,
                                                CytosineMutationFrequencyOfTotal = NA_real_,
                                                CytosineMutationFrequencyOfMotif = NA_real_))
  expect_true(validObject(obj.mut))
  expect_true(isEmpty(obj.mut))
  obj.mut@AllMutations <- data.frame(Position = 1, MutationFrequency = 15)
  expect_true(validObject(obj.mut))
  expect_false(isEmpty(obj.mut))

  # tas.sequences
  obj.seq <- new("tas.sequences", Table = data.frame(Sequences = "",
                                                     Count = NA_real_,
                                                     Percent = NA_real_,
                                                     Indels = "",
                                                     BasesChanged = NA_real_,
                                                     AA = "",
                                                     ProteinMutation = ""),
                                  Supplemental = data.table(Sequence = "",
                                                            UMIs = list(),
                                                            IDs = list()))
  expect_true(validObject(obj.seq))
  expect_true(isEmpty(obj.seq))
  obj.seq@Table <- data.frame(Sequences = "ATGCAG",
                              Count = 100,
                              Percent = 50,
                              Indels = "WT",
                              BasesChanged = 0,
                              AA = "MQ",
                              ProteinMutation = "WT")
  expect_true(validObject(obj.seq))
  expect_false(isEmpty(obj.seq))

  # tas.aid.tables
  obj.aid <- new("tas.aid.tables", WRCH = data.frame(Motif = "",
                                                     Start = NA_integer_,
                                                     End = NA_integer_,
                                                     Cytosine = NA_integer_,
                                                     MotifMutagenesis = NA_real_,
                                                     CytosineMutagenesis = NA_real_),
                                   WRCY = data.frame(Motif = "",
                                                     Start = NA_integer_,
                                                     End = NA_integer_,
                                                     Cytosine = NA_integer_,
                                                     MotifMutagenesis = NA_real_,
                                                     CytosineMutagenesis = NA_real_))
  expect_true(validObject(obj.aid))
  expect_true(isEmpty(obj.aid))
  obj.aid@WRCH <- data.frame(Motif = "WRCH",
                             Start = 1L,
                             End = 4L,
                             Cytosine = 3L,
                             MotifMutagenesis = 11.3,
                             CytosineMutagenesis = 10.7)
  expect_true(validObject(obj.aid))
  expect_false(isEmpty(obj.aid))

  # tas.dna.repair
  obj.dna <- new("tas.dna.repair", WT = NA_real_,
                                   NHEJ = NA_real_,
                                   MMEJ = NA_real_,
                                   BaseChange = NA_real_,
                                   IndelBaseChange = NA_real_,
                                   Other = NA_real_)
  expect_true(validObject(obj.dna))
  expect_true(isEmpty(obj.dna))
  obj.dna@WT <- 97.4
  obj.dna@NHEJ <- 0.5
  obj.dna@MMEJ <- 0.1
  obj.dna@BaseChange <- 1.8
  obj.dna@IndelBaseChange <- 0
  obj.dna@Other <- 0.2
  expect_true(validObject(obj.dna))
  expect_false(isEmpty(obj.dna))

  # tas.alignment
  obj.align <- new("tas.alignment", DNA = empty.pass(),
                                    AA = empty.pass())
  expect_true(validObject(obj.align))
  expect_true(isEmpty(obj.align))
  obj.align@DNA <- pairwiseAlignment(DNAStringSet(c("ATGCAG", "ATGGAG")), DNAString("ATGCAG"))
  expect_error(validObject(obj.align), "pairwise alignments are not equal")
  obj.align@AA <- pairwiseAlignment(AAStringSet(c("MQ", "ME")), AAString("MQ"))
  expect_true(validObject(obj.align))
  expect_false(isEmpty(obj.align))

  # AmpliconSequencing
  obj.as <- new("AmpliconSequencing", Alignment = new("tas.alignment"),
                                      Sequences = new("tas.sequences"),
                                      Mutations = new("tas.mutations"),
                                      MutationTypes = new("tas.dna.repair"),
                                      AIDTables = new("tas.aid.tables"),
                                      Settings = new("tas.object.settings"))
  expect_true(validObject(obj.as))
  expect_true(isEmpty(obj.as))
  obj.as@Mutations@AllMutations <- data.frame(Position = 1, MutationFrequency = 15)
  expect_true(validObject(obj.as))
  expect_false(isEmpty(obj.as))
})
