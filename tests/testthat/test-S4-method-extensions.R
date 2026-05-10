# -------
# isEmpty
# -------

test_that("Validate tasAnalyzer extensions for function isEmpty.", {

  #tas.object.settings
  obj.set <- new("tas.object.settings", Name = "unknown",
                 IsAntibody = FALSE,
                 MergedFASTQPath = "",
                 ReferenceSequence = list(""),
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
  obj.mut <- new("tas.mutations", DNA = list(AllMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                                                        CytosineMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                                                        NonCytosineMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                                                        MotifSums = c(AverageAllMutations = NA_real_,
                                                                      AverageCytosineMutations = NA_real_,
                                                                      AverageNonCytosineMutations = NA_real_,
                                                                      FrequencyOfAllMutationsAtCytosines = NA_real_,
                                                                      FrequencyOfAllMutationsAtNonCytosines = NA_real_)),
                                             AA = list(AllMutations = data.frame(Position = numeric(), MutationFrequency = numeric()),
                                                       MutationMatrix = matrix()),
                                             AIDTables = list(WRCH = data.frame(Motif = character(), Start = integer(), End = integer(), Cytosine = integer(), CytosineMutationFrequency = numeric()),
                                                              WRCY = data.frame(Motif = character(), Start = integer(), End = integer(), Cytosine = integer(), CytosineMutationFrequency = numeric()))
                  )
  expect_true(validObject(obj.mut))
  expect_true(isEmpty(obj.mut))
  obj.mut@DNA$AllMutations <- data.frame(Position = 1, MutationFrequency = 15)
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
                 Supplemental = data.table::data.table(Index = NA_integer_,
                                                       UMIs = list(),
                                                       IDs = list(),
                                                       IndelStart = NA_real_,
                                                       IndelType = NA_character_,
                                                       RefIdx = NA_real_),
                 Alignments = list(DNA = list(.empty.pass()),
                                   AA = list(.empty.pass())),
                 ReadCounts = NA_integer_)
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

  # AmpliconSequencing
  obj.as <- new("AmpliconSequencing", Sequences = new("tas.sequences"),
                                      Mutations = new("tas.mutations"),
                                      MutationTypes = new("tas.dna.repair"),
                                      Settings = new("tas.object.settings"))
  expect_true(validObject(obj.as))
  expect_true(isEmpty(obj.as))
  obj.as@Mutations@DNA$AllMutations <- data.frame(Position = 1, MutationFrequency = 15)
  expect_true(validObject(obj.as))
  expect_false(isEmpty(obj.as))

  # PairedAmpliconSequencing
  obj.pas <- new("PairedAmpliconSequencing", Control = new("AmpliconSequencing"),
                 Experimental = new("AmpliconSequencing"))
  expect_true(validObject(obj.pas))
  expect_true(isEmpty(obj.pas))
  obj.pas@Control@Mutations@DNA$AllMutations <- data.frame(Position = 1, MutationFrequency = 15)
  expect_true(validObject(obj.pas))
  expect_false(isEmpty(obj.pas))
})


test_that("Extended show() functions print to console", {
  expect_output(show(methods::new("tas.object.settings")))
  expect_output(show(methods::new("tas.sequences")))
  expect_output(show(methods::new("tas.mutations")))
  expect_output(show(methods::new("tas.dna.repair")))
  expect_output(show(methods::new("AmpliconSequencing")))
  expect_output(test.results <- analyzeAmplicon(test.settings))
  expect_output(show(test.results))
  expect_output(test.paired.results <- pairedAnalyzeAmplicon(file.path(tempdir(), "pair-homo-ctrl-merged.fastq.gz"),
                                                             file.path(tempdir(), "pair-homo-expt-merged.fastq.gz"),
                                                             full.settings.ctrl))
  expect_output(show(test.paired.results))
})
