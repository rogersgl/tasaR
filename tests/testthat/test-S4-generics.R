###############################################################################
# Testing helper functions
# -------------

fakeIlluminaId <- function() {
  str_c("A90033:281:WB24762987th.Miseq:1:", sample(100:50000, 1), ":", sample(100:50000, 1), ":AAGAGGCA+CGGAGAGA")
}

makeTestSeqTable <- function() {
  t <- data.frame(Sequences = c("ATGCAG", "ATGGAG", "ATGCG"),
                  Count = c(5, 2, 1),
                  Percent = c(5, 2, 1)/sum(55, 17, 3)*100,
                  Indels = c("", "", "-1"),
                  BasesChanged = (c(0, 1, 0)),
                  AA = c("MQ", "ME", "M"),
                  ProteinMutation = c("WT", "Q2E", "Q2-"))
  s <- data.table(Index = c(1L, 2L, 3L),
                  UMIs = list(list("TGACG", "ACATA", "CGGTA", "GAACA", "TGGCA"), list("GGGCC", "ACCTC"), list("TCGTA")),
                  IDs = list(list(list(fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId()),
                               list(fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId()),
                               list(fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId()),
                               list(fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId()),
                               list(fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId())
                               ),
                          list(list(fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId()),
                               list(fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId())),
                          list(list(fakeIlluminaId(), fakeIlluminaId(), fakeIlluminaId()))
                  ),
                  IndelStart = NA_real_,
                  IndelType = NA_character_)
  da <- pairwiseAlignment(DNAStringSet(c("ATGCAG", "ATGGAG")), DNAString("ATGCAG"))
  pa <- pairwiseAlignment(AAStringSet(c("MQ", "ME")), AAString("MQ"))
  ReadCounts <- c(Merged = 53L, Filtered = 44L, UMIs = 8L, UniqueSequences = 3L)

  new("tas.sequences", Table = t, Supplemental = s, Alignments = list(DNA = da, AA = pa), ReadCounts = ReadCounts)
}


makeTestMutation <- function() {
  a <- data.frame(Position = c(1, 2, 3, 4, 5), MutationFrequency = c(0.1, 8.93, 57.02, 1.21, 2.76))
  c <- data.frame(Position = 3, MutationFrequency = 57.02)
  nc <- data.frame(Position = c(1, 2, 4, 5), MutationFrequency = c(0.1, 8.93, 1.21, 2.76))
  ms <-  c(AverageAllMutations = 0.91,
                AverageCytosineMutations = 7.72,
                AverageNonCytosineMutations = 0.63,
                FrequencyOfAllMutationsAtCytosines = 77.3,
                FrequencyOfAllMutationsAtNonCytosines = 22.7)
  ap <- data.frame(Position = c(1, 2), MutationFrequency = c(4.3, 48.7))
  mm <- consensusMatrix(pairwiseAlignment(AAStringSet(c("MQ", "ME")), AAStringSet("MQ")))
  mm <- mm[1:(nrow(mm)-3),]

  aidt <- list(WRCH = data.frame(Motif = "WRCH", Start = 1, End = 4, Cytosine = 3, CytosineMutationFrequency = 57.2),
               WRCY = data.frame(Motif = "WRCY", Start = 1, End = 4, Cytosine = 3, CytosineMutationFrequency = 57.2))

  new("tas.mutations", DNA = list(AllMutations = a, CytosineMutations = c, NonCytosineMutations = nc, MotifSums = ms),
      AA = list(AllMutations = ap,
                MutationMatrix = mm),
      AIDTables = aidt
      )
}


makeTestMutationTypes <- function() {
  new("tas.dna.repair", WT = 100, NHEJ = 0, MMEJ = 0, BaseChange = 0, IndelBaseChange = 0, Other = 0)
}


makeTestSettings <- function() {
  new("tas.object.settings", Name = "test",
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
      AntibodyRegions = c(FR1Start = 1L,
                          CDR1Start = 73L,
                          FR2Start = 97L,
                          CDR2Start = 148L,
                          FR3Start = 172L,
                          CDR3Start = 286L,
                          FR4Start = 310L,
                          FR4End = 342L))
}

###############################################################################
# Getter methods
# -------------

### getAlignments ###

test_that("getAlignments, getDNAalign, and getAAalign return the expected outputs.", {

  # class tas.alignment
  obj <- makeTestSeqTable()
  aln <- getAlignments(obj)

  expect_true(class(aln) == "list")
  expect_all_true(names(aln) == c("DNA", "AA"))

  expect_true(class(aln$DNA) == "PairwiseAlignmentsSingleSubject")
  expect_true(class(aln$AA) == "PairwiseAlignmentsSingleSubject")

  expect_all_true(pattern(aln$DNA) == DNAStringSet(c("ATGCAG", "ATGGAG")))
  expect_all_true(subject(aln$DNA) == DNAStringSet(c("ATGCAG", "ATGCAG")))

  expect_all_true(pattern(aln$AA) == AAStringSet(c("MQ", "ME")))
  expect_all_true(subject(aln$AA) == AAStringSet(c("MQ", "MQ")))

  # class AmpliconSequencing
  obj.as <- new("AmpliconSequencing", Sequences = obj)
  aln.as <- getAlignments(obj.as)
  expect_true(identical(aln.as, aln))

  #getDNAalign
  aln.dna <- getDNAalign(obj)
  expect_true(identical(aln.dna, aln$DNA))
  aln.dna.as <- getDNAalign(obj.as)
  expect_true(identical(aln.dna.as, aln$DNA))

  #getAAalign
  aln.aa <- getAAalign(obj)
  expect_true(identical(aln.aa, aln$AA))
  aln.aa.as <- getAAalign(obj.as)
  expect_true(identical(aln.aa.as, aln$AA))
})



### getSequenceTable ###

test_that("Sequence Table getter methods return the expected outputs.", {
  obj <- makeTestSeqTable()
  seq <- getSequenceTable(obj)

  expect_true(class(seq) == "data.frame")
  expect_all_true(colnames(seq) == c("Sequences", "Count", "Percent", "Indels", "BasesChanged", "AA", "ProteinMutation"))
  expect_all_true(seq$Sequences == c("ATGCAG", "ATGGAG", "ATGCG"))
  expect_all_true(seq$Count == c(5, 2, 1))
  expect_all_true(seq$Percent == c(5, 2, 1)/sum(55, 17, 3)*100)
  expect_all_true(seq$Indels == c("", "", "-1"))
  expect_all_true(seq$BasesChanged == (c(0, 1, 0)))
  expect_all_true(seq$AA == c("MQ", "ME", "M"))
  expect_all_true(seq$ProteinMutation == c("WT", "Q2E", "Q2-"))

  dna.seq <- getSequencesDNA(obj)
  expect_equal(dna.seq, seq$Sequences)
  aa.seq <- getSequencesAA(obj)
  expect_equal(aa.seq, seq$AA)

  obj.as <- new("AmpliconSequencing", Sequences = obj)
  seq.as <- getSequenceTable(obj.as)
  expect_true(identical(seq.as, seq))
  dna.seq.as <- getSequencesDNA(obj.as)
  expect_equal(dna.seq.as, seq.as$Sequences)
  aa.seq.as <- getSequencesAA(obj.as)
  expect_equal(aa.seq.as, seq.as$AA)

})

test_that("getSequenceSupplemental returns the expected outputs.", {
  obj <- makeTestSeqTable()
  sup <- getSequenceSupplemental(obj)

  expect_all_true(class(sup) == c("data.table", "data.frame"))
  expect_all_true(colnames(sup) == c("Sequence", "UMIs", "IDs", "IndelStart", "IndelType"))
  expect_all_true(sup$Sequence == c("ATGCAG", "ATGGAG", "ATGCG"))
  expect_true(identical(sup$UMIs, list(list("TGACG", "ACATA", "CGGTA", "GAACA", "TGGCA"), list("GGGCC", "ACCTC"), list("TCGTA"))))
  expect_true(length(sup$IDs) == 3)
  expect_true(length(unlist(sup$IDs)) == 42)
  expect_all_true(str_detect(unlist(sup$IDs), "A90033:281:WB24762987th.Miseq:1:"))

  obj.as <- new("AmpliconSequencing", Sequences = obj)
  sup.as <- getSequenceSupplemental(obj.as)
  expect_true(identical(sup.as, sup))
})



### getMutations ###

test_that("Getter functions for mutations return the expected outputs.", {

  ### class tas.mutations ###
  obj <- makeTestMutation()

  ### DNA ###
  mut.dna <- getDNAMutations(obj)

  expect_true(class(mut.dna) == "list")
  expect_true(length(mut.dna) == 4)
  expect_all_true(names(mut.dna) == c("AllMutations", "CytosineMutations", "NonCytosineMutations", "MotifSums"))
  expect_true(identical(mut.dna$AllMutations, data.frame(Position = c(1, 2, 3, 4, 5), MutationFrequency = c(0.1, 8.93, 57.02, 1.21, 2.76))))
  expect_true(identical(mut.dna$CytosineMutations, data.frame(Position = 3, MutationFrequency = 57.02)))
  expect_true(identical(mut.dna$NonCytosineMutations, data.frame(Position = c(1, 2, 4, 5), MutationFrequency = c(0.1, 8.93, 1.21, 2.76))))
  expect_true(identical(mut.dna$MotifSums, c(AverageAllMutations = 0.91,
                                         AverageCytosineMutations = 7.72,
                                         AverageNonCytosineMutations = 0.63,
                                         FrequencyOfAllMutationsAtCytosines = 77.3,
                                         FrequencyOfAllMutationsAtNonCytosines = 22.7)))
  # getMutationDistribution
  mut.a <- getMutationDistributionDNA(obj)
  expect_true(identical(mut.a, mut.dna$AllMutations))

  # getMutationDistributionCytosine
  mut.c <- getMutationDistributionCytosine(obj)
  expect_true(identical(mut.c, mut.dna$CytosineMutations))

   # getMutationDistributionNonCytosine
  mut.nc <- getMutationDistributionNonCytosine(obj)
  expect_true(identical(mut.nc, mut.dna$NonCytosineMutations))

  # getMutationMotifSums
  mut.ms <- getMutationMotifSums(obj)
  expect_true(identical(mut.ms, mut.dna$MotifSums))


  # class AmpliconSequencing #
  obj.as <- new("AmpliconSequencing", Mutations = obj)
  mut.as <- getDNAMutations(obj.as)
  expect_true(identical(mut.as, mut.dna))
  mut.as.a <- getMutationDistributionDNA(obj.as)
  expect_true(identical(mut.as.a, mut.as$AllMutations))
  mut.as.c <- getMutationDistributionCytosine(obj.as)
  expect_true(identical(mut.as.c, mut.as$CytosineMutations))
  mut.as.nc <- getMutationDistributionNonCytosine(obj.as)
  expect_true(identical(mut.as.nc, mut.as$NonCytosineMutations))
  mut.as.ms <- getMutationMotifSums(obj.as)
  expect_true(identical(mut.as.ms, mut.as$MotifSums))


  ### AA ###
  mut.aa <- getMutationDistributionAA(obj)
  expect_true(class(mut.aa) == "data.frame")
  expect_all_true(colnames(mut.aa) == c("Position", "MutationFrequency"))
  expect_true(identical(mut.aa, data.frame(Position = c(1, 2), MutationFrequency = c(4.3, 48.7))))

  mut.mm <- getMutationMatrixAA(obj)
  expect_all_true(class(mut.mm) == c("matrix", "array"))

  mm <- consensusMatrix(pairwiseAlignment(AAStringSet(c("MQ", "ME")), AAStringSet("MQ")))
  mm <- mm[1:(nrow(mm)-3),]
  expect_true(identical(mut.mm, mm))

  mut.aa.as <- getMutationDistributionAA(obj.as)
  expect_true(identical(mut.aa.as, mut.aa))
  mut.mm.as <- getMutationMatrixAA(obj.as)
  expect_true(identical(mut.mm.as, mut.mm))


  ### AID Tables ###

  aidt <- getAIDTables(obj)
  expect_true(class(aidt) == "list")
  expect_all_true(names(aidt) == c("WRCH", "WRCY"))
  expect_true(identical(aidt$WRCH, data.frame(Motif = "WRCH", Start = 1, End = 4, Cytosine = 3, CytosineMutationFrequency = 57.2)))
  expect_true(identical(aidt$WRCY, data.frame(Motif = "WRCY", Start = 1, End = 4, Cytosine = 3, CytosineMutationFrequency = 57.2)))

  aidt.h <- getWRCHTable(obj)
  expect_true(identical(aidt.h, aidt$WRCH))
  aidt.y <- getWRCYTable(obj)
  expect_true(identical(aidt.y, aidt$WRCY))

  aidt.as <- getAIDTables(obj.as)
  expect_true(identical(aidt.as, aidt))
  aidt.h.as <- getWRCHTable(obj.as)
  expect_true(identical(aidt.h.as, aidt.h))
  aidt.y.as <- getWRCYTable(obj.as)
  expect_true(identical(aidt.y.as, aidt.y))

})



### getMutationTypes ###


test_that("getMutationTypes returns the expected output.", {

  ### class tas.mutations ###
  obj <- makeTestMutationTypes()
  mt <- getMutationTypes(obj)

  expect_true(class(mt) == "data.frame")
  expect_all_true(colnames(mt) == c("WT", "NHEJ", "MMEJ", "BaseChange", "IndelBaseChange", "Other"))
  expect_all_true(unlist(mt) == c(100, 0, 0, 0, 0, 0))

  ### class AmpliconSequencing ###

  obj.as <- new("AmpliconSequencing", MutationTypes = obj)
  mt.as <- getMutationTypes(obj.as)
  expect_true(identical(mt.as, mt))
})



### getSettings ###

test_that("getSettings returns the expected output.", {

  ### class tas.settings ###
  obj <- makeTestSettings()
  s <- getSettings(obj)

  expect_true(class(s) == "list")
  l <- as.list(make.set.df())
  l <- c(head(l, -8), list(AntibodyRegions = unlist(tail(l,8))))
  expect_true(identical(s, l))

  ### class AmpliconSequencing ###
  obj.as <- new("AmpliconSequencing", Settings = obj)
  s.as <- getSettings(obj.as)
  expect_true(identical(s.as, s))

})


