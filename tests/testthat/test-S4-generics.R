###############################################################################
# Testing helper functions
# -------------

makeTestAlignment <- function() {
  da <- pairwiseAlignment(DNAStringSet(c("ATGCAG", "ATGGAG")), DNAString("ATGCAG"))
  pa <- pairwiseAlignment(AAStringSet(c("MQ", "ME")), AAString("MQ"))
  new("tas.alignment", DNA = da, AA = pa)
}

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
  s <- data.table(Sequence = c("ATGCAG", "ATGGAG", "ATGCG"),
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
                  ))

  new("tas.sequences", Table = t, Supplemental = s)
}


makeTestMutation <- function() {
  a <- data.frame(Position = c(1, 2, 3, 4, 5), MutationFrequency = c(0.1, 8.93, 57.02, 1.21, 2.76))
  c <- data.frame(Position = 3, MutationFrequency = 57.02)
  nc <- data.frame(Position = c(1, 2, 4, 5), MutationFrequency = c(0.1, 8.93, 1.21, 2.76))
  ms <- c(0.1, 0.2, 4, 11.2, 93, 2.1, 1.3, 7, 57.02, 2.3, 5, 87)
  names(ms) <- c("TotalPercentMutated", "TotalPercentMutatedInMotif", "TotalPercentMutatedInCytosine", "MotifMutationAverage", "MotifMutationFrequencyofTotal", "MotifMutationNonCAverage", "MotifMutationNonCFrequencyOfTotal", "CytosineMutationAverage", "NonCytosineMutationAverage", "CytosineMutationFrequencyOfTotal", "CytosineMutationFrequencyOfMotif")
  new("tas.mutations", AllMutations = a, CytosineMutations = c, NonCytosineMutations = nc, MotifSums = ms)
}


makeTestMutationTypes <- function() {
  new("tas.dna.repair", WT = 100, NHEJ = 0, MMEJ = 0, BaseChange = 0, IndelBaseChange = 0, Other = 0)
}


makeTestAIDTables <- function() {
  h <- data.frame(Motif = "WRCH", Start = 1, End = 4, Cytosine = 3, MotifMutagenesis = 64.3, CytosineMutagenesis = 57.2)
  y <- data.frame(Motif = "WRCY", Start = 1, End = 4, Cytosine = 3, MotifMutagenesis = 64.3, CytosineMutagenesis = 57.2)
  new("tas.aid.tables", WRCH = h, WRCY = y)
}


makeTestSettings <- function() {
  new("tas.object.settings", Name = "test",
                      IsAntibody = TRUE,
                      MeasureSHM = TRUE,
                      MeasureDNARepairTypes = FALSE,
                      MergedFASTQPath = file.path(tempdir(), "test-merged.fastq"),
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
  obj <- makeTestAlignment()
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
  obj.as <- new("AmpliconSequencing", Alignment = makeTestAlignment())
  aln.as <- getAlignments(obj)
  expect_true(identical(aln, aln.as))

  #getDNAalign
  aln.dna <- getDNAalign(obj)
  expect_true(identical(aln.dna, aln$DNA))
  aln.dna.as <- getDNAalign(aln.as)
  expect_true(identical(aln.dna.as, aln$DNA))

  #getAAalign
  aln.aa <- getAAalign(obj)
  expect_true(identical(aln.aa, aln$AA))
  aln.aa.as <- getDNAalign(aln.as)
  expect_true(identical(aln.aa.as, aln$AA))
})



### getSequenceTable ###

test_that("getSequenceTable returns the expected outputs.", {
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

  obj.as <- new("AmpliconSequencing", Sequences = obj)
  seq.as <- getSequenceTable(obj.as)
  expect_true(identical(seq.as, seq))
})

test_that("getSequenceSupplemental returns the expected outputs.", {
  obj <- makeTestSeqTable()
  sup <- getSequenceSupplemental(obj)

  expect_all_true(class(sup) == c("data.table", "data.frame"))
  expect_all_true(colnames(sup) == c("Sequence", "UMIs", "IDs"))
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

test_that("getMutations and sub-functions return the expected outputs.", {

  ### class tas.mutations ###
  obj <- makeTestMutation()
  mut <- getMutations(obj)

  expect_true(class(mut) == "list")
  expect_true(length(mut) == 4)
  expect_all_true(names(mut) == c("AllMutations", "CytosineMutations", "NonCytosineMutations", "MotifSums"))
  expect_true(identical(mut$AllMutations, data.frame(Position = c(1, 2, 3, 4, 5), MutationFrequency = c(0.1, 8.93, 57.02, 1.21, 2.76))))
  expect_true(identical(mut$CytosineMutations, data.frame(Position = 3, MutationFrequency = 57.02)))
  expect_true(identical(mut$NonCytosineMutations, data.frame(Position = c(1, 2, 4, 5), MutationFrequency = c(0.1, 8.93, 1.21, 2.76))))
  expect_true(identical(mut$MotifSums, setNames(c(0.1, 0.2, 4, 11.2, 93, 2.1, 1.3, 7, 57.02, 2.3, 5, 87),
                                                c("TotalPercentMutated", "TotalPercentMutatedInMotif", "TotalPercentMutatedInCytosine", "MotifMutationAverage", "MotifMutationFrequencyofTotal", "MotifMutationNonCAverage", "MotifMutationNonCFrequencyOfTotal", "CytosineMutationAverage", "NonCytosineMutationAverage", "CytosineMutationFrequencyOfTotal", "CytosineMutationFrequencyOfMotif")
                                                )))
  # getMutationDistribution
  mut.a <- getMutationDistribution(obj)
  expect_true(identical(mut.a, mut$AllMutations))

  # getMutationDistributionCytosine
  mut.c <- getMutationDistributionCytosine(obj)
  expect_true(identical(mut.c, mut$CytosineMutations))

   # getMutationDistributionNonCytosine
  mut.nc <- getMutationDistributionNonCytosine(obj)
  expect_true(identical(mut.nc, mut$NonCytosineMutations))

  # getMutationMotifSums
  mut.ms <- getMutationMotifSums(obj)
  expect_true(identical(mut.ms, mut$MotifSums))


  ### class AmpliconSequencing ###
  obj.as <- new("AmpliconSequencing", Mutations = obj)
  mut.as <- getMutations(obj.as)
  expect_true(identical(mut.as, mut))
  mut.as.a <- getMutationDistribution(obj.as)
  expect_true(identical(mut.as.a, mut.as$AllMutations))
  mut.as.c <- getMutationDistributionCytosine(obj.as)
  expect_true(identical(mut.as.c, mut.as$CytosineMutations))
  mut.as.nc <- getMutationDistributionNonCytosine(obj.as)
  expect_true(identical(mut.as.nc, mut.as$NonCytosineMutations))
  mut.as.ms <- getMutationMotifSums(obj.as)
  expect_true(identical(mut.as.ms, mut.as$MotifSums))
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




### getAIDTables ###

test_that("getAIDTables and associated functions return the expected output.", {

  ### clas tas.aid.tables ###
  obj <- makeTestAIDTables()
  aid <- getAIDTables(obj)

  expect_true(class(aid) == "list")
  expect_all_true(names(aid) == c("WRCH", "WRCY"))
  expect_true(identical(aid$WRCH, data.frame(Motif = "WRCH", Start = 1, End = 4, Cytosine = 3, MotifMutagenesis = 64.3, CytosineMutagenesis = 57.2)))
  expect_true(identical(aid$WRCY, data.frame(Motif = "WRCY", Start = 1, End = 4, Cytosine = 3, MotifMutagenesis = 64.3, CytosineMutagenesis = 57.2)))

  # getWRCHTables
  aid.h <- getWRCHTable(obj)
  expect_true(identical(aid.h, aid$WRCH))

  # getWRCYTables
  aid.y <- getWRCYTable(obj)
  expect_true(identical(aid.y, aid$WRCY))


  ### class AmpliconSequencing ###

  obj.as <- new("AmpliconSequencing", AIDTables = obj)
  aid.as <- getAIDTables(obj.as)
  expect_true(identical(aid.as, aid))

  # getWRCHTables
  aid.as.h <- getWRCHTable(obj.as)
  expect_true(identical(aid.as.h, aid$WRCH))

  # getWRCYTables
  aid.as.y <- getWRCYTable(obj.as)
  expect_true(identical(aid.as.y, aid$WRCY))

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


