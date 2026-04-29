
# --------------
# findAIDtargets
# --------------

test_that("findAIDtargets locates AID hotspot motifs and returns correctly formatted outputs.", {
  test.aid <- findAIDtargets(test.settings)

  expect_true(class(test.aid) == "list")
  expect_all_true(names(test.aid) == c("WRCH", "WRCY"))
  expect_all_true(lapply(test.aid, class) == "data.frame")
  expect_all_true(colnames(test.aid$WRCH) == c("Motif", "Start", "End", "Cytosine"))
  expect_true(nrow(test.aid$WRCH) == 36)

  expect_all_true(test.aid$WRCH$Motif == c('WRCH', 'WRCH', 'DGYW', 'DGYW', 'WRCH', 'WRCH', 'DGYW', 'WRCH', 'DGYW', 'WRCH', 'DGYW', 'WRCH', 'WRCH', 'DGYW', 'DGYW', 'DGYW', 'WRCH', 'DGYW', 'WRCH', 'DGYW', 'DGYW', 'DGYW', 'WRCH', 'WRCH', 'WRCH', 'DGYW', 'WRCH', 'DGYW', 'WRCH', 'WRCH', 'DGYW', 'WRCH', 'DGYW', 'DGYW', 'WRCH', 'DGYW'))
  expect_all_true(test.aid$WRCH$Start == as.integer(c(5, 24, 24, 30, 32, 35, 43, 58, 58, 63, 73, 85, 91, 104, 121, 138, 141, 145, 151, 151, 172, 186, 201, 208, 228, 233, 239, 239, 247, 257, 273, 283, 283, 299, 337, 337)))
  expect_all_true(test.aid$WRCH$End == as.integer(c(8, 27, 27, 33, 35, 38, 46, 61, 61, 66, 76, 88, 94, 107, 124, 141, 144, 148, 154, 154, 175, 189, 204, 211, 231, 236, 242, 242, 250, 260, 276, 286, 286, 302, 340, 340)))
  expect_all_true(test.aid$WRCH$Cytosine == as.integer(c(7, 26, 25, 31, 34, 37, 44, 60, 59, 65, 74, 87, 93, 105, 122, 139, 143, 146, 153, 152, 173, 187, 203, 210, 230, 234, 241, 240, 249, 259, 274, 285, 284, 300, 339, 338)))

  expect_true(nrow(test.aid$WRCY) == 25)
  expect_all_true(test.aid$WRCY$Motif == c('WRCY', 'WRCY', 'RGYW', 'WRCY', 'RGYW', 'WRCY', 'RGYW', 'WRCY', 'RGYW', 'WRCY', 'WRCY', 'RGYW', 'RGYW', 'RGYW', 'WRCY', 'RGYW', 'RGYW', 'RGYW', 'WRCY', 'WRCY', 'WRCY', 'RGYW', 'RGYW', 'WRCY', 'RGYW'))
  expect_all_true(test.aid$WRCY$Start == as.integer(c(5, 24, 30, 35, 43, 58, 58, 63, 73, 85, 91, 104, 121, 138, 141, 145, 151, 172, 208, 228, 257, 273, 299, 337, 337)))
  expect_all_true(test.aid$WRCY$End == as.integer(c(8, 27, 33, 38, 46, 61, 61, 66, 76, 88, 94, 107, 124, 141, 144, 148, 154, 175, 211, 231, 260, 276, 302, 340, 340)))
  expect_all_true(test.aid$WRCY$Cytosine == as.integer(c(7, 26, 31, 37, 44, 60, 59, 65, 74, 87, 93, 105, 122, 139, 143, 146, 152, 173, 210, 230, 259, 274, 300, 339, 338)))

})


# ----------------
# measureMutations
# ----------------

test_that("measureMutations and associated getter functions work.", {
  test.seq <- makeTestSequences(test.settings)
  expect_output(obj <- measureMutations(test.seq, test.settings))

  expect_true(class(obj) == "tas.mutations")
  expect_all_true(slotNames(obj) == slotNames("tas.mutations"))
  expect_all_true(slotNames(obj) == c("DNA", "AA", "AIDTables", "Antibody"))

  ### DNA ###

  dna <- getDNAMutations(obj)
  expect_true(class(dna) == "list")
  expect_all_true(names(dna) == c("AllMutations", "CytosineMutations", "NonCytosineMutations", "MotifSums"))
  expect_all_true(colnames(dna$AllMutations) == c("Position", "MutationFrequency"))
  expect_true(nrow(dna$AllMutations) == 344)
  expect_equal(dna$AllMutations[dna$AllMutations$MutationFrequency!=0,], data.frame(Position = c(93, 153, 255.5, 279.5, 285, 297, 298, 299, 300, 301),
                                                                                    MutationFrequency = c(12.5, 12.5, 12.5, 25, 12.5, 12.5, 25, 25, 25, 12.5)),
               ignore_attr = c("row.names"))
  expect_all_true(colnames(dna$CytosineMutations) == c("Position", "MutationFrequency"))
  expect_true(nrow(dna$CytosineMutations) == 36)
  expect_equal(dna$CytosineMutations[dna$CytosineMutations$MutationFrequency!=0,], data.frame(Position = c(93, 153, 285, 300),
                                                                                    MutationFrequency = c(12.5, 12.5, 12.5, 25)),
               ignore_attr = c("row.names"))
  expect_all_true(colnames(dna$NonCytosineMutations) == c("Position", "MutationFrequency"))
  expect_true(nrow(dna$NonCytosineMutations) == 344-36)
  expect_equal(dna$NonCytosineMutations[dna$NonCytosineMutations$MutationFrequency!=0,], data.frame(Position = c(255.5, 279.5, 297, 298, 299, 301),
                                                                                              MutationFrequency = c(12.5, 25, 12.5, 25, 25, 12.5)),
               ignore_attr = c("row.names"))
  expect_equal(dna$MotifSums, c(AverageAllMutations = 0.5087209,
                                AverageCytosineMutations = 1.7361111,
                                AverageNonCytosineMutations = 0.3652597,
                                FrequencyOfAllMutationsAtCytosines = 35.7142857,
                                FrequencyOfAllMutationsAtNonCytosines = 64.2857143))

  expect_equal(getMutationDistributionDNA(obj), dna$AllMutations)
  expect_equal(getMutationDistributionCytosine(obj), dna$CytosineMutations)
  expect_equal(getMutationDistributionNonCytosine(obj), dna$NonCytosineMutations)


  ### AA ###

  aa.pos <- getMutationDistributionAA(obj)
  expect_equal(nrow(aa.pos), 114)
  expect_equal(which(aa.pos$MutationFrequency != 0), c(51, 85, 93, 99, 113))
  expect_equal(aa.pos$MutationFrequency[which(aa.pos$MutationFrequency != 0)], c(12.5, 12.5, 25, 12.5, 12.5))

  aa.mm <- getMutationMatrixAA(obj)
  expect_equal(which(colSums(aa.mm) != 0), which(aa.pos$MutationFrequency != 0))
  expect_equal(colSums(aa.mm)[which(colSums(aa.mm) != 0)], aa.pos$MutationFrequency[which(aa.pos$MutationFrequency != 0)])
  expect_equal(rowSums(aa.mm)[rowSums(aa.mm) != 0], setNames(c(12.5, 25, 37.5), c("R", "-", "+")))


  ### AID Tables ###

  aidt <- getAIDTables(obj)
  expect_equal(class(aidt), "list")
  expect_equal(names(aidt), c("WRCH", "WRCY"))
  expect_equal(sapply(aidt, class), c(WRCH = "data.frame", WRCY = "data.frame"))
  expect_equal(sapply(aidt, nrow), c(WRCH = 36, WRCY = 25))
  expect_equal(aidt$WRCH$CytosineMutationFrequency, c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12.5, 0, 0, 0, 0, 0, 0, 12.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12.5, 25, 0, 0))
  expect_equal(aidt$WRCY$CytosineMutationFrequency, c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 0, 0))
  expect_equal(getWRCHTable(obj), aidt$WRCH)
  expect_equal(getWRCYTable(obj), aidt$WRCY)

  ### Antibody ###
  abml <- getAbMutations(obj)
  expect_true(class(abml) == "list")
  expect_true(names(abml) == "RegionMutations")
  expect_true(class(abml$RegionMutations) == "numeric")
  expect_all_true(names(abml$RegionMutations) == c("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"))
  expect_equal(abml$RegionMutations, c(FR1 = -100.000000, CDR1 = 1.785714, FR2 = -100.000000, CDR2 = 1.785714, FR3 = -78.571429, CDR3 = 714.285714, FR4 = -100.000000))


  ####################
  # AmpliconSequencing
  # ------------------

  obj.as <- new("AmpliconSequencing", Mutations = obj)
  expect_equal(getDNAMutations(obj), getDNAMutations(obj.as))
  expect_equal(getMutationDistributionDNA(obj), getMutationDistributionDNA(obj.as))
  expect_equal(getMutationDistributionCytosine(obj), getMutationDistributionCytosine(obj.as))
  expect_equal(getMutationDistributionNonCytosine(obj), getMutationDistributionNonCytosine(obj.as))
  expect_equal(getMutationMotifSums(obj), getMutationMotifSums(obj.as))
  expect_equal(getMutationDistributionAA(obj), getMutationDistributionAA(obj.as))
  expect_equal(getMutationMatrixAA(obj), getMutationMatrixAA(obj.as))
  expect_equal(getAIDTables(obj), getAIDTables(obj.as))
  expect_equal(getWRCHTable(obj), getWRCHTable(obj.as))
  expect_equal(getWRCYTable(obj), getWRCYTable(obj.as))
  expect_equal(getAbMutations(obj), getAbMutations(obj.as))
})




