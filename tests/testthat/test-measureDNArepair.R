# ----------------
# measureDNArepair
# ----------------

test_that("measureDNArepair works", {
  test.seq <- makeTestSequences(test.settings)
  obj <- measureDNArepair(test.seq, test.settings)
  expect_true(class(obj) == "tas.dna.repair")

  dr <- getMutationTypes(obj)
  expect_equal(colnames(dr), c("WT", "NHEJ", "MMEJ", "BaseChange", "IndelBaseChange", "Other"))
  expect_equal(class(dr), "data.frame")
  expect_equal(dr$WT, 25)
  expect_equal(dr$NHEJ, 12.5)
  expect_equal(dr$MMEJ, 25)
  expect_equal(dr$BaseChange, 12.5)
  expect_equal(dr$IndelBaseChange, 12.5)
  expect_equal(dr$Other, 12.5)

  obj.as <- new("AmpliconSequencing", MutationTypes = obj)
  expect_equal(getMutationTypes(obj), getMutationTypes(obj.as))
})
