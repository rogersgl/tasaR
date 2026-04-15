################################################################################
########################## TAS_LOAD_DEPENDENCIES ###############################
################################################################################

test_that("Required packages are loaded", {
  packages <- c("BiocManager",
                "R.utils",
                "zip",
                "data.table",
                "openxlsx",
                "stringr",
                "ggplot2",
                "ggseqlogo",
                "parallel",
                "parallelly",
                "seqinr",
                "ape",
                "Biostrings",
                "ShortRead",
                "pwalign",
                "msa",
                "ggmsa")

  #check all
  expect_all_true(packages %in% .packages())

  #check 1 by 1
  for (i in packages){
    expect_true(i %in% .packages())
  }
})
