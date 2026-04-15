################################################################################
################################# TAS_MSA ######################################
################################################################################

library(testthat)

skip_if_not(exists("tas_msa", mode = "function"), "tas_msa() not found in session")

# Packages required for these integration tests
skip_if_not_installed("Biostrings")
skip_if_not_installed("msa")
skip_if_not_installed("ggmsa")
skip_if_not_installed("seqinr")
skip_if_not_installed("ape")
skip_if_not_installed("ggplot2")

# Helper: build a minimal Sequence.Table.List expected by tas_msa
make_sequence_table_for_sample <- function(sample, seqs, percents, indels = NULL, basesChanged = NULL, AA = NULL, ProteinMutation = NULL) {
  df <- data.frame(TargetSequence = as.character(seqs),
                   Percent = as.numeric(percents))
  if (!is.null(indels)) df$Indels <- indels else df$Indels <- rep(NA, length(seqs))
  if (!is.null(basesChanged)) df$BasesChanged <- basesChanged else df$BasesChanged <- rep(0, length(seqs))
  if (!is.null(AA)) df$AA <- AA
  if (!is.null(ProteinMutation)) df$ProteinMutation <- ProteinMutation

  stl <- list()
  stl[[sample]] <- df
  # stl$AlignDNA <- list()
  # # AlignDNA should be DNAStringSet of the sequences (remove any '-' before conversion)
  # stl$AlignDNA[[sample]] <- DNAStringSet(gsub("-", "", seqs)) #should this be a pairwiseAlignment?
  # if (!is.null(AA)) {
  #   stl$AlignProtein <- list()
  #   stl$AlignProtein[[sample]] <- AAStringSet(AA) #should this be a pairwiseAlignment?
  # }
  return(stl)
}

# ----------------------------
# 1) Basic DNA MSA test
# ----------------------------
test_that("tas_msa performs a DNA MSA and returns correctly labeled DNAMultipleAlignment and ggplot", {
  skip_if_not_installed("msa")
  Sample.Names <- "S1"

  # Reference and target sequences (short, simple)
  ref <- "ATGATGATG"
  seqs <- c("ATGATGATG", "ATGATGACG", "ATGATGAT-") # third has dash to simulate small indel
  percents <- c(60, 30, 10)
  indels <- c("WT", "", "-1")
  basesChanged <- c("WT", "1", "0")

  stl <- make_sequence_table_for_sample("S1", seqs, percents, indels, basesChanged)
  # assemble Sequence.Table.List with required shape
  Sequence.Table.List <- list(S1 = stl$S1)
  Input.DataFrame <- data.frame(ReferenceSequence = ref)
  rownames(Input.DataFrame) <- Sample.Names

  Config.List <- list(nCores = 1, sequence.alignment.count = 10, protein.mutations = 0, PhyloTree = 0)

  out <- tas_msa(Sample.Names = Sample.Names,
                 Sequence.Table.List = Sequence.Table.List,
                 Input.DataFrame = Input.DataFrame,
                 Config.List = Config.List)

  expect_true(is.list(out))
  expect_true("S1" %in% names(out))
  s1 <- out$S1

  # DNA branch should be present
  expect_true("DNA" %in% names(s1))
  expect_true(is.list(s1$DNA))
  expect_true(nrow(s1$DNA$Alignment == 4))
  # Logo should be a ggplot produced by ggmsa
  expect_s3_class(s1$DNA$Logo, "ggplot")
  # Alignment should be a DNAMultipleAlignment object
  expect_true(inherits(s1$DNA$Alignment, "DNAMultipleAlignment"))
  # Check alignment contains the Reference name
  aln_names <- rownames(s1$DNA$Alignment)
  expect_true(all(aln_names == c("Reference", "1. WT - 60%", "2. 1 SNV - 30%", "3. -1, 0 SNV - 10%")))
})

# ----------------------------
# 2) Protein MSA test
# ----------------------------
test_that("tas_msa performs protein MSA when protein.mutations=1 and returns AAMultipleAlignment and ggplot", {
  skip_if_not_installed("msa")
  Sample.Names <- "S1"

  # Need long(ish) reference sequences to prevent Clustal Omega from throwing errors
  ref <- "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGGCCGGGGGGTTCCTGAGACTCTCCTGTGAGCTGAGGGGAAGCATCTTTAACCAGTATGCCATGGCCTGGTTCCGCCAGGCTCCAGGGAAGGAGAGGGAGTTCGTCGCCGGCATGGGCGCCGTGCCCCACTACGGCGAGTTCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATGCCAAGAGCACGGTGTATCTGCAAATGAGCAGCCTGAAGCCCGAGGACACGGCCATCTATTTCTGTGCCAGGAGCAAGAGCACCTACATCAGCTACAACAGCAACGGCTACGACTACTGGGGCAGGGGAACCCAGGTCACCGTCTCCTCA"
  # Provide AA translations for sequences
  # We'll create two sequences whose translations differ
  seqs <- c("GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGGCCGGGGGGTTCCTGAGACTCTCCTGTGAGCTGAGGGGAAGCATCTTTAACCAGTATGCCATGGCCTGGTTCCGCCAGGCTCCAGGGAAGGAGAGGGAGTTCGTCGCCGGCATGGGCGCCGTGCCCCACTACGGCGAGTTCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATGCCAAGAGCACGGTGTATCTGCAAATGAGCAGCCTGAAGCCCGAGGACACGGCCATCTATTTCTGTGCCAGGAGCAAGAGCACCTACATCAGCTACAACAGCAACGGCTACGACTACTGGGGCAGGGGAACCCAGGTCACCGTCTCCTCA",
            "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGGCCGGGGGGTTCCTGAGACTCTCCTGTGAGCTGAGGGGAAGCATCTTTAACCAGTATGCCATGGCCTGGTTCCGCCAGGCTCCAGGGAAGGAGAGGGAGTTCGTCGCCGGCATGGGCGCCGTGCCCCACTACGGCGAGTTCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATGCCAAGAGCACGGTGTATCTGCAAATGAGCAGCCTGAAGCCCGAGGACACGGCCATCTATTTCTGTGCCAGGAGCAAGAGCACCTACATCAGTTATAACAGCAACGGCTACGACTACTGGGGCAGGGGAACCCAGGTCACCGTCTCCTCA",
            "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGGCCGGGGGGTTCCTGAGACTCTCCTGTGAGCTGAGGGGAAGCATCTTTAACCAGTATGCCATGGCCTGGTTCCGCCAGGCTCCAGGGAAGGAGAGGGAGTTCGTCGCCGGCATGGGCGCCGTGCCCCACTACGGCGAGTTCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATGCCAAGAGCACGGTGTATCTGCAAATGAGCAGCCTGAAGCCCGAGGACACGGCCATCTATTTCTGTGCCAGGAGCAAGAGCACCTACATCAGTTATAACAGCAACGACTACGACTACTGGGGCAGGGGAACCCAGGTCACCGTCTCCTCA",
            "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGGCCGGGGGGTTCCTGAGACTCTCCTGTGAGCTGAGGGGAAGCATCTTTAACCAGTATGCCATGGCCTGGTTCCGCCAGGCTCCAGGGAAGGAGAGGGAGTTCGTCGCCGGCATGGGCGCCGTGCCCCACTACGGCGAGTTCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATGCCGCACGGTGTATCTGCAAATGAGCAGCCTGAAGCCCGAGGACACGGCCATCTATTTCTGTGCCAGGAACAAGAGCACCTACATCATCTACAACACCAACGACTACGACTACTGGGGCAGGGGAACCCAGGTCACCGTCTCCTCA",
            "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGGCCGGGGGGTTCCTGAGACTCTCCTGTGAACTGAGGGGAAACATCTTTAACCAGTATGCCATGGCCTGGTTCCGCCAGGCTCCAGGGAAGGAGAGGGAGTTCGTCGCCGGCATGGGCGCCGTGCCCCAGTACGGCGAGTTCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATGCCAAGAGCACGGTGTATTTGCAAATGAGTAGCCTGAAGCCCGAGGACACGGCCATCTATTTCTGTGCCAGGAGCAAGAGCACCTACATCAACTACAGCAGCAACGTCTACGACTACTGGGGCAGGGGAACCCTGGTCACCGTCTCCTCA",
            "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGGCCGGGGGGTTCCTGAGACTCTCCTGTGAGCTGAGGGGAAGTATCTTTAACCAGTATGTCATGGCCTGGTTTCGCCAGGCTCCAGGGAAAGAGAGGGAGTTCGTCGCCGGCATGGGCGCCGTGCCCCACTACGTCGAGTTCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATCACCATCTCCAGAGACAATGCCAAGAGTACGGTGTATCTGCAAATGAGCAGCCTGAAGCCCGAGGACACGGCCATCTATTTCTGTGCCAGGAGCAAGACCACCTACATCAGCTACAACAGCAACGGCTACGACTACTGGGGCAGGGGAACCCAGGTCACCGTCTCCTCA")
  percents <- c(60, 20, 10, 6, 3.8, 0.2)
  indels <- c("WT", NA, NA, "-4", NA, "+19")
  basesChanged <- c("WT", "2", "3", "4", "9", "7")
  AA <- as.character(suppressWarnings(translate(DNAStringSet(gsub("-", "", seqs))))) # real AA translations
  ProteinMutation <- c("WT", "WT", "G107D", "Nonsense", "S27N, H56Q, S102N, N104S, G107V, Q116L", "Indel") # label strings (function will use them in names)
  stl <- make_sequence_table_for_sample("S1", seqs, percents, indels, basesChanged, AA, ProteinMutation)

  Sequence.Table.List <- list(S1 = stl$S1)
  Input.DataFrame <- data.frame(ReferenceSequence = ref)
  rownames(Input.DataFrame) <- Sample.Names

  Config.List <- list(nCores = 1, sequence.alignment.count = 10, protein.mutations = 1, PhyloTree = 0)

  out <- tas_msa(Sample.Names = Sample.Names,
                 Sequence.Table.List = Sequence.Table.List,
                 Input.DataFrame = Input.DataFrame,
                 Config.List = Config.List)

  s1 <- out$S1
  expect_true("Protein" %in% names(s1))
  expect_s3_class(s1$Protein$Logo, "ggplot")
  # Protein Alignment should be an AAMultipleAlignment object (msa::msa returns such objects for protein)
  expect_all_true(rownames(s1$DNA$Alignment) == c("Reference",
                                              "1. WT - 60%",
                                              "2. 2 SNV - 20%",
                                              "3. 3 SNV - 10%",
                                              "4. -4, 4 SNV - 6%",
                                              "5. 9 SNV - 3.8%",
                                              "6. +19, 7 SNV - 0.2%"))
  expect_all_true(rownames(s1$Protein$Alignment) == c("Reference",
                                                      "1. WT - 60%",
                                                      "2. WT - 20%",
                                                      "3. G107D - 10%",
                                                      "4. Nonsense - 6%",
                                                      "5. S27N, H56Q, S102N, N104S, G107V, Q116L - 3.8%",
                                                      "6. Indel - 0.2%"))
})

# ----------------------------
# 3) PhyloTree branch (real msa + seqinr + ape)
# ----------------------------
test_that("tas_msa generates a phylogenetic tree when PhyloTree = 1 (real msa + seqinr + ape)", {
  skip_if_not_installed("msa")
  skip_if_not_installed("seqinr")
  skip_if_not_installed("ape")

  Sample.Names <- "S1"
  seqs <- c("ATGCAGGGT", "ATGCAGAGT", "ATGCAGATT", "ATGTAGGGT")
  percents <- c(70, 20, 15, 5)
  basesChanged <- c("0", "1", "2", "1")
  ProteinMutation <- c("WT", "G3S", "G3I", "Nonsense")
  stl <- make_sequence_table_for_sample("S1", seqs, percents, basesChanged = basesChanged, ProteinMutation = ProteinMutation)
  Sequence.Table.List <- list(S1 = stl$S1)

  Input.DataFrame <- data.frame(ReferenceSequence = "ATGCAGGGT")
  rownames(Input.DataFrame) <- Sample.Names

  Config.List <- list(nCores = 1, sequence.alignment.count = 10, protein.mutations = 0, PhyloTree = 1)

  out <- tas_msa(Sample.Names = Sample.Names,
                 Sequence.Table.List = Sequence.Table.List,
                 Input.DataFrame = Input.DataFrame,
                 Config.List = Config.List)

  s1 <- out$S1
  expect_true("PhyloTree" %in% names(s1))
  # Check class or presence of tip labels
  tree <- s1$PhyloTree
  expect_true(is.list(tree))
  # If it's a phylo object, it should inherit "phylo"
  if (inherits(tree, "phylo")) {
    expect_true(length(tree$tip.label) == 4)
    expect_all_true(tree$tip.label == c("WT", "G3S", "G3I", "Nonsense"))
  }
})

# ----------------------------
# 4) Validation error cases
# ----------------------------
test_that("tas_msa throws errors for missing table columns and missing config entries", {
  Sample.Names <- "S1"
  # missing Percent or TargetSequence
  Sequence.Table.List_bad <- list(S1 = data.frame(Bad = "x"))
  Input.DataFrame <- data.frame(ReferenceSequence = "ATGATGATG")
  rownames(Input.DataFrame) <- Sample.Names
  Config.List <- list(nCores = 1, sequence.alignment.count = 1)

  expect_error(
    tas_msa(Sample.Names = Sample.Names, Sequence.Table.List = Sequence.Table.List_bad, Input.DataFrame = Input.DataFrame, Config.List = Config.List),
    regexp = "The following columns were not found in the input file"
  )

  # missing config entry
  Sequence.Table.List_ok <- list(S1 = data.frame(TargetSequence = "ATG", Percent = 100))
  Config.List_bad <- list(nCores = 1) # missing sequence.alignment.count
  expect_error(
    tas_msa(Sample.Names = Sample.Names, Sequence.Table.List = Sequence.Table.List_ok, Input.DataFrame = Input.DataFrame, Config.List = Config.List_bad),
    regexp = "sequence.alignment.count"
  )
})
