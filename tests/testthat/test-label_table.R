################################################################################
############################# TAS_LABEL_TABLE ##################################
################################################################################

skip_if_not(exists("tas_label_table", mode = "function"), "tas_label_table() not found in this session")

# These tests require Bioconductor packages
skip_if_not_installed("Biostrings")
skip_if_not_installed("S4Vectors")
skip_if_not_installed("BiocGenerics")

# Helper: create a copy of the function that uses lapply instead of mclapply
make_fcopy <- function(fn) {
  fcopy <- fn
  env_mock <- new.env(parent = environment(fn))
  env_mock$mclapply <- function(X, FUN, mc.cores = 1, ...) lapply(X, FUN, ...)
  environment(fcopy) <- env_mock
  return(fcopy)
}

# Helper: wrapper to create simple Sequence.Table.List structure from character vectors
# Sequence.Table.List[[sample]] should be a data.frame with column TargetSequence
make_sequence_table_list <- function(sample_names, seqs_list) {
  # seqs_list is a named list of character vectors, each vector are TargetSequence entries
  out <- list()
  for (s in sample_names) {
    out[[s]] <- data.frame(TargetSequence = as.character(seqs_list[[s]]), stringsAsFactors = FALSE)
  }
  return(out)
}

# Helper: create Reference.Sequences.DNA as list of DNAString objects named by samples
make_ref_dna_list <- function(sample_names, ref_seqs) {
  l <- lapply(ref_seqs, DNAString)
  names(l) <- sample_names
  return(l)
}

# ------------------------
# Test 1: Sample labeling
# ------------------------
test_that("tas_label_table labels sequences correctly", {
  fcopy <- make_fcopy(tas_label_table)

  Sample.Names <- c("S1")
  seqs_list <- list(S1 = c("ACGTACGT", # WT
                           "ACGACGT", # -1
                           "ACGTTACGT", # +1
                           "ACGGTCGG", # 3 SNVs
                           "ACTTGTACT", # +2 and -1
                           "AGGAGGT", # -1 and 2 SNVs
                           "ACGTACG", # -1 at 3' end (missed del)
                           "GACGTACGT" # +1 at 5' end (missed ins)
                           ))
  Sequence.Table.List <- make_sequence_table_list(Sample.Names, seqs_list)
  Reference.Sequences.DNA <- make_ref_dna_list(Sample.Names, list(S1 = "ACGTACGT"))

  Config.List <- list(nCores = 1, protein.mutations = 0)

  out <- fcopy(Sample.Names = Sample.Names,
               Sequence.Table.List = Sequence.Table.List,
               Reference.Sequences.DNA = Reference.Sequences.DNA,
               Config.List = Config.List)

  # Expect list with S1 (table) and AlignDNA appended, no AlignProtein with protein.mutations = 0
  expect_true(is.list(out))
  expect_true(all(c("S1","AlignDNA") %in% names(out)))
  expect_true(length(out$AlignDNA$S1) == 8)
  expect_false(all(c("S1","AlignProtein") %in% names(out)))

  # Expect column labels of table
  s1_table <- out$S1
  expect_true("TargetSequence" %in% colnames(s1_table))
  expect_true("Indels" %in% colnames(s1_table))
  expect_true("BasesChanged" %in% colnames(s1_table))

  # Expected results and comparison
  anticipated_indels <- c("WT", "-1", "+1", NA, "+2, -1", "-1", "-1", "+1")
  anticipated_bc <- c("WT", NA, NA, 3, NA, 2, NA, NA)

  expect_equal(s1_table$Indels, anticipated_indels)
  expect_equal(s1_table$BasesChanged, anticipated_bc)
})


# ------------------------
# Test 2: multiple WT error
# ------------------------

test_that("tas_label_table produces error if multiple WT sequences are detected", {
  fcopy <- make_fcopy(tas_label_table)

  Sample.Names <- c("S1")
  seqs_list <- list(S1 = c("ACGTACGT", # WT
                           "ACGTACGT")) # WT 2

  Sequence.Table.List <- make_sequence_table_list(Sample.Names, seqs_list)
  Reference.Sequences.DNA <- make_ref_dna_list(Sample.Names, list(S1 = "ACGTACGT"))

  Config.List <- list(nCores = 1, protein.mutations = 0)

  expect_error(
    fcopy(Sample.Names = Sample.Names,
          Sequence.Table.List = Sequence.Table.List,
          Reference.Sequences.DNA = Reference.Sequences.DNA,
          Config.List = Config.List),
    regexp = "More than 1 WT sequence detected.",
    fixed = FALSE
  )
})


# ------------------------
# Test 3: protein mutation branch
# ------------------------
test_that("tas_label_table protein branch produces AA and ProteinMutation columns", {
  fcopy <- make_fcopy(tas_label_table)

  Sample.Names <- c("S1")
  # Choose a reference coding sequence (multiple of 3). Example: ATG GCT = M A
  # We'll create one WT and one with a single nucleotide substitution that changes amino acid
  ref_seq <- "ATGGCTGAG"  # translates to "MAE"
  seqs_list <- list(S1 = c("ATGGCTGAG",
                           "ATGACTGAG", # mutation G->A at pos 4 causing A2T
                           "ATCGCTGAT", # mutation G->C at post 3 and G->T at pos 9 causing M1I, E3D
                           "ATGGCT", # indel multiple of 3 causing E3-
                           "ATGTAGGAG", # mutation GCT->TAG at 3-6 nonsense mutation
                           "ATGCTGAG" # indel -1 at pos 4 frameshift
                           ))
  Sequence.Table.List <- make_sequence_table_list(Sample.Names, seqs_list)
  Reference.Sequences.DNA <- make_ref_dna_list(Sample.Names, list(S1 = ref_seq))

  cfg <- list(nCores = 1, protein.mutations = 1)

  out <- fcopy(Sample.Names = Sample.Names,
               Sequence.Table.List = Sequence.Table.List,
               Reference.Sequences.DNA = Reference.Sequences.DNA,
               Config.List = cfg)

  # Expect protein alignments appended to output list
  expect_true(all(c("S1","AlignProtein") %in% names(out)))
  expect_true(length(out$AlignProtein$S1) == 6)

  # Expect AA and ProteinMutation columns appended to Sequence.Table[[S1]]
  s1_table <- out$S1
  expect_true("AA" %in% colnames(s1_table))
  expect_true("ProteinMutation" %in% colnames(s1_table))

  # Sequences and labeling as expected
  anticipated_aa <- c("MAE", "MTE", "IAD", "MA", "M*", "ML")
  anticipated_pmut <- c("WT", "A2T", "M1I, E3D", "Indel", "Nonsense", "Indel")

  expect_equal(s1_table$AA, anticipated_aa)
  expect_equal(s1_table$ProteinMutation, anticipated_pmut)
})



