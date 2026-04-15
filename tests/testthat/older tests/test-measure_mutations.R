################################################################################
########################## TAS_MEASURE_MUTATIONS ###############################
################################################################################

skip_if_not(exists("tas_measure_mutations", mode = "function"), "tas_measure_mutations() not found in this session")

# The function relies on Bioconductor and openxlsx for workbook creation.
skip_if_not_installed("Biostrings")
skip_if_not_installed("S4Vectors")
skip_if_not_installed("BiocGenerics")
skip_if_not_installed("openxlsx")


# Helper to return a copy of the function with deterministic behavior (if needed)
make_fcopy <- function(fn) {
  # function uses lapply internally, so we can use fn directly
  fn
}

# Helper to build Sequence.Table.List in the exact shape the function expects:
# - Sequence.Table.List[[sample]] is a data.frame with at least 2 columns (TargetSequence, Reads)
# - Sequence.Table.List$AlignDNA[[sample]] is a DNAStringSet corresponding to each TargetSequence
# - Optionally AlignProtein list for protein branch
make_sequence_table_list_for_sample <- function(sample_name, seqs, counts, ref_dna, add_protein = FALSE) {
  # seqs: character vector of sequences (TargetSequence)
  # counts: numeric vector same length as seqs (reads per sequence)
  if(length(seqs) != length(counts)){
    stop("seqs and counts must match length")
    }
  df <- data.frame(TargetSequence = as.character(seqs), Reads = as.integer(counts))
  out <- list()
  out[sample_name] <- list(df)

  # AlignDNA must be a list within Sequence.Table.List named AlignDNA, later used as Sequence.Table.List$AlignDNA[[x]]
  seqs_dna <- DNAStringSet(seqs)
  out$AlignDNA <- list()
  out$AlignDNA[sample_name] <- list(pairwiseAlignment(seqs_dna, ref_dna))

  if (add_protein) {
    # translate each DNA sequence (remove gaps/hyphens if present)
    prots <- suppressWarnings(translate(DNAStringSet(gsub("-", "", seqs))))
    out$AlignProtein <- list()
    ref_prot <- suppressWarnings(translate(DNAStringSet(ref_dna)))
    out$AlignProtein[sample_name] <- list(pairwiseAlignment(prots, ref_prot))
  }
  return(out)
}

# Test 1: basic DNA mutagenesis computed correctly
test_that("tas_measure_mutations computes DNA mutagenesis and returns workbooks (no protein, no AID)", {
  skip_if_not_installed("Biostrings")
  fcopy <- make_fcopy(tas_measure_mutations)

  Sample.Names <- "S1"

  # Reference sequence of length 6
  ref_dna <- "AAAAAA"

  # Two observed unique sequences: one WT, one mutated at position 3 (A->C)
  seqs <- c("AAAAAA", "AACAAA")  # note second has 'C' at position 3 (indexing 1-based)
  counts <- c(5L, 1L)            # majority WT, one mutated

  # Build Sequence.Table.List with AlignDNA and table entry
  Sequence.Table.List <- make_sequence_table_list_for_sample("S1", seqs, counts, ref_dna, add_protein = FALSE)

  # Input.DataFrame row with required columns; Antibody = "no" to skip antibody region code
  Input.DataFrame <- data.frame(ReferenceSequence = "AAAAAA",
                                Antibody = "no",
                                FR1Start = NA,
                                CDR1Start = NA,
                                FR2Start = NA,
                                CDR2Start = NA,
                                FR3Start = NA,
                                CDR3Start = NA,
                                FR4Start = NA,
                                stringsAsFactors = FALSE)
  rownames(Input.DataFrame) <- Sample.Names

  Config.List <- list(nCores = 1,
                      protein.mutations = 0,
                      measure.shm = 0,
                      dna.repair.pathways = 0)

  out <- fcopy(Sample.Names = Sample.Names,
               Sequence.Table.List = Sequence.Table.List,
               AID.Targets = list(),
               Input.DataFrame = Input.DataFrame,
               Config.List = Config.List)

  # Output is a list per sample plus Workbooks
  expect_true(is.list(out))
  expect_true("S1" %in% names(out))
  expect_true(is.numeric(out$S1$MutagenesisDNA))
  # length of mutagenesis vector should equal reference length
  expect_equal(length(out$S1$MutagenesisDNA), nchar(ref_dna))

  # mutated position is the 3rd base (we expect >0 mutagenesis there)
  expect_true(out$S1$MutagenesisDNA[3] > 0)
  # other positions should be zero (or near zero) given constructed data
  expect_true(all(out$S1$MutagenesisDNA[-3] == 0))

  # Workbooks present in returned list
  expect_true(is.list(out$Workbooks))
  expect_true("Mut_pos_wb" %in% names(out$Workbooks))
  expect_true("Sequences_wb" %in% names(out$Workbooks))
})

# Test 2: protein mutation branch
test_that("tas_measure_mutations computes protein mutagenesis when protein.mutations=1", {
  skip_if_not_installed("Biostrings")
  fcopy <- make_fcopy(tas_measure_mutations)

  Sample.Names <- "S1"
  # Choose reference coding sequence length multiple of 3
  ref_dna <- "ATGCGC"  # translates to "MR" (two amino acids)
  # Create sequences: WT and one with a single nt change causing amino acid change at codon 2
  seqs <- c("ATGCGC", "ATGCAC")  # second sequence pos 5 G->A cuasing R2H mutation
  counts <- c(3L, 1L)
  Sequence.Table.List <- make_sequence_table_list_for_sample("S1", seqs, counts, ref_dna, add_protein = TRUE)

  Input.DataFrame <- data.frame(ReferenceSequence = "ATGCGC",
                                Antibody = "no",
                                FR1Start = NA,
                                CDR1Start = NA,
                                FR2Start = NA,
                                CDR2Start = NA,
                                FR3Start = NA,
                                CDR3Start = NA,
                                FR4Start = NA,
                                stringsAsFactors = FALSE)
  rownames(Input.DataFrame) <- Sample.Names

  Config.List <- list(nCores = 1,
                      protein.mutations = 1,
                      measure.shm = 0,
                      dna.repair.pathways = 0)

  out <- fcopy(Sample.Names = Sample.Names,
               Sequence.Table.List = Sequence.Table.List,
               AID.Targets = NULL,
               Input.DataFrame = Input.DataFrame,
               Config.List = Config.List)

  # MutagenesisProtein should be present in sample output
  expect_true(is.list(out$S1))
  expect_true("MutagenesisDNA" %in% names(out$S1))
  expect_true("MutagenesisProtein" %in% names(out$S1))

  # the MutagenesisProtein element should be numeric (percentage per AA position)
  if (is.list(out$S1$MutagenesisProtein)) {
    # some code paths return list(MutagenesisProtein = numeric_vector)
    mp <- out$S1$MutagenesisProtein[[1]]
  } else {
    mp <- out$S1$MutagenesisProtein
  }
  expect_true(is.numeric(mp))
  # length of protein mutagenesis vector equals number of translated AAs
  expect_equal(length(mp), nchar(translate(DNAString(ref_dna))))

  # expect mutation at protein position 2 frequency >0
  expect_true(out$S1$MutagenesisProtein[2] > 0)
  # other positions should be zero (or near zero) given constructed data
  expect_true(all(out$S1$MutagenesisProtein[-2] == 0))
})



# Test 3: AID motif measurements (measure.shm = 1) including protein
test_that("tas_measure_mutations summarizes AID motif mutagenesis when measure.shm=1", {
  skip_if_not_installed("Biostrings")
  fcopy <- make_fcopy(tas_measure_mutations)

  Sample.Names <- "S1"
  ref_dna <- "AAGCTA" # AGCT = WRCY/WRCH
  # Two sequences, one mutated at position 4 (AID cytosine)
  seqs <- c("AAGCTA", "AAGTTA")
  counts <- c(7L, 3L)
  Sequence.Table.List <- make_sequence_table_list_for_sample("S1", seqs, counts, ref_dna, add_protein = TRUE)

  # Construct AID.Targets for S1: WRCH group should include a motif that covers position 4
  # Provide a small WRCH data.frame with Start, End, Cytosine to cover position 4
  wrc_df <- data.frame(Motif = character(0), Start = integer(0), End = integer(0), Cytosine = integer(0))
  wrch_df <- data.frame(Motif = "WRCH", Start = 2L, End = 5L, Cytosine = 4L)
  wrcy_df <- data.frame(Motif = "WRCY", Start = 2L, End = 5L, Cytosine = 4L)

  # The function expects AID.Targets[[x]] to be a list with $WRCY and $WRCH data.frames.
  AID.Targets <- list()
  AID.Targets$S1 <- list(WRCY = wrcy_df, WRCH = wrch_df)

  Input.DataFrame <- data.frame(ReferenceSequence = "AAGCTA",
                                Antibody = "no",
                                FR1Start = NA,
                                CDR1Start = NA,
                                FR2Start = NA,
                                CDR2Start = NA,
                                FR3Start = NA,
                                CDR3Start = NA,
                                FR4Start = NA)
  rownames(Input.DataFrame) <- Sample.Names

  Config.List <- list(nCores = 1,
                      protein.mutations = 1,
                      measure.shm = 1,
                      dna.repair.pathways = 0,
                      read.frequency.limit = 0)

  out <- fcopy(Sample.Names = Sample.Names,
               Sequence.Table.List = Sequence.Table.List,
               AID.Targets = AID.Targets,
               Input.DataFrame = Input.DataFrame,
               Config.List = Config.List)

  # Output contains AID in sample result
  expect_true(is.list(out$S1))
  expect_true("AID" %in% names(out$S1))

  # AID should be a list with elements for WRCY and WRCH
  aid_res <- out$S1$AID
  expect_true(is.list(aid_res))
  expect_true(all(c("WRCY","WRCH") %in% names(aid_res)))

  # WRCH positions returned should include Cytosine and MotifMutagenesis columns
  wrch_positions <- aid_res$WRCH$Positions
  expect_true(is.data.frame(wrch_positions))
  expect_true(all(c("CytosineMutagenesis", "MotifMutagenesis") %in% colnames(wrch_positions)))

  wrcy_positions <- aid_res$WRCY$Positions
  expect_true(is.data.frame(wrcy_positions))
  expect_true(all(c("CytosineMutagenesis", "MotifMutagenesis") %in% colnames(wrcy_positions)))

  # MotifSums should be a data.frame present
  expect_true(is.data.frame(aid_res$WRCH$MotifSums))
  anticipated_sums <- data.frame(TotalPercentMutated = 5,
                                 TotalPercentMutatedInMotif = 5,
                                 TotalPercentMutatedInCytosine = 5,
                                 MotifMutationAverage = 30,
                                 MotifMutationFrequencyofTotal = 100,
                                 MotifMutationNonCAverage = 0,
                                 MotifMutationNonCFrequencyOfTotal = 0,
                                 CytosineMutationAverage = 30,
                                 NonCytosineMutationAverage = 0,
                                 CytosineMutationFrequencyOfTotal = 100,
                                 CytosineMutationFrequencyOfMotif = 100)
  expect_equal(aid_res$WRCH$MotifSums, anticipated_sums)

  # MutCyt and MutNonCyt vectors should be numeric
  expect_true(is.numeric(aid_res$WRCH$MutCyt))
  expect_true(is.numeric(aid_res$WRCH$MutNonCyt))

  # WRCY results should be the same as WRCH
  expect_equal(aid_res$WRCH$Positions[,-(colnames(aid_res$WRCH$Positions) %in% "Motif")],
               aid_res$WRCY$Positions[,-(colnames(aid_res$WRCY$Positions) %in% "Motif")])
  expect_equal(aid_res$WRCH[-(names(aid_res$WRCH) %in% "Positions")],
               aid_res$WRCY[-(names(aid_res$WRCY) %in% "Positions")])

  # Workbooks should include sheets created for Mut_pos and Sequences
  expect_true(is.list(out$Workbooks))
  expect_true("Mut_pos_wb" %in% names(out$Workbooks))
  expect_true("Sequences_wb" %in% names(out$Workbooks))
})


# Helper to quickly build Sequence.Table.List in the shape the function expects:
# - Sequence.Table.List[[sample]] is a data.frame with at least two columns: TargetSequence, Reads (counts)
# - Sequence.Table.List$AlignDNA[[sample]] is a DNAStringSet with one element per TargetSequence
# - Optionally AlignProtein for protein tests
# make_stl <- function(sample, seqs, counts, add_protein = FALSE) {
#   if (length(seqs) != length(counts)) stop("seqs and counts must match")
#   df <- data.frame(TargetSequence = as.character(seqs),
#                    Reads = as.integer(counts),
#                    stringsAsFactors = FALSE)
#   stl <- list()
#   stl[[sample]] <- df
#   stl$AlignDNA <- list()
#   stl$AlignDNA[[sample]] <- DNAStringSet(seqs)
#   if (add_protein) {
#     stl$AlignProtein <- list()
#     stl$AlignProtein[[sample]] <- suppressWarnings(translate(DNAStringSet(gsub("-", "", seqs))))
#   }
#   return(stl)
# }

# ------------------------------
# Test: antibody-region branch (DNA + protein)
# ------------------------------
test_that("tas_measure_mutations antibody-region branch and protein-region sub-branch", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ggseqlogo")

  fcopy <- make_fcopy(tas_measure_mutations)

  # set up
  Sample.Names <- "S1"

  # Reference of length 21 (for easy region partitioning); will be used as WT
  ref_dna <- "ATGAAACCCGGGAAGCTACAGATGAAACCCGGGAAGCTACAG"  # length 42 -> 14 AA codons for 7 antibody regions
  # create TargetSequences: majority WT and one mutated at pos 5
  seqs <- c(ref_dna, paste0(substr(ref_dna,1,4),"T",substr(ref_dna,6,nchar(ref_dna)))) # mutate position 5
  counts <- c(10L, 2L)

  Sequence.Table.List <- make_sequence_table_list_for_sample("S1", seqs, counts, ref_dna, add_protein = TRUE)

  # Build Input.DataFrame with Antibody = "yes" and strictly increasing region starts
  # We'll choose coordinates such that regions are in frame and 1 amino acid each
  Input.DataFrame <- data.frame(
    ReferenceSequence = ref_dna,
    Antibody = "yes",
    FR1Start = 1,
    CDR1Start = 7,
    FR2Start = 13,
    CDR2Start = 19,
    FR3Start = 25,
    CDR3Start = 31,
    FR4Start = 37
  )
  rownames(Input.DataFrame) <- Sample.Names

  # Config: enable protein.mutations so the protein-region branch runs
  Config.List <- list(nCores = 1, protein.mutations = 1, measure.shm = 0, dna.repair.pathways = 0)

  # Run function (should not error)
  out <- tas_measure_mutations(Sample.Names = Sample.Names,
                               Sequence.Table.List = Sequence.Table.List,
                               AID.Targets = NULL,
                               Input.DataFrame = Input.DataFrame,
                               Config.List = Config.List)

  # Assert output structure
  expect_true(is.list(out))
  expect_true("S1" %in% names(out))
  sample_out <- out$S1

  # RegionCoordinates and RegionMutagenesisDNA should be present for antibody case
  expect_true("RegionCoordinates" %in% names(sample_out))
  expect_true("RegionMutagenesisDNA" %in% names(sample_out))

  # RegionCoordinates should be numeric and contain the named regions and END1
  coords <- sample_out$RegionCoordinates
  expect_true(all(c("FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4","END1") %in% names(coords)))
  expect_true(is.numeric(coords["FR1"]))  # sanity

  # Each region entry in RegionMutagenesisDNA should be a list element with Mut_Pos and Mut_Freq etc.
  region_list <- sample_out$RegionMutagenesisDNA
  expect_true(is.list(region_list))
  expect_true(length(region_list) == 7) # FR1..FR4 except END1 (function sets Ab_regions length)
  # Check first region contains expected fields
  first_region <- region_list[[1]]
  expect_true(all(c("Mut_Pos","Mut_Freq","Mut_Avg","Mut_Expected","Mut_Norm","Logo") %in% names(first_region)))

  # Because protein.mutations = 1, RegionMutagenesisProtein should be present
  expect_true("RegionMutagenesisProtein" %in% names(sample_out))
  prot_regions <- sample_out$RegionMutagenesisProtein
  expect_true(is.list(prot_regions))
  expect_true(length(prot_regions) == 7)
  expect_true(all(sapply(prot_regions, function(x) "Mut_Pos" %in% names(x))))
})





# ------------------------------
# Test: dna.repair.pathways = 1 branch
# ------------------------------
test_that("tas_measure_mutations dna.repair.pathways=1 produces mut_types workbook and sums to ~100", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("openxlsx")

  fcopy <- make_fcopy(tas_measure_mutations)

  # prepare sample and Sequence.Table.List entries that include Indels, BasesChanged, Percent
  Sample.Names <- c("S1","S2")

  # We'll make small reference sequences
  ref1 <- "AAGCTA"
  ref2 <- "ATGAAACCCGGGAAGCTACAGATGAAACCCGGGAAGCTACAG"

  # Build Sequence.Table.List with required columns for each sample
  Sequence.Table.List <- list()
  Sequence.Table.List$AlignDNA <- list()
  # S1 rows: WT, NHEJ(-1), MMEJ(-5), Change only (WT indel but BasesChanged>0), Indel+BaseChange
  seqs_s1 <- c("AAGCTA","AACTA","ATA","TAGCCA","AACTG")
  reads_s1 <- c(3,2,2,2,1)
  indels_s1 <- c("WT","-1","-3","","-1")
  baseschanged_s1 <- c(0,0,0,2,1)
  percent_s1 <- c(30,20,20,20,10) # sums to 100
  Sequence.Table.List[["S1"]] <- data.frame(TargetSequence = seqs_s1,
                            Reads = reads_s1,
                            Indels = indels_s1,
                            BasesChanged = baseschanged_s1,
                            Percent = percent_s1,
                            stringsAsFactors = FALSE)
  Sequence.Table.List$AlignDNA[["S1"]] <- pairwiseAlignment(DNAStringSet(seqs_s1), DNAString(ref1))

  # S2 simple WT-only sample
  seqs_s2 <- c("ATGAAACCCGGGAAGCTACAGATGAAACCCGGGAAGCTACAG")
  reads_s2 <- c(10)
  indels_s2 <- c("WT")
  baseschanged_s2 <- c(0)
  percent_s2 <- c(100)
  Sequence.Table.List[["S2"]] <- data.frame(TargetSequence = seqs_s2,
                            Reads = reads_s2,
                            Indels = indels_s2,
                            BasesChanged = baseschanged_s2,
                            Percent = percent_s2)
  Sequence.Table.List$AlignDNA[["S2"]] <- pairwiseAlignment(DNAStringSet(seqs_s2), DNAString(ref2))

  # Input.DataFrame rows
  Input.DataFrame <- data.frame(
    ReferenceSequence = c(ref1, ref2),
    Antibody = c("no","no"),
    FR1Start = c(NA, NA), CDR1Start = c(NA, NA), FR2Start = c(NA, NA),
    CDR2Start = c(NA, NA), FR3Start = c(NA, NA), CDR3Start = c(NA, NA),
    FR4Start = c(NA, NA))
  rownames(Input.DataFrame) <- Sample.Names

  # Config with dna.repair.pathways enabled
  Config.List <- list(nCores = 1, protein.mutations = 0, measure.shm = 0, dna.repair.pathways = 1)

  # Run function
  out <- fcopy(Sample.Names = Sample.Names,
               Sequence.Table.List = Sequence.Table.List,
               AID.Targets = list(),
               Input.DataFrame = Input.DataFrame,
               Config.List = Config.List)

  # The function should have returned Output.Mutagenesis list
  expect_true(is.list(out))

  # Should include Workbooks and specifically mut_types_wb because dna.repair.pathways=1
  expect_true("Workbooks" %in% names(out))
  expect_true("mut_types_wb" %in% names(out$Workbooks))

  # check results of classification analysis
  anticipated_sums <- data.frame(c("WT", "NHEJ", "MMEJ", "Base Change", "Indel + Base Change", "Other"),
                                 S1 = c(30, 20, 20, 20, 10, 0),
                                 S2 = c(100, 0, 0, 0, 0, 0))
  colnames(anticipated_sums) <- c("", Sample.Names)
  expect_equal(readWorkbook(out$Workbooks$mut_types_wb), anticipated_sums)

})
