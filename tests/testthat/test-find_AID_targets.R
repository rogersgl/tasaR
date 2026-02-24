################################################################################
########################### TAS_FIND_AID_TARGETS ###############################
################################################################################

skip_if_not(exists("tas_find_AID_Targets", mode = "function"), "tas_find_AID_Targets() not found in this session")
skip_if_not_installed("Biostrings")

# helper: copy function and make mclapply deterministic
make_fcopy <- function(fn) {
  fcopy <- fn
  env_mock <- new.env(parent = environment(fn))
  env_mock$mclapply <- function(X, FUN, mc.cores = 1, ...) lapply(X, FUN, ...)
  environment(fcopy) <- env_mock
  fcopy
}

test_that("tas_find_AID_Targets finds WRCY/RGYW and overlapping WRCH/DGYW motifs correctly (single sample)", {
  fcopy <- make_fcopy(tas_find_AID_Targets)

  Sample.Names <- "S1"
  # Construct a reference sequence containing:
  # - "AACT" (WRCY) at positions 3-6 (TT A A C T GG...) => cyt = e-1 = 5
  # - "AGCA" (DGYW) at positions 9-12 => cyt = e-2 = 10
  ref_seq <- DNAString("TTAACTGGAGCA")
  Reference.Sequences.DNA <- list(S1 = ref_seq)

  cfg <- list(nCores = 1)

  out <- fcopy(Sample.Names = Sample.Names,
               Reference.Sequences.DNA = Reference.Sequences.DNA,
               Config.List = cfg)

  # Output should be a list with sample name
  expect_true(is.list(out))
  expect_true("S1" %in% names(out))

  s1 <- out$S1
  # For S1, s1 should be a list with components WRCY and WRCH
  expect_true(is.list(s1))
  expect_true(all(c("WRCY","WRCH") %in% names(s1)))

  # WRCY table should contain the WRCY/RGYW matches (AACT/AGCA)
  wrcy_df <- s1$WRCY
  expect_true(is.data.frame(wrcy_df))

  # check # of rows
  expect_true(nrow(wrcy_df) == 2)

  # find the row corresponding to motif "WRCY"
  row_wrcy <- wrcy_df[wrcy_df$Motif == "WRCY", , drop = FALSE]
  expect_true(nrow(row_wrcy) == 1)

  #find the row corresponding to motif "RGYW"
  row_rgyw <- wrcy_df[wrcy_df$Motif == "RGYW", , drop = FALSE]
  expect_true(nrow(row_rgyw) == 1)

  # check Start/End and Cytosine calculation: for AACT at positions 3-6, cyt = 6-1 = 5
  expect_equal(as.integer(row_wrcy$Start), 3L)
  expect_equal(as.integer(row_wrcy$End), 6L)
  expect_equal(as.integer(row_wrcy$Cytosine), 5L)

  # check Start/End and Cytosine calculation: for AGCA at positions 9-12, cyt = 12-2 = 10
  expect_equal(as.integer(row_wrcy$Start), 9L)
  expect_equal(as.integer(row_wrcy$End), 12L)
  expect_equal(as.integer(row_wrcy$Cytosine), 10L)

  # WRCH table should include DGYW matches (function groups DGYW into WRCH output)
  wrch_df <- s1$WRCH
  expect_true(is.data.frame(wrch_df))

  # find WRCH rows (motif column should equal "WRCH")
  row_wrch <- wrch_df[wrch_df$Motif == "WRCH", , drop = FALSE]
  expect_true(nrow(row_wrch) == 2)

  # find DGYW row (motif column should equal "DGYW")
  row_dgyw <- wrch_df[wrch_df$Motif == "DGYW", , drop = FALSE]
  expect_true(nrow(row_dgyw) == 1)

  # For AACT at positions 3-6, cyt = 6-1 = 5
  expect_equal(as.integer(row_wrch$Start[1]), 3L)
  expect_equal(as.integer(row_wrch$End[1]), 6L)
  expect_equal(as.integer(row_wrch$Cytosine[1]), 5L)

  # For AGCA (WRCH) at positions 9-12, cyt = e-1 = 12-1 = 11
  expect_equal(as.integer(row_wrch$Start[2]), 9L)
  expect_equal(as.integer(row_wrch$End[2]), 12L)
  expect_equal(as.integer(row_wrch$Cytosine[2]), 11L)

  # For AGCA (overlapping DGYW) at positions 9-12, cyt = e-2 = 12-2 = 10
  expect_equal(as.integer(row_dgyw$Start), 9L)
  expect_equal(as.integer(row_dgyw$End), 12L)
  expect_equal(as.integer(row_dgyw$Cytosine), 10L)

  # Rows in each returned df should be ordered by Start (increasing)
  expect_equal(order(wrcy_df$Start), seq_len(nrow(wrcy_df)))
  expect_equal(order(wrch_df$Start), seq_len(nrow(wrch_df)))
})

test_that("tas_find_AID_Targets supports returns empty data frames for no matches", {
  fcopy <- make_fcopy(tas_find_AID_Targets)

  Sample.Names <- c("S1","S2")
  # S1 contains a WRCY motif (AACT)
  ref1 <- DNAString("TTAACTGG")
  # S2 contains no AID motifs (use homopolymer of G)
  ref2 <- DNAString("GGGGGGGGGG")

  Reference.Sequences.DNA <- list(S1 = ref1, S2 = ref2)
  cfg <- list(nCores = 1)

  out <- fcopy(Sample.Names = Sample.Names,
               Reference.Sequences.DNA = Reference.Sequences.DNA,
               Config.List = cfg)

  expect_true(is.list(out))
  expect_true(all(Sample.Names %in% names(out)))

  # S1 should have one WRCY match
  w1 <- out$S1$WRCY
  expect_true(is.data.frame(w1))
  expect_true(nrow(w1) >= 1)
  expect_true(any(w1$Motif == "WRCY"))

  # S2 should return empty data.frames (0 rows) for both WRCY and WRCH
  w2_wrcy <- out$S2$WRCY
  w2_wrch <- out$S2$WRCH
  expect_true(is.data.frame(w2_wrcy))
  expect_true(is.data.frame(w2_wrch))
  expect_equal(nrow(w2_wrcy), 0)
  expect_equal(nrow(w2_wrch), 0)
})
