
################################################################################
################################ TAS_FILTER ####################################
################################################################################

skip_if_not(exists("tas_filter", mode = "function"), "tas_filter() not found in this session")

# Helper: minimal valid Input.DataFrame for tas_filter
make_valid_input_df_filter <- function(sample_name = "S1",
                                       fwd_ext = "AAA", fwd_primer = "AC",
                                       rev_ext = "GGG", rev_primer = "TT",
                                       amp_len = 150, max_del = 0, max_ins = 0,
                                       fwd_type = "Barcode", rev_type = "Barcode") {
  data.frame(
    SampleName = sample_name,
    ForwardExtensionType = fwd_type,
    ForwardExtension = fwd_ext,
    ForwardPrimer = fwd_primer,
    ReverseExtensionType = rev_type,
    ReverseExtension = rev_ext,
    ReversePrimer = rev_primer,
    AmpliconLength = amp_len,
    MaxDeletion = max_del,
    MaxInsertion = max_ins,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

make_valid_config_filter <- function(wd = tempdir(), nCores = 1) {
  list(WorkingDirectory = wd, nCores = nCores)
}

# 1) Input column validation
test_that("tas_filter errors when required input columns are missing", {
  df <- make_valid_input_df_filter()
  df$ForwardPrimer <- NULL
  cfg <- make_valid_config_filter()

  expect_error(
    tas_filter(Sample.Names = "S1",
               Reads.List = list(), # empty - won't get that far
               Input.DataFrame = df,
               Config.List = cfg),
    regexp = "The following columns were not found in the input file",
    ignore.case = TRUE
  )
})

# 2) Config entries validation
test_that("tas_filter errors when required config entries are missing", {
  df <- make_valid_input_df_filter()
  cfg <- list(nCores = 1) # missing WorkingDirectory

  expect_error(
    tas_filter(Sample.Names = "S1",
               Reads.List = list(),
               Input.DataFrame = df,
               Config.List = cfg),
    regexp = "The following columns were not found in the input file",
    ignore.case = TRUE
  )
})

# The following tests require Biostrings and ShortRead - skip if not installed
skip_if_not_installed("Biostrings")
skip_if_not_installed("ShortRead")

# 3) Non-DNA primer check triggers error
test_that("tas_filter rejects non-DNA characters in ForwardPrimer / ReversePrimer", {
  df_bad_fwd <- make_valid_input_df_filter(fwd_primer = "A1")  # contains '1' -> invalid
  cfg <- make_valid_config_filter()

  # Provide a minimal Reads.List object to reach the primer check
  # Create a ShortReadQ with one dummy read
  s <- DNAStringSet("ACGTACGT")
  q <- BStringSet(rep("IIIIIIII", length(s)))
  id <- BStringSet("r1")
  srq <- ShortReadQ(sread = s, quality = q, id = id)
  reads_list <- list(S1 = srq)

  expect_error(
    tas_filter(Sample.Names = "S1",
               Reads.List = reads_list,
               Input.DataFrame = df_bad_fwd,
               Config.List = cfg),
    regexp = "ForwardPrimer must be a DNA sequence",
    ignore.case = TRUE
  )

  df_bad_rev <- make_valid_input_df_filter(rev_primer = "N*") # '*' invalid
  expect_error(
    tas_filter(Sample.Names = "S1",
               Reads.List = reads_list,
               Input.DataFrame = df_bad_rev,
               Config.List = cfg),
    regexp = "ReversePrimer must be a DNA sequence",
    ignore.case = TRUE
  )
})


test_that("tas_filter multi-sample behavior: filters each sample independently", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ShortRead")
  skip_if_not(exists("tas_filter", mode = "function"))

  # Copy function and set up a mock environment to control mclapply and DNA_ALPHABET
  fcopy <- tas_filter
  env_mock <- new.env(parent = environment(tas_filter))

  # PROVIDE DNA_ALPHABET as a VECTOR of single characters (exact 15-letter alphabet requested)
  env_mock$DNA_ALPHABET <- strsplit("ACGTMRWSYKVHDBN", "")[[1]]

  # Force mclapply to behave like lapply for predictable test behavior
  env_mock$mclapply <- function(X, FUN, mc.cores = 1, ...) lapply(X, FUN, ...)

  # Keep original parent for other lookups (stringr, Biostrings, etc.)
  environment(fcopy) <- env_mock

  # Define Sample.Names
  Sample.Names <- c("S1", "S2", "S3")

  # Prepare Input.DataFrame with rownames matching Sample.Names (function indexes by name)
  df_rows <- data.frame(
    SampleName = Sample.Names,
    ForwardExtensionType = c("Barcode","Barcode","Barcode"),
    ForwardExtension = c("AAA", "AAA", "AT"),
    ForwardPrimer = c("AC", "AC", "G"),
    ReverseExtensionType = c("Barcode","Barcode","Barcode"),
    ReverseExtension = c("GGG", "CCC", "TT"),
    ReversePrimer = c("TT", "GG", "AA"),
    AmpliconLength = c(150, 150, 150),
    MaxDeletion = c(0, 0, 0),
    MaxInsertion = c(0, 0, 0),
    stringsAsFactors = FALSE
  )
  # ensure rows are addressable by sample name
  rownames(df_rows) <- Sample.Names

  cfg <- list(WorkingDirectory = tempdir(), nCores = 1)

  # Build Reads.List: named list of ShortReadQ objects
  # S1: 2 reads, both contain both patterns -> expect 2 returned
  seqs_s1 <- DNAStringSet(c(
    paste0("AAAAC", "GGG", "GGG", "TTT", "AACCC"),  # contains both
    paste0("AAAAC", "AACCC")                       # contains both
  ))
  q_s1 <- BStringSet(sapply(seqs_s1,function(x){
    str_flatten(rep("I", nchar(x)))
  }))
  id_s1 <- BStringSet(paste0("s1_read", seq_along(seqs_s1)))
  srq_s1 <- ShortReadQ(sread = seqs_s1, quality = q_s1, id = id_s1)

  # S2: 2 reads, contain forward pattern only (no reverse complement) -> expect 0 returned
  seqs_s2 <- DNAStringSet(c(
    paste0("AAAAC", "GGG", "GGG", "TTT"),
    "AAAACGGG"
  ))
  q_s2 <- BStringSet(sapply(seqs_s2,function(x){
    str_flatten(rep("I", nchar(x)))
  }))
  id_s2 <- BStringSet(paste0("s2_read", seq_along(seqs_s2)))
  srq_s2 <- ShortReadQ(sread = seqs_s2, quality = q_s2, id = id_s2)

  # S3: 3 reads, only one contains both patterns -> expect 1 returned
  seqs_s3 <- DNAStringSet(c(
    "ATGGGGATGCCC",       # contains forward only
    "ATGCCCGGTTAA",         # contains both (has ATG and TTAA)
    "CCCCGGGG"         # contains none
  ))
  q_s3 <- BStringSet(sapply(seqs_s3,function(x){
    str_flatten(rep("I", nchar(x)))
  }))
  id_s3 <- BStringSet(paste0("s3_read", seq_along(seqs_s3)))
  srq_s3 <- ShortReadQ(sread = seqs_s3, quality = q_s3, id = id_s3)

  Reads.List <- list(S1 = srq_s1, S2 = srq_s2, S3 = srq_s3)

  # Run filter
  out <- fcopy(Sample.Names = Sample.Names,
               Reads.List = Reads.List,
               Input.DataFrame = df_rows,
               Config.List = cfg)

  # Basic structure checks
  expect_type(out, "list")
  expect_named(out, Sample.Names)

  # Check per-sample counts
  expect_true(length(out$S1) == 2)  # both reads matched
  expect_true(length(out$S2) == 0)  # none matched
  expect_true(length(out$S3) == 1)  # only one matched

  # Inspect returned sequences to ensure they include expected patterns
  s1_returned_seqs <- as.character(out$S1@sread)
  expect_true(all(grepl("AAAAC", s1_returned_seqs)))
  expect_true(all(grepl("AACCC", s1_returned_seqs)))

  s3_returned_seq <- as.character(out$S3@sread[1])
  expect_true(grepl("ATG", s3_returned_seq))
  expect_true(grepl("TTAA", s3_returned_seq) || grepl("TTAA", as.character(reverseComplement(DNAStringSet(s3_returned_seq)))),
              info = "S3 returned sequence must contain forward pattern and reverse complement pattern")
})

