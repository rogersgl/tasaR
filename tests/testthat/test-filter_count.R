################################################################################
############################# TAS_FILTER_COUNT #################################
################################################################################

skip_if_not(exists("tas_filter_count", mode = "function"), "tas_filter_count() not found in this session")

# These tests require ShortRead and Bioconductor helpers
skip_if_not_installed("ShortRead")
skip_if_not_installed("S4Vectors")
skip_if_not_installed("BiocGenerics")

# Helper: write a valid FASTQ file with nRecords at path (plain text)
write_fastq_n <- function(path, nRecords = 1) {
  rec <- "@SEQ_ID\nACGT\n+\n!!!!\n"
  # Ensure directory exists
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # Write nRecords consecutively
  con <- file(path, "w")
  for (i in seq_len(nRecords)) cat(rec, file = con)
  close(con)
  invisible(TRUE)
}

# Helper: write a gzipped FASTQ with nRecords
write_gz_fastq_n <- function(path, nRecords = 1) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "w")
  rec <- "@SEQ_ID\nACGT\n+\n!!!!\n"
  for (i in seq_len(nRecords)) cat(rec, file = con)
  close(con)
  invisible(TRUE)
}

# Helper to create a ShortReadQ of given length (so that length() works on it)
make_srq_n <- function(n) {
  seqs <- DNAStringSet(rep("ACGT", n))
  qs <- BStringSet(vapply(width(seqs), function(len) paste(rep("I", len), collapse = ""), character(1)))
  ids <- BStringSet(paste0("r", seq_len(n)))
  ShortReadQ(sread = seqs, quality = qs, id = ids)
}

test_that("tas_filter_count produces expected columns and results with mock data", {
  skip_if_not_installed("ShortRead")
  # Setup WD
  tmpwd <- file.path(tempdir(), "tas_filter_count_test1")
  if (dir.exists(tmpwd)){
    unlink(tmpwd, recursive = TRUE)
    }
  dir.create(tmpwd, recursive = TRUE)
  WD <- paste0(normalizePath(tmpwd, mustWork = TRUE), "/")

  # make unpaired and merged directories
  dir.create(file.path(WD, "unpaired"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(WD, "merged"), recursive = TRUE, showWarnings = FALSE)

  # S1, S2 = no UMIs
  # S3, S4 = UMIs
  # S5 = 0 merged reads
  Sample.Names <- c("S1", "S2", "S3", "S4", "S5")

  # Construct Input.DataFrame with ForwardFASTQFileName and SampleName
  Input.DataFrame <- data.frame(
    SampleName = Sample.Names,
    ForwardFASTQFileName = c("s1_R1.fastq", "s2_R1.fastq", "s3_R1.fastq", "s4_R1.fastq", "s5_R1.fastq"),
    stringsAsFactors = FALSE
  )

  # Create raw (unpaired) FASTQ files (nRecords = # of reads)
  write_fastq_n(file.path(WD, "unpaired", Input.DataFrame$ForwardFASTQFileName[1]), nRecords = 8)
  write_fastq_n(file.path(WD, "unpaired", Input.DataFrame$ForwardFASTQFileName[2]), nRecords = 5)
  write_fastq_n(file.path(WD, "unpaired", Input.DataFrame$ForwardFASTQFileName[3]), nRecords = 9)
  write_fastq_n(file.path(WD, "unpaired", Input.DataFrame$ForwardFASTQFileName[4]), nRecords = 12)
  write_fastq_n(file.path(WD, "unpaired", Input.DataFrame$ForwardFASTQFileName[5]), nRecords = 2)

  # Create merged gz FASTQ (nRecords = # of reads)
  write_gz_fastq_n(file.path(WD, "merged", paste0(Sample.Names[1], "-merged.fastq.gz")), nRecords = 7)
  write_gz_fastq_n(file.path(WD, "merged", paste0(Sample.Names[2], "-merged.fastq.gz")), nRecords = 4)
  write_gz_fastq_n(file.path(WD, "merged", paste0(Sample.Names[3], "-merged.fastq.gz")), nRecords = 7)
  write_gz_fastq_n(file.path(WD, "merged", paste0(Sample.Names[4], "-merged.fastq.gz")), nRecords = 11)
  write_gz_fastq_n(file.path(WD, "merged", paste0(Sample.Names[5], "-merged.fastq.gz")), nRecords = 1)

  # Reads.Filtered.List: lengths (Filtered) are numbers fed to make_srq_n() function
  Reads.Filtered.List <- list(
    S1 = make_srq_n(6),
    S2 = make_srq_n(3),
    S3 = make_srq_n(6),
    S4 = make_srq_n(10),
    S5 = make_srq_n(1)
  )

  # Sequence.Table.List: per-sample tables (data.frames) with UniqueSequences counts
  Sequence.Table.List <- list(
    S1 = data.frame(TargetSequence = c("ATG","GCA"),
                    Reads = c(3, 3)),
    S2 = data.frame(TargetSequence = c("ATG"),
                    Reads = 3),
    S3 = data.frame(TargetSequence = c("ATG"),
                    `UMI Count` = 2, check.names = FALSE),
    S4 = data.frame(TargetSequence = c("ATG","GCA"),
                    `UMI Count` = c(2, 1), check.names = FALSE),
    S5 = data.frame(TargetSequence = character(),
                    Reads = character())
  )
  # For All, set entries empty (no UMIs)
  Sequence.Table.List$All <- list(S1 = data.frame(),
                                  S2 = data.frame(),
                                  S3 = data.frame(TargetSequence = (c("ATG", "ATG")),
                                                  Reads = c(3, 3),
                                                  UMI = c("C","G")),
                                  S4 = data.frame(TargetSequence = (c("ATG", "ATG", "GCA")),
                                                  Reads = c(4, 3, 3),
                                                  UMI = c("C","G", "T")),
                                  S5 = data.frame())

  cfg <- list()  # not used by function except for passing

  # Call function
  res <- tas_filter_count(Sample.Names = Sample.Names,
                          Reads.Filtered.List = Reads_Filtered_List <- Reads.Filtered.List,
                          Sequence.Table.List = Sequence.Table.List,
                          Input.DataFrame = Input.DataFrame,
                          Config.List = cfg,
                          WD = WD)

  # Expected columns: Raw, Merged, Filtered, UniqueSequences
  expect_true(all(c("Raw", "Merged", "Filtered", "UMIs", "UniqueSequences") %in% colnames(res)))

  # Rownames equal sample names
  expect_equal(rownames(res), Sample.Names)

  # Validate numeric counts
  anticipated_raw <- c(8, 5, 9, 12, 2)
  anticipated_merged <- c(7, 4, 7, 11, 1)
  anticipated_filtered <- c(6, 3, 6, 10, 1)
  anticipated_UMIs <- c(NA, NA, 2, 3, NA)
  anticipated_unique <- c(2, 1, 1, 2, 0)

  expect_equal(res$Raw, anticipated_raw)
  expect_equal(res$Merged, anticipated_merged)
  expect_equal(res$Filtered, anticipated_filtered)
  expect_equal(res$UMIs, anticipated_UMIs)
  expect_equal(res$UniqueSequences, anticipated_unique)

  # cleanup
  try(unlink(tmpwd, recursive = TRUE), silent = TRUE)
})
