
################################################################################
############################## TAS_SEQ_TABLE ###################################
################################################################################

skip_if_not(exists("tas_seq_table", mode = "function"), "tas_seq_table() not found in this session")
skip_if_not_installed("Biostrings")
skip_if_not_installed("ShortRead")
skip_if_not_installed("data.table")
skip_if_not_installed("S4Vectors")
skip_if_not_installed("BiocGenerics")

# -------------------------
# Helpers
# -------------------------

# Return a copy of `fn` whose environment has deterministic mclapply = lapply
# This keeps other lookups resolving to the original function environment.
make_fcopy <- function(fn) {
  fcopy <- fn
  env_mock <- new.env(parent = environment(fn))
  env_mock$mclapply <- function(X, FUN, mc.cores = 1, ...) lapply(X, FUN, ...)
  # Provide DNA_ALPHABET vector where required by other functions (some functions expect it)
  env_mock$DNA_ALPHABET <- strsplit("ACGTMRWSYKVHDBN", "")[[1]]
  environment(fcopy) <- env_mock
  return(fcopy)
}

# Construct a ShortReadQ from a character vector of DNA sequences
make_srq <- function(seq_chars, id_prefix = "r") {
  seqs <- DNAStringSet(seq_chars)
  # Construct a simple quality string of 'I' for each sequence
  qs <- BStringSet(vapply(width(seqs), function(len) paste(rep("I", len), collapse = ""), character(1)))
  ids <- BStringSet(paste0(id_prefix, seq_along(seqs)))
  srq <- ShortReadQ(sread = seqs, quality = qs, id = ids)
  return(srq)
}

# -------------------------
# Test (uses helpers)
# -------------------------

test_that("tas_seq_table UMI branch: groups UMIs, removes UMIs with <3 reads, returns tables;
          non-UMI branch: counts unique sequences and filters by read.frequency.limit", {

            fcopy <- make_fcopy(tas_seq_table)

            Sample.Names <- c("S1", "S2")
            df <- data.frame(InsertStart = c(7, 7),
                             InsertEnd = c(9, 9),
                             ForwardExtensionType = c("UMI", ""),
                             ForwardExtension = c("NNY", ""),
                             ForwardPrimer = c("AAA", "AAA"),
                             ReverseExtensionType = c("", ""),
                             ReverseExtension = c("", ""),
                             ReversePrimer = c("", ""),
                             ReferenceSequence = c("REF","REF")
            )
            rownames(df) <- Sample.Names

            # Build reads as in previous test design
            reads <- c(
              paste0("GGC", "AAA", "CGC", "TTTT"),
              paste0("GGC", "AAA", "CGC", "TTTT"),
              paste0("GGC", "AAA", "CGC", "TTTT"),
              paste0("CAT", "AAA", "CCC", "GGGG"),
              paste0("CAT", "AAA", "CCC", "GGGG"),
              paste0("GGG", "AAA", "CCC", "CCCC"),
              paste0("GGG", "AAA", "CCC", "CCCC"),
              paste0("GGG", "AAA", "GGG", "CCCC"),
              paste0("GGG", "AAA", "GGG", "AAAA")
            )
            srq <- make_srq(reads)
            Reads.Filtered.List <- list(S1 = srq, S2 = srq)

            cfg <- list(nCores = 1, protein.mutations = FALSE, read.frequency.limit = 30)

            out <- fcopy(Sample.Names = Sample.Names,
                         Reads.Filtered.List = Reads.Filtered.List,
                         Input.DataFrame = df,
                         Config.List = cfg)

            # Basic checks
            expect_type(out, "list")
            expect_named(out, c("S1", "S2", "All"))
            expect_named(out$All, c("S1", "S2"))

            seq_all <- out$All$S1
            expect_true(all(c("UMI", "TargetSequence", "Reads") %in% colnames(seq_all)))
            expect_true(isEmpty(out$All$S2))

            # AAA should be present with Reads == 3 and TargetSequence == "TTTT"
            umi_row <- seq_all[seq_all$UMI == "GGC", , drop = FALSE]
            expect_true(nrow(umi_row) == 1)
            expect_equal(as.integer(umi_row$Reads), 3)
            expect_equal(as.character(umi_row$TargetSequence), "CGC")

            # CAT UMI filtered out (<3)
            expect_false(any(seq_all$UMI == "CAT"))

            # No "Rejected" present
            expect_false(any(seq_all$TargetSequence == "Rejected"))

            # Aggregated Sequence.Table should contain counts for "CGC"
            seq_table <- out$S1
            expect_true(all(c("TargetSequence", "UMI Count", "Percent") %in% colnames(seq_table)))
            expect_true(any(seq_table$TargetSequence == "CGC"))
            cgc_row <- seq_table[seq_table$TargetSequence == "CGC", , drop = FALSE]
            expect_true(as.integer(cgc_row$`UMI Count`) >= 1)

            #check for S2 (no UMIs)
            seq_table <- out$S2
            expect_true(all(c("TargetSequence", "Reads", "Percent") %in% colnames(seq_table)))

            # Only "CCC" and "CGC" should remain with threshold 30 (%)
            expect_true(nrow(seq_table) == 2)
            expect_equal(as.character(seq_table$TargetSequence[1]), "CCC")
            expect_equal(as.character(seq_table$TargetSequence[2]), "CGC")
            expect_equal(sum(as.numeric(seq_table$Reads)), 7)
          })

