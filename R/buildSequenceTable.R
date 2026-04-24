###############################################################################
# buildSequenceTable()
# --------------------

#' Make a labeled table of sequences
#'
#' @description
#' Uses the settings specified by a tas.settings object, including a file path to a merged .fastq file, to generate a labeled table of sequences, with counts, frequencies, DNA mutations, and AA mutations annotated.
#'
#' Will use barcodes for sequence filtering/demultiplexing if provided. Will bin sequences and count by UMIs if provided (only supports 1 UMI currently). Otherwise, will bin sequences and filter based on the global setting tasGlobalSettings$read.frequency.limit (the minimum % to accept a sequence).
#'
#' @param settings An object of S4 class tas.settings
#'
#' @returns An S4 object of class tas.sequences
#' @export
#'
#' @examples
#' buildSequenceTable(settings)
buildSequenceTable <- function(settings) {
  filter.counts <- numeric()
  if (!file.exists(settings@MergedFASTQPath)) {stop("FASTQ file not found.")}
  cat("Reading .fastq file...\n")
  reads <- ShortRead::readFastq(settings@MergedFASTQPath)
  filter.counts <- c(filter.counts, Merged = length(reads))
  cat("Filtering reads...\n")
  reads.filtered <- filterSequences(reads, settings)
  if (length(reads.filtered) == 0) {stop("No reads after filtering. Check that your primer and extension sequences match your sample and that primers were not trimmed during merging.")}
  filter.counts <- c(filter.counts, Filtered = length(reads.filtered))
  reads <- NULL
  #cat("Building sequence table...\n")
  sequenceTable(reads.filtered, settings, filter.counts)
}





###############################################################################
# Helper functions
# ------

### filter reads by primers + extensions ###

filterSequences <- function(reads, settings) {
  fwd <- DNAString(str_c(settings@ForwardExtension, settings@ForwardPrimer))
  rev <- DNAString(str_c(settings@ReverseExtension, settings@ReversePrimer))
  temp <- ShortRead::sread(reads)
  temp_rc <- reverseComplement(temp)
  fwd_IR <- vmatchPattern(fwd,temp)
  fwd_idx <- which((elementNROWS(fwd_IR)==1))
  dff <- as.data.frame(fwd_IR)
  fwd_idx <- fwd_idx[dff$start[dff$group %in% fwd_idx]==1]
  rev_IR <- vmatchPattern(rev,temp_rc,fixed = FALSE)
  rev_idx <- which((elementNROWS(rev_IR)==1))
  dfr <- as.data.frame(rev_IR)
  rev_idx <- rev_idx[dfr$start[dfr$group %in% rev_idx]==1]
  both_idx <- intersect(fwd_idx,rev_idx)
  return(reads[both_idx])
}

### make sequence table ###

sequenceTable <- function(reads.filtered, settings, filter.counts = NA) {


  if (tolower(settings@ForwardExtensionType) == "umi" && tolower(settings@ReverseExtensionType) == "umi") {
    stop("tasaR does not currently support dual UMIs on both ends of the amplicon.")
  } else if (tolower(settings@ForwardExtensionType) == "umi") {
    umi.pos <- data.frame(start = 1, end = Biostrings::nchar(settings@ForwardExtension))
  } else if (tolower(settings@ReverseExtensionType) == "umi") {
    umi.pos <- data.frame(start = -Biostrings::nchar(settings@ReverseExtension))
  } else {
    umi.pos <- NA
  }

  if (any(!is.na(umi.pos))) {
    if (is.null(umi.pos$end)) {
      umi_temp <- ShortRead::sread(narrow(reads.filtered, start = umi.pos$start))
    } else {
      umi_temp <- ShortRead::sread(narrow(reads.filtered, start = umi.pos$start, end = umi.pos$end))
    }
    cat("Binning UMIs...\n")
    # extract id, seq, umi from Reads.Filtered.List and remove any sequences containing N called nt's
    id_temp <- as.character(ShortRead::id(reads.filtered))
    seq_temp <- ShortRead::sread(narrow(reads.filtered, start = settings@InsertStart, end = settings@InsertEnd))
    n_idx <- which(elementNROWS(Biostrings::vmatchPattern("n", seq_temp)) == 0)
    umi_temp <- as.character(umi_temp[n_idx])
    seq_temp <- as.character(seq_temp[n_idx])
    id_temp <- id_temp[n_idx]

    #use data.table to bin by UMI, remove UMIs with <3 reads
    dt_temp <- data.table::data.table(umi=umi_temp, seq=seq_temp, id=id_temp)
    dt_id_umi <- dt_temp[, .(id = list(id), count=length(id)), by="umi"]
    dt_id_umi <- dt_id_umi[dt_id_umi$count>=3,]

    #group sequences
    dt_id_seq <- dt_temp[, .(id = list(id),count=length(id)), by="seq"]

    #ungroup umis and sequences while retaining grouping information in umi or seq.group columns
    dt_id_umi_match <- data.table::data.table(umi = rep(dt_id_umi$umi,dt_id_umi$count), umi.id = unlist(dt_id_umi$id))
    dt_id_seq_match <- data.table::data.table(seq.group = rep(1:nrow(dt_id_seq),dt_id_seq$count), seq.id = unlist(dt_id_seq$id), seq = rep(dt_id_seq$seq,dt_id_seq$count))

    #merge umi and seq data tables based on FASTQ IDs
    dt_id_merge <- data.table::merge.data.table(dt_id_umi_match, dt_id_seq_match, by.x="umi.id", by.y="seq.id")

    #re-group merged data.tables by umi
    dt_merge_group <- dt_id_merge[, .(seq = umi_pileup(seq), count = length(seq.group), id = list(umi.id)), by="umi"]
    setorder(dt_merge_group,-count)
    dt_merge_group <- dt_merge_group[seq != "Rejected"]
    dt <- dt_merge_group[, .(N = .N, umis = list(umi), ids = list(id)), by = seq]
    setorder(dt, -N)
    umi.count <- length(unlist(dt$umis))
    unique.count <- nrow(dt)

    ##########

  } else { #make dt if no UMIs
    cat("Binning sequences...\n")
    id_temp <- as.character(ShortRead::id(reads.filtered))
    seq_temp <- ShortRead::sread(narrow(reads.filtered ,start = settings@InsertStart, end = settings@InsertEnd))
    n_idx <- which(elementNROWS(Biostrings::vmatchPattern("n", seq_temp)) == 0)
    seq_temp <- as.character(seq_temp[n_idx])
    id_temp <- id_temp[n_idx]
    dt_all <- data.table(seq = seq_temp, id = id_temp)
    dt <- dt_all[, .(N = .N, ids = list(id)), by = seq]
    min_freq <- sum(dt$N)*tasGlobalSettings$read.frequency.limit/100
    dt <- dt[N > min_freq]
    setorder(dt, -N)
    dt <- cbind(dt, data.table(umis = list(NA)))
  }

  ##########

  cat("Labeling DNA mutations...\n")
  Reference.Sequence.DNA <- Biostrings::DNAString(settings@ReferenceSequence)
  Reads.Unique.DNA <- Biostrings::DNAStringSet(dt$seq)
  Pairwise.Aligned.DNA <- pwalign::pairwiseAlignment(Reads.Unique.DNA, Reference.Sequence.DNA)

  # write sequences to table with deletions marked by -
  dt$seq <- as.character(pattern(Pairwise.Aligned.DNA))

  # table of indel counts and sizes
  dt.indel <- data.table(iNum = insertion(nindel(Pairwise.Aligned.DNA))[,"Length"],
                         iSize = str_c("+", insertion(nindel(Pairwise.Aligned.DNA))[,"WidthSum"]),
                         dNum = deletion(nindel(Pairwise.Aligned.DNA))[,"Length"],
                         dSize = str_c("-", deletion(nindel(Pairwise.Aligned.DNA))[,"WidthSum"]))

  # find sequences with terminal deletions missed by indel()
  len.check <- data.table(pWidth = width(pattern(Pairwise.Aligned.DNA)) - pwalign::nchar(Reference.Sequence.DNA),
                          idSize = insertion(nindel(Pairwise.Aligned.DNA))[,"WidthSum"] - deletion(nindel(Pairwise.Aligned.DNA))[,"WidthSum"])
  idx.end.del <- which(len.check$pWidth != len.check$idSize)
  end.del.size <- str_c("-", str_count(as.character(alignedPattern(Pairwise.Aligned.DNA[idx.end.del])), "-") - as.numeric(dt.indel$dSize[idx.end.del]))
  end.del.start <- str_locate(as.character(alignedPattern(Pairwise.Aligned.DNA[idx.end.del])), "-")[,"start"]
  dt.indel$dNum[idx.end.del] <- as.numeric(dt.indel$dNum[idx.end.del]) + 1

  # concatenate multiple insertions
  ins.temp <- insertion(Pairwise.Aligned.DNA)
  ins.midx <- which(dt.indel$iNum>1)
  dt.indel$iSize[ins.midx] <- lapply(ins.midx, function(x){
    str_flatten(str_c("+", ins.temp[[x]]@width), collapse = ", ")
  })

  # concatenate multilple deletions, including from ends
  del.temp <- deletion(Pairwise.Aligned.DNA)
  del.midx <- c(which(dt.indel$dNum>1), idx.end.del)
  dt.indel$dSize[del.midx] <- lapply(del.midx, function(x){
    s <- str_flatten(str_c("-", del.temp[[x]]@width), collapse = ", ")
    if (x %in% idx.end.del) {
      s <- str_c(s, end.del.size[x == idx.end.del], collapse = ", ")
    }
    return(s)
  })

  # empty strings w/o indels
  dt.indel$iSize[dt.indel$iNum == 0] <- ""
  dt.indel$dSize[dt.indel$dNum == 0] <- ""

  # concatenate insertions & deletions together into a single output string
  str.indel <- unlist(sapply(1:(nrow(dt.indel)), function(x){
    if (dt.indel$iNum[x] == 0 && dt.indel$dNum[x] == 0){
      return("")
    } else if (dt.indel$iNum[x] == 0 && dt.indel$dNum[x] != 0){
      return(dt.indel$dSize[x])
    } else if (dt.indel$iNum[x] != 0 && dt.indel$dNum[x] == 0) {
      return(dt.indel$iSize[x])
    } else {
      str_c(dt.indel$iSize[x], dt.indel$dSize[x], sep = ", ")
    }
  }))

  # table of mismatches
  dt.mm <- data.table(pwalign::mismatchTable(Pairwise.Aligned.DNA))[, .N, by = PatternId]
  vec.mm <- unlist(sapply(1:length(Pairwise.Aligned.DNA), function(x){
    if (x %in% dt.mm$PatternId){
      return(dt.mm$N[dt.mm$PatternId == x])
    } else {
      return(0)
    }
  }))

  idx.wt <- which(str.indel == "" & vec.mm == 0)
  if (length(idx.wt) > 1){
    stop("More than 1 WT sequence detected.")
  }
  str.indel[idx.wt] <- "WT"
  dt <- cbind(dt, data.table(Indels = str.indel, BasesChanged = vec.mm))




  ############################################################################################
  # Notes about annotating protein mutations:
  #
  # Insertions and deletions (indels) pose a significant challenge for quantifying the frequency
  # of mutations in the protein sequence at each position. Because they can cause frameshifts,
  # their effects can be propogated through the rest of the downstream sequence. To avoid this
  # overestimation of downstream mutations, protein sequences are truncated at the site of the
  # first indel within the Reference Sequence, and marked with a "-" for deletion and "+" for
  # insertion.
  #
  # This allows pairwise alignments to ignore confounding downstream mutations and focus on
  # where the mutation actually occurred. This also ensures that detection of protein mismatches
  # depends on specific mutations rather than frameshifts.
  ############################################################################################

  cat("Labeling protien mutations...\n")

  Reference.Sequence.Protein <- AAStringSet(suppressWarnings(translate(Reference.Sequence.DNA)))

  indel.start.dt <- rbindlist(list(cbind(as.data.table(deletion(Pairwise.Aligned.DNA))[,c("group", "start")], data.table(type = "deletion")),
                                   cbind(as.data.table(insertion(Pairwise.Aligned.DNA))[,c("group", "start")], data.table(type = "insertion")),
                                   cbind(data.table(group = idx.end.del, start = end.del.start), data.table(type = "deletion")),
                                   data.table(group = 1:length(Pairwise.Aligned.DNA))
  ), fill = TRUE)[, .(minStart = if (all(is.na(start))) {NA_integer_} else {min(start, na.rm = TRUE)}, type = type[which.min(start)]), by = group]
  indel.start.dt <- indel.start.dt[!is.na(indel.start.dt$group),]
  setorder(indel.start.dt, group)

  prot.reads <- as.character(suppressWarnings(translate(Reads.Unique.DNA)))

  # truncate at first indel
  prot.trunc <- sapply(1:nrow(indel.start.dt), function(x) {
    if (is.na(indel.start.dt$minStart[x])) {
      return(prot.reads[x])
    } else {
      trunc <- ceiling(indel.start.dt$minStart[x]/3) - 1
      st <- str_trunc(prot.reads[x], trunc, ellipsis = "")
      if (indel.start.dt$type[x] == "deletion") {
        return(str_c(st, "-"))
      } else if (indel.start.dt$type[x] == "insertion") {
        return(str_c(st, "+"))
      }
    }
  })

  # detect any non-indel nonsense mutations
  idx_ns <- which(str_detect(prot.trunc, "[*]"))
  ns.pos <- na.omit(str_locate(prot.trunc, "[*]")[,"start"])
  prot.trunc[idx_ns] <- unlist(sapply(seq_along(idx_ns), function(x) {
    str_trunc(prot.trunc[idx_ns[x]], (ns.pos[x]), ellipsis = "")
  }))

  prot.align <- pairwiseAlignment(AAStringSet(prot.trunc), Reference.Sequence.Protein)

  #Label each sequence with protein mutations
  aa.mut <- data.frame(rep("",length(prot.align)))
  colnames(aa.mut) <- "ProteinMutation"

  mmT <- mismatchTable(prot.align)
  idx_indel <- !Biostrings::nchar(prot.trunc)==Biostrings::nchar(Reference.Sequence.Protein)
  idx_WT <- !seq_along(prot.align) %in% mmT$PatternId # indel multiples of 3 identified incorrectly
  idx_WT[which((idx_indel+idx_WT)==2)] <- FALSE # remove false +ve from idx_WT
  idx_mm <- which((idx_indel+idx_WT)==0)
  mmT <- mmT[which(mmT$PatternId %in% idx_mm),]
  mmChar <- cbind(data.frame(PatternId=mmT$PatternId),data.frame(Mutation=str_c(mmT[,"SubjectSubstring"],mmT[,"SubjectStart"],mmT[,"PatternSubstring"])))

  mutCollapse <- function(x){
    str_flatten(x,", ")
  }

  # error protection, aggregate throws error if mmChar is empty
  if (!isEmpty(mmChar)){
    aa.mut$ProteinMutation[idx_mm] <- aggregate(Mutation ~ PatternId, data = mmChar, FUN = mutCollapse)[,"Mutation"]
  }
  aa.mut$ProteinMutation[idx_indel] <- "Indel"
  aa.mut$ProteinMutation[idx_WT] <- "WT"
  aa.mut$ProteinMutation[str_detect(prot.trunc,"\\*")] <- "Nonsense"

  dt <- cbind(dt, data.table(AA = prot.trunc, aa.mut))

  if (exists("filter.counts")) {
    if (all(is.na(filter.counts))) {
      ReadCounts <- NA_integer_
    } else {
      ReadCounts <-  setNames(as.integer(c(filter.counts, umi.count, unique.count)), c(names(filter.counts), "UMIs", "UniqueSequences"))
    }
  } else {ReadCounts <- NA_integer_}

  ##########

  new("tas.sequences", Table = data.frame(Sequences = dt$seq,
                                          Count = dt$N,
                                          Percent = (dt$N / sum(dt$N) * 100),
                                          Indels = dt$Indels,
                                          BasesChanged = dt$BasesChanged,
                                          AA = dt$AA,
                                          ProteinMutation = dt$ProteinMutation),
      Supplemental = data.table(data.table(Index = 1:nrow(dt),
                                           UMIs = dt$umis,
                                           IDs = dt$ids,
                                           IndelStart = indel.start.dt$minStart,
                                           IndelType = indel.start.dt$type)),
      Alignments = list(DNA = Pairwise.Aligned.DNA,
                        AA = prot.align),
      ReadCounts = ReadCounts)
}






# possible framework for fastqstreamer implementation
# slower for amplicon-ez sized fastq files and seems not needed

# infile <- "/Users/geoff/Google Drive/Amplicon EZ/test5/analyzed 2/merged/E1-J3-T3-merged.fastq.gz"
# outfile <- file.path(tempdir(), "filtered.fastq.gz")
#
# strm <- FastqStreamer(infile, n = 1e6)
#
# first <- TRUE
#
# repeat {
#   fq <- yield(strm)
#   if (length(fq) == 0) break
#
#   fq.filt <- filterSequences(fq, settings)
#
#   if (length(fq.filt) > 0) {
#     writeFastq(
#       fq.filt,
#       outfile,
#       mode = if (first) "w" else "a",
#       compress = TRUE
#     )
#     first <- FALSE
#   }
# }
#
# close(strm)
