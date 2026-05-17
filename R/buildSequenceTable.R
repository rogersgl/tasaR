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
buildSequenceTable <- function(settings,
                               # min.read.frequency = 0.1,
                               # paired.analysis.ctrl = FALSE,
                               # diploid = TRUE,
                               ...
                               # with.nuclease = FALSE,
                               # gRNA.seq = NULL,
                               # manual.cut.site = NULL
                               ) {
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
  sequenceTable(reads.filtered,
                settings,
                filter.counts,
                # min.read.frequency = 0.1,
                ...
                # paired.analysis.ctrl,
                # diploid,
                # with.nuclease = FALSE,
                # gRNA.seq = NULL,
                # manual.cut.site = NULL
                )
}





###############################################################################
# Helper functions
# ------

### filter reads by primers + extensions ###

filterSequences <- function(reads, settings) {
  fwd <- Biostrings::DNAString(stringr::str_c(settings@ForwardExtension, settings@ForwardPrimer))
  rev <- Biostrings::DNAString(stringr::str_c(settings@ReverseExtension, settings@ReversePrimer))
  temp <- ShortRead::sread(reads)
  temp_rc <- Biostrings::reverseComplement(temp)
  fwd_IR <- Biostrings::vmatchPattern(fwd, temp)
  fwd_idx <- which((S4Vectors::elementNROWS(fwd_IR)==1))
  dff <- as.data.frame(fwd_IR)
  fwd_idx <- fwd_idx[dff$start[dff$group %in% fwd_idx]==1]
  rev_IR <- Biostrings::vmatchPattern(rev,temp_rc,fixed = FALSE)
  rev_idx <- which((S4Vectors::elementNROWS(rev_IR)==1))
  dfr <- as.data.frame(rev_IR)
  rev_idx <- rev_idx[dfr$start[dfr$group %in% rev_idx]==1]
  both_idx <- S4Vectors::intersect(fwd_idx,rev_idx)
  return(reads[both_idx])
}


### make sequence table ###

sequenceTable <- function(reads.filtered,
                          settings,
                          filter.counts = NA,
                          min.read.frequency = 0.1,
                          paired.analysis.ctrl = FALSE,
                          diploid = TRUE,
                          with.nuclease = FALSE,
                          gRNA.seq = NULL,
                          manual.cut.site = NULL,
                          ...) {


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
      umi_temp <- ShortRead::sread(ShortRead::narrow(reads.filtered, start = umi.pos$start))
    } else {
      umi_temp <- ShortRead::sread(ShortRead::narrow(reads.filtered, start = umi.pos$start, end = umi.pos$end))
    }
    cat("Binning UMIs...\n")
    # extract id, seq, umi from Reads.Filtered.List and remove any sequences containing N called nt's
    id_temp <- as.character(ShortRead::id(reads.filtered))
    seq_temp <- ShortRead::sread(ShortRead::narrow(reads.filtered, start = settings@InsertStart, end = settings@InsertEnd))
    n_idx <- which(S4Vectors::elementNROWS(Biostrings::vmatchPattern("n", seq_temp)) == 0)
    umi_temp <- as.character(umi_temp[n_idx])
    seq_temp <- as.character(seq_temp[n_idx])
    id_temp <- id_temp[n_idx]

    #use data.table to bin by UMI, remove UMIs with <3 reads
    dt_temp <- data.table(umi=umi_temp, seq=seq_temp, id=id_temp)
    dt_id_umi <- dt_temp[, .(id = list(id), count=length(id)), by="umi"]
    dt_id_umi <- dt_id_umi[dt_id_umi$count>=3,]

    #group sequences
    dt_id_seq <- dt_temp[, .(id = list(id),count=length(id)), by="seq"]

    #ungroup umis and sequences while retaining grouping information in umi or seq.group columns
    dt_id_umi_match <- data.table(umi = rep(dt_id_umi$umi,dt_id_umi$count), umi.id = BiocGenerics::unlist(dt_id_umi$id))
    dt_id_seq_match <- data.table(seq.group = rep(1:nrow(dt_id_seq), dt_id_seq$count), seq.id = BiocGenerics::unlist(dt_id_seq$id), seq = rep(dt_id_seq$seq, dt_id_seq$count))

    #merge umi and seq data tables based on FASTQ IDs
    dt_id_merge <- merge.data.table(dt_id_umi_match, dt_id_seq_match, by.x="umi.id", by.y="seq.id")

    #re-group merged data.tables by umi
    dt_merge_group <- dt_id_merge[, .(seq = .umi_pileup(seq), count = length(seq.group), id = list(umi.id)), by="umi"]
    setorder(dt_merge_group,-count)
    dt_merge_group <- dt_merge_group[seq != "Rejected"]
    dt <- dt_merge_group[, .(N = .N, umis = list(umi), ids = list(id)), by = seq]
    setorder(dt, -N)
    umi.count <- length(BiocGenerics::unlist(dt$umis))
    unique.count <- nrow(dt)

    ##########

  } else { #make dt if no UMIs
    cat("Binning sequences...\n")
    id_temp <- as.character(ShortRead::id(reads.filtered))
    seq_temp <- ShortRead::sread(ShortRead::narrow(reads.filtered ,start = settings@InsertStart, end = settings@InsertEnd))
    n_idx <- which(S4Vectors::elementNROWS(Biostrings::vmatchPattern("n", seq_temp)) == 0)
    seq_temp <- as.character(seq_temp[n_idx])
    id_temp <- id_temp[n_idx]
    dt_all <- data.table(seq = seq_temp, id = id_temp)
    dt <- dt_all[, .(N = .N, ids = list(id)), by = seq]
    min_freq <- sum(dt$N)*min.read.frequency/100
    dt <- dt[N > min_freq]
    setorder(dt, -N)
    dt <- cbind(dt, data.table(umis = list(NA)))
    umi.count <- NA
    unique.count <- nrow(dt)
  }

  ##########


  if (paired.analysis.ctrl) {
    t <- dt[(dt$N > 0.3*dt$N[1]),] # get reads with at least 30% variant allele frequency
    if (diploid && nrow(t)>2) {
      warning("Detected more than 2 potential reference sequences (allele frequency > 30%).")
    }
    if (length(unique(Biostrings::nchar(t$seq))) > 1) {
      stop("Potential wild-type sequences are of different lengths. Analysis is not currently supported by tasaR.")
    }
    Reference.Sequence.DNA <- as.list(Biostrings::DNAStringSet(t$seq))
  } else {
    Reference.Sequence.DNA <- lapply(settings@ReferenceSequence, Biostrings::DNAString)
  }

  cat("Labeling mutations...\n")

  if (length(Reference.Sequence.DNA) > 1) {
    seqs <- Biostrings::DNAStringSet(dt$seq)
    aln <- lapply(Reference.Sequence.DNA, function(x) {
      pwalign::pairwiseAlignment(seqs, x)
      })
    aln.idx <- apply(as.data.table(lapply(aln, BiocGenerics::score)), 1, which.max)
    Reads.Unique.DNA <- lapply(seq_along(Reference.Sequence.DNA), function(x) {
      Biostrings::DNAStringSet(dt$seq[aln.idx==x])
    })
    Pairwise.Aligned.DNA <- sapply(seq_along(Reference.Sequence.DNA), function(x) {
      aln[[x]][aln.idx==x]
    })
  } else {
    aln.idx <- rep(1, nrow(dt))
    Reads.Unique.DNA <- list(Biostrings::DNAStringSet(dt$seq))
    Pairwise.Aligned.DNA <- lapply(Reference.Sequence.DNA, function(x) {
      pwalign::pairwiseAlignment(Reads.Unique.DNA[[1]], x)
    })
  }

  Reference.Sequence.Protein <- suppressWarnings(lapply(Reference.Sequence.DNA, Biostrings::translate))
  Reads.Unique.Protein <- suppressWarnings(lapply(Reads.Unique.DNA, Biostrings::translate))

  # TODO: thinking about putting the .nucleaseAlignmentCorrection here for pairedAnalyzeAmplicon
  # add type = c("single", "paired") to arguments and flag upstream
  #
  # Pros: would fix issues with misaligned sequences around the cut site, which lead to incorrect results
  # in the measureMutations by position function.
  #
  # Cons: would break the entire workflow since the object becomes a DNAMultipleAlignment instead of a
  # PairwiseAlignmentsSingleSubject object. So everything downstream in measureMutations that uses
  # the PairwiseAlignment would need to be adjusted.
  #
  # Seems like roughtly a wash - could insert it here and fix everything else downstream with
  # multiple pathways, or could do it afterwards and redo/replace the analyses. That's wasted CPU cycles,
  # but it might be easier.
  #
  # Alternatively, instead of diverging paths, could make an entirely new measureMutations function for the
  # nuclease Paired alignments.

  # if (type == "paired") {
  #   Pairwise.Aligned.DNA <- lapply(Pairwise.Aligned.DNA, function(x) {
  #
  #   })
  # }


  #############

  label.list <- lapply(seq_along(Pairwise.Aligned.DNA), function(x) {

    dna.align <- Pairwise.Aligned.DNA[[x]]
    dt <- dt[aln.idx == x,]

    # write sequences to table with deletions marked by -
    dt$seq <- as.character(pwalign::pattern(dna.align))

    # table of indel counts and sizes
    dt.indel <- data.table(iNum = pwalign::insertion(pwalign::nindel(dna.align))[,"Length"],
                           iSize = stringr::str_c("+", pwalign::insertion(pwalign::nindel(dna.align))[,"WidthSum"]),
                           dNum = pwalign::deletion(pwalign::nindel(dna.align))[,"Length"],
                           dSize = stringr::str_c("-", pwalign::deletion(pwalign::nindel(dna.align))[,"WidthSum"]))

    # find sequences with terminal deletions missed by indel()
    len.check <- data.table(pWidth = BiocGenerics::width(pwalign::pattern(dna.align)) - pwalign::nchar(Reference.Sequence.DNA[[x]]),
                            idSize = pwalign::insertion(pwalign::nindel(dna.align))[,"WidthSum"] - pwalign::deletion(pwalign::nindel(dna.align))[,"WidthSum"])
    idx.end.del <- which(len.check$pWidth != len.check$idSize)
    end.del.size <- stringr::str_c("-", stringr::str_count(as.character(pwalign::alignedPattern(dna.align[idx.end.del])), "-") - as.numeric(dt.indel$dSize[idx.end.del]))
    end.del.start <- stringr::str_locate(as.character(pwalign::alignedPattern(dna.align[idx.end.del])), "-")[,"start"]
    dt.indel$dNum[idx.end.del] <- as.numeric(dt.indel$dNum[idx.end.del]) + 1

    # concatenate multiple insertions
    ins.temp <- pwalign::insertion(dna.align)
    ins.midx <- which(dt.indel$iNum>1)
    dt.indel$iSize[ins.midx] <- S4Vectors::lapply(ins.midx, function(x){
      stringr::str_flatten(stringr::str_c("+", ins.temp[[x]]@width), collapse = ", ")
    })

    # concatenate multilple deletions, including from ends
    del.temp <- pwalign::deletion(dna.align)
    del.midx <- c(which(dt.indel$dNum>1), idx.end.del)
    dt.indel$dSize[del.midx] <- S4Vectors::lapply(del.midx, function(x){
      s <- stringr::str_flatten(stringr::str_c("-", del.temp[[x]]@width), collapse = ", ")
      if (x %in% idx.end.del) {
        s <- stringr::str_c(s, end.del.size[x == idx.end.del], collapse = ", ")
      }
      return(s)
    })

    # empty strings w/o indels
    dt.indel$iSize[dt.indel$iNum == 0] <- ""
    dt.indel$dSize[dt.indel$dNum == 0] <- ""

    # concatenate insertions & deletions together into a single output string
    str.indel <- BiocGenerics::unlist(sapply(1:(nrow(dt.indel)), function(x){
      if (dt.indel$iNum[x] == 0 && dt.indel$dNum[x] == 0){
        return("")
      } else if (dt.indel$iNum[x] == 0 && dt.indel$dNum[x] != 0){
        return(dt.indel$dSize[x])
      } else if (dt.indel$iNum[x] != 0 && dt.indel$dNum[x] == 0) {
        return(dt.indel$iSize[x])
      } else {
        stringr::str_c(dt.indel$iSize[x], dt.indel$dSize[x], sep = ", ")
      }
    }))

    # table of mismatches
    dt.mm <- data.table(pwalign::mismatchTable(dna.align))[, .N, by = PatternId]
    vec.mm <- BiocGenerics::unlist(sapply(1:length(dna.align), function(x){
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



  ################################################################################################
  # Notes about annotating protein mutations:                                                    #
  #                                                                                              #
  # Insertions and deletions (indels) pose a significant challenge for quantifying the frequency #
  # of mutations in the protein sequence at each position. Because they can cause frameshifts,   #
  # their effects can be propogated through the rest of the downstream sequence. To avoid this   #
  # overestimation of downstream mutations, protein sequences are truncated at the site of the   #
  # first indel within the Reference Sequence, and marked with a "-" for deletion and "+" for    #
  # insertion.                                                                                   #
  #                                                                                              #
  # This allows pairwise alignments to ignore confounding downstream mutations and focus on      #
  # where the mutation actually occurred. This also ensures that detection of protein mismatches #
  # depends on specific mutations rather than frameshifts.                                       #
  ################################################################################################

    prot.wt <- Reference.Sequence.Protein[[x]]

    indel.start.dt <- rbindlist(list(cbind(as.data.table(pwalign::deletion(dna.align))[,c("group", "start")], data.table(type = "deletion")),
                                     cbind(as.data.table(pwalign::insertion(dna.align))[,c("group", "start")], data.table(type = "insertion")),
                                     cbind(data.table(group = idx.end.del, start = end.del.start), data.table(type = "deletion")),
                                     data.table(group = 1:length(dna.align))
    ), fill = TRUE)[, .(minStart = if (all(is.na(start))) {NA_integer_} else {min(start, na.rm = TRUE)}, type = type[which.min(start)]), by = group]
    indel.start.dt <- indel.start.dt[!is.na(indel.start.dt$group),]
    setorder(indel.start.dt, group)

    prot.reads <- as.character(Reads.Unique.Protein[[x]])

    # truncate at first indel
    prot.trunc <- sapply(1:nrow(indel.start.dt), function(x) {
      if (is.na(indel.start.dt$minStart[x])) {
        return(prot.reads[x])
      } else {
        trunc <- ceiling(indel.start.dt$minStart[x]/3) - 1
        st <- stringr::str_trunc(prot.reads[x], trunc, ellipsis = "")
        if (indel.start.dt$type[x] == "deletion") {
          return(stringr::str_c(st, "-"))
        } else if (indel.start.dt$type[x] == "insertion") {
          return(stringr::str_c(st, "+"))
        }
      }
    })

    # detect any non-indel nonsense mutations
    idx_ns <- which(stringr::str_detect(prot.trunc, "[*]"))
    ns.pos <- stats::na.omit(stringr::str_locate(prot.trunc, "[*]")[,"start"])
    prot.trunc[idx_ns] <- BiocGenerics::unlist(sapply(seq_along(idx_ns), function(x) {
      stringr::str_trunc(prot.trunc[idx_ns[x]], (ns.pos[x]), ellipsis = "")
    }))

    prot.align <- pwalign::pairwiseAlignment(Biostrings::AAStringSet(prot.trunc), prot.wt)

    #Label each sequence with protein mutations
    aa.mut <- data.frame(rep("",length(prot.align)))
    colnames(aa.mut) <- "ProteinMutation"

    mmT <- pwalign::mismatchTable(prot.align)
    idx_indel <- !Biostrings::nchar(prot.trunc)==Biostrings::nchar(prot.wt)
    idx_WT <- !seq_along(prot.align) %in% mmT$PatternId # indel multiples of 3 identified incorrectly
    idx_WT[which((idx_indel+idx_WT)==2)] <- FALSE # remove false +ve from idx_WT
    idx_mm <- which((idx_indel+idx_WT)==0)
    mmT <- mmT[which(mmT$PatternId %in% idx_mm),]
    mmChar <- cbind(data.frame(PatternId=mmT$PatternId),data.frame(Mutation=stringr::str_c(mmT[,"SubjectSubstring"],mmT[,"SubjectStart"],mmT[,"PatternSubstring"])))

    # error protection, aggregate throws error if mmChar is empty
    if (!isEmpty(mmChar)){
      aa.mut$ProteinMutation[idx_mm] <- stats::aggregate(Mutation ~ PatternId, data = mmChar, FUN = stringr::str_flatten_comma)[,"Mutation"]
    }
    aa.mut$ProteinMutation[idx_indel] <- "Indel"
    aa.mut$ProteinMutation[idx_WT] <- "WT"
    aa.mut$ProteinMutation[stringr::str_detect(prot.trunc,"\\*")] <- "Nonsense"

    dt.out <- cbind(dt, data.table(AA = prot.trunc, aa.mut), indel.start.dt, data.table(RefIdx = aln.idx[aln.idx == x]))
    list(DT = dt.out, pAln = prot.align)
  })

  # re-integrate separate data tables and order
  dt.list <- purrr::map(label.list, "DT")
  pAln.list <- purrr::map(label.list, "pAln")
  names(pAln.list) <- NULL
  dt.dt <- rbindlist(dt.list)
  setorder(dt.dt, -N)

  if (with.nuclease) {
    msa.dna <- lapply(seq_along(Pairwise.Aligned.DNA), function(x) {
      ref <- unique(as.character(pwalign::unaligned(pwalign::subject(Pairwise.Aligned.DNA[[x]]))))
      seqs <- as.character(pwalign::pattern(Pairwise.Aligned.DNA[[x]]))
      seqs <- .alignment.seqs(seqs,
                              ref,
                              seq.indels = dt.list[[x]]$Indels,
                              seq.basesChanged = dt.list[[x]]$BasesChanged,
                              seq.percent = as.numeric(dt.list[[x]]$N/sum(dt.dt$N)*100),
                              gRNA.window = 1:nchar(ref),
                              seq.count = 10)
      msa.aln <- msa::msaClustalW(seqs, order = "input")
      if (!is.null(gRNA.seq)){
        cut.site <- .findCutSitesFromSequence(gRNA.seq, ref)
      } else if (!is.null(manual.cut.site)){
        cut.site <- manual.cut.site
      } else {stop("Could not find nuclease cut site. Too confused to continue.")}
      .nucleaseAlignmentCorrection(msa.aln, cut.site, reference = 1L)
    })
  } else {
    msa.dna <- list(Biostrings::DNAMultipleAlignment())
  }


  if (exists("filter.counts")) {
    if (all(is.na(filter.counts))) {
      ReadCounts <- NA_integer_
    } else {
      ReadCounts <-  stats::setNames(as.integer(c(filter.counts, umi.count, unique.count)), c(names(filter.counts), "UMIs", "UniqueSequences"))
    }
  } else {ReadCounts <- NA_integer_}

  ##########

  new("tas.sequences", Table = data.frame(Sequences = dt.dt$seq,
                                          Count = dt.dt$N,
                                          Percent = (dt.dt$N / sum(dt.dt$N) * 100),
                                          Indels = dt.dt$Indels,
                                          BasesChanged = dt.dt$BasesChanged,
                                          AA = dt.dt$AA,
                                          ProteinMutation = dt.dt$ProteinMutation),
      Supplemental = data.table(data.table(Index = 1:nrow(dt.dt),
                                           UMIs = dt.dt$umis,
                                           IDs = dt.dt$ids,
                                           IndelStart = dt.dt$minStart,
                                           IndelType = dt.dt$type,
                                           RefIdx = dt.dt$RefIdx)),
      Alignments = list(DNA = Pairwise.Aligned.DNA, # now a list of pairwiseAlignments
                        AA = pAln.list, # now a list of pairwiseAlignments
                        msaDNA = msa.dna),
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
