
tas_filter2 <- function(reads, settings) {
  fwd <- DNAString(str_c(settings@ForwardExtension, settings@ForwardPrimer))
  rev <- DNAString(str_c(settings@ReverseExtension, settings@ReversePrimer))
  temp <- reads@sread
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


tas_sequence_table <- function(reads.filtered, settings) {
  if (tolower(settings@ForwardExtensionType) == "umi" && tolower(settings@ReverseExtensionType) == "umi") {
    stop("tasAnalyzer does not currently support dual UMIs on both ends of the amplicon.")
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
  }

  # extract id, seq, umi from Reads.Filtered.List and remove any sequences containing N called nt's
  id_temp <- as.character(ShortRead::id(reads.filtered))
  seq_temp <- ShortRead::sread(narrow(reads.filtered ,start = settings@InsertStart, end = settings@InsertEnd))
  n_idx <- which(elementNROWS(Biostrings::vmatchPattern("n", seq_temp)) == 0)
  umi_temp <- as.character(umi_temp[n_idx])
  seq_temp <- as.character(seq_temp[n_idx])
  id_temp <- id_temp[n_idx]

  #use data.table to bin by UMI, remove UMIs with <3 reads
  dt_temp <- data.table::data.table(umi=umi_temp, seq=seq_temp, id=id_temp)
  dt_id_umi <- dt_temp[, .(id = list(id),count=length(id)), by="umi"]
  dt_id_umi <- dt_id_umi[dt_id_umi$count>=3,]

  #group sequences
  dt_id_seq <- dt_temp[, .(id = list(id),count=length(id)), by="seq"]

  #ungroup umis and sequences while retaining grouping information in umi or seq.group columns
  dt_id_umi_match <- data.table::data.table(umi = rep(dt_id_umi$umi,dt_id_umi$count), umi.id = unlist(dt_id_umi$id))
  dt_id_seq_match <- data.table::data.table(seq.group = rep(1:nrow(dt_id_seq),dt_id_seq$count), seq.id = unlist(dt_id_seq$id), seq = rep(dt_id_seq$seq,dt_id_seq$count))

  #merge umi and seq data tables based on FASTQ IDs
  dt_id_merge <- data.table::merge.data.table(dt_id_umi_match, dt_id_seq_match, by.x="umi.id", by.y="seq.id")



  #re-group merged data.tables by umi
  dt_merge_group <- dt_id_merge[, .(seq = umi_pileup(seq), count = length(seq.group)), by="umi"]
  setorder(dt_merge_group,-count)
  dt_merge_group <- dt_merge_group[seq != "Rejected"]
  dt <- dt_merge_group[, .N, by=seq]

}

Sequence.Table.All <- list()
Sequence.Table.All[Sample.Names] <- mclapply(Sample.Names,function(x){

  # print(str_c("all , ",x))
  if (toupper(Input.DataFrame[x,"ForwardExtensionType"])=="UMI"){
    umi.pos <- data.frame(start=1,end=Biostrings::nchar(Input.DataFrame[x,"ForwardExtension"]))
  }else if (toupper(Input.DataFrame[x,"ReverseExtensionType"])=="UMI"){
    umi.pos <- data.frame(start=-Biostrings::nchar(Input.DataFrame[x,"ReverseExtension"]))
  }else{
    umi.pos <- NA
  }
  if (any(!is.na(umi.pos))){
    if (is.null(umi.pos$end)){
      umi_temp <- narrow(Reads.Filtered.List[[x]], start = umi.pos$start)@sread
    }else{
      umi_temp <- narrow(Reads.Filtered.List[[x]], start = umi.pos$start, end = umi.pos$end)@sread
    }
    #extract id, seq, umi from Reads.Filtered.List and remove any sequences containing N called nt's
    id_temp <- as.character(Reads.Filtered.List[[x]]@id)
    seq_temp <- narrow(Reads.Filtered.List[[x]],start = Input.DataFrame[x,"InsertStart"],end = Input.DataFrame[x,"InsertEnd"])@sread
    n_idx <- which(elementNROWS(Biostrings::vmatchPattern("n",seq_temp))==0)
    umi_temp <- as.character(umi_temp[n_idx])
    seq_temp <- as.character(seq_temp[n_idx])
    id_temp <- id_temp[n_idx]

    #use data.table to bin by UMI, remove UMIs with <3 reads
    dt_temp <- data.table::data.table(umi=umi_temp, seq=seq_temp, id=id_temp)
    dt_id_umi <- dt_temp[, .(id = list(id),count=length(id)), by="umi"]
    dt_id_umi <- dt_id_umi[dt_id_umi$count>=3,]

    #group sequences
    dt_id_seq <- dt_temp[, .(id = list(id),count=length(id)), by="seq"]

    #ungroup umis and sequences while retaining grouping information in umi or seq.group columns
    dt_id_umi_match <- data.table::data.table(umi = rep(dt_id_umi$umi,dt_id_umi$count), umi.id = unlist(dt_id_umi$id))
    dt_id_seq_match <- data.table::data.table(seq.group = rep(1:nrow(dt_id_seq),dt_id_seq$count), seq.id = unlist(dt_id_seq$id), seq = rep(dt_id_seq$seq,dt_id_seq$count))

    #merge umi and seq data tables based on FASTQ IDs
    dt_id_merge <- data.table::merge.data.table(dt_id_umi_match,dt_id_seq_match,by.x="umi.id",by.y="seq.id")

    #function to determine optimal sequence for each umi and filter non-passing umis as "Rejected"
    umi_pileup <- function(i){
      t <- table(i)
      d <- data.table(seq = names(t),count=as.numeric(t))
      setorder(d,-count)
      if((length(d$count)==1) || (d$count[1]>=(d$count[2]+2) && (d$count[2]<d$count[1]*0.7))){
        return(d$seq[1])
      }else{
        return("Rejected")
      }
    }

    #re-group merged data.tables by umi
    dt_merge_group <- dt_id_merge[, .(seq = umi_pileup(seq), count = length(seq.group)), by="umi"]
    setorder(dt_merge_group,-count)
    dt_merge_group <- dt_merge_group[seq != "Rejected"]
    dt <- dt_merge_group[, .N, by=seq]

    Reference.Sequence.DNA <- DNAString(Input.DataFrame[x, "ReferenceSequence"])
    Reads.Unique.DNA <- DNAStringSet(dt$seq)
    Pairwise.Aligned.DNA <- pairwiseAlignment(Reads.Unique.DNA, Reference.Sequence.DNA)

    ins <-  indel(Pairwise.Aligned.DNA)@insertion
    idx.ins <- which(!sapply(ins,S4Vectors::isEmpty))
    ins.num <- rep("",length(ins))
    ins.num[idx.ins] <- lapply(ins[idx.ins],function(y){
      ins.widths <- y@width
      widths.split <- unlist(str_split(ins.widths,"-"))
      output <- str_c("+",widths.split)
      return(output)
    })

    del <- indel(Pairwise.Aligned.DNA)@deletion
    idx.del <- which(!sapply(del,S4Vectors::isEmpty))
    del.num <- rep("",length(del))
    del.num[idx.del] <- lapply(del[idx.del],function(y){
      del.widths <- y@width
      widths.split <- unlist(str_split(del.widths,"-"))
      output <- str_c("-",widths.split)
      return(output)
    })

    Indel <- lapply(1:length(del.num),function(y){
      a <- c(ins.num[[y]],del.num[[y]])
      b <- a[Biostrings::nchar(a)>0]
      if (!S4Vectors::isEmpty(b)){
        output <- str_flatten(b, collapse = ", ")
      }else{
        output <- NA
      }
      return(output)
    })
    Indel.df <- data.frame(unlist(Indel))
    colnames(Indel.df) <- "Indels"

    seq.mm <- data.frame(BiocGenerics::table(mismatchTable(Pairwise.Aligned.DNA)$PatternId))
    if (!S4Vectors::isEmpty(seq.mm)){
      colnames(seq.mm) <- c("Idx","Muts")
    }
    seq.mm$Idx <- as.numeric(as.character(seq.mm$Idx))
    bc <- rep(NA,nrow(Indel.df))
    bc[seq.mm$Idx] <- seq.mm$Muts
    nMut.df <- data.frame(BasesChanged = bc)

    #correct for missing indels (at ends of the sequence are not picked up by the indel() function)
    idx.indel <- which(Biostrings::nchar(dt$seq) != Biostrings::nchar(Reference.Sequence.DNA))
    idx.missed.indel <- idx.indel[!idx.indel %in% c(idx.ins,idx.del)]
    missed.indel.size <- Biostrings::nchar(dt$seq)[idx.missed.indel]-Biostrings::nchar(Reference.Sequence.DNA)
    missed.indel.size[missed.indel.size>0] <- str_c("+",missed.indel.size[missed.indel.size>0])
    Indel.df$Indels[idx.missed.indel] <- missed.indel.size

    idx.wt.temp <- which(!1:nrow(Indel.df) %in% seq.mm$Idx)
    idx.wt <- idx.wt.temp[is.na(Indel.df$Indels[idx.wt.temp])]
    if (length(idx.wt) > 1){
      stop(str_c("More than 1 WT sequence detected for sample ", x))
    }else if (length(idx.wt) < 1){
      stop(str_c("No WT sequence detected for sample ", x))
    }
