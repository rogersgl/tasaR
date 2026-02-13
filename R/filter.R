#' @name tas_filter
#' @title Filtering on paired end reads.
#' @description
#' Filters merged paired end read files. Supports 5' or 3' indexing barcodes or unique molecular identifiers (UMIs).
#' @param Sample.Names A character vector with names of all samples.
#' @param Reads.List A list of merged reads in ShortReads format. Names of each list entry must correspond to Sample.Names.
#' @param Input.DataFrame A data frame imported from the .xlsx file specifying details of each sample.
#' @param Config.List A list containing configuration parameters for tasAnalyzer.
#' @returns A list of filtered sequences, in class ShortRead.
#' @export
tas_filter <- function(Sample.Names,Reads.List,Input.DataFrame,Config.List){
  Input.Columns <- c("SampleName",
                     "ForwardExtensionType",
                     "ForwardExtension",
                     "ForwardPrimer",
                     "ReverseExtensionType",
                     "ReverseExtension",
                     "ReversePrimer",
                     "AmpliconLength",
                     "MaxDeletion",
                     "MaxInsertion")

  Config.Entries <- c("WorkingDirectory",
                      "nCores")

  if (FALSE %in% (Input.Columns %in% colnames(Input.DataFrame))){
    stop(c("The following columns were not found in the input file: ",
           str_flatten(Input.Columns[!Input.Columns %in% colnames(Input.DataFrame)],collapse = ", ")))
  }

  if (FALSE %in% (Config.Entries %in% names(Config.List))){
    stop(c("The following columns were not found in the input file: ",
           str_flatten(Config.Entries[!Config.Entries %in% names(Config.List)],collapse = ", ")))
  }

  regex.dna <- str_flatten(DNA_ALPHABET[1:15])

  if (FALSE %in% (nzchar(Input.DataFrame$ForwardExtensionType)) & TRUE %in% (str_detect(Input.DataFrame$ForwardExtension, paste0("[^", regex.dna, "]")))){
    stop("ForwardExtension must be a DNA sequence if ForwardExtensionType is provided.")
  }

  if (FALSE %in% (nzchar(Input.DataFrame$ReverseExtensionType)) & TRUE %in% (str_detect(Input.DataFrame$ReverseExtension, paste0("[^", regex.dna, "]")))){
    stop("ReverseExtension must be a DNA sequence if ReverseExtensionType is provided.")
  }

  if (TRUE %in% (str_detect(Input.DataFrame$ForwardPrimer, paste0("[^", regex.dna, "]")))){
    stop("ForwardPrimer must be a DNA sequence.")
  }

  if (TRUE %in% (str_detect(Input.DataFrame$ReversePrimer, paste0("[^", regex.dna, "]")))){
    stop("ReversePrimer must be a DNA sequence.")
  }

  Reads.Filtered <- list()
  Reads.Filtered[Sample.Names] <- mclapply(Sample.Names,function(x){
    fwd <- DNAString(str_c(Input.DataFrame[x,"ForwardExtension"],Input.DataFrame[x,"ForwardPrimer"]))
    rev <- DNAString(str_c(Input.DataFrame[x,"ReverseExtension"],Input.DataFrame[x,"ReversePrimer"]))
    temp <- Reads.List[[x]]@sread
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
    Reads.List[[x]][both_idx]
  },mc.cores = Config.List$nCores)

  return(Reads.Filtered)
}



#' @name tas_seq_table
#' @title Makes data frames with details of all filtered sequences for each sample.
#' @description
#' Processes filtered reads to generate data frames with details of all filtered sequences for each sample. Supports samples with or without unique moleculear identifiers (UMIs). If UMIs are included, sequences will be processed based on UMIs. Otherwise, unique sequences are merged based on read count.
#' @param Sample.Names A character vector with names of all samples.
#' @param Reads.Filtered.List A list of filtered reads in ShortReads format (output from tas_filter). Names of each list entry must correspond to Sample.Names.
#' @param Input.DataFrame A data frame imported from the .xlsx file specifying details of each sample.
#' @param Config.List A list containing configuration parameters for tasAnalyzer.
#' @returns A list of sequence tables (data frames) for each sample. If UMIs were included, also returns a similar table with all sequences for counting filter statistics.
#' @import data.table
#' @export
tas_seq_table <- function(Sample.Names,Reads.Filtered.List,Input.DataFrame,Config.List){

  Input.Columns <- c("InsertStart",
                     "InsertEnd",
                     "ReferenceSequence",
                     "ForwardExtensionType",
                     "ForwardExtension",
                     "ForwardPrimer",
                     "ReverseExtensionType",
                     "ReverseExtension",
                     "ReversePrimer")

  Config.Entries <- c("nCores",
                      "protein.mutations")

  if(FALSE %in% (Input.Columns %in% colnames(Input.DataFrame))){
    stop(c("The following columns were not found in the input file: ",
           str_flatten(Input.Columns[!Input.Columns %in% colnames(Input.DataFrame)],collapse = ", ")))
  }

  if(FALSE %in% (Config.Entries %in% names(Config.List))){
    stop(c("The following columns were not found in the input file: ",
           str_flatten(Config.Entries[!Config.Entries %in% names(Config.List)],collapse = ", ")))
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
    if (!is.na(umi.pos)){
      if (is.null(umi.pos$end)){
        umi_temp <- narrow(Reads.Filtered.List[[x]], start = umi.pos$start)@sread
      }else{
        umi_temp <- narrow(Reads.Filtered.List[[x]], start = umi.pos$start, end = umi.pos$end)@sread
      }
      seq_temp <- narrow(Reads.Filtered.List[[x]],start = Input.DataFrame[x,"InsertStart"],end = Input.DataFrame[x,"InsertEnd"])@sread
      n_idx <- which(elementNROWS(Biostrings::vmatchPattern("n",seq_temp))==0)
      umi_temp <- umi_temp[n_idx]
      seq_temp <- seq_temp[n_idx]

      umi_pileup <- function(i){
        t <- ShortRead::tables(i, n = length(i))
        d <- t$distribution$nOccurrences[order(t$distribution$nOccurrences,decreasing = TRUE)]
        if((length(d)==1) | (d[1]>=(d[2]+2) & (d[2]<d[1]*0.7))){
          return(names(t$top[1]))
        }else{
          return("Rejected")
        }
      }
      a <- ShortRead::tables(umi_temp,n=length(umi_temp))$top[ShortRead::tables(umi_temp,n=length(umi_temp))$top>=3]
      umi_df <- data.frame(UMI=names(a),Reads=a,row.names=NULL)
      b <- BiocGenerics::match(as.character(umi_temp),umi_df$UMI)
      dt <- data.table::data.table(value = b, idx = seq_along(b),seq = as.character(seq_temp))
      result_dt <- dt[, list(sequences = umi_pileup(DNAStringSet(seq))), by = "value"]
      df <- as.data.frame(result_dt)
      df <- na.omit(df[order(df$value),])
      df <- data.frame(df[,2])
      colnames(df) <- "TargetSequence"
      umi_df <- cbind(umi_df,df)
      Output <- umi_df[!df$TargetSequence=="Rejected",]
      Output <- Output[,c("TargetSequence","Reads","UMI")]
      return(Output)
    }else{
      Output <- list()
      return(Output)
    }
  }, mc.cores = Config.List$nCores)

  Sequence.Table <- list()
  Sequence.Table[Sample.Names] <- mclapply(Sample.Names,function(x){
    # print(str_c("table , ",x))
    if (!S4Vectors::isEmpty(Sequence.Table.All[[x]])){
      uCounts <- as.data.frame(BiocGenerics::table(Sequence.Table.All[[x]]$TargetSequence))
      colnames(uCounts) <- c("TargetSequence","UMI Count")
      uCounts <- uCounts[order(uCounts$`UMI Count`,decreasing = TRUE),]
      uCounts <- cbind(uCounts, data.frame(Percent=(uCounts$`UMI Count`/sum(uCounts$`UMI Count`))*100))
      uCounts$TargetSequence <- as.character(uCounts$TargetSequence)
      return(uCounts)
    }else if (S4Vectors::isEmpty(Sequence.Table.All[[x]])){
      seq_temp <- narrow(Reads.Filtered.List[[x]],start = Input.DataFrame[x,"InsertStart"],end = Input.DataFrame[x,"InsertEnd"])@sread
      seq_temp <- seq_temp[which(elementNROWS(vmatchPattern("n",seq_temp))==0)]
      seq_temp_t <- BiocGenerics::table(seq_temp)
      uCounts <- data.frame(Reads=data.frame(seq_temp_t,row.names = NULL))
      colnames(uCounts) <- c("TargetSequence","Reads")
      uCounts$TargetSequence <- as.character(uCounts$TargetSequence)
      uCounts <- uCounts[order(uCounts$Reads,decreasing = TRUE),]
      Percent <- uCounts$Reads/sum(uCounts$Reads)*100
      uCounts <- uCounts[Percent>=as.numeric(Config.List$read.frequency.limit),]
      uCounts <- cbind(uCounts,(uCounts$Reads/sum(uCounts$Reads)*100))
      colnames(uCounts) <- c("TargetSequence","Reads","Percent")
      return(uCounts)
    }
  },mc.cores = Config.List$nCores)

  Sequence.Table <- c(Sequence.Table,list(All=Sequence.Table.All))

  return(Sequence.Table)
}




#' @name tas_label_table
#' @title Labels sequence tables with mutations.
#' @description
#' Classifies each unique sequence for each sample based on indel type and number of mutated bases. If the amino acid (AA) module is enabled, also determines protein sequences and mutations for each DNA sequence.
#' @param Sample.Names A character vector with names of all samples.
#' @param Sequence.Table.List A list of data frames containing each unique DNA sequence for all samples. Names of each list entry must correspond to Sample.Names.
#' @param Reference.Sequences.DNA A list of DNAString objects containing expected DNA sequences for all samples. Names of each list entry must correspond to Sample.Names.
#' @param Input.DataFrame A data frame imported from the .xlsx file specifying details of each sample.
#' @param Config.List A list containing configuration parameters for tasAnaly <- er.
#' @returns A list of labeled sequence tables (data frames) for each sample.
#' @export
tas_label_table <- function(Sample.Names,Sequence.Table.List,Reference.Sequences.DNA,Input.DataFrame,Config.List){
  Sequence.Table <- Sequence.Table.List
  #pairwise alignments of DNA

  Reads.Unique.DNA <- list()
  Reads.Unique.DNA[Sample.Names] <- mclapply(Sample.Names,function(x){
    DNAStringSet(Sequence.Table[[x]]$TargetSequence)
  },mc.cores = Config.List$nCores)

  Pairwise.Aligned.DNA <- list()
  Pairwise.Aligned.DNA[Sample.Names] <- mclapply(Sample.Names,function(x){
    pairwiseAlignment(Reads.Unique.DNA[[x]],Reference.Sequences.DNA[[x]])
  },mc.cores = Config.List$nCores)

  Sequence.Table[Sample.Names] <- mclapply(Sample.Names,function(x){
    # print(x)
    ins <-  indel(Pairwise.Aligned.DNA[[x]])@insertion
    idx.ins <- which(!sapply(ins,S4Vectors::isEmpty))
    ins.num <- lapply(1:length(ins),function(y){
      ins.widths <- ins[y]@unlistData@width
      if (!S4Vectors::isEmpty(ins.widths)){
        widths.split <- unlist(str_split(ins.widths,"-"))
        output <- str_c("+",widths.split)
      }else{
        output <- ""
      }
      return(output)
    })

    del <- indel(Pairwise.Aligned.DNA[[x]])@deletion
    idx.del <- which(!sapply(del,S4Vectors::isEmpty))
    del.num <- lapply(1:length(del),function(y){
      del.widths <- del[y]@unlistData@width
      if (!S4Vectors::isEmpty(del.widths)){
        widths.split <- unlist(str_split(del.widths,"-"))
        output <- str_c("-",widths.split)
      }else{
        output <- ""
      }
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

    seq.mm <- data.frame(BiocGenerics::table(mismatchTable(Pairwise.Aligned.DNA[[x]])$PatternId))
    if (!S4Vectors::isEmpty(seq.mm)){
      colnames(seq.mm) <- c("Idx","Muts")
    }
    seq.mm$Idx <- as.numeric(as.character(seq.mm$Idx))
    bc <- rep(NA,nrow(Indel.df))
    bc[seq.mm$Idx] <- seq.mm$Muts
    nMut.df <- data.frame(BasesChanged = bc)

    #correct for missing indels (at ends of the sequence are not picked up by the indel() function)
    idx.indel <- which(Biostrings::nchar(Sequence.Table[[x]]$TargetSequence)!=Biostrings::nchar(Reference.Sequences.DNA[[x]]))
    idx.missed.indel <- idx.indel[!idx.indel %in% c(idx.ins,idx.del)]
    missed.indel.size <- Biostrings::nchar(Sequence.Table[[x]]$TargetSequence)[idx.missed.indel]-Biostrings::nchar(Reference.Sequences.DNA[[x]])
    missed.indel.size[missed.indel.size>0] <- str_c("+",missed.indel.size[missed.indel.size>0])
    Indel.df$Indels[idx.missed.indel] <- missed.indel.size

    idx.wt.temp <- which(!1:nrow(Indel.df) %in% seq.mm$Idx)
    idx.wt <- idx.wt.temp[is.na(Indel.df$Indels[idx.wt.temp])]
    if (length(idx.wt)!=1){
      stop("More than 1 WT sequence detected.")
    }

    Indel.df$Indels[idx.wt] <- "WT"
    nMut.df$BasesChanged[idx.wt] <- "WT"

    cbind(Sequence.Table[[x]],Indel.df,nMut.df)

  }, mc.cores = Config.List$nCores)

  Sequence.Table <- c(Sequence.Table,list(AlignDNA=Pairwise.Aligned.DNA))

  #check if protein analyses should be performed (Config.List$protein.mutations = 1)
  if (Config.List$protein.mutations==1){

    Reference.Sequences.Protein <- list()
    Reference.Sequences.Protein <- lapply(Reference.Sequences.DNA,function(x){
      suppressWarnings(translate(x))
    })

    Reads.Unique.Protein <- list()
    Reads.Unique.Protein[Sample.Names] <- mclapply(Sample.Names,function(x){
      suppressWarnings(translate(DNAStringSet(gsub("-","",DNAStringSet(Reads.Unique.DNA[[x]])))))
    },mc.cores = Config.List$nCores)

    Pairwise.Aligned.Protein <- list()
    Pairwise.Aligned.Protein[Sample.Names] <- mclapply(Sample.Names,function(x){
      pairwiseAlignment(Reads.Unique.Protein[[x]],Reference.Sequences.Protein[[x]])
    },mc.cores = Config.List$nCores)

    #label sequences with protein sequence and mutations
    Sequence.Table[Sample.Names] <- mclapply(Sample.Names,function(x){
      # print(x)
      #amino acid sequences
      aa <- as.character(Reads.Unique.Protein[[x]])

      if (!S4Vectors::isEmpty(which(Biostrings::nchar(aa)!=Biostrings::nchar(Reference.Sequences.Protein[[x]])))){
        ns <- str_detect(aa,"\\*")
        s <- str_locate(aa[ns],"\\*")[,"start"]
        names(s) <- NULL
        aa[ns] <-  sapply(1:length(s),function(y){
          i <- which(ns)[y]
          str_trunc(aa[i],s[y],ellipsis = "")
        })
      }

      #Label each sequence with protein mutations
      PAP <- Pairwise.Aligned.Protein[[x]]
      WT_Protein <- Reference.Sequences.Protein[[x]]
      aa.mut <- data.frame(rep("",length(PAP)))
      colnames(aa.mut) <- "ProteinMutation"

      mmT <- mismatchTable(PAP)
      idx_indel <- !Biostrings::nchar(aa)==Biostrings::nchar(WT_Protein)
      idx_WT <- !seq_along(PAP) %in% mmT$PatternId #indel multiples of 3 identified incorrectly
      idx_WT[which((idx_indel+idx_WT)==2)] <- FALSE #remove false +ve from idx_WT
      idx_mm <- which((idx_indel+idx_WT)==0)
      mmT <- mmT[which(mmT$PatternId %in% idx_mm),]
      mmChar <- cbind(data.frame(PatternId=mmT$PatternId),data.frame(Mutation=str_c(mmT[,"SubjectSubstring"],mmT[,"SubjectStart"],mmT[,"PatternSubstring"])))

      mutCollapse <- function(x){
        str_flatten(x,", ")
      }

      aa.mut$ProteinMutation[idx_mm] <- aggregate(Mutation ~ PatternId,data = mmChar,FUN = mutCollapse)[,"Mutation"]
      aa.mut$ProteinMutation[idx_indel] <- "Indel"
      aa.mut$ProteinMutation[idx_WT] <- "WT"
      aa.mut$ProteinMutation[str_detect(aa.mut$ProteinMutation,"\\*")] <- "Nonsense"

      cbind(Sequence.Table[[x]],data.frame(AA=aa),aa.mut)

    }, mc.cores = Config.List$nCores)

    Sequence.Table <- c(Sequence.Table,list(AlignProtein=Pairwise.Aligned.Protein))
  }

  return(Sequence.Table)
}



#' @name tas_write_filtered
#' @title Writes and compresses .fastq files after filtering.
#' @description
#' After sequences are filtered, writes .fastq files of filtered reads, and compresses files as .gz files.
#' @param Sample.Names A character vector with names of all samples.
#' @param Reads.Filtered.List A list of filtered reads in ShortReads format (output from tas_filter). Names of each list entry must correspond to Sample.Names.
#' @param Config.List A list containing configuration parameters for tasAnalyzer.
#' @param WD A character vector indicating the file path to the working directory containing input files.
#' @returns No returns within R. Writes files to disk.
#' @export
tas_write_filtered <- function(Sample.Names,Reads.Filtered.List,Config.List,WD){
  if (dir.exists(str_c(WD,"filtered"))==FALSE){
    dir.create(str_c(WD,"filtered"))
  }

  invisible(mclapply(seq_along(Reads.Filtered.List),function(x){
    if (file.exists(str_c(WD,"filtered/",Sample.Names[x],"-filtered.fastq.gz"))==TRUE){
      print(str_c("The following file was replaced: ",str_c(WD,"filtered/",Sample.Names[x],"-filtered.fastq.gz")))
      file.remove(str_c(WD,"filtered/",Sample.Names[x],"-filtered.fastq.gz"))
    }
    writeFastq(Reads.Filtered.List[[x]],str_c(WD,"filtered/",Sample.Names[x],"-filtered.fastq.gz"))
  },mc.cores = Config.List$nCores))
}




#' @name tas_filter_count
#' @title Counts number of filtered reads at each step of analysis.
#' @description
#' After filtration and sequence processing, counts the number of filtered reads after each step. Writes the output to logs/Filtering.xlsx.
#' @param Sample.Names A character vector with names of all samples.
#' @param Reads.Filtered.List A list of filtered reads in ShortReads format (output from tas_filter). Names of each list entry must correspond to Sample.Names.
#' @param Sequence.Table.List A list of data frames containing each unique DNA sequence for all samples. Names of each list entry must correspond to Sample.Names.
#' @param Input.DataFrame The data frame imported from the .xlsx file specifying details of each sample.
#' @param Config.List A list containing configuration parameters for tasAnalyzer.
#' @param WD A character vector indicating the file path to the working directory containing input files.
#' @returns A data frame with read filtration counts.
#' @export
tas_filter_count <- function(Sample.Names,Reads.Filtered.List,Sequence.Table.List,Input.DataFrame,Config.List,WD){
  Read.Filter.Count <- data.frame(t(sapply(seq_along(Sample.Names),function(x){
    a <- length(ShortRead::readFastq(str_c(WD,"unpaired/",Input.DataFrame$ForwardFASTQFileName[[x]])))
    b <- length(ShortRead::readFastq(str_c(WD,"merged/",Input.DataFrame$SampleName[[x]],"-merged.fastq.gz")))
    c(a,b)
  })))
  rownames(Read.Filter.Count) <- Sample.Names
  colnames(Read.Filter.Count) <- c("Raw","Merged")

  Read.Filter.Count <- cbind(Read.Filter.Count,data.frame(Filtered=unlist(lapply(Reads.Filtered.List,length))))

  df <- data.frame(unlist(sapply(Sample.Names,function(x){
    if (!S4Vectors::isEmpty(Sequence.Table.List$All[[x]])){
      BiocGenerics::nrow(Sequence.Table.List$All[[x]])
    }else{
      NA
    }
  })))
  colnames(df) <- "UMIs"
  if (FALSE %in% is.na(df)){
    Read.Filter.Count <- cbind(Read.Filter.Count,df)
  }

  df2 <- data.frame(sapply(Sequence.Table.List[Sample.Names],BiocGenerics::nrow))
  colnames(df2) <- "UniqueSequences"
  Read.Filter.Count <- cbind(Read.Filter.Count,df2)

  return(Read.Filter.Count)
}
