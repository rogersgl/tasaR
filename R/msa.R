#' @name tas_msa
#' @title Perform multiple sequence alignments of filtered and annotated sequences
#' @description
#' Uses labeled sequence tables to perform sequence alignments of top common sequences. The number of sequences to be aligned is specified in Config.List by sequence.alignment.count. Can be resource intensive, and increases greatly with more sequences. A default of 10 sequences is recommended as reasonable to run locally and also for readability of output figures. Will automatically trim sequence numbers if fewer sequence species are included for certain samples. Also supports optional amino acid analysis with the protein.mutations parameter or generation of phylogenetic trees with the PhyloTree parameter in Config.List.
#' @param Sample.Names A character vector with names of all samples.
#' @param Sequence.Table.List A list of data frames containing each unique DNA sequence for all samples. Names of each list entry must correspond to Sample.Names.
#' @param Input.DataFrame A data frame imported from the .xlsx file specifying details of each sample.
#' @param Config.List A list containing configuration parameters for tasAnalyzer.
#'
#' @returns A list containing ggplot objects of the specified analyses for each sample.
#' @export
tas_msa <- function(Sample.Names,Sequence.Table.List,Input.DataFrame,Config.List){

  Input.Columns <- "ReferenceSequence"

  Table.Columns <- c("TargetSequence",
                     "Percent")


  Config.Entries <- c("nCores",
                      "sequence.alignment.count")

  if(FALSE %in% (Input.Columns %in% colnames(Input.DataFrame))){
    stop(c("The following columns were not found in the input dataframe: ",
           str_flatten(Input.Columns[!Input.Columns %in% colnames(Input.DataFrame)],collapse = ", ")))
  }

  if(FALSE %in% (Config.Entries %in% names(Config.List))){
    stop(c("The following columns were not found in the config list: ",
           str_flatten(Config.Entries[!Config.Entries %in% names(Config.List)],collapse = ", ")))
  }

  Output_MSA <- list()
  Output_MSA[Sample.Names] <- lapply(Sample.Names,function(x){

    if(FALSE %in% (Table.Columns %in% colnames(Sequence.Table.List[[x]]))){
      stop(c("The following columns were not found in the input file: ",
             str_flatten(Table.Columns[!Table.Columns %in% colnames(Sequence.Table.List[[x]])],collapse = ", ")))
    }

    dna_WT <- DNAString(Input.DataFrame[x,"ReferenceSequence"])
    msa.length <- min(c(length(Sequence.Table.List[[x]]$TargetSequence)),as.numeric(Config.List$sequence.alignment.count))
    top <- c(DNAStringSet(dna_WT),DNAStringSet(Sequence.Table.List[[x]]$TargetSequence[1:msa.length]))
    ilabel <- Sequence.Table.List[[x]]$Indels[1:msa.length]
    ilabel[is.na(ilabel)] <- ""

    blabel <- Sequence.Table.List[[x]]$BasesChanged[1:msa.length]
    blabel[is.na(blabel)] <- ""
    blabel[Biostrings::nchar(blabel)>0] <- str_c(blabel[Biostrings::nchar(blabel)>0]," SNV")
    blabel[str_detect(blabel,"WT")] <- "WT"

    tlabel <- rep("",msa.length)
    tlabel[Biostrings::nchar(ilabel)>0 & Biostrings::nchar(blabel)>0] <- str_c(ilabel[Biostrings::nchar(ilabel)>0 & Biostrings::nchar(blabel)>0],
                                                       ", ",
                                                       blabel[Biostrings::nchar(ilabel)>0 & Biostrings::nchar(blabel)>0])
    tlabel[Biostrings::nchar(ilabel)>0 & Biostrings::nchar(blabel)==0] <- ilabel[Biostrings::nchar(ilabel)>0 & Biostrings::nchar(blabel)==0]
    tlabel[Biostrings::nchar(ilabel)==0 & Biostrings::nchar(blabel)>0] <- blabel[Biostrings::nchar(ilabel)==0 & Biostrings::nchar(blabel)>0]

    tlabel[str_detect(tlabel,"WT")] <- "WT"

    names(top) <- c("Reference",str_c(1:msa.length,". ",tlabel," - ",
                                      round(Sequence.Table.List[[x]]$Percent[1:msa.length],digits=2),
                                      "%"))
    sink(tempfile())
    top_align <- msa(top,order = "input",method = "ClustalOmega",type = "dna")
    sink()
    top_align <- DNAMultipleAlignment(as(top_align,"BStringSet"))
    msaDNA <- suppressMessages(
      ggmsa(top_align,
            consensus_views = TRUE,
            ref = "Reference",
            color = "Taylor_NT",
            char_width = 0.7,
            border = NA,
            seq_name = TRUE)+
        coord_cartesian()+
        facet_msa(field = 100)+
        theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = "in"))+
        theme(axis.text = element_text(size = 6))
    )
    Output <- list(DNA=list(Logo=msaDNA,Alignment=top_align))

    if (Config.List$protein.mutations==1){
      Table.Columns.Prot <- c("ProteinMutation","AA")
      if(FALSE %in% (Table.Columns.Prot %in% colnames(Sequence.Table.List[[x]]))){
        stop(c("The following columns were not found in the input file: ",
               str_flatten(Table.Columns.Prot[!Table.Columns.Prot %in% colnames(Sequence.Table.List[[x]])],collapse = ", ")))
      }

      prot_WT <- suppressWarnings(translate(dna_WT))

      y <- c(as.character(prot_WT),Sequence.Table.List[[x]]$AA[1:msa.length])
      ns <- str_detect(y,"\\*")
      s <- str_locate(y[ns],"\\*")[,"start"]
      names(s) <- NULL
      if(isEmpty(y[which(width(y)!=Biostrings::nchar(prot_WT))])==FALSE){
        y[ns] <- sapply(seq_along(y[ns]),function(z){
          str_trunc(y[ns][z],s[z]-1,ellipsis = "")
        })
      }

      topP <- AAStringSet(y)
      names(topP) <- c("Reference",str_c(1:msa.length,". ",
                                         Sequence.Table.List[[x]]$ProteinMutation[1:msa.length],
                                         " - ",
                                         round(Sequence.Table.List[[x]]$Percent[1:msa.length],digits=2),
                                         "%"))
      sink(tempfile())
      topP_align <- msa(topP,order = "input",method = "ClustalOmega",type = "protein")
      sink()
      topP_align <- AAMultipleAlignment(as(topP_align,"BStringSet"))

      msaProt <- suppressMessages(
        ggmsa(topP_align,
              consensus_views = TRUE,
              ref = "Reference",
              color = "Chemistry_AA",
              char_width = 0.7,
              border = NA,
              seq_name = TRUE)+
          coord_cartesian()+
          facet_msa(field = 50)+
          theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = "in"))+
          theme(axis.text = element_text(size = 10))
      )

      Output <- c(Output,list(Protein=list(Logo=msaProt,Alignment=topP)))
    }

    if (Config.List$PhyloTree==1){
      sink(tempfile())
      m <- msa(DNAStringSet(Sequence.Table.List[[x]]$TargetSequence[1:msa.length]),type = "dna",verbose = FALSE,method = "ClustalOmega")
      sink()
      ms <- msaConvert(m,type="seqinr::alignment")
      msd <- dist.alignment(ms)
      tree <- bionj(msd)
      tree$tip.label <- Sequence.Table.List[[x]]$ProteinMutation[1:msa.length]
      Output <- c(Output,list(PhyloTree=tree))
    }

    return(Output)
  }) #Single-threaded because this loop is very memory-intensive

  return(Output_MSA)
}





