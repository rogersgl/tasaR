#' @name tas_export_msa
#' @title Export multiple sequence alignment graphics
#' @description
#' Exports the multiple sequence alignment graphs created by tas_msa as .png (alignments) or .pdf (phylogenetic tree) objects.
#'
#' @param Sample.Names A character vector with names of all samples.
#' @param Output.MSA.List A list of multiple sequence alignment graphs for each sample. Names of each list entry must correspond to Sample.Names.
#' @param Config.List A list containing configuration parameters for tasAnalyzer.
#' @param WD A character vector indicating the file path to the working directory containing input files.
#'
#' @returns No returns within R. Writes image files to disk.
#' @export
tas_export_msa <- function(Sample.Names,Output.MSA.List,Config.List,WD){

  if (dir.exists(str_c(WD,"Export/"))==FALSE){
    dir.create(str_c(WD,"Export/"))
  }

  if (dir.exists(str_c(WD,"Export/MSA/"))==FALSE){
    dir.create(str_c(WD,"Export/MSA/"))
  }

  if (dir.exists(str_c(WD,"Export/MSA/DNA/"))==FALSE){
    dir.create(str_c(WD,"Export/MSA/DNA/"))
  }

  if (Config.List$protein.mutations==1){
    if (dir.exists(str_c(WD,"Export/MSA/Protein/"))==FALSE){
      dir.create(str_c(WD,"Export/MSA/Protein/"))
    }
  }

  if (Config.List$PhyloTree==1){
    if (dir.exists(str_c(WD,"Export/MSA/tree/"))==FALSE){
      dir.create(str_c(WD,"Export/MSA/tree/"))
    }
  }


  for (x in Sample.Names){
    ggsave(str_c(WD,"Export/MSA/DNA/",x,"-MSA.png"),plot = Output.MSA.List[[x]]$DNA$Logo,width = 11,height = 8.5*(nrow(Output.MSA.List[[x]]$DNA$Alignment)/(as.numeric(Config.List$sequence.alignment.count)+1)),units = "in")
    if (Config.List$protein.mutations==1){
      ggsave(str_c(WD,"Export/MSA/Protein/",x,"-MSA.png"),plot = Output.MSA.List[[x]]$Protein$Logo,width = 11,height = 8.5*(nrow(Output.MSA.List[[x]]$DNA$Alignment)/(as.numeric(Config.List$sequence.alignment.count)+1)),units = "in")
    }
    if (Config.List$PhyloTree==1){
      invisible(pdf(str_c(WD,"Export/MSA/tree/",x,"-tree.pdf"), width = 6, height = 6))
      invisible(plot(Output.MSA.List[[x]]$PhyloTree,cex=0.7))
      invisible(dev.off())
    }
  }
}

#' @name tas_export_mutations
#' @title Export graphs and spreadsheets
#' @description
#' Generates graphs and workbooks showing the results of sequence labeling and mutation analysis. Writes image files and .xlsx spreadsheets to disk for archive and compatibility with other software analyses.
#'
#' @param Sample.Names A character vector with names of all samples.
#' @param Output.Mutations.List A list of the analyzed mutation outputs for each sample. Names of each list entry must correspond to Sample.Names.
#' @param Sequence.Table.List A list of data frames containing each unique DNA sequence for all samples. Names of each list entry must correspond to Sample.Names.
#' @param Reference.Sequences.DNA A list containing the expected DNA sequences for each sample.
#' @param Input.DataFrame A data frame imported from the .xlsx file specifying details of each sample.
#' @param Config.List A list containing configuration parameters for tasAnalyzer.
#' @param WD A character vector indicating the file path to the working directory containing input files.
#'
#' @returns A list adding new tables and graphs for each sample to the original Output.Mutations.List.
#' @export
tas_export_mutations <- function(Sample.Names, Output.Mutations.List, Sequence.Table.List, Reference.Sequences.DNA, Input.DataFrame, Config.List, WD){

  if (dir.exists(str_c(WD,"Export/"))==FALSE){
    dir.create(str_c(WD,"Export/"))
  }

  if (dir.exists(str_c(WD,"Export/seqlogo/"))==FALSE){
    dir.create(str_c(WD,"Export/seqlogo/"))
  }

  for (i in Sample.Names){
    if (dir.exists(str_c(WD,"Export/seqlogo/",i,"/"))==FALSE){
      dir.create(str_c(WD,"Export/seqlogo/",i,"/"))
    }
  }

  if (dir.exists(str_c(WD,"Export/graphs/"))==FALSE){
    dir.create(str_c(WD,"Export/graphs/"))
  }

  if (dir.exists(str_c(WD,"Export/graphs/Mutagenesis/"))==FALSE){
    dir.create(str_c(WD,"Export/graphs/Mutagenesis/"))
  }

  if (dir.exists(str_c(WD,"Export/graphs/Cytosines labeled/"))==FALSE){
    dir.create(str_c(WD,"Export/graphs/Cytosines labeled/"))
  }

  if (dir.exists(str_c(WD,"Export/graphs/Box plots/"))==FALSE){
    dir.create(str_c(WD,"Export/graphs/Box plots/"))
  }

  if (dir.exists(str_c(WD,"Export/results/"))==FALSE){
    dir.create(str_c(WD,"Export/results/"))
  }

  if (dir.exists(str_c(WD,"Export/graphs/histograms/"))==FALSE){
    dir.create(str_c(WD,"Export/graphs/histograms/"))
  }

  Ab_regions <- c("FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4")

  Output.Mutations.List$Graphs <-  lapply(Sample.Names,function(i){
    Graphs <- list()
    #save sequence logo plots
    if(Input.DataFrame[i,"Antibody"]=="Yes"){
      Ab_pos <- c(as.numeric(Input.DataFrame[i,str_c(Ab_regions,"Start")]),Biostrings::nchar(Reference.Sequences.DNA[[i]])+1)
      names(Ab_pos) <- c(Ab_regions,"END1")
      for (r in 1:length(Ab_regions)){
        seq.len <- Ab_pos[r+1]-Ab_pos[r]+1
        ggsave(str_c(WD,"Export/seqlogo/",i,"/",Ab_regions[r],".pdf"),
               plot = Output.Mutations.List[[i]]$RegionMutagenesisDNA[[Ab_regions[r]]]$Logo,
               width=seq.len*0.2,
               height=1,
               units="in")
      }
    }

    #seq logo plots of whole sequence, length 100 bp
    dna.align <- rep(Sequence.Table.List$AlignDNA[[i]],Sequence.Table.List[[i]][,2])
    consensus_DNA <- consensusMatrix(dna.align)[c("A","C","G","T","-"),]
    n <- ncol(consensus_DNA)
    r <- ceiling(n/100)
    nr <- n/r
    s <- split(1:n,ceiling((1:n)/100))
    consensus_DNA_s <- list()
    consensus_DNA_s <- lapply(s,function(y){
      consensus_DNA[,y]
    })


    logo.len <- lapply(consensus_DNA_s,function(y){
      suppressMessages(ggseqlogo(y,method='prob')+
                         theme_logo(base_size=2.25)+
                         theme(axis.text.x = element_blank(),
                               axis.text.y = element_text(size = 2.5),
                               axis.title.y = element_text(size = 2.5))+
                         scale_x_continuous(expand = expansion(mult = 0, add = 0)))
    })


    for (y in 1:length(logo.len)){
      ggsave(str_c(WD,"Export/seqlogo/",i,"/",names(logo.len[y]),".pdf"),
             plot = logo.len[y],
             width=n*0.05+0.5,
             height=0.25,
             units="in")
    }

    ###########
    ### DNA ###
    ###########

    #Mutation bar graphs
    temp <- na.omit(readWorkbook(Output.Mutations.List$Workbooks$Mut_pos_wb, sheet = "MutAll")[i])
    temp <- cbind(data.frame(seq_along(1:nrow(temp))),temp)
    colnames(temp) <- c("nt","Mutated")

    #bars only
    MutBar <- ggplot(data = temp,
                     aes(x = nt,y = Mutated))+
      geom_col(stat = "identity",fill = "black",color = "black",linewidth = 0.1)+
      labs(y = "% Mutation")+
      theme_classic()+
      theme(axis.line = element_line(linewidth = 0.3, linetype = "solid",
                                     colour = "black"))
    ggsave(str_c(WD,"Export/graphs/Mutagenesis/",i,".pdf"),plot = MutBar,width = 6,height = 3,units = "in")

    if (Config.List$measure.shm==1){
      #with WRCH cytosines labeled
      MutBarC <- MutBar+
        geom_text(aes(label = ifelse((nt %in% Output.Mutations.List[[i]]$AID$WRCH$Positions$Cytosine), "C", "")),
                  position = position_dodge(width = .9), vjust = 0,
                  color = "red")
      ggsave(str_c(WD,"Export/graphs/Cytosines labeled/",i,"-C.pdf"),plot = MutBarC,width = 6,height = 3,units = "in")

      #dot plot with mutation frequencies of WRCH-C and nonC base pairs for each sample
      temp2 <- data.frame(Mutation=c(Output.Mutations.List[[i]]$AID$WRCH$MutCyt,Output.Mutations.List[[i]]$AID$WRCH$MutNonCyt),
                          Type=c(rep("AID Cytosines",length(Output.Mutations.List[[i]]$AID$WRCH$MutCyt)),rep("Other Bases",length(Output.Mutations.List[[i]]$AID$WRCH$MutNonCyt))))
      MutCnC <- ggplot(data = temp2,aes(x = Type,y = Mutation))+
        geom_boxplot(width = 0.5,outlier.size = 1)+
        labs(x = "Base Type",y = "% Mutation")+
        theme_classic()+
        theme(axis.line = element_line(linewidth = 0.3, linetype = "solid",
                                       colour = "black"))+
        scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(str_c(WD,"Export/graphs/Box plots/",i,".pdf"),plot = MutCnC,width = 1.5,height = 3,units = "in")
    }

    h <- suppressWarnings(as.numeric(Sequence.Table.List[[i]]$BasesChanged))
    h[Sequence.Table.List[[i]]$BasesChanged=="WT"] <- 0
    h <- rep(h,Sequence.Table.List[[i]][,2])
    h <- data.frame(NumberOfMutations=na.omit(h))
    MutH <- ggplot(data = h, aes(x = NumberOfMutations))+
      geom_histogram(stat = "bin",binwidth = 1)+
      labs(x = "# of mutations (indels omitted)",y = colnames(Sequence.Table.List[[i]])[2])+
      theme_classic()+
      theme(axis.line = element_line(linewidth = 0.3, linetype = "solid",
                                     colour = "black"))+
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    ggsave(str_c(WD,"Export/graphs/histograms/",i,"-hist.pdf"),plot = MutH,width = 4.5,height = 3,units = "in")

    Graphs[[i]] <- list(MutBar_DNA=MutBar,MutH=MutH)
    if (Config.List$measure.shm==1){
      Graphs[[i]] <- c(Graphs[[i]],list(MutBarC_DNA=MutBarC,MutCnC_DNA=MutCnC))
    }

    ###############
    ### PROTEIN ###
    ###############

    if (Config.List$protein.mutations==1){
      tempP <- Output.Mutations.List[[i]]$MutagenesisProtein
      tempP <- cbind(data.frame(seq_along(1:length(tempP))),tempP)
      colnames(tempP) <- c("AA","Mutated")

      #bars only
      MutBarP <- ggplot(data = tempP,
                        aes(x = AA,y = Mutated))+
        geom_col(stat = "identity",fill = "black",color = "black",linewidth = 0.1)+
        labs(y = "% Mutation")+
        theme_classic()+
        theme(axis.line = element_line(linewidth = 0.3, linetype = "solid",
                                       colour = "black"))
      ggsave(str_c(WD,"Export/graphs/Mutagenesis/",i,"-Protein.pdf"),plot = MutBarP,width = 6,height = 3,units = "in")

      #with WRCH cytosines labeled
      if (Config.List$measure.shm==1){
        AA_Cyt <- unique(floor((Output.Mutations.List[[i]]$AID$WRCH$Positions$Cytosine+3)/3))

        MutBarCP <- MutBarP+
          geom_text(aes(label = ifelse((AA %in% AA_Cyt), "C", "")),
                    position = position_dodge(width = .9), vjust = 0,
                    color = "red")
        ggsave(str_c(WD,"Export/graphs/Cytosines labeled/",i,"-C-Protein.pdf"),plot = MutBarCP,width = 6,height = 3,units = "in")


        #dot plot with mutation frequencies of WRCH-C and nonC base pairs for each sample
        temp2P <- data.frame(Mutation=tempP$Mutated,Type=rep("Other Amino Acids",length(tempP$Mutated)))
        temp2P[AA_Cyt,2] <- "Contains AID Cytosine"
        MutCnCP <- ggplot(data = temp2P,aes(x = Type,y = Mutation))+
          geom_boxplot(width = 0.5,outlier.size = 1)+
          labs(x = "Amino Acid",y = "% Mutation")+
          theme_classic()+
          theme(axis.line = element_line(linewidth = 0.3, linetype = "solid",
                                         colour = "black"))+
          scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(str_c(WD,"Export/graphs/Box plots/",i,"-Protein.pdf"),plot = MutCnCP,width = 1.5,height = 3,units = "in")
        Graphs[[i]] <- c(Graphs[[i]],list(MutBar_Protein=MutBarP,MutBarC_Protein=MutBarCP,MutCnC_Protein=MutCnCP))
      }
    }
    return(Graphs)
  }) #single threaded due to memory usage

  #summary plot for total % mutation across sequence for each sample
  if (dir.exists(str_c(WD,"Export/graphs/Summary"))==FALSE){
    dir.create(str_c(WD,"Export/graphs/Summary"))
  }

  tmut <- sapply(Sample.Names,function(y){
    sum(Output.Mutations.List[[y]]$MutagenesisDNA)/length(Output.Mutations.List[[y]]$MutagenesisDNA)
  })
  temp4 <- data.frame(Sample = Sample.Names,
                      Mutation = tmut,
                      row.names = NULL)
  colnames(temp4) <- c("Sample","Mutation")
  temp4$Sample <- factor(temp4$Sample,levels = temp4$Sample)
  MutAll <- ggplot(data = temp4,
                   aes(x = Sample, y = Mutation))+
    geom_col(stat = "identity",fill = "black",color = "black")+
    labs(y = "% Mutation across all bases")+
    theme_classic()+
    theme(axis.line = element_line(linewidth = 0.3, linetype = "solid",
                                   colour = "black"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(str_c(WD,"Export/graphs/Summary/All nt mutation summary.pdf"),plot = MutAll,width = 0.75*length(Sample.Names),height = 3,units = "in")

  Output.Mutations.List$SummaryGraphs <- list(MutAll=MutAll)

  #overview plot of showing mutagenesis of AID cytosines for all samples
  if (Config.List$measure.shm==1){

    MotifSums <- readWorkbook(Output.Mutations.List$Workbooks$Mut_pos_wb, sheet = "MotifSums")

    temp3 <- data.frame(Samples = colnames(MotifSums[,2:ncol(MotifSums),drop=FALSE]),
                        Mutation = as.numeric(MotifSums[MotifSums[,1]=="CytosineMutationAverage",2:ncol(MotifSums)]),
                        row.names = NULL)
    colnames(temp3) <- c("Sample","Mutation")
    temp3$Sample <- factor(temp3$Sample,levels = temp3$Sample)

    MutC_All <- ggplot(data = temp3,
                       aes(x = Sample, y = Mutation))+
      geom_col(stat = "identity",fill = "black",color = "black")+
      labs(y = "% Mutation at AID Cytosines")+
      theme_classic()+
      theme(axis.line = element_line(linewidth = 0.3, linetype = "solid",
                                     colour = "black"))+
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(str_c(WD,"Export/graphs/Summary/Cytosine mutation summary.pdf"),plot = MutC_All,width = 0.75*length(Sample.Names),height = 3,units = "in")

    mut_mut <- readWorkbook(Output.Mutations.List$Workbooks$Mut_pos_wb, sheet = "MutAll")
    mut_samp <- character()
    mut_pos <- numeric()
    for (i in seq_along(colnames(mut_mut))){
      mut_samp <- c(mut_samp,rep(colnames(mut_mut)[i],length(mut_mut[,i])))
      mut_pos <- c(mut_pos,1:length(mut_mut[,i]))
    }
    hm <- data.frame(Samples=mut_samp,nt=mut_pos,Mutation=unlist(mut_mut,use.names = FALSE))
    mut_hm <- ggplot(data = hm,
                     aes(x = nt, y = Samples, fill = Mutation))+
      geom_tile()+
      scale_fill_distiller(palette = "YlGnBu")+
      theme_classic()+
      theme(axis.line = element_line(linewidth = 0, linetype = "solid",
                                     colour = "black"))+
      scale_x_continuous(expand = expansion(mult = c(0.02, 0)))
    ggsave(str_c(WD,"Export/graphs/Summary/Mutation heatmap.pdf"),plot = MutC_All,width = 3,height = 0.75*length(Sample.Names),units = "in")

    nmps <- sapply(Sequence.Table.List,function(x){
      x$Indels[x$BasesChanged=="WT"] <- NA
      x$BasesChanged[x$BasesChanged=="WT"] <- 0
      as.numeric(x$BasesChanged[is.na(x$Indels)])
    })
    nmps.maxlen <- max(sapply(nmps,length))
    nmps <- unlist(lapply(nmps,function(x){
      c(x,rep(NA,nmps.maxlen-length(x)))
    }),use.names=FALSE)
    nmps.samp <- unlist(lapply(Sample.Names,function(x){
      rep(x,nmps.maxlen)
    }),use.names=FALSE)

    nmps_df <- data.frame(Sample = nmps.samp, Mutations = nmps)
    nMutPerSeq <- suppressWarnings(ggplot(nmps_df,
                                          aes(x = Sample, y = Mutations))+
                                     geom_jitter(alpha = 0.12)+
                                     labs(y = "# of mutations / sequence")+
                                     theme_classic()+
                                     theme(axis.line = element_line(linewidth = 0.3, linetype = "solid",
                                     colour = "black"))+
                                     scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
                                     theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    ggsave(str_c(WD,"Export/graphs/Summary/Mutations per sequence dotplot.pdf"),plot = nMutPerSeq,width = 0.75*length(Sample.Names),height = 3,units = "in")
    Output.Mutations.List$SummaryGraphs <- c(Output.Mutations.List$SummaryGraphs,list(MutC_All=MutC_All,MutHM=mut_hm,nMutPerSeq=nMutPerSeq))

    if (Config.List$measure.shm==1){
      mutc_mut <- readWorkbook(Output.Mutations.List$Workbooks$Mut_pos_wb, sheet = "MutCyt")
      mutc_samp <- character()
      mutc_pos <- numeric()
      for (i in seq_along(colnames(mutc_mut))){
        mutc_samp <- c(mutc_samp,rep(colnames(mutc_mut)[i],length(mutc_mut[,i])))
        mutc_pos <- c(mutc_pos,1:length(mutc_mut[,i]))
      }
      hmc <- data.frame(Samples=mutc_samp,nt=mutc_pos,Mutation=unlist(mutc_mut,use.names = FALSE))
      mutc_hm <- ggplot(data = hmc,
                       aes(x = nt, y = Samples, fill = Mutation))+
        geom_tile()+
        scale_fill_distiller(palette = "YlGnBu")+
        theme_classic()+
        theme(axis.line = element_line(linewidth = 0, linetype = "solid",
                                       colour = "black"))+
        scale_x_continuous(expand = expansion(mult = c(0.02, 0)))+
        labs(x = "AID cytosine #")
      ggsave(str_c(WD,"Export/graphs/Summary/Mutation cytosines heatmap.pdf"),plot = MutC_All,width = 0.75*length(Sample.Names),height = 3,units = "in")
      Output.Mutations.List$SummaryGraphs <- c(Output.Mutations.List$SummaryGraphs,list(MutC_HM=mutc_hm))
    }

  }

  if (Config.List$dna.repair.pathways==1){
    #stacked bar graph of mutation types
    mut_label <- c("WT","NHEJ","MMEJ","Base Change","Indel + Base Change","Other")
    mut.types <- data.frame(sapply(Sample.Names,function(x){
      idx_wt <- which(Sequence.Table.List[[x]]$Indels=="WT")
      idx_nhej <- which(suppressWarnings(as.numeric(Sequence.Table.List[[x]]$Indels) >= -2))
      idx_mmej <- which(suppressWarnings(as.numeric(Sequence.Table.List[[x]]$Indels) < -2))
      idx_change <- which(suppressWarnings(as.numeric(Sequence.Table.List[[x]]$BasesChanged) > 0))
      idx_indel_bc <- intersect(c(idx_nhej,idx_mmej),idx_change)
      idx_nhej <- idx_nhej[!idx_nhej %in% idx_indel_bc]
      idx_mmej <- idx_mmej[!idx_mmej %in% idx_indel_bc]
      idx_change <- idx_change[!idx_change %in% idx_indel_bc]
      idx_other <- which(!1:nrow(Sequence.Table.List[[x]]) %in% c(idx_wt,idx_nhej,idx_mmej,idx_change,idx_indel_bc))

      sums <- c(sum(Sequence.Table.List[[x]]$Percent[idx_wt]),
                sum(Sequence.Table.List[[x]]$Percent[idx_nhej]),
                sum(Sequence.Table.List[[x]]$Percent[idx_mmej]),
                sum(Sequence.Table.List[[x]]$Percent[idx_change]),
                sum(Sequence.Table.List[[x]]$Percent[idx_indel_bc]),
                sum(Sequence.Table.List[[x]]$Percent[idx_other]))
      names(sums) <- mut_label
      if (abs(sum(sums)-100)>1e-7){
        stop(str_c("Frequency summation error detected in indel types for sample: "),x)
      }
      return(sums)
    }))
    colnames(mut.types) <- Sample.Names
    df <- data.frame(Samples = unlist(lapply(Sample.Names,function(x){
                        rep(x,length(mut_label))
                      }),use.names = FALSE),
                     MutationTypes = rep(rownames(mut.types),length(Sample.Names)),
                     Percentage = unlist(mut.types,use.names = FALSE))
    df$Samples <- factor(df$Samples,levels = Sample.Names)
    df$MutationTypes <- factor(df$MutationTypes, levels = mut_label)
    mut.types.graph <- ggplot(df,aes(x = Samples, y = Percentage, fill = MutationTypes)) +
      geom_col(position = position_stack(reverse = TRUE),stat = "identity")+
      #labs(y = "% Mutation at AID Cytosines")+
      theme_classic()+
      theme(axis.line = element_line(linewidth = 0.3, linetype = "solid",
                                     colour = "black"))+
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.text = element_text(size = 9),
            axis.text = element_text(size = 9))+
      scale_fill_brewer(type = "qual",palette = "Set1")+
      labs(x = "Sample",
           fill = "Mutation Type")
    ggsave(str_c(WD,"Export/graphs/Summary/Mutation Types.pdf"), plot = mut.types.graph, width = 4.5, height = 3, units = "in")
    write.xlsx(mut.types,str_c(WD,"Export/results/Mutation Types.xlsx"),overwrite = TRUE,colNames = TRUE, rowNames = TRUE)
    Output.Mutations.List$SummaryGraphs <- c(Output.Mutations.List$SummaryGraphs,list(MutTypes=mut.types.graph))
  }

  saveWorkbook(Output.Mutations.List$Workbooks$Mut_pos_wb,str_c(WD,"Export/results/Mutations.xlsx"),overwrite = TRUE)
  saveWorkbook(Output.Mutations.List$Workbooks$Sequences_wb,str_c(WD,"Export/results/Sequences.xlsx"),overwrite = TRUE)
  if (Config.List$measure.shm==1){
    saveWorkbook(Output.Mutations.List$Workbooks$WRCY_wb,str_c(WD,"Export/results/WRCY Tables.xlsx"),overwrite = TRUE)
    saveWorkbook(Output.Mutations.List$Workbooks$WRCH_wb,str_c(WD,"Export/results/WRCH Tables.xlsx"),overwrite = TRUE)
  }

  return(Output.Mutations.List)

}
