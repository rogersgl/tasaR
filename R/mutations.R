#' @name tas_find_AID_targets
#' @title Finds and annotates AID hotspot motifs.
#' @description
#' Uses the reference sequences of each sample being analyzed to identify AID hotspot motifs. Identification is based on both the common WRCY and the slightly more inclusive WRCH motif. Evaluates both strands.
#'
#' @param Sample.Names A character vector with names of all samples.
#' @param Reference.Sequences.DNA A list containing the expected DNA sequences for each sample.
#' @param Config.List A list containing configuration parameters for tasAnalyzer.
#'
#' @returns A list of AID hotspot coordinates for each sample.
#' @export
tas_find_AID_Targets <- function(Sample.Names, Reference.Sequences.DNA, Config.List){
  AID.Motifs <- c("WRCY","RGYW","WRCH","DGYW")

  AID.Targets <- list()
  AID.Targets[Sample.Names] <- lapply(Reference.Sequences.DNA,function(x){
    a <- list()
    a[AID.Motifs] <- lapply(AID.Motifs,function(y){
      b <- matchPattern(y,x,fixed = FALSE)
      s <- start(b)
      e <- end(b)
      if (y=="WRCY" | y=="WRCH"){
        cyt <- e-1
      }
      if (y=="RGYW" | y=="DGYW"){
        cyt <- e-2
      }
      data.frame(Motif=rep(y,length(b)),Start=s,End=e,Cytosine=cyt)
    })
    f <- rbind(a[["WRCY"]],a[["RGYW"]])
    f <- f[order(f$Start),]
    g <- rbind(a[["WRCH"]],a[["DGYW"]])
    g <- g[order(g$Start),]
    list(WRCY=f,WRCH=g)
  })

  return(AID.Targets)
}


#' @name tas_measure_mutations
#' @title Measure mutations in sequence tables
#' @description
#' Uses provided reference DNA sequences and the sequence tables generated during filtration to measure mutations in the sequence population. Optional modules support analysis of AID hotspot motifs, breakdown by antibody regions if applicable, and the translated amino acid sequence.
#'
#' @param Sample.Names A character vector with names of all samples.
#' @param Sequence.Table.List A list of data frames containing each unique DNA sequence for all samples. Names of each list entry must correspond to Sample.Names.
#' @param AID.Targets A list of AID hotspot coordinates for each sample.
#' @param Input.DataFrame A data frame imported from the .xlsx file specifying details of each sample.
#' @param Config.List A list containing configuration parameters for tasAnalyzer.
#'
#' @returns A list providing the analyzed outputs for each sample.
#' @export
tas_measure_mutations <- function(Sample.Names, Sequence.Table.List, AID.Targets, Input.DataFrame, Config.List){

   Input.Columns <- c("ReferenceSequence",
                      "Antibody",
                      "FR1Start",
                      "CDR1Start",
                      "FR2Start",
                      "CDR2Start",
                      "FR3Start",
                      "CDR3Start",
                      "FR4Start")

  if(FALSE %in% (Input.Columns %in% colnames(Input.DataFrame))){
    stop(c("The following columns were not found in the input file: ",
           str_flatten(Input.Columns[!Input.Columns %in% colnames(Input.DataFrame)],collapse = ", ")))
  }

  Output.Mutagenesis <- list()
  Output.Mutagenesis[Sample.Names] <- lapply(Sample.Names, function(x){

    WT_DNA <- DNAStringSet(Input.DataFrame[x,"ReferenceSequence"])
    seqs <- DNAStringSet(Sequence.Table.List[[x]]$TargetSequence)
    dna.align <- pairwiseAlignment(seqs, WT_DNA)
    dna.align.all <- rep(dna.align, Sequence.Table.List[[x]][,2])

    #calculate consensus matrixes and mutagenesis frequency by position for DNA
    WT_DNA <- DNAStringSet(Input.DataFrame[x,"ReferenceSequence"])
    consensus_DNA_WT <- consensusMatrix(WT_DNA)[c("A","C","G","T","-"),]
    consensus_DNA_WT_flip <- 1-consensus_DNA_WT
    consensus_DNA <- consensusMatrix(dna.align)[c("A","C","G","T","-"),]
    mutagenesis_DNA <- (colSums(consensus_DNA*consensus_DNA_WT_flip)/length(dna.align))*100
    Output <- list(MutagenesisDNA=mutagenesis_DNA)

    #calculate consensus matrixes and mutagenesis frequency by position for Protein
    if (Config.List$protein.mutations==1){
      WT_Protein <- suppressWarnings(AAStringSet(translate(WT_DNA)))
      p.seqs <-  AAStringSet(Sequence.Table.List[[x]]$AA)
      #p.seqs <- AAStringSet(rep(Sequence.Table.List[[x]]$AA, Sequence.Table.List[[x]][,2]))
      protein.align <- pairwiseAlignment(p.seqs, WT_Protein)
      protein.align <- rep(protein.align, Sequence.Table.List[[x]][,2])
      #protein.align <- rep(Sequence.Table.List$AlignProtein[[x]], Sequence.Table.List[[x]]$Reads)

      AA <- c(AA_STANDARD,"*","-")
      consensus_Protein_WT <- pwalign::consensusMatrix(WT_Protein)[AA,]
      consensus_Protein_WT_flip <- 1-consensus_Protein_WT

      # observed crashing here without the additional alignedPattern() call inside consensusMatrix
      # Saw hard crashes of R, as well as error messages suggesting a misassigned pointer issue (perhaps during callout to C?)
      # Example: Error in h(simpleError(msg, call)) : error in evaluating the argument 'x' in selecting a method for function 'consensusMatrix': R_ExternalPtrTag: argument of type LISTSXP is not an external pointer
      # Suspect there may be a bug in the method dispatch in the pwalign package for an amino acid PairwiseAlignedSingleSubject? Odd though because DNA worked fine.
      consensus_Protein <- pwalign::consensusMatrix(alignedPattern(protein.align))[AA,]


      mutagenesis_Protein <- list(MutagenesisProtein=(colSums(consensus_Protein*consensus_Protein_WT_flip)/length(protein.align))*100)
      Output <- c(Output,mutagenesis_Protein)
    }

    #AID mutagenesis
    if (Config.List$measure.shm==1){
      AID_Mutagenesis <- lapply(AID.Targets[[x]],function(y){
        cyt_mut <- mutagenesis_DNA[y$Cytosine]
        motif_mut <- sapply(1:nrow(y),function(z){sum(mutagenesis_DNA[y$Start[z]:y$End[z]])})
        motif_mut_nonc <- sapply(1:nrow(y),function(z){
          i <- y$Start[z]:y$End[z]
          j <- i[!y$Start[z]:y$End[z] %in% y$Cytosine[z]]
          sum(mutagenesis_DNA[j])
        })
        a <- data.frame(MotifMutagenesis=motif_mut,CytosineMutagenesis=cyt_mut)
        Positions <- cbind(y,a)

        motif_pos <- unique(as.numeric(sapply(1:nrow(y),function(z){
          y$Start[z]:y$End[z]
        })))
        cyt_pos <- unique(y$Cytosine)
        nonc_pos <- which((!1:Biostrings::nchar(WT_DNA) %in% cyt_pos)==TRUE)
        motif_pos_nonc <- motif_pos[!motif_pos %in% cyt_pos]

        total_mut <- sum(mutagenesis_DNA)
        total_motif_mut <- sum(mutagenesis_DNA[motif_pos])
        motif_mut_avg <- mean(motif_mut)
        motif_mut_freq <- sum(mutagenesis_DNA[motif_pos])/total_mut*100
        motif_mut_nonc_avg <- mean(motif_mut_nonc)
        motif_mut_nonc_freq <- sum(mutagenesis_DNA[motif_pos_nonc])/total_mut*100
        cyt_mut_avg <- mean(cyt_mut)
        nonc_mut_avg <- mean(mutagenesis_DNA[nonc_pos])
        cyt_mut_freq <- (sum(cyt_mut)/total_mut)*100
        cyt_mut_freq_of_motif <- (sum(cyt_mut)/sum(motif_mut))*100

        Motif_Sums <- data.frame(TotalPercentMutated=total_mut/(length(mutagenesis_DNA)*100)*100,
                                 TotalPercentMutatedInMotif=total_motif_mut/(length(mutagenesis_DNA)*100)*100,
                                 TotalPercentMutatedInCytosine=sum(cyt_mut)/(length(mutagenesis_DNA)*100)*100,
                                 MotifMutationAverage=motif_mut_avg,
                                 MotifMutationFrequencyofTotal=motif_mut_freq,
                                 MotifMutationNonCAverage=motif_mut_nonc_avg,
                                 MotifMutationNonCFrequencyOfTotal=motif_mut_nonc_freq,
                                 CytosineMutationAverage=cyt_mut_avg,
                                 NonCytosineMutationAverage=nonc_mut_avg,
                                 CytosineMutationFrequencyOfTotal=cyt_mut_freq,
                                 CytosineMutationFrequencyOfMotif=cyt_mut_freq_of_motif)

        Mut_nonc <- mutagenesis_DNA[nonc_pos]
        Mut_c <- Positions$CytosineMutagenesis

        list(Positions=Positions,MotifSums=Motif_Sums,MutCyt=Mut_c,MutNonCyt=Mut_nonc)
      })
      Output <- c(Output,list(AID=AID_Mutagenesis))
    }

    #mutagenesis by region
    Ab_regions <- c("FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4")

    if ((tolower(Input.DataFrame[x,"Antibody"])=="yes")==TRUE){
      Ab_coord <- as.numeric(Input.DataFrame[x,str_c(Ab_regions,rep("Start",length(Ab_regions)))])
      Ab_coord <- c(Ab_coord,Biostrings::nchar(WT_DNA)+1) #END1 is length +1 to mark the "start" of the "end" segment so that it can be used similar to all the other start fragments such that END1-1=the length of the antibody
      names(Ab_coord) <- c(Ab_regions,"END1")
      #check if coordinates supplied are reasonable
      if (TRUE %in% is.na(Ab_coord) | FALSE %in% (sort(as.matrix(Ab_coord),na.last=TRUE)==as.matrix(Ab_coord)) | TRUE %in% (Ab_coord>Ab_coord[["END1"]])){
        print(str_c("Sample ",x))
        print(str_c("Regions: ",names(Input.DataFrame[x,str_c(Ab_regions,rep("Start",length(Ab_regions)))])[is.na(Input.DataFrame[x,str_c(Ab_regions,rep("Start",length(Ab_regions)))])]))
        stop("The above sample was labeled as an antibody, but the listed coordinates were not provided or valid")
      }
      #DNA mutations by region
      Region_Mutations_DNA <- list()
      Region_Mutations_DNA[Ab_regions] <- lapply(seq_along(Ab_regions),function(y){
        total_mut <- sum(mutagenesis_DNA)
        region_start <- Ab_coord[[y]]
        region_end <- Ab_coord[[y+1]]-1
        if (region_start>=region_end){
          print(str_c("Sample ",x))
          print(str_c("Region: ",Ab_regions[y]))
          stop("Invalid coordinates for the above region.")
        }
        v <- mutagenesis_DNA[region_start:region_end]
        f <- (sum(v)/total_mut)*100
        m <- mean(v)
        e <- length(region_start:region_end)/(Ab_coord[["END1"]]-1)*100
        n <- f-e
        logo <- ggseqlogo(consensus_DNA[,region_start:region_end],method='prob') +
          theme_logo(base_size=9) +
          theme(axis.text.x = element_blank())
        list(Mut_Pos=v,Mut_Freq=f,Mut_Avg=m,Mut_Expected=e,Mut_Norm=n,Logo=logo)
      })

      Output <- c(Output,list(RegionCoordinates=Ab_coord,RegionMutagenesisDNA=Region_Mutations_DNA))

      #protein mutations by region
      if (Config.List$protein.mutations==1){
        Region_Mutations_Protein <- list()
        Region_Mutations_Protein[Ab_regions] <- lapply(seq_along(Ab_regions),function(y){
          region_start_prot <- ceiling(Ab_coord[[y]]/3)
          region_end_prot <- floor((Ab_coord[[y+1]]-1)/3)
          if (region_start_prot!=round(region_start_prot) | region_end_prot!=round(region_end_prot)){
            print(str_c("Sample ",x))
            print(str_c("Region: ",Ab_regions[y]))
            stop("Coordinate error for above region: not a multiple of 3.")
          }
          total_mut_prot <- sum(mutagenesis_Protein[[1]])
          vp <- mutagenesis_Protein[[1]][region_start_prot:region_end_prot]
          fp <- (sum(vp)/total_mut_prot)*100
          mp <- mean(vp)
          ep <- length(region_end_prot:region_start_prot)/floor((Ab_coord[["END1"]]-1)/3)*100
          np <- fp-ep

          # ggseqlogo doesn't work on regions < length 2
          if(length(region_start_prot:region_end_prot) > 1){
            logo <- ggseqlogo(consensus_Protein[,region_start_prot:region_end_prot], method='prob') +
              theme_logo(base_size=9) +
              theme(axis.text.x = element_blank()) +
              theme(legend.position="none")
          }else{
            logo <- "Could not generate sequence logo plot"
          }
          list(Mut_Pos=vp, Mut_Freq=fp, Mut_Avg=mp, Mut_Expected=ep, Mut_Norm=np, Logo=logo)
        })
        Output <- c(Output,list(RegionMutagenesisProtein=Region_Mutations_Protein))
      }
    }

    return(Output)
  }) #single threaded due to potential for high memory usage

  #measure mutation types as likely DNA repair pathways
  if (Config.List$dna.repair.pathways==1){
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
      if (abs(sum(sums) - 100) > 1e-7){
        stop(str_c("Frequency summation error detected in indel types for sample: "),x)
      }
      return(sums)
    }))
    colnames(mut.types) <- Sample.Names
    rownames(mut.types) <- mut_label
    mut_types_wb <- createWorkbook("Mutation Types.xlsx")
    addWorksheet(mut_types_wb, "Sheet1")
    writeData(mut_types_wb, "Sheet1", mut.types, rowNames = TRUE)
  }

  #summary mutation data compilation

  if (Config.List$measure.shm==1){
    maxlen <- max(Biostrings::nchar(Input.DataFrame$ReferenceSequence))
    Mut_pos_all <- data.frame(matrix(nrow=maxlen,ncol=length(Sample.Names)))
    colnames(Mut_pos_all) <- Sample.Names

    lencyt <- numeric()
    lennonc <- numeric()
    for(i in Sample.Names){
      lencyt <- c(lencyt,nrow(Output.Mutagenesis[[i]]$AID$WRCH$Positions))
      lennonc <- c(lennonc,length(Output.Mutagenesis[[i]]$AID$WRCH$MutNonCyt))
    }
    maxlencyt <- max(lencyt)
    Mut_pos_cyt <- data.frame(matrix(nrow=maxlencyt,ncol=length(Sample.Names)))
    colnames(Mut_pos_cyt) <- Sample.Names
    maxlennonc <- max(lennonc)
    Mut_pos_nonc <- data.frame(matrix(nrow=maxlennonc,ncol=length(Sample.Names)))
    colnames(Mut_pos_nonc) <- Sample.Names

    Motif_Sums <- data.frame(matrix(nrow=ncol(Output.Mutagenesis[[1]]$AID$WRCH$MotifSums),ncol=length(Sample.Names)))
    colnames(Motif_Sums) <- Sample.Names
    rownames(Motif_Sums) <- colnames(Output.Mutagenesis[[1]]$AID$WRCH$MotifSums)

    for(i in Sample.Names){
      a <- Output.Mutagenesis[[i]]$MutagenesisDNA
      b <- c(a,rep(NA,maxlen-length(a)))
      Mut_pos_all[i] <- b

      d <- Output.Mutagenesis[[i]]$AID$WRCH$MutCyt
      e <- c(d,rep(NA,maxlencyt-length(d)))
      Mut_pos_cyt[i] <- e

      f <- Output.Mutagenesis[[i]]$AID$WRCH$MutNonCyt
      g <- c(f,rep(NA,maxlennonc-length(f)))
      Mut_pos_nonc[i] <- g

      Motif_Sums[i] <- as.numeric(Output.Mutagenesis[[i]]$AID$WRCH$MotifSums)
    }
  }else{
    maxlen <- max(Biostrings::nchar(Input.DataFrame$ReferenceSequence))
    Mut_pos_all <- data.frame(matrix(nrow=maxlen,ncol=length(Sample.Names)))
    colnames(Mut_pos_all) <- Sample.Names
    for(i in Sample.Names){
      a <- Output.Mutagenesis[[i]]$MutagenesisDNA
      b <- c(a,rep(NA,maxlen-length(a)))
      Mut_pos_all[i] <- b
    }
  }

  Mut_pos_wb <- createWorkbook("Mut_pos.xlsx")
  addWorksheet(Mut_pos_wb,"MutAll")
  writeData(Mut_pos_wb,"MutAll",Mut_pos_all)
  if (Config.List$measure.shm==1){
    addWorksheet(Mut_pos_wb,"MutCyt")
    writeData(Mut_pos_wb,"MutCyt",Mut_pos_cyt)
    addWorksheet(Mut_pos_wb,"MutNonC")
    writeData(Mut_pos_wb,"MutNonC",Mut_pos_nonc)
    addWorksheet(Mut_pos_wb,"MotifSums")
    writeData(Mut_pos_wb,"MotifSums",Motif_Sums,rowNames = TRUE)
    WRCY_wb <- createWorkbook("WRCY.xlsx")
    WRCH_wb <- createWorkbook("WRCH.xlsx")
  }

  Sequences_wb <- createWorkbook("Sequences.xlsx")

  for(i in Sample.Names){
    #output sequence tables
    addWorksheet(Sequences_wb,i)
    writeData(Sequences_wb,i,Sequence.Table.List[[i]])

    #output AID tables
    if (Config.List$measure.shm==1){
      addWorksheet(WRCY_wb,i)
      writeData(WRCY_wb,i,Output.Mutagenesis[[i]]$AID$WRCY$Positions)
      addWorksheet(WRCH_wb,i)
      writeData(WRCH_wb,i,Output.Mutagenesis[[i]]$AID$WRCH$Positions)
    }
  }



  Output.Mutagenesis$Workbooks <- list()
  Output.Mutagenesis$Workbooks <- list(Mut_pos_wb = Mut_pos_wb, Sequences_wb = Sequences_wb)
  if (Config.List$measure.shm==1){
    Output.Mutagenesis$Workbooks <- c(Output.Mutagenesis$Workbooks, list(WRCY_wb = WRCY_wb, WRCH_wb = WRCH_wb))
  }
  if (Config.List$dna.repair.pathways==1){
    Output.Mutagenesis$Workbooks <- c(Output.Mutagenesis$Workbooks, list(mut_types_wb = mut_types_wb))
  }
  return(Output.Mutagenesis)
}


