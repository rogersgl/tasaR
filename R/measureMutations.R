###############################################################################
# measureMutations()
# ------------------

measureMutations <- function(sequence.table, settings) {
    config <- getSettings(settings)
    t <- getSequenceTable(sequence.table)


    AID.targets <- findAIDtargets(config)


    cat("Measuring DNA mutations...\n")

    wt.dna <- DNAString(config$ReferenceSequence)
    seq.len <- Biostrings::nchar(wt.dna)
    seqs <- DNAStringSet(t$Sequences)
    # pw <- pairwiseAlignment(seqs, wt.dna)
    # dna.align <- rep(pw, t$Count)
    dna.align <- rep(getDNAalign(sequence.table), t$Count)


    ### All Mutations ###

    mm.df <- mismatchSummary(dna.align)$subject
    dt.mm <- data.table(Position = mm.df$SubjectPosition, N = mm.df$Count)

    dt.i <- data.table(Position = (start(unlist(insertion(dna.align))) + 0.5))
    dt.i <- dt.i[, .N, by = Position]
    setorder(dt.i, -N)

    dt.d <- data.table(start = start(unlist(deletion(dna.align))), width = width(unlist(deletion(dna.align))))
    dt.d <- dt.d[, .(start = start, end = (start + width -1))]
    dt.d <- dt.d[, .(Position = unlist(Map(seq, start, end)))]
    dt.d <- dt.d[, .N, by = Position]

    dt.merge <- rbindlist(list(dt.mm, dt.i, dt.d))[, .(MutationFrequency = sum(N)/sum(t$Count)*100), by = Position]
    idx.nz <- dt.merge$Position
    idx.z <- (1:seq.len)[!(1:seq.len %in% idx.nz)]
    dt.merge <- rbind(dt.merge, data.table(Position = idx.z, MutationFrequency = rep(0, length(idx.z))))
    setorder(dt.merge, Position)

    AllMutations <- as.data.frame(dt.merge)


    ### (Non)Cytosine Mutations ###

    pos.c <- AID.targets$WRCH$Cytosine
    pos.nonc <- (dt.merge$Position)[!(dt.merge$Position %in% pos.c)]

    CytosineMutations <- as.data.frame(dt.merge[(dt.merge$Position %in% pos.c),])
    NonCytosineMutations <- as.data.frame(dt.merge[(dt.merge$Position %in% pos.nonc),])


    ### Motif Sums ###
    AverageAllMutations <-  mean(AllMutations$MutationFrequency)
    AverageCytosineMutations <-  mean(CytosineMutations$MutationFrequency)
    AverageNonCytosineMutations <-  mean(NonCytosineMutations$MutationFrequency)
    FrequencyOfAllMutationsAtCytosines <- sum(CytosineMutations$MutationFrequency)/sum(AllMutations$MutationFrequency)*100
    FrequencyOfAllMutationsAtNonCytosines <- sum(NonCytosineMutations$MutationFrequency)/sum(AllMutations$MutationFrequency)*100
    ms <- setNames(c(AverageAllMutations, AverageCytosineMutations, AverageNonCytosineMutations, FrequencyOfAllMutationsAtCytosines, FrequencyOfAllMutationsAtNonCytosines),
                   c("AverageAllMutations", "AverageCytosineMutations", "AverageNonCytosineMutations", "FrequencyOfAllMutationsAtCytosines", "FrequencyOfAllMutationsAtNonCytosines"))

    ### AID Tables ###
    AID.tables <- list(WRCH = cbind(AID.targets$WRCH, data.frame(CytosineMutationFrequency = dt.merge$MutationFrequency[dt.merge$Position %in% AID.targets$WRCH$Cytosine])),
                       WRCY = cbind(AID.targets$WRCY, data.frame(CytosineMutationFrequency = dt.merge$MutationFrequency[dt.merge$Position %in% AID.targets$WRCY$Cytosine])))
    row.names(AID.tables$WRCH) <- NULL
    row.names(AID.tables$WRCY) <- NULL





    ### Antibody Regions ###
    if (config$IsAntibody) {
      region.coords <- getAntibodyCoordinates(config$AntibodyRegions)
      total.mut <- sum(dt.merge$MutationFrequency)
      region.mut.dev <- lapply(region.coords, function(x) {
        region.mut <- dt.merge[(dt.merge$Position %in% x),"MutationFrequency"]
        (sum(region.mut)/total.mut-length(x)/seq.len)/(length(x)/seq.len)*100 # percent deviation from expected mutation rate if distribution was even across the sequence
      })
    }

    ### Protein Mutations ###

    wt.prot <- suppressWarnings(AAStringSet(translate(wt.dna)))
    prot.len <- Biostrings::nchar(wt.prot)

    prot.align <- rep(getAAalign(sequence.table), t$Count)
    ids <- getSequenceSupplemental(sequence.table)$IndelStart
    idt <- getSequenceSupplemental(sequence.table)$IndelType

    dt.p <- as.data.table(mismatchSummary(prot.align)$subject)[!Pattern %in% c("-", "+")]
    dt.p <- dt.p[, .(.N, MutationFrequency = sum(Count)/sum(t$Count)*100, Position = SubjectPosition), by = SubjectPosition][,c("Position", "MutationFrequency")]
    dt.idp <- data.table(Position = ceiling(ids/3), MutationFrequency = t$Percent)[!is.na(Position)]
    dt.pl <- rbindlist(list(dt.p, dt.idp,
                           data.table(Position = which(!1:prot.len %in% c(dt.p$Position, dt.idp$Position)), MutationFrequency = rep(0, length(which(!1:prot.len %in% c(dt.p$Position, dt.idp$Position))))))
                           , fill = TRUE)[, .(MutationFrequency = sum(MutationFrequency)), by = Position]

    setorder(dt.pl, Position)

    all.mut.prot <- as.data.frame(dt.pl)

    cm.wt <- consensusMatrix(wt.prot)
    cm.prot <- consensusMatrix(AAStringSet(getSequencesAA(sequence.table)))
    cm.prot <- cm.prot/(sum(t$Count))*100
    cm.prot[cm.wt==1] <- 0

    new("tas.mutations", DNA = list(AllMutations = AllMutations,
                                    CytosineMutations = CytosineMutations,
                                    NonCytosineMutations = NonCytosineMutations,
                                    MotifSums = ms),
                         AA = list(AllMutations = all.mut.prot, # data.frame
                                   MutationMatrix = cm.prot), # matrix
                         AIDTables = AID.tables) # list(WRCH, WRCY)

}



###############################################################################
# Helper functions
# ------

findAIDtargets <- function(settings) {
  if (class(settings) == "tas.object.settings") {
    settings <- getSettings(settings)
  }
  if (class(settings) != "list") {
    stop("settings must be a tas.object.settings object or a list")
  }

  AID.Motifs <- c("WRCY","RGYW","WRCH","DGYW")

  aid <- list()
  aid[AID.Motifs] <- lapply(AID.Motifs, function(x){
    pat.views <- matchPattern(x, DNAString(settings$ReferenceSequence), fixed = FALSE)
    s <- start(pat.views)
    e <- end(pat.views)

    if (x=="WRCY" | x=="WRCH"){
      cyt <- e-1
    }
    if (x=="RGYW" | x=="DGYW"){
      cyt <- e-2
    }
    data.frame(Motif=rep(x,length(pat.views)), Start=s, End=e, Cytosine=cyt)
  })

  y.df <- rbind(aid[["WRCY"]], aid[["RGYW"]])
  y.df <- y.df[order(y.df$Start),]
  h.df <- rbind(aid[["WRCH"]], aid[["DGYW"]])
  h.df <- h.df[order(h.df$Start),]

  list(WRCH = h.df, WRCY = y.df)
}


getAntibodyCoordinates <- function(abr) {
  out <- list(length(abr))
  for (i in seq_along(abr)[-length(abr)]) {
    out[[i]] <- abr[i]:(abr[i+1] - 1)
  }
  names(out) <- str_remove(names(abr), "Start")[-length(abr)]
  return(out)
}


