#' @include utils-color-palettes.R
#'
# helpers

gg.tasar.defaults <- list(
  y_scale = list(expand = ggplot2::expansion(mult = c(0, .05))),
  x_scale = list(expand = ggplot2::expansion(mult = c(0, 0))),
  theme = list(text = ggplot2::element_text(size = 8),
               axis.line = ggplot2::element_line(linewidth = 0.5/.pt),
               axis.ticks = ggplot2::element_line(linewidth = 0.5/.pt),
               axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)),
               axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
               legend.title = ggplot2::element_blank(),
               legend.box.spacing = ggplot2::unit(-0.1, "in"),
               legend.key.height = ggplot2::unit(0.1, "in"),
               legend.key.width = ggplot2::unit(0.2, "in"),
               legend.key.spacing.y = ggplot2::unit(0.05, "in"))
)



label.cdrs <- function(plot, abr) {
  if (!is(plot, "ggplot")) {
    stop("Object 'plot' must be a ggplot object.")
  }
  p <- plot + ggplot2::annotate("rect", xmin = c(abr["CDR1Start"], abr["CDR2Start"], abr["CDR3Start"]),
                       xmax = c(abr["FR2Start"]-1, abr["FR3Start"]-1, abr["FR4Start"]-1),
                       ymin = c(0, 0, 0),
                       ymax = c(Inf, Inf, Inf),
                       color = "gray90",
                       fill = "gray90") +
    ggplot2::annotate("text", x = c(mean(c(abr["CDR1Start"], abr["FR2Start"] - 1)), mean(c(abr["CDR2Start"], abr["FR3Start"] - 1)), mean(c(abr["CDR3Start"], abr["FR4Start"] - 1))),
             y = Inf, label = c("CDR1", "CDR2", "CDR3"), size = 7/.pt, vjust = -1) +
    ggplot2::coord_cartesian(clip = "off")
}



# graph dna mut.pos

gg.dna.pos <- function(input) {

  sett <- getSettings(input)

  p <- ggplot2::ggplot(data = getMutationDistributionDNA(input), ggplot2::aes(x = Position, y = MutationFrequency))

  if (sett$IsAntibody) {
    abr <- sett$AntibodyRegions
    p <- label.cdrs(p, abr)
  }
  p <- p + ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::scale_x_continuous(expand =  gg.tasar.defaults$x_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
                 axis.line = gg.tasar.defaults$theme$axis.line,
                 axis.ticks = gg.tasar.defaults$theme$axis.ticks,
                 axis.title.y = gg.tasar.defaults$theme$axis.title.y,
                 axis.title.x = gg.tasar.defaults$theme$axis.title.x) +
    ggplot2::ylab("Mutation (%)") +
    ggplot2::xlab("Position (nt)") +
    ggh4x::force_panelsizes(rows = ggplot2::unit(1.5, "in"), cols = ggplot2::unit(4, "in"))
  return(p)
}


# graph dna mut.pos w/ cytosine labels

gg.dna.pos.cyt <- function(input) {

  sett <- getSettings(input)

  cyt <- getMutationDistributionCytosine(input)
  cyt <- cbind(cyt, data.frame(Type = rep("AID Cytosines", nrow(cyt))))
  ncyt <- getMutationDistributionNonCytosine(input)
  ncyt <- cbind(ncyt, data.frame(Type = rep("Other Bases", nrow(ncyt))))
  df <- rbind(cyt, ncyt)

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = Position, y = MutationFrequency, fill = Type))
  if (sett$IsAntibody) {
    abr <- sett$AntibodyRegions
    p <- label.cdrs(p, abr)
  }
  p <- p +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::scale_x_continuous(expand =  gg.tasar.defaults$x_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
             axis.line = gg.tasar.defaults$theme$axis.line,
             axis.ticks = gg.tasar.defaults$theme$axis.ticks,
             axis.title.y = gg.tasar.defaults$theme$axis.title.y,
             axis.title.x = gg.tasar.defaults$theme$axis.title.x,
             legend.title = gg.tasar.defaults$theme$legend.title,
             legend.box.spacing = gg.tasar.defaults$theme$legend.box.spacing,
             legend.key.height = gg.tasar.defaults$theme$legend.key.height,
             legend.key.width = gg.tasar.defaults$theme$legend.key.width,
             legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    ggplot2::ylab("Mutation (%)") +
    ggplot2::xlab("Position (nt)") +
    ggplot2::scale_fill_manual(values = c("red", "black")) +
       ggh4x::force_panelsizes(rows = ggplot2::unit(1.5, "in"), cols = ggplot2::unit(4, "in"))
  return(p)
}


# graph box plot of %mutations at cytosine and non-cytosines

gg.aid.box <- function(input) {

  cyt <- getMutationDistributionCytosine(input)
  cyt <- cbind(cyt, data.frame(Type = rep("AID Cytosines", nrow(cyt))))
  ncyt <- getMutationDistributionNonCytosine(input)
  ncyt <- cbind(ncyt, data.frame(Type = rep("Other Bases", nrow(ncyt))))
  df <- rbind(cyt, ncyt)

  ggplot2::ggplot(data = df, ggplot2::aes(x = Type, y = MutationFrequency)) +
    ggplot2::geom_boxplot(width = 0.5, outlier.size = 1) +
    ggplot2::labs(x = "Base Type",y = "% Mutation") +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x) +
    ggh4x::force_panelsizes(rows = ggplot2::unit(1.5, "in"), cols = ggplot2::unit(1.5, "in"))
}


# graph aa mut.pos

gg.aa.pos <- function(input) {
  sett <- getSettings(input)
  p <- ggplot2::ggplot(data = getMutationDistributionAA(input), ggplot2::aes(x = Position, y = MutationFrequency))
  if (sett$IsAntibody) {
    abr <- ceiling(sett$AntibodyRegions/3)
    p <- label.cdrs(p, abr)
  }
  p <- p +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::scale_x_continuous(expand =  gg.tasar.defaults$x_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x) +
    ggplot2::ylab("Mutation (%)") +
    ggplot2::xlab("Position (AA)") +
    ggh4x::force_panelsizes(rows = ggplot2::unit(1.5, "in"), cols = ggplot2::unit(4, "in"))
  return(p)
}


gg.aa.pos.cyt <- function(input) {
  sett <- getSettings(input)
  df <- getMutationDistributionAA(input)
  df <- cbind(df, data.frame(Type = rep("Other Bases", nrow(df))))
  df$Type[df$Position %in% unique(ceiling(getMutationDistributionCytosine(input)$Position/3))] <- "AID Cytosines"

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = Position, y = MutationFrequency, fill = Type))

  if (sett$IsAntibody) {
    abr <- ceiling(sett$AntibodyRegions/3)
    p <- label.cdrs(p, abr)
  }
  p <- p +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::scale_x_continuous(expand =  gg.tasar.defaults$x_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.box.spacing = gg.tasar.defaults$theme$legend.box.spacing,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    ggplot2::ylab("Mutation (%)") +
    ggplot2::xlab("Position (AA)") +
    ggplot2::scale_fill_manual(values = c("red", "black")) +
    ggh4x::force_panelsizes(rows = ggplot2::unit(1.5, "in"), cols = ggplot2::unit(4, "in"))
  return(p)
}




# graph aa stacked bar mutation graph

gg.aa.muts.stacked <- function(input) {
  mm <- getMutationMatrixAA(input)[AAvec,]
  df <- data.frame(Position = c(sapply(1:ncol(mm), rep, nrow(mm))),
                   AA = rep(AAvec, ncol(mm)),
                   Frequency = c(mm)
                   )

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = Position, y = Frequency, fill = AA))
  p <- p +
    ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE), stat = "identity", width = 1) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::scale_x_continuous(expand =  gg.tasar.defaults$x_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.box.spacing = gg.tasar.defaults$theme$legend.box.spacing,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    ggplot2::ylab("% Mutated Residues") +
    ggplot2::xlab("Position (AA)") +
    ggplot2::scale_fill_manual(values = aa_letter3) +
    ggh4x::force_panelsizes(rows = ggplot2::unit(2, "in"), cols = ggplot2::unit(5, "in"))
  return(p)
}


gg.dna.mut.count.hist <- function(input) {

  stl <- getSequenceTable(input)

  idx.id <- which(nzchar(stl$Indels))
  idx.wt <- which(stl$Indels == "WT")
  idx.id <- idx.id[idx.id != idx.wt]

  nid <- stringr::str_count(stl$Indels[idx.id], ",") + 1

  nmut <- stl$BasesChanged
  nmut[idx.id] <- nmut[idx.id] + nid

  p <- ggplot2::ggplot(data = data.frame(Mutations = nmut), ggplot2::aes(x = Mutations)) +
    ggplot2::geom_histogram(binwidth = 1) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::scale_x_continuous(expand =  gg.tasar.defaults$x_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
              axis.line = gg.tasar.defaults$theme$axis.line,
              axis.ticks = gg.tasar.defaults$theme$axis.ticks,
              axis.title.y = gg.tasar.defaults$theme$axis.title.y,
              axis.title.x = gg.tasar.defaults$theme$axis.title.x) +
    ggplot2::ylab("# of sequences") +
    ggplot2::xlab("# of mutations") +
    ggh4x::force_panelsizes(rows = ggplot2::unit(1, "in"), cols = ggplot2::unit(1.5, "in"))

  return(p)
}





# graph mutation types stacked bar graph

gg.dna.mut.types <- function(input, type = "bar") {

  mt <- getMutationTypes(input)

  mt.df <- data.frame(Names = c("WT", "NHEJ", "MMEJ", "Bases Changed", "Indel + Bases Changed", "Other"),
                      Frequency = as.numeric(mt),
                      Sample = rep("", nrow(mt)))

  mt.df$Names <- factor(mt.df$Names, unique(mt.df$Names))

  # freq.ymax <- cumsum(mt.df$Frequency)
  # freq.ymin <- c(0, head(freq.ymax, n=-1))
  # freq.ypos <- (freq.ymax + freq.ymin) / 2
  # geom_text(x = 1, y = freq.ypos, label = colnames(mt))
  # text overlaps too much to be useful with this

  p <- ggplot2::ggplot(data = mt.df, ggplot2::aes(x = 2, y = Frequency, fill = Names)) +
    ggplot2::geom_bar(position = ggplot2::position_stack(reverse = TRUE), stat = "identity", width = 1.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
                axis.line = gg.tasar.defaults$theme$axis.line,
                axis.ticks.y = gg.tasar.defaults$theme$axis.ticks,
                axis.title.y = gg.tasar.defaults$theme$axis.title.y,
                legend.title = gg.tasar.defaults$theme$legend.title,
                axis.title.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank()) +
    ggplot2::scale_fill_brewer(type = "qual", palette = "Set1") +
    ggplot2::ylab("Frequency (%)")

  if (type == "bar") {
    p <- p +
      ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
      ggplot2::theme(legend.box.spacing = ggplot2::unit(0, "in"),
            legend.key.height = gg.tasar.defaults$theme$legend.key.height,
            legend.key.width = gg.tasar.defaults$theme$legend.key.width,
            legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
      ggplot2::scale_x_continuous(limits = c(1, 3)) +
      ggh4x::force_panelsizes(rows = ggplot2::unit(1, "in"), cols = ggplot2::unit(0.3, "in"))
  }

  if (type == "donut") {
    p <- p +
      ggplot2::theme_void() +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
            legend.box.spacing = ggplot2::unit(0, "in"),
            legend.key.height = gg.tasar.defaults$theme$legend.key.height,
            legend.key.width = gg.tasar.defaults$theme$legend.key.width,
            legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::xlim(0.1, 3.5) +
      ggh4x::force_panelsizes(rows = ggplot2::unit(1.5, "in"), cols = ggplot2::unit(1.5, "in"))
  }

  if (type == "pie") {
    p <- p +
      ggplot2::theme_void() +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
            legend.box.spacing = ggplot2::unit(0, "in"),
            legend.key.height = gg.tasar.defaults$theme$legend.key.height,
            legend.key.width = gg.tasar.defaults$theme$legend.key.width,
            legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
      ggplot2::coord_polar(theta = "y") +
      ggh4x::force_panelsizes(rows = ggplot2::unit(1.5, "in"), cols = ggplot2::unit(1.5, "in"))
  }
  return(p)
}


# ggmsa alignments - depreciated, replaced by TeX-based workflow

# gg.dna.align <- function(input, seq.number = 10) {
#
#   sett <- getSettings(input)
#   len <- length(getSequencesDNA(input))
#   if (len < seq.number) {seq.number <- len}
#   stl <- getSequenceTable(input)[(1:seq.number),]
#
#   indel.label <- stl$Indels
#   indel.idx <- which(nzchar(indel.label))
#
#   n <- stringr::str_c(ifelse(nzchar(stl$Indels), stl$Indels, ""),
#              ifelse(stl$BasesChanged != 0,
#                 ifelse(nzchar(stl$Indels), ", ", "") |>
#                     stringr::str_c(stl$BasesChanged, " SNV"), ""),
#              stringr::str_c(" - ", round(stl$Percent, 2), "%")
#             )
#   n <- c("Reference", n)
#   seqs <- Biostrings::DNAStringSet(c(sett$ReferenceSequence[[1]], # a bit of a crude bypass, just take the 1st ref seq. If >1 will still have WT marking on other alleles. However, choice of ref #1 is not guaranteed consistent.
#     stringr::str_remove_all(stl$Sequences, "[+\\-]")))
#
#   names(seqs) <- n
#   sink(tempfile())
#   align <- msa::msaClustalW(seqs, order = "input")
#   sink()
#
#   p <- suppressMessages(
#     ggmsa::ggmsa(Biostrings::DNAMultipleAlignment(as(align,"BStringSet")),
#           consensus_views = TRUE,
#           ref = "Reference",
#           color = "Taylor_NT",
#           char_width = 0.7,
#           border = NA,
#           seq_name = TRUE)+
#       ggplot2::coord_cartesian()+
#       ggmsa::facet_msa(field = 100)+
#       ggplot2::theme(plot.margin = ggplot2::margin(0.5,0.5,0.5,0.5,unit = "in"))+
#       ggplot2::theme(axis.text = ggplot2::element_text(size = 6))
#   )
#
#   return(p)
# }
#
#
# # aa alignment
#
# gg.aa.align <- function(input, seq.number = 10) {
#
#   sett <- getSettings(input)
#   len <- length(getSequencesAA(input))
#   if (len < seq.number) {seq.number <- len}
#   stl <- getSequenceTable(input)[(1:seq.number),]
#
#   n <- stringr::str_c(1:seq.number, ". ", stl$ProteinMutation,
#              stringr::str_c("- ", round(stl$Percent, 2), "%")
#   )
#   n <- c("Reference", n)
#   ref <- as.character(suppressWarnings(Biostrings::translate(Biostrings::DNAString(sett$ReferenceSequence[[1]])))) # a bit of a crude bypass, just take the 1st ref seq. If >1 will still have WT marking on other alleles. However, choice of ref #1 is not guaranteed consistent.
#   seqs <- Biostrings::AAStringSet(stringr::str_remove_all(c(ref, stl$AA), "[+\\-]"))
#
#   # names(seqs) <- as.character(seq_along(seqs))
#   names(seqs) <- n
#
#   sink(tempfile())
#   align <- msa::msaClustalW(seqs, order = "input", type = "protein")
#   sink()
#
#   p <- suppressMessages(
#     ggmsa::ggmsa(Biostrings::AAMultipleAlignment(as(align,"BStringSet")),
#           consensus_views = TRUE,
#           ref = "Reference",
#           color = "Taylor_NT",
#           char_width = 0.7,
#           border = NA,
#           seq_name = TRUE)+
#       ggplot2::coord_cartesian()+
#       ggmsa::facet_msa(field = 100)+
#       ggplot2::theme(plot.margin = ggplot2::margin(0.5,0.5,0.5,0.5,unit = "in"))+
#       ggplot2::theme(axis.text = ggplot2::element_text(size = 6))
#   )
#
#   return(p)
# }

# graph histogram of # of mutations/seq
