


gg.summary.all.ms <- function(sample.names, dt.ms) {
  df.ams <- data.frame(Sample = sample.names,
                       mut = BiocGenerics::unlist(dt.ms["AverageAllMutations", -1], use.names = FALSE),
                       type = rep("All Bases", length(sample.names)))
  ggplot2::ggplot(data = df.ams, ggplot2::aes(x = Sample, y = mut)) +
    ggplot2::geom_bar(stat = "identity", color = "black", fill = "black", width = 0.75) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    ggplot2::ylab("% Mutation\n (All Bases)") +
    ggplot2::xlab("Sample") +
    ggh4x::force_panelsizes(rows = ggplot2::unit(1, "in"), cols = ggplot2::unit(0.3 * length(sample.names), "in"))
}


gg.summary.cyt.ms <- function(sample.names, dt.ms) {
  df.cms <- data.frame(Sample = sample.names,
                       mut = BiocGenerics::unlist(dt.ms["AverageCytosineMutations", -1], use.names = FALSE),
                       type = rep("AID Cytosines", length(sample.names)))

  ggplot2::ggplot(data = df.cms, ggplot2::aes(x = Sample, y = mut, fill = type)) +
    ggplot2::geom_bar(stat = "identity", color = "black", fill = "black", width = 0.75) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    ggplot2::ylab("% Mutation \n(AID Cytosines)") +
    ggplot2::xlab("Sample") +
    ggh4x::force_panelsizes(rows = ggplot2::unit(1, "in"), cols = ggplot2::unit(0.3 * length(sample.names), "in"))
}



gg.dna.pos.facet <- function(sample.names, merged.pos.list) {

  df.amd <- as.data.frame(rbindlist(S4Vectors::lapply(sample.names, function(x) {
    data.table(pos = BiocGenerics::unlist(merged.pos.list$AllMutationsDNA[,"Position"]),
               freq = merged.pos.list$AllMutationsDNA[, ..x],
               sample = rep(x, nrow(merged.pos.list$AllMutationsDNA))
    )
  }), use.names = FALSE))
  colnames(df.amd) <- c("Position", "freq", "sample")

  ggplot2::ggplot(data = df.amd, ggplot2::aes(x = Position, y = freq, fill = sample)) +
    ggplot2::geom_col(position = "identity", color = "black", fill = "black") +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::scale_x_continuous(expand = gg.tasar.defaults$x_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          panel.border = ggplot2::element_rect(colour = "black", fill = NA, linetype = "solid"),
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    ggplot2::facet_grid(scales = "fixed",
               rows = ggplot2::vars(sample), switch = "y") +
    ggplot2::theme(strip.text.y.left = ggplot2::element_text(size = 8),
          strip.background = ggh4x::element_part_rect(color = "black", linewidth = 0.25, side = "tlbr"),
          strip.placement = "outside") +
    ggplot2::ylab("Mutation Frequency (%)") +
    ggh4x::force_panelsizes(rows = ggplot2::unit(0.75 * length(sample.names), "in"), cols = 0.05 * nrow(df.amd)/length(sample.names))
}


gg.aa.pos.facet <- function(sample.names, merged.pos.list) {

  df.amp <- as.data.frame(rbindlist(S4Vectors::lapply(sample.names, function(x) {
    data.table(pos = BiocGenerics::unlist(merged.pos.list$AllMutationsAA[,"Position"]),
               freq = merged.pos.list$AllMutationsAA[, ..x],
               sample = rep(x, nrow(merged.pos.list$AllMutationsAA))
    )
  }), use.names = FALSE))
  colnames(df.amp) <- c("Position", "freq", "sample")

  ggplot2::ggplot(data = df.amp, ggplot2::aes(x = Position, y = freq, fill = sample)) +
    ggplot2::geom_col(position = "identity", color = "black", fill = "black") +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::scale_x_continuous(expand = gg.tasar.defaults$x_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          panel.border = ggplot2::element_rect(colour = "black", fill = NA, linetype = "solid"),
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    ggplot2::facet_grid(scales = "fixed",
               rows = ggplot2::vars(sample), switch = "y") +
    ggplot2::theme(strip.text.y.left = ggplot2::element_text(size = 8),
          strip.background = ggh4x::element_part_rect(color = "black", linewidth = 0.25, side = "tlbr"),
          strip.placement = "outside") +
    ggplot2::ylab("Mutation Frequency (%)") +
    ggh4x::force_panelsizes(rows = ggplot2::unit(0.75 * length(sample.names), "in"), cols = 0.05 * nrow(df.amp)/length(sample.names))
}




gg.dna.cyt.hm <- function(sample.names, merged.pos.list) {

  df.cmd <- as.data.frame(rbindlist(S4Vectors::lapply(sample.names, function(x) {
    data.table(cyt = BiocGenerics::unlist(merged.pos.list$CytosineMutations[,"Position"]),
               freq = merged.pos.list$CytosineMutations[, ..x],
               sample = rep(x, nrow(merged.pos.list$CytosineMutations)),
               type = rep("AID Cytosine", nrow(merged.pos.list$CytosineMutations))
    )
  }), use.names = FALSE))
  colnames(df.cmd) <- c("Position", "freq", "sample", "type")
  df.cmd$Position <- factor(df.cmd$Position, levels = unique(df.cmd$Position))

  ggplot2::ggplot(data = df.cmd, ggplot2::aes(x = Position, y = sample, fill = freq)) +
    ggplot2::geom_tile(color = "NA", height = 0.9) +
    ggplot2::scale_fill_distiller(palette = "YlGnBu", labels = function(fill) paste0(fill, "%"))+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0, linetype = "solid",
                                   colour = "black"),
          axis.text.x = ggplot2::element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major = ggplot2::element_blank(),
          legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0), limits = rev) +
    ggplot2::ylab("Sample") +
    ggh4x::force_panelsizes(rows = ggplot2::unit(0.6 * length(sample.names), "in"), cols = 0.12 * nrow(df.cmd)/length(sample.names))
}


gg.aid.box.batch <- function(sample.names, merged.pos.list) {

  df.cmd <- as.data.frame(rbindlist(S4Vectors::lapply(sample.names, function(x) {
    data.table(cyt = BiocGenerics::unlist(merged.pos.list$CytosineMutations[,"Position"]),
               freq = merged.pos.list$CytosineMutations[, ..x],
               sample = rep(x, nrow(merged.pos.list$CytosineMutations)),
               type = rep("AID Cytosine", nrow(merged.pos.list$CytosineMutations))
    )
  }), use.names = FALSE))
  colnames(df.cmd) <- c("Position", "freq", "sample", "type")

  df.ncmd <- as.data.frame(rbindlist(S4Vectors::lapply(sample.names, function(x) {
    data.table(ncyt = BiocGenerics::unlist(merged.pos.list$NonCytosineMutations[,"Position"]),
               freq = merged.pos.list$NonCytosineMutations[, ..x],
               sample = rep(x, nrow(merged.pos.list$NonCytosineMutations)),
               type = rep("Other Bases", nrow(merged.pos.list$NonCytosineMutations))
    )
  }), use.names = FALSE))
  colnames(df.ncmd) <- c("Position", "freq", "sample", "type")

  df.aid <- rbind(df.ncmd, df.cmd)
  df.aid$type <- factor(df.aid$type, levels = c("Other Bases", "AID Cytosine"))

  ggplot2::ggplot(data = df.aid, ggplot2::aes(x = sample, y = freq, fill = type)) +
    ggplot2::geom_boxplot(color = "black", outlier.shape = 21, outlier.size = 1, linewidth = 0.25, outlier.stroke = 0.25) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggplot2::scale_x_discrete(expand = gg.tasar.defaults$x_scale$expand) +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    ggplot2::labs(x = "Sample", y = "Mutation Frequency (%)") +
    ggplot2::scale_fill_manual(values = c("white", "red")) +
    ggh4x::force_panelsizes(rows = ggplot2::unit(1, "in"), cols = ggplot2::unit(0.4 * length(sample.names), "in"))
}


gg.dna.mut.types.batch <- function(sample.names, merged.pos.list, dt.mt){

  df.mt <- as.data.frame(rbindlist(S4Vectors::lapply(sample.names, function(x) {
    data.table(rn = factor(dt.mt$rn,
                           levels = c("WT", "NHEJ", "MMEJ", "BaseChange", "IndelBaseChange", "Other")),
               freq = dt.mt[, ..x],
               sample = rep(x, nrow(dt.mt))
    )
  }), use.names = FALSE))
  colnames(df.mt) <- c("type", "freq", "sample")

  ggplot2::ggplot(data = df.mt, ggplot2::aes(x = sample, y = freq, fill = type)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_stack(reverse = TRUE), width = 0.8, color = "black", linewidth = 0.2) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks.y = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    ggplot2::scale_fill_brewer(type = "qual", palette = "Set1") +
    ggplot2::ylab("Frequency (%)") +
    ggplot2::xlab("Sample") +
    ggplot2::scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    ggh4x::force_panelsizes(rows = ggplot2::unit(1, "in"), cols = ggplot2::unit(0.3 * length(sample.names), "in"))
}

