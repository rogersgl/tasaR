

gg.summary.all.ms <- function(sample.names, dt.ms) {
  df.ams <- data.frame(Sample = sample.names,
                       mut = unlist(dt.ms["AverageAllMutations", -1], use.names = FALSE),
                       type = rep("All Bases", length(sample.names)))
  ggplot(data = df.ams, aes(x = Sample, y = mut)) +
    geom_bar(stat = "identity", color = "black", fill = "black", width = 0.75) +
    theme_classic() +
    scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    ylab("% Mutation\n (All Bases)") +
    xlab("Sample") +
    force_panelsizes(rows = unit(1, "in"), cols = unit(0.3 * length(sample.names), "in"))
}

gg.summary.cyt.ms <- function(sample.names, dt.ms) {
  df.cms <- data.frame(Sample = sample.names,
                       mut = unlist(dt.ms["AverageCytosineMutations", -1], use.names = FALSE),
                       type = rep("AID Cytosines", length(sample.names)))

  ggplot(data = df.cms, aes(x = Sample, y = mut, fill = type)) +
    geom_bar(stat = "identity", color = "black", fill = "black", width = 0.75) +
    theme_classic() +
    scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    ylab("% Mutation \n(AID Cytosines)") +
    xlab("Sample") +
    force_panelsizes(rows = unit(1, "in"), cols = unit(0.3 * length(sample.names), "in"))
}



gg.dna.pos.facet <- function(sample.names, merged.pos.list) {

  df.amd <- as.data.frame(rbindlist(lapply(sample.names, function(x) {
    data.table(pos = unlist(merged.pos.list$AllMutationsDNA[,"Position"]),
               freq = merged.pos.list$AllMutationsDNA[, ..x],
               sample = rep(x, nrow(merged.pos.list$AllMutationsDNA))
    )
  }), use.names = FALSE))
  colnames(df.amd) <- c("Position", "freq", "sample")

  ggplot(data = df.amd, aes(x = Position, y = freq, fill = sample)) +
    geom_col(position = "identity", color = "black", fill = "black") +
    theme_classic() +
    scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    scale_x_continuous(expand = gg.tasar.defaults$x_scale$expand) +
    theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          panel.border = element_rect(colour = "black", fill = NA, linetype = "solid"),
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    facet_grid(scales = "fixed",
               rows = vars(sample), switch = "y") +
    theme(strip.text.y.left = element_text(size = 8),
          strip.background = element_part_rect(color = "black", linewidth = 0.25, side = "tlbr"),
          strip.placement = "outside") +
    ylab("Mutation Frequency (%)") +
    force_panelsizes(rows = unit(0.75 * length(sample.names), "in"), cols = 0.05 * nrow(df.amd)/length(sample.names))
}


gg.aa.pos.facet <- function(sample.names, merged.pos.list) {

  df.amp <- as.data.frame(rbindlist(lapply(sample.names, function(x) {
    data.table(pos = unlist(merged.pos.list$AllMutationsAA[,"Position"]),
               freq = merged.pos.list$AllMutationsAA[, ..x],
               sample = rep(x, nrow(merged.pos.list$AllMutationsAA))
    )
  }), use.names = FALSE))
  colnames(df.amp) <- c("Position", "freq", "sample")

  ggplot(data = df.amp, aes(x = Position, y = freq, fill = sample)) +
    geom_col(position = "identity", color = "black", fill = "black") +
    theme_classic() +
    scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    scale_x_continuous(expand = gg.tasar.defaults$x_scale$expand) +
    theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          panel.border = element_rect(colour = "black", fill = NA, linetype = "solid"),
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    facet_grid(scales = "fixed",
               rows = vars(sample), switch = "y") +
    theme(strip.text.y.left = element_text(size = 8),
          strip.background = element_part_rect(color = "black", linewidth = 0.25, side = "tlbr"),
          strip.placement = "outside") +
    ylab("Mutation Frequency (%)") +
    force_panelsizes(rows = unit(0.75 * length(sample.names), "in"), cols = 0.05 * nrow(df.amp)/length(sample.names))
}




gg.dna.cyt.hm <- function(sample.names, merged.pos.list) {

  df.cmd <- as.data.frame(rbindlist(lapply(sample.names, function(x) {
    data.table(cyt = unlist(merged.pos.list$CytosineMutations[,"Position"]),
               freq = merged.pos.list$CytosineMutations[, ..x],
               sample = rep(x, nrow(merged.pos.list$CytosineMutations)),
               type = rep("AID Cytosine", nrow(merged.pos.list$CytosineMutations))
    )
  }), use.names = FALSE))
  colnames(df.cmd) <- c("Position", "freq", "sample", "type")
  df.cmd$Position <- factor(df.cmd$Position, levels = unique(df.cmd$Position))

  ggplot(data = df.cmd, aes(x = Position, y = sample, fill = freq)) +
    geom_tile(color = "NA", height = 0.9) +
    scale_fill_distiller(palette = "YlGnBu", labels = function(fill) paste0(fill, "%"))+
    theme_classic()+
    theme(axis.line = element_line(linewidth = 0, linetype = "solid",
                                   colour = "black"),
          axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major = element_blank(),
          legend.title = element_blank()
    ) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0), limits = rev) +
    ylab("Sample") +
    force_panelsizes(rows = unit(0.6 * length(sample.names), "in"), cols = 0.12 * nrow(df.cmd)/length(sample.names))
}


gg.aid.box.batch <- function(sample.names, merged.pos.list) {

  df.cmd <- as.data.frame(rbindlist(lapply(sample.names, function(x) {
    data.table(cyt = unlist(merged.pos.list$CytosineMutations[,"Position"]),
               freq = merged.pos.list$CytosineMutations[, ..x],
               sample = rep(x, nrow(merged.pos.list$CytosineMutations)),
               type = rep("AID Cytosine", nrow(merged.pos.list$CytosineMutations))
    )
  }), use.names = FALSE))
  colnames(df.cmd) <- c("Position", "freq", "sample", "type")

  df.ncmd <- as.data.frame(rbindlist(lapply(sample.names, function(x) {
    data.table(ncyt = unlist(merged.pos.list$NonCytosineMutations[,"Position"]),
               freq = merged.pos.list$NonCytosineMutations[, ..x],
               sample = rep(x, nrow(merged.pos.list$NonCytosineMutations)),
               type = rep("Other Bases", nrow(merged.pos.list$NonCytosineMutations))
    )
  }), use.names = FALSE))
  colnames(df.ncmd) <- c("Position", "freq", "sample", "type")

  df.aid <- rbind(df.ncmd, df.cmd)
  df.aid$type <- factor(df.aid$type, levels = c("Other Bases", "AID Cytosine"))

  ggplot(data = df.aid, aes(x = sample, y = freq, fill = type)) +
    geom_boxplot(color = "black", outlier.shape = 21, outlier.size = 1, linewidth = 0.25, outlier.stroke = 0.25) +
    theme_classic() +
    scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    scale_x_discrete(expand = gg.tasar.defaults$x_scale$expand) +
    theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          axis.title.x = gg.tasar.defaults$theme$axis.title.x,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    labs(x = "Sample", y = "Mutation Frequency (%)") +
    scale_fill_manual(values = c("white", "red")) +
    force_panelsizes(rows = unit(1, "in"), cols = unit(0.4 * length(sample.names), "in"))
}


gg.dna.mut.types.batch <- function(sample.names, merged.pos.list){

  df.mt <- as.data.frame(rbindlist(lapply(sample.names, function(x) {
    data.table(rn = factor(dt.mt$rn,
                           levels = c("WT", "NHEJ", "MMEJ", "BaseChange", "IndelBaseChange", "Other")),
               freq = dt.mt[, ..x],
               sample = rep(x, nrow(dt.mt))
    )
  }), use.names = FALSE))
  colnames(df.mt) <- c("type", "freq", "sample")

  ggplot(data = df.mt, aes(x = sample, y = freq, fill = type)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.8, color = "black", linewidth = 0.2) +
    theme_classic() +
    theme(text = gg.tasar.defaults$theme$text,
          axis.line = gg.tasar.defaults$theme$axis.line,
          axis.ticks.y = gg.tasar.defaults$theme$axis.ticks,
          axis.title.y = gg.tasar.defaults$theme$axis.title.y,
          legend.title = gg.tasar.defaults$theme$legend.title,
          legend.key.height = gg.tasar.defaults$theme$legend.key.height,
          legend.key.width = gg.tasar.defaults$theme$legend.key.width,
          legend.key.spacing.y = gg.tasar.defaults$theme$legend.key.spacing.y) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    ylab("Frequency (%)") +
    xlab("Sample") +
    scale_y_continuous(expand = gg.tasar.defaults$y_scale$expand) +
    force_panelsizes(rows = unit(1, "in"), cols = unit(0.3 * length(sample.names), "in"))
}

