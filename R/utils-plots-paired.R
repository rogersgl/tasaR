
gg.dna.pos.paired <- function(ctrl.aln, expt.aln, gRNA.window) {

  ctrl <- getMutationDistributionDNA(paired.sequencing.results@Control)[gRNA.window,]
  expt <- getMutationDistributionDNA(paired.sequencing.results@Experimental)[gRNA.window,]

  df <- rbind(cbind(expt, data.frame(Sample = rep("expt", nrow(expt)))),
              cbind(ctrl, data.frame(Sample = rep("ctrl", nrow(ctrl))))
              )

  # distribute insertions with 0.5 position indexes to whole indexes on either side
  # needed because ggplot2 stacked
  ins.idx <- which((df$Position %% 1) != 0)
  df.ins <- df[ins.idx,]
  df$MutationFrequency[ins.idx - 1] <- df$MutationFrequency[ins.idx - 1] + df.ins$MutationFrequency/2
  df$MutationFrequency[ins.idx + 1] <- df$MutationFrequency[ins.idx + 1] + df.ins$MutationFrequency/2
  df <- df[-ins.idx,]

  p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = Position, y = MutationFrequency, fill = Sample)) +
    ggplot2::geom_col(stat = "identity", width = 1, position = "identity") +
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
    ggplot2::scale_fill_manual(labels = c("Control", "Experimental"), values = c("gray70", "red")) +
    ggplot2::geom_vline(xintercept = median(gRNA.window + 0.5),
                        linetype = "dashed",
                        linewidth = 0.4) +
    ggh4x::force_panelsizes(rows = ggplot2::unit(1.5, "in"), cols = ggplot2::unit(length(gRNA.window)/25, "in"))
}
