




#' Multi-sample analysis
#'
#' @description
#' Function providing for batch analysis of samples. Sample settings must be
#' provided, either as a .csv file (1 sample per row), or as an R object (a list
#' or data.frame) containing the same information. See \code{\link{readSettings}} for
#' more information.
#'
#'
#' @param input.settings.list A character vector describing a path to a .csv file to read settings from. Alternatively, an R data.frame or list containing those settings.
#'
#' @returns A list of AmpliconSequencing objects
#' @export
#'
#' @examples
#' batchAnalyseAmplicon('/.../folder/settings.csv')
#' batchAnalyzeAmplicon(data.frame)
#' batchAnalyzeAmplicon(list)
batchAnalyzeAmplicon <- function(input.settings.list) {
  settings.list <- list()

  # file path to .csv
  if (class(input.settings.list) == "character") {

    if (!str_ends(input.settings.list, ".csv")) {stop("A character vector input.settings.list must specify the file path to a .csv file.")}
    if (!file.exists(input.settings.list)) {stop("File not found: ", input.settings.list)}

    input.df <- read.csv(input.settings.list)

    settings.list[input.df$Name] <- lapply(1:nrow(input.df), function(x) {
      temp <- readSettings(input.df[x,])
      validObject(temp)
      return(temp)
    })
  }

  # R list of settings
  if (class(input.settings.list) == "list") {
    settings.list <- lapply(1:length(input.settings.list), function(x) {
      if (class(x) == "tas.object.settings"){
        validObject(x)
        return(x)
      } else {
          temp <- readSettings(x)
          validObject(temp)
          return(temp)
        }
      for (i in 1:length(settings.list)) {
        names(settings.list[i]) <- getSettings(settings.list[i])$Name
      }
    })
  }

  sample.names <- names(settings.list)
  results.list <- list()
  results.list[sample.names] <- lapply(settings.list, function(x) {
    analyzeAmplicon(x)
  })
  return(results.list)
}



#' Summary analyses for batch processing
#'
#' @description
#' Performs comparison analyses of the results from \code{\link{batchAnalyzeAmlicon}}. Results are returned
#' as an R list, and can be saved as a variable. Also includes the option to directly export the results
#' to a folder specified by the user. Treats all samples equally, so is unaware of pairings of control/
#' experimental samples to direct analysis.
#'
#'
#' @param results.list A list of AmpliconSequencing objects from \code{\link{batchAnalyzeAmplicon}}
#' @param export A logical specifying whether the resulting tables and graphs should also be written to disk
#' @param path (Optional) If \code{export = TRUE}, the file path to save the exported information to
#' @param suppressConsoleOutput Logical vector specifying whether to print results to console or not, if applicable
#'
#' @returns A list of results, with sub-lists for tables and graphs
#' @export
#'
#' @examples
#' batchSummarize(results.list)
#' batchSummarize(results.list, export = TRUE, path = './save_in_folder/')
batchSummarize <- function(results.list, export = FALSE, path = NULL, suppressConsoleOutput = FALSE) {

  # extract and collate data from list of S4 objects
  sample.names <- names(results.list)

  dt.funcs <- list("getMutationDistributionDNA",
                   "getMutationDistributionCytosine",
                   "getMutationDistributionNonCytosine",
                   "getMutationDistributionAA"
                   )

  merged.pos.list <- lapply(dt.funcs, function(func){
    dt <- as.data.table(Reduce(function(x, y) merge(x, y, by = "Position", all = TRUE), lapply(results.list, func)))
    dt[is.na(dt)] <- 0
    colnames(dt)[colnames(dt) != "Position"] <- sample.names
    return(dt)
  })
  names(merged.pos.list) <- c("AllMutationsDNA", "CytosineMutations", "NonCytosineMutations", "AllMutationsAA")

  dt.ms <- as.data.table(Reduce(function(x, y) merge(x, y, by = "rn", all = TRUE), lapply(results.list, function(x){
    data.table(getMutationMotifSums(x), keep.rownames = TRUE)
    })))
  colnames(dt.ms)[colnames(dt.ms) != "rn"] <- sample.names

  dt.mt <- as.data.table(Reduce(function(x, y) merge(x, y, by = "rn", all = TRUE), lapply(results.list, function(x){
    data.table(t(getMutationTypes(x)), keep.rownames = TRUE)
  })))
  colnames(dt.mt)[colnames(dt.mt) != "rn"] <- sample.names

  # compare average mut freq at all bases
  p.ams <- gg.summary.all.ms(sample.names, dt.ms)

  # compare average mut freq at AID cytosines
  p.cms <- gg.summary.cyt.ms(sample.names, dt.ms)

  # compare all dna mutations by pos
  p.amd <- gg.dna.pos.facet(sample.names, merged.pos.list)

  # compare all aa mutations by pos
  p.amp <- gg.aa.pos.facet(sample.names, merged.pos.list)

  # compare all cytosine mutations by pos as heatmap
  p.cmd <- gg.dna.cyt.hm(sample.names, merged.pos.list)

  # compare aid cytosines vs other bases boxplots
  p.aid <- gg.aid.box.batch(sample.names, merged.pos.list)

  # compare dna repair types stacked bar
  p.mt <- gg.dna.mut.types.batch(sample.names, merged.pos.list)


  batsum <- list(Graphs = list(AvgMutAll = p.ams,
                            AvgMutCyt = p.cms,
                            MutPosDNA = p.amd,
                            MutPosAA = p.amp,
                            MutCytHM = p.cmd,
                            AIDBox = p.aid,
                            MutTypes = p.mt),
              Tables = c(list(MotifSums = as.data.frame(dt.ms)),
                              lapply(merged.pos.list, as.data.frame),
                         list(MutationTypes = as.data.frame(dt.mt))
                         )
              )

  if (export && !is.null(path)) {
    if (!dir.exists(file.path(tempdir(), "export"))) {dir.create(file.path(tempdir(), "export"))}
    write.csv(batsum$Tables$MotifSums, file = file.path(tempdir(), "export", "Motif Sums (batch).csv"))
    write.csv(batsum$Tables$AllMutationsDNA, file = file.path(tempdir(), "export", "All DNA Mutations (batch).csv"))
    write.csv(batsum$Tables$CytosineMutations, file = file.path(tempdir(), "export", "Cytosine DNA Mutations (batch).csv"))
    write.csv(batsum$Tables$NonCytosineMutations, file = file.path(tempdir(), "export", "Non-Cytosine DNA Mutations (batch).csv"))
    write.csv(batsum$Tables$AllMutationsAA, file = file.path(tempdir(), "export", "All Protein Mutations (batch).csv"))
    write.csv(batsum$Tables$MutationTypes, file = file.path(tempdir(), "export", "DNA Repair Types (batch).csv"))

    ggsave("Average DNA Mutation Rate.pdf", batsum$Graphs$AvgMutAll, path = file.path(tempdir(), "export"))
    ggsave("Average AID Cytosine Mutation Rate.pdf", batsum$Graphs$AvgMutCyt, path = file.path(tempdir(), "export"))
    ggsave("DNA Mutation Distributions.pdf", batsum$Graphs$MutPosDNA, path = file.path(tempdir(), "export"))
    ggsave("Protein Mutation Distributions.pdf", batsum$Graphs$MutPosAA, path = file.path(tempdir(), "export"))
    ggsave("AID Cytosine Mutation Heatmap.pdf", batsum$Graphs$MutCytHM, path = file.path(tempdir(), "export"))
    ggsave("AID Mutation Boxplot.pdf", batsum$Graphs$AIDBox, path = file.path(tempdir(), "export"))
    ggsave("DNA Repair Types.pdf", batsum$Graphs$MutTypes, path = file.path(tempdir(), "export"))

    if (!dir.exists(path)) {dir.create(path)}
    file.copy(file.path(tempdir(), "export"), path, recursive = TRUE)
  } else if (export && is.null(path)) {
    warning("Could not export results, path not provided.")
  }

  if (suppressConsoleOutput) {
    return(invisible(batsum))
  } else {
    return(batsum)
  }
}


