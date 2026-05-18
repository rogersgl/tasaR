#' Analyze a pair of sequences
#'
#' @description
#' Analyzes a matched pair of control and experimental samples for mutations in the latter.
#' Does not require a reference sequence, and can automatically use the control sequence to
#' determine reference sequences, including if the DNA is heterozygous. Sequences will be aligned
#' and labeled based on their closest reference sequence.
#'
#' Parameters can be specified by either a tas.object.settings object as the settings parameter,
#' or the required settings can be inputted into the parameters of the function.
#'
#' @param file.control A character vector specifying the file path to the control merged fastq file.
#' @param file.experimental A character vector specifying the file path to the experimental merged fastq file(s).
#' @param settings (Optional) A tas.object.settings object specifying settings for the sequences. File path will be overridden by the function variables above.
#' @param fwd.primer (Optional) A character vector of the DNA sequence of the forward primer.
#' @param rev.primer (Optional) A character vector of the DNA sequence of the reverse primer.
#' @param fwd.exten.type (Optional) A character vector of the type of forward extension. Can be "Barcode", "UMI", or "".
#' @param fwd.exten.seq (Optional) A character vector of the DNA sequence of the forward extension.
#' @param rev.exten.type (Optional) A character vector of the type of reverse extension. Can be "Barcode", "UMI", or "".
#' @param rev.exten.seq (Optional) A character vector of the DNA sequence of the reverse extension.
#' @param diploid (Optional) A logical vector specifying whether the sequence is diploid or not. Default is TRUE, affects whether a warning message is displayed if more than 2 potential reference sequences are detected in the control sample.
#' @param with.nuclease A logical vector specifying whether this is a nuclease experiment to analyze DNA scarring. Modifies workflow of basic analysis pipeline to perform gap correction around nuclease cut site.
#' @param gRNA.seq DNA sequence of Cas9 gRNA 20 bp protospacer. Aligned to reference sequence to determine nuclease cut site.
#' @param manual.cut.site A numeric vector to specify the nuclease cut site manually in the reference sequence. Currently only supports a single cut site.
#' @param ... (Optional) Any additional parameters to pass to child functions
#'
#' @returns An S4 object of class PairedAmpliconSequencing (if length(file.experimental) == 1).
#' A list of AmpliconSequencing objects if length(file.experimental) > 1.
#' @export
#'
#' @examples
#' \dontrun{
#'   pairedAnalyzeAmplicon(file.control, file.experimental, settings = tas.object.settings)
#' }
pairedAnalyzeAmplicon <- function(file.control,
                                  file.experimental,
                                  settings = NULL,
                                  fwd.primer = "",
                                  rev.primer = "",
                                  fwd.exten.type = "",
                                  fwd.exten.seq = "",
                                  rev.exten.type = "",
                                  rev.exten.seq = "",
                                  diploid = TRUE,
                                  with.nuclease = TRUE,
                                  gRNA.seq = NULL,
                                  manual.cut.site = NULL,
                                  ...) {
  if (is.null(settings) && !all(nzchar(c(fwd.primer, rev.primer, fwd.exten.type, fwd.exten.seq, rev.exten.type, rev.exten.seq)))) {
    stop("Must provide settings as either a tas.object.settings object or via the required variables in this function. See 'readSettings' help file.")
  } else if (is(settings, "tas.object.settings")) {
    config <- settings
  } else if (is(settings, "data.frame") || is(settings, "list")) {
    config <- readSettings(settings)
  } else if (is.null(settings) && any(nzchar(c(fwd.primer, rev.primer, fwd.exten.type, fwd.exten.seq, rev.exten.type, rev.exten.seq)))){
    config <- methods::new("tas.object.settings",
                           ForwardPrimer = fwd.primer,
                           ReversePrimer = rev.primer,
                           ForwardExtensionType = fwd.exten.type,
                           ForwardExtension = fwd.exten.seq,
                           ReverseExtensionType = rev.exten.type,
                           ReverseExtension = rev.exten.seq
                           )
  }

  if (!file.exists(file.control)) {stop("Control fastq file not found.")}
  if (!file.exists(file.experimental)) {stop("Experimental fastq file not found.")}

  config.ctrl <- config
  config.ctrl@MergedFASTQPath <- file.control
  if (is.na(config.ctrl@InsertStart)) {
    config.ctrl@InsertStart <- as.integer(1 + nchar(stringr::str_c(config.ctrl@ForwardExtension, config.ctrl@ForwardPrimer)))
  }
  if (is.na(config.ctrl@InsertEnd)) {
    config.ctrl@InsertEnd <- as.integer(-1 - nchar(stringr::str_c(config.ctrl@ReverseExtension, config.ctrl@ReversePrimer)))
  }

  seq.ctrl <- buildSequenceTable(config.ctrl, paired.analysis.ctrl = TRUE, diploid, with.nuclease = with.nuclease, gRNA.seq = gRNA.seq, manual.cut.site = manual.cut.site)
  config.ctrl@ReferenceSequence <- lapply(seq.ctrl@Alignments$DNA, function(x) {
    as.character(pwalign::unaligned(pwalign::subject(x)))
  })
  config.ctrl@AmpliconLength <- as.integer(nchar(config.ctrl@ReferenceSequence[[1]]) + config.ctrl@InsertStart + abs(config.ctrl@InsertEnd) - 2)
  config.ctrl@Name <- "Control"

  as.ctrl <- new("AmpliconSequencing", Sequences = seq.ctrl,
                  Mutations = measureMutationsNuclease(seq.ctrl, config.ctrl),
                  MutationTypes = measureDNArepair(seq.ctrl, config.ctrl),
                  Settings = config.ctrl)

  config.expt <- config.ctrl

  as.expt.list <- list(length(file.experimental))
  for (i in 1:length(file.experimental)) {
    config.expt@MergedFASTQPath <- file.experimental[i]
    config.expt@Name <- ifelse(!is.null(names(file.experimental)),
                               names(file.experimental)[i],
                               stringr::str_c("Experimental_", as.character(i)))
    seq.expt <- buildSequenceTable(config.expt, paired.analysis.ctrl = FALSE, with.nuclease = with.nuclease, gRNA.seq = gRNA.seq, manual.cut.site = manual.cut.site)
    as.expt.list[[i]] <- new("AmpliconSequencing", Sequences = seq.expt,
                             Mutations = measureMutationsNuclease(seq.expt, config.expt),
                             MutationTypes = measureDNArepair(seq.expt, config.expt),
                             Settings = config.expt)
    names(as.expt.list)[i] <- getSettings(as.expt.list[[i]])$Name
  }

  if (length(file.experimental) == 1) {
    pas <- methods::new("PairedAmpliconSequencing", Control = as.ctrl, Experimental = as.expt.list[[1]])
    return(pas)
  } else {
    pas.list <- c(list(Control = as.ctrl), as.expt.list)
    return(pas.list)
  }
}





