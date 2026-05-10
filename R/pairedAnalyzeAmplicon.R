


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
#' @param file.control A character vector specifying the file path to the control merged fastq file
#' @param file.experimental A character vector specifying the file path to the experimental merged fastq file
#' @param settings (Optional) A tas.object.settings object specifying settings for the sequences. File path will be overridden by the function variables above.
#' @param fwd.primer (Optional) A character vector of the DNA sequence of the forward primer
#' @param rev.primer (Optional) A character vector of the DNA sequence of the reverse primer
#' @param fwd.exten.type (Optional) A character vector of the type of forward extension. Can be "Barcode", "UMI", or ""
#' @param fwd.exten.seq (Optional) A character vector of the DNA sequence of the forward extension
#' @param rev.exten.type (Optional) A character vector of the type of reverse extension. Can be "Barcode", "UMI", or ""
#' @param rev.exten.seq (Optional) A character vector of the DNA sequence of the reverse extension
#' @param diploid (Optional) A logical vector specifying whether the sequence is diploid or not. Default is TRUE, affects whether a warning message is displayed if more than 2 potential reference sequences are detected in the control sample.
#' @param ... (Optional) Any additional parameters to pass to child functions
#'
#' @returns An S4 object of class PairedAmpliconSequencing
#' @export
#'
#' @examples
#' pairedAnalyzeAmplicon(file.control, file.experimental, settings = tas.object.settings)
pairedAnalyzeAmplicon <- function(file.control,
                                  file.experimental,
                                  settings = NULL,
                                  fwd.primer = "",
                                  rev.primer = "",
                                  fwd.exten.type = "",
                                  fwd.exten.seq = "",
                                  rev.exten.type = "",
                                  rev.exten.seq = "",
                                  diploid = TRUE, ...) {
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

  seq.ctrl <- buildSequenceTable(config.ctrl, paired.analysis.ctrl = TRUE, diploid)
  config.ctrl@ReferenceSequence <- lapply(seq.ctrl@Alignments$DNA, function(x) {
    as.character(pwalign::unaligned(pwalign::subject(x)))
  })
  config.ctrl@AmpliconLength <- as.integer(nchar(config.ctrl@ReferenceSequence[[1]]) + config.ctrl@InsertStart + abs(config.ctrl@InsertEnd) - 2)
  config.ctrl@Name <- "Control"

  as.ctrl <- new("AmpliconSequencing", Sequences = seq.ctrl,
                  Mutations = measureMutations(seq.ctrl, config.ctrl),
                  MutationTypes = measureDNArepair(seq.ctrl, config.ctrl),
                  Settings = config.ctrl)

  config.expt <- config.ctrl
  config.expt@MergedFASTQPath <- file.experimental
  config.expt@Name <- "Experimental"

  as.expt <- analyzeAmplicon(config.expt)

  methods::new("PairedAmpliconSequencing", Control = as.ctrl, Experimental = as.expt)
}


exportNucleaseAnalysis <- function(paired.seq.results,
                                   nuclease.type = "Cas9",
                                   gRNA.seq = NULL,
                                   manual.cut.sites = NULL) {
  if (!is(paired.seq.results, "PairedAmpliconSequencing")){stop("Parameter paired.seq.results must be a PairedAmpliconSequencing object.")}
  match.arg(nuclease.type, choices = c("Cas9", "Cas12a", "Custom"))
  if (is.null(gRNA.seq) && is.null(manual.cut.sites)){stop("Cut site(s) must be specified either by providing the gRNA sequence or manually as a numeric vector. See help at ?exportNucleaseAnalysis.")}

  ref.seq <- lapply(getSettings(paired.seq.results@Control)$ReferenceSequence, Biostrings::DNAStringSet)
  if (!is.null(gRNA.seq) && is.null(manual.cut.sites)){
    manual.cut.sites <- .findCutSitesFromSequence(gRNA.seq, ref.seq)
  }

  # TODO: finish this function
  # have cut sites, now what? brain fried.
}


