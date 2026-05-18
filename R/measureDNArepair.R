

#' Estimate DNA repair pathway usage
#'
#' @description
#' Estimates the frequency of DNA repair pathway usage in a given sequence table of DNA sequences. Classification as follows:
#'
#' _WT_: unmodified wild-type reference sequences
#'
#' _NHEJ_: non-homologous end joining, insertions and deletions of -1 and -2
#'
#' _MMEJ_: microhomology-mediated end joining, deletions larger than -2
#'
#' _BaseChange_: base mutations without an indel scar
#'
#' _IndelBaseChange_: sequences with both an indel and base changes
#'
#' _Other_: other sequences (containing multiple indels, etc.)
#'
#' @param sequence.table S4 object of class tas.sequences
#' @param settings S4 object of class tas.object.settings
#'
#' @returns S4 object of class tas.dna.repair with estimates of DNA repair types
#' @export
#'
#' @examples
#' \dontrun{
#'   measureDNArepair(sequence.table, settings)
#' }
measureDNArepair <- function(sequence.table, settings) {
  config <- getSettings(settings)
  t <- getSequenceTable(sequence.table)

  idx.other <- BiocGenerics::grep("[,]", t$Indels) # multiple indels
  idx.ibc <- which((nzchar(t$Indels) + (t$BasesChanged != 0)) == 2)[!which((nzchar(t$Indels) + (t$BasesChanged != 0)) == 2) %in% idx.other]
  idx.nhej <- BiocGenerics::grep("\\+|-(?:1|2)\\b", t$Indels)[!BiocGenerics::grep("\\+|-(?:1|2)\\b", t$Indels) %in% c(idx.other, idx.ibc)]
  idx.mmej <- BiocGenerics::grep("[-]", t$Indels)[!BiocGenerics::grep("[-]", t$Indels) %in% c(idx.other, idx.ibc, idx.nhej)]
  idx.mm <- which(((!nzchar(t$Indels)) + (t$BasesChanged != 0) == 2))[!(which(((!nzchar(t$Indels)) + (t$BasesChanged != 0) == 2)) %in% c(idx.other, idx.ibc))]
  idx.wt <- which(t$Indels == "WT")

  idx.all <- c(idx.other, idx.ibc, idx.nhej, idx.mmej, idx.mm, idx.wt)
  if (length(idx.all) != length(unique(idx.all))) {stop("Duplicate indexes detected.")}
  if (length(idx.all) != nrow(t)) {stop("Invalid number of indexes")}
  if (any(!1:nrow(t) %in% idx.all)) {stop("Missing indexes detected: ", which(!1:nrow(t) %in% idx.all))}
  wt <- sum(t$Percent[idx.wt])
  nhej <- sum(t$Percent[idx.nhej])
  mmej <- sum(t$Percent[idx.mmej])
  mm <- sum(t$Percent[idx.mm])
  ibc <- sum(t$Percent[idx.ibc])
  other <- sum(t$Percent[idx.other])

  # TODO: add heterozygous SNP handling
  # TODO: additional support for Cas9/nuclease stuff


  new("tas.dna.repair", WT = wt,
                        NHEJ = nhej,
                        MMEJ = mmej,
                        BaseChange = mm,
                        IndelBaseChange = ibc,
                        Other = other)

}
