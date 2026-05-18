# TODO: implement helper function .isTinytexOk to check before export
#       establish export functions
#       build tests

#' @include utils-color-palettes.R
NULL

#' Write a TeXshade file from a multiple sequence alignment
#'
#' Generates a standalone LaTeX (\code{.tex}) document containing an embedded
#' FASTA alignment (via the \code{filecontents*} environment) and a single
#' \code{texshade} block. The output highlights:
#' \itemize{
#'   \item the reference sequence with residue-specific background colors
#'   \item mismatched residues in other sequences with residue-specific colors
#'   \item matching residues in other sequences with a configurable text color
#'     on a white background
#'   \item gap characters rendered consistently as black text on a white
#'     background
#' }
#'
#' The function supports both nucleotide and amino acid multiple alignments.
#' Set \code{seq_type = "N"} for DNA/RNA-style alignments or \code{seq_type = "P"}
#' for protein alignments.
#'
#' Residue selections are computed in gap-free TeXshade residue coordinates, and
#' shading regions are constructed so they do not span alignment gaps.
#'
#' The resulting document is intended for compilation with pdfLaTeX (for example
#' via TinyTeX).
#'
#' @param aln A multiple alignment object inheriting from \code{MultipleAlignment}
#'   (for example, from \code{msa::msaClustalW()}).
#' @param outfile Character. Path to the output \code{.tex} file.
#' @param reference Reference sequence, either a numeric index or a sequence name.
#'   This sequence is placed first in the exported alignment and fully color-shaded
#'   by residue.
#' @param color_scheme Character. Built-in palette selector. Supported values
#'   are \code{"auto"} (default). For DNA sequences, \code{"nt_default"},
#'   \code{"nt_okabe_ito"}, \code{"nt_shapely"}, and \code{"nt_alt"}. For protein
#'  sequences, \code{"aa_chemistry"}, \code{"aa_shapely"}, \code{"aa_letter1"},
#'  \code{"aa_letter2"}, and \code{"aa_letter3"}.
#'   With \code{"auto"}, the palette is chosen from \code{seq_type}.
#' @param color_map Optional named character vector of residue-to-color mappings.
#'   If supplied, this overrides \code{color_scheme}. Values may be LaTeX/xcolor
#'   names or 6-digit hex colors.
#' @param match_text_color Character. Color used for matching residues in non-reference
#'   rows when displayed on a white background. May be a LaTeX color name or hex value.
#' @param residues_per_line Integer. Target number of residues per line in the
#'   TeXshade output. The actual rendered number may be lower depending on page
#'   width, font, and labels.
#' @param cut.site Numeric. Points to the nt immediately to the left of the nuclease
#'   cut site within the windowed sub-sequence used for alignment printing.
#' @param include_indels Logical. If \code{TRUE}, residues aligned opposite gaps in
#'   the reference are treated as mismatches and shaded accordingly.
#' @param sanitize_names Logical. If \code{TRUE}, sequence names are converted to
#'   LaTeX-safe identifiers by replacing non-alphanumeric characters with
#'   underscores.
#' @param font_family Character. Monospaced font family to use under pdfLaTeX.
#'   Supported values are \code{"inconsolata"} (default), \code{"beramono"},
#'   \code{"courier"}, and \code{"lmodern"}.
#' @param seq_type Character. TeXshade sequence type: \code{"N"} for nucleotides
#'   or \code{"P"} for proteins.
#'
#' @return Invisibly returns the normalized path to \code{outfile}.
#'
#' @details
#' Colors specified as hex values are converted internally to LaTeX color
#' definitions using \code{\\definecolor} with the \code{HTML} model. Named colors
#' are passed through directly.
#'
#' Gap characters (\code{-}) are styled independently of residue shading and are
#' forced to render as black text on a white background for consistency.
#'
#' Page layout is configured using the \code{geometry} package (landscape mode
#' with reduced margins) to maximize the number of residues per line.
#'
#' @examples
#' \dontrun{
#' library(msa)
#' seqs <- DNAStringSet(c(
#'   ref = "ACGTACGTACGT",
#'   s1  = "ACGTTCGTAC-T",
#'   s2  = "ACG--CGTACGT"
#' ))
#'
#' aln <- msa::msaClustalW(seqs)
#' write_texshade_msa(
#'   aln,
#'   outfile = "alignment.tex",
#'   reference = "ref",
#'   color_scheme = "nt_default",
#'   match_text_color = "#A6A6A6",
#'   font_family = "inconsolata",
#'   seq_type = "N"
#' )
#' }
#'
#' @seealso \code{\link{make_texshade_from_dna}},
#'   \code{\link{make_texshade_from_aa}}, \code{\link[msa]{msaClustalW}}
#'
#' @export
write_texshade_msa <- function(aln,
                               outfile,
                               reference = 1L,
                               seq_type = c("N", "P"),
                               color_scheme = c("auto", "nt_default", "nt_okabe_ito", "nt_shapely","nt_alt",
                                                "aa_chemistry", "aa_shapely", "aa_letter1", "aa_letter2", "aa_letter3"),
                               color_map = NULL,
                               match_text_color = "A6A6A6",
                               residues_per_line = 100L,
                               cut.site = NULL,
                               include_indels = TRUE,
                               sanitize_names = TRUE,
                               font_family = "inconsolata") {
  if (!inherits(aln, "MultipleAlignment")) {
    stop("aln must inherit from a MultipleAlignment class.")
  }

  seq_type <- match.arg(seq_type)
  color_scheme <- match.arg(color_scheme)

  # Manually set sequence names
  aln_names <- BiocGenerics::rownames(aln)
  # BiocGenerics::rownames(aln) <- 1:nrow(aln)
  aln_names <- stringr::str_replace_all(aln_names, "%", "\\\\%")
  tex_aln_names <- sapply(seq_along(aln_names), function(x) {
    stringr::str_c("  \\nameseq{",x,"}{", aln_names[x], "}")
  })

  # cut site
  if (!is.null(cut.site)){
    tex_cut_sites <- sapply(cut.site, function(x) {
      stringr::str_c("  \\feature{top}{1}{", x, "..", x, "}{restriction[Black]}{Cut Site}")
    })
    tex_cut_sites <- c("\\featuresfootnotesize", "\\featurestt", tex_cut_sites)
  } else {
    tex_cut_sites <- ""
  }

  aln <- .nucleaseAlignmentCorrection(aln, cut.site, reference)

  # Extract aligned rows (gaps preserved)
  seq_class <- switch(
    seq_type,
    N = "DNAStringSet",
    P = "AAStringSet"
  )

  dss <- as(aln, seq_class)
  seqs0 <- as.character(dss)
  ids0  <- names(dss)

  if (is.null(ids0)) {
    ids0 <- rep("", length(seqs0))
  }
  blank <- !nzchar(ids0)
  if (any(blank)) {
    ids0[blank] <- paste0("seq", which(blank))
  }

  if (sanitize_names) {
    ids0 <- gsub("[^A-Za-z0-9_.-]", "_", ids0)
  }

  # Resolve reference row in the ORIGINAL alignment
  ref_idx <- if (is.numeric(reference)) {
    as.integer(reference)[1L]
  } else {
    match(reference, ids0)
  }
  if (is.na(ref_idx) || ref_idx < 1L || ref_idx > length(seqs0)) {
    stop("reference does not resolve to a valid row.")
  }

  # Reorder so the reference is row 1 in the exported FASTA / TeXshade output
  ord <- c(ref_idx, setdiff(seq_along(seqs0), ref_idx))
  seqs <- seqs0[ord]
  ids  <- ids0[ord]

  # Map original row index -> TeXshade row index after reordering
  tex_idx_map <- match(seq_along(seqs0), ord)

  # normalize color definitions
  palettes <- list(
    nt_default = nt_default,
    nt_okabe_ito = nt_okabe_ito,
    nt_shapely = nt_shapely,
    nt_alt = nt_alt,
    aa_chemistry = aa_chemistry,
    aa_shapely = aa_shapely,
    aa_letter1 = aa_letter1,
    aa_letter2 = aa_letter2,
    aa_letter3 = aa_letter3
  )

  if (!is.null(color_map)) {
    color_scheme <- color_map
  } else if (color_scheme == "auto") {
    color_scheme <- switch(
      seq_type,
      N = palettes$nt_default,
      P = palettes$aa_chemistry
    )
  } else {
    color_scheme <- palettes[[color_scheme]]
  }

  if (is.null(color_scheme) || length(color_scheme) == 0L) {
    stop("color_scheme resolved to an empty palette.")
  }

  if (is.null(names(color_scheme)) || all(!nzchar(names(color_scheme)))) {
    stop("color_scheme must be a named vector of residue colors.")
  }

  color_scheme <- color_scheme[!is.na(names(color_scheme))]
  names(color_scheme) <- toupper(names(color_scheme))

  is_hex <- function(x) grepl("^#?[0-9A-Fa-f]{6}$", x)
  resolve_tex_color <- function(x, prefix = "c") {
    if (is_hex(x)) {
      list(name = paste0(prefix, "_", gsub("^#", "", x)),
           def  = sprintf("\\definecolor{%s}{HTML}{%s}",
                          paste0(prefix, "_", gsub("^#", "", x)),
                          sub("^#", "", x)))
    } else {
      list(name = x, def = NULL)
    }
  }

  tex_color_names <- paste0("base", names(color_scheme))
  names(tex_color_names) <- names(color_scheme)
  color_defs <- character(0)
  for (b in names(color_scheme)) {
    col <- resolve_tex_color(color_scheme[[b]], prefix = paste0("base", b))
    tex_color_names[[b]] <- col$name
    if (!is.null(col$def)) color_defs <- c(color_defs, col$def)
  }

  # Color for matching residues on white background
  mt <- resolve_tex_color(match_text_color, prefix = "match")
  match_text_color_tex <- mt$name
  if (!is.null(mt$def)) color_defs <- c(color_defs, mt$def)

  match_ranges <- function(ref_seq, qry_seq) {
    r <- strsplit(toupper(ref_seq), "", fixed = TRUE)[[1L]]
    q <- strsplit(toupper(qry_seq), "", fixed = TRUE)[[1L]]
    if (length(r) != length(q)) {
      stop("Aligned sequences must have the same length.")
    }

    qry_pos <- 0L
    runs <- character(0)
    in_run <- FALSE
    start <- end <- NA_integer_

    flush_run <- function() {
      if (in_run) {
        runs <<- c(runs, sprintf("%d..%d", start, end))
        in_run <<- FALSE
      }
    }

    for (k in seq_along(r)) {
      rc <- r[[k]]
      qc <- q[[k]]

      if (qc != "-") qry_pos <- qry_pos + 1L

      keep <- (qc != "-") && (rc == qc)

      if (keep) {
        if (!in_run) {
          start <- qry_pos
          in_run <- TRUE
        }
        end <- qry_pos
      } else {
        flush_run()
      }
    }

    flush_run()
    paste(runs, collapse = ",")
  }

  collapse_ranges <- function(pos) {
    pos <- sort(unique(as.integer(pos)))
    if (!length(pos)) return("")
    grp <- cumsum(c(TRUE, diff(pos) != 1L))
    pieces <- vapply(split(pos, grp), function(v) {
      sprintf("%d..%d", v[1L], v[length(v)])
    }, character(1L))
    paste(pieces, collapse = ",")
  }

  # Position lists are in TeXshade residue numbering (gap-free residue counts).
  # This is important because TeXshade's manual selection syntax is residue-based.
  base_positions <- function(aln_seq) {
    chars <- strsplit(toupper(aln_seq), "", fixed = TRUE)[[1L]]
    pos <- 0L
    out <- setNames(vector("list", length(color_scheme)), names(color_scheme))
    for (nm in names(out)) out[[nm]] <- integer(0)

    for (ch in chars) {
      if (ch != "-") {
        pos <- pos + 1L
        if (ch %in% names(out)) out[[ch]] <- c(out[[ch]], pos)
      }
    }
    out
  }

  mismatch_positions <- function(ref_seq, qry_seq) {
    r <- strsplit(toupper(ref_seq), "", fixed = TRUE)[[1L]]
    q <- strsplit(toupper(qry_seq), "", fixed = TRUE)[[1L]]
    if (length(r) != length(q)) {
      stop("Aligned sequences must have the same length.")
    }

    ref_pos <- 0L
    qry_pos <- 0L
    out <- setNames(vector("list", length(color_scheme)), names(color_scheme))
    for (nm in names(out)) out[[nm]] <- integer(0)

    for (k in seq_along(r)) {
      rc <- r[[k]]
      qc <- q[[k]]

      if (rc != "-") ref_pos <- ref_pos + 1L
      if (qc != "-") qry_pos <- qry_pos + 1L

      if (qc == "-") next

      mismatch <- (rc == "-") || (rc != qc)
      if (mismatch) {
        #if (rc == "-" && !include_indels) next
        if (qc %in% names(out)) out[[qc]] <- c(out[[qc]], qry_pos)
      }
    }
    out
  }

  tex_file <- normalizePath(outfile, mustWork = FALSE)
  fasta_name <- paste0(tools::file_path_sans_ext(basename(outfile)), ".fasta")

  # Build embedded FASTA content
  fasta_lines <- c(
    sprintf("\\begin{filecontents*}{%s}", fasta_name),
    unlist(Map(function(id, seq) c(paste0(">", id), seq), ids, seqs), use.names = FALSE),
    "\\end{filecontents*}"
  )

  # Build shading commands
  cmd_lines <- c(
    #sprintf("\\orderseqs{%s}", paste(seq_along(seqs), collapse = "-")),
    sprintf("\\donotshade{%s}", paste(seq_along(seqs), collapse = ",")),
    "\\hideconsensus",
    sprintf("\\residuesperline{%d}", as.integer(residues_per_line)),
    "\\showruler{top}{1}",
    "\\rulersteps{10}",
    "\\gapchar{-}",
    "\\gapcolors{Black}{White}"
  )

  # Reference row is row 1 after reordering
  ref_base_pos <- base_positions(seqs[[1L]])
  for (b in names(color_scheme)) {
    rng <- collapse_ranges(ref_base_pos[[b]])
    if (nzchar(rng)) {
      cmd_lines <- c(cmd_lines,
                     sprintf("\\shaderegion{1}{%s}{Black}{%s}", rng, tex_color_names[[b]]))
    }
  }

  # Mismatch shading for all other rows
  for (orig_i in setdiff(seq_along(seqs0), ref_idx)) {
    tex_i <- tex_idx_map[[orig_i]]
    mm <- mismatch_positions(seqs0[[ref_idx]], seqs0[[orig_i]])
    for (b in names(color_scheme)) {
      rng <- collapse_ranges(mm[[b]])
      if (nzchar(rng)) {
        cmd_lines <- c(cmd_lines,
                       sprintf("\\shaderegion{%d}{%s}{Black}{%s}", tex_i, rng, tex_color_names[[b]]))
      }
    }

    # Matching residues get white background + custom foreground
    rng <- match_ranges(seqs0[[ref_idx]], seqs0[[orig_i]])
    if (nzchar(rng)) {
      cmd_lines <- c(cmd_lines,
                     sprintf("\\shaderegion{%d}{%s}{%s}{White}", tex_i, rng, match_text_color_tex))
    }
  }

  # Set alignment font
  font_family <- tolower(font_family)

  font_lines <- switch(font_family,
                       "inconsolata" = c(
                         "\\usepackage[T1]{fontenc}",
                         "\\usepackage{inconsolata}",
                         "\\renewcommand*\\familydefault{\\ttdefault}"
                       ),
                       "beramono" = c(
                         "\\usepackage[T1]{fontenc}",
                         "\\usepackage[scaled]{beramono}",
                         "\\renewcommand*\\familydefault{\\ttdefault}"
                       ),
                       "courier" = c(
                         "\\usepackage[T1]{fontenc}",
                         "\\usepackage{courier}",
                         "\\renewcommand*\\familydefault{\\ttdefault}"
                       ),
                       "lmodern" = c(
                         "\\usepackage[T1]{fontenc}",
                         "\\usepackage{lmodern}",
                         "\\renewcommand*\\familydefault{\\ttdefault}"
                       ),
                       stop("Unsupported font_family. Choose from: inconsolata, beramono, courier, lmodern")
  )

  tex <- c(
    "\\documentclass{article}",
    "\\usepackage[landscape,margin=0.5in]{geometry}",
    "\\usepackage[dvipsnames,svgnames,HTML]{xcolor}",
    font_lines,
    "\\usepackage{texshade}",
    "\\pagestyle{empty}",
    "\\begin{document}",
    fasta_lines,
    color_defs,
    sprintf("\\begin{texshade}{%s}", fasta_name),
    sprintf("  \\seqtype{%s}", seq_type),
    sprintf("  \\shownames{left}"),
    tex_aln_names,
    paste0("  ", cmd_lines),
    tex_cut_sites,
    "\\end{texshade}",
    "\\end{document}",
    ""
  )

  writeLines(tex, tex_file, useBytes = TRUE)
  invisible(tex_file)
}

#' Align DNA sequences and write a TeXshade file
#'
#' Convenience wrapper that passes a sequence alignment to \code{write_texshade_msa()}
#' to generate a standalone TeXshade document. Then calls the tinytex package
#' to output a formatted pdf of the generated .tex file.
#'
#' This wrapper is intended for aligned DNA sequences. The resulting TeX file
#' highlights the reference sequence by nucleotide, shades mismatches in other
#' sequences, renders matching residues with a configurable text color on a white
#' background, and styles gap characters consistently.
#'
#' @param aln A DNA sequence alignment of class \code{DNAMultipleAlignment}.
#' @param outfile Character. Path to the output \code{.tex} file.
#' @param reference Reference sequence, either a numeric index or a sequence name.
#'   Passed through to \code{\link{write_texshade_msa}} after alignment.
#' @param color_scheme A character vector. Specifies the color scheme to be used by
#'   \code{\link{write_texshade_msa}}. Default is 'auto'.
#' @param ... Additional arguments passed to \code{\link{write_texshade_msa}}, such as
#'   \code{match_text_color}, \code{residues_per_line},
#'   \code{include_indels}, \code{sanitize_names}, and \code{font_family}.
#'
#' @return Invisibly returns the normalized path to \code{outfile}.
#'
#' @details
#' The wrapper performs export only; all TeXshade formatting and color logic
#' are handled by \code{write_texshade_msa()}. See that manual page for the full
#' set of output customization options.
#'
#' @seealso \code{\link{write_texshade_msa}}, \code{\link[msa]{msaClustalW}}
#'
#' @export
make_texshade_from_dna <- function(aln,
                                   outfile,
                                   reference = 1L,
                                   color_scheme = "auto",
                                   ...) {
  tex_file <- write_texshade_msa(aln, outfile = outfile, reference = reference, seq_type = "N", color_scheme = color_scheme, ...)
  tinytex::pdflatex(tex_file)
  return(tex_file)
}


#' Align protein sequences and write a TeXshade file
#'
#' Convenience wrapper that passes a protein sequence alignment to
#' \code{write_texshade_msa()} to generate a standalone TeXshade document.
#' Then calls the tinytex package to output a formatted pdf of the generated .tex file.
#'
#' This wrapper is intended for aligned amino acid sequences. The resulting TeX file
#' highlights the reference sequence by residue, shades mismatches in other
#' sequences, renders matching residues with a configurable text color on a white
#' background, and styles gap characters consistently.
#'
#' @param aln A protein sequence alignment of class \code{AAMultipleAlignment}.
#' @param outfile Character. Path to the output \code{.tex} file.
#' @param reference Reference sequence, either a numeric index or a sequence name.
#'   Passed through to \code{\link{write_texshade_msa}} after alignment.
#' @param color_scheme A character vector. Specifies the color scheme to be used by
#'   \code{\link{write_texshade_msa}}. Default is 'auto'.
#' @param ... Additional arguments passed to \code{\link{write_texshade_msa}}, such as
#'   \code{match_text_color}, \code{residues_per_line},
#'   \code{include_indels}, \code{sanitize_names}, and \code{font_family}.
#'
#' @return Invisibly returns the normalized path to \code{outfile}.
#'
#' @details
#' The wrapper performs export only; all TeXshade formatting and color logic
#' are handled by \code{write_texshade_msa()}. See that manual page for the full
#' set of output customization options.
#'
#' @seealso \code{\link{write_texshade_msa}}, \code{\link[msa]{msaClustalW}}
#'
#' @export
make_texshade_from_aa <- function(aln,
                                  outfile,
                                  reference = 1L,
                                  color_scheme = "auto",
                                  ...) {
  tex_file <- write_texshade_msa(aln, outfile = outfile, reference = reference, seq_type = "P", color_scheme = color_scheme, ...)
  tinytex::pdflatex(tex_file)
  return(tex_file)
}




# library(msa)
# library(Biostrings)
#
# seqs <- DNAStringSet(c(
#   ref = "ACGTTGCAACGT",
#   q1  = "ACGTCGCAACGT",
#   q2  = "ACGTTGGAACGT",
#   q3  = "ATGTTGCAATGT"
# ))
#
# # Option 1: create the MSA first, then export
# aln <- msaClustalW(seqs, order = "input")
# write_texshade_msa(aln, "alignment_report.tex", reference = "ref")
#
# # Option 2: do both steps at once
# make_texshade_from_dna(seqs, "alignment_report.tex", reference = "ref")







