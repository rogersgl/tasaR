.read_batch_manifest <- function(path) {
  df <- data.table::fread(path, check.names = FALSE)
  .validate_batch_manifest(df)
  if (nrow(df) == 0) stop("Batch manifest is empty.")
  df
}

.remove_queued_sample <- function(queue_df, queue_id) {

  idx <- match(queue_id, queue_df[["Queue ID"]])

  if (is.na(idx)) {
    return(list(
      queued = queue_df,
      removed = NULL,
      success = FALSE
    ))
  }

  row <- queue_df[idx, , drop = FALSE]

  updated_queue <- queue_df[-idx, , drop = FALSE]

  list(
    queued = updated_queue,
    removed = row,
    success = TRUE
  )
}

.check_fastq_pair <- function(sn, r1, r2, other_samples = character(), other_files = character()) {
  sn_ok <- identical(make.names(sn), sn)
  r1_ok <- grepl("\\.fastq(\\.gz)?$", r1, ignore.case = TRUE)
  r2_ok <- grepl("\\.fastq(\\.gz)?$", r2, ignore.case = TRUE)

  sn_unique <- !sn %in% other_samples
  r1_unique <- !r1 %in% other_files
  r2_unique <- !r2 %in% other_files

  issues <- c(
    if (!sn_ok) "Invalid name",
    if (!r1_ok) "Invalid R1 FASTQ",
    if (!r2_ok) "Invalid R2 FASTQ",
    if (!sn_unique) "Duplicate name",
    if (!r1_unique) "Duplicate R1 file",
    if (!r2_unique) "Duplicate R2 file"
  )

  if (length(issues) == 0) "Ready" else paste(issues, collapse = "<br>")
}

.validate_batch_manifest <- function(df) {
  required_cols <- c("Sample Name", "Left (R1) File", "Right (R2) File")
  missing <- setdiff(required_cols, names(df))

  if (length(missing) > 0) {
    stop(
      paste("Missing required columns:", paste(missing, collapse = ", "))
    )
  }

  invisible(TRUE)
}

.sample_table_with_status <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  .validate_batch_manifest(df)

  if (nrow(df) == 0) {
    df$Status <- character(0)
    return(df)
  }

  df$Status <- vapply(seq_len(nrow(df)), function(i) {
    other_samples <- if (nrow(df) > 1) df[["Sample Name"]][-i] else character(0)
    other_files <- if (nrow(df) > 1) {
      unlist(df[-i, c("Left (R1) File", "Right (R2) File")], use.names = FALSE)
    } else {
      character(0)
    }

    .check_fastq_pair(
      df[i, "Sample Name", drop = TRUE],
      df[i, "Left (R1) File", drop = TRUE],
      df[i, "Right (R2) File", drop = TRUE],
      other_samples = other_samples,
      other_files = other_files
    )
  }, character(1))

  df
}

.get_merge_settings <- function(input) {
  list(algorithm = input$merge_engine,
       min_overlap = input$min_overlap,
       min_length = input$min_length,
       max_length = input$max_length,
       merge_similarity = input$merge_similarity,
       forward_primer = input$forward_primer,
       reverse_primer = input$reverse_primer,
       merge_filetype = input$merge_filetype,
       merge_discard_ambig = input$merge_discard_ambig,
       merge_seq_build_log = input$merge_seq_build_log, #b/B
       merge_file_errors = input$merge_file_errors,  #f/F
       merge_statistics = input$merge_statistics, #s/S
       merge_seq_reconstruction_log = input$merge_seq_reconstruction_log, #r/R
       merge_kmer = input$merge_kmer, #k/K
       merge_prefix = input$merge_prefix
       )
}

# pandaseq_merge_files <- function(forward_fastq,
#                                  reverse_fastq,
#                                  output_fastq,
#                                  log_file = NULL,
#                                  min_length = NULL,
#                                  max_length = NULL,
#                                  extra_args = c("-F", "-d", "bFSrk"),
#                                  verbose = FALSE)

# merge_inputs
# Sample_Name
# R1_Path
# R2_Path
# Source
.run_merge_pipeline <- function(merge_inputs, settings){
  log.path <- file.path(tempdir(), "merge", stringr::str_c(merge_inputs$Sample_Name,"-log.txt"))
  output.path <- file.path(tempdir(), "merge", stringr::str_c(if(nzchar(settings$merge_prefix)) stringr::str_c(settings$merge_prefix, "-"),
                                                                         merge_inputs$Sample_Name,
                                                                         if(settings$merge_filetype == "FASTQ"){".fastq"
                                                                           } else if(settings$merge_filetype == "FASTA"){".fasta"
                                                                             } else(stop("Merge file type unclear. Too confused to continue."))))

  extra_args <- .pandaseq_extra_args(settings)
  cat(extra_args)
  for (i in 1:nrow(merge_inputs)) {
    this_row <- merge_inputs[i,]
    pandaseq_merge_files(forward_fastq = this_row$R1_Path,
                         reverse_fastq = this_row$R2_Path,
                         output_fastq = output.path,
                         log_file = log.path,
                         min_length = settings$min_length,
                         max_length = settings$max_length,
                         extra_args = extra_args,
                         verbose = FALSE
    )
  }

}


.pandaseq_extra_args <- function(settings) {
  # -o minoverlap
  # -t threshold
  # -p forwardprimer
  # -q reverseprimer
  # -N
  # -d flags
  c("-F",
    "-o", settings$min_overlap,
    "-t", settings$merge_similarity,
    if(any(settings$merge_seq_build_log,
               settings$merge_file_errors,
               settings$merge_statistics,
               settings$merge_seq_reconstruction_log,
               settings$merge_kmer)){
      c("-d", stringr::str_c(ifelse(settings$merge_seq_build_log, "B", "b"),
                             ifelse(settings$merge_file_errors, "F", "f"),
                             ifelse(settings$merge_statistics, "S", "s"),
                             ifelse(settings$merge_seq_reconstruction_log, "R", "r"),
                             ifelse(settings$merge_kmer, "K", "k"))

      )
    } else {character()},
    if(nzchar(settings$forward_primer)) c("-p", "settings$forward_primer") else character(),
    if(nzchar(settings$reverse_primer)) c("-q", "settings$reverse_primer") else character(),
    if(settings$merge_discard_ambig) "-N" else character()
    )
}

