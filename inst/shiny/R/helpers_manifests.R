.write_mutation_working_manifest <- function(df, path, order = character(), id_col = "Batch Queue ID") {
  if (is.null(path) || is.null(df) || nrow(df) == 0 || "Note" %in% names(df)) {
    return(invisible(FALSE))
  }

  out <- .order_sample_df(df, order, id_col)
  if (id_col %in% names(out)) out[[id_col]] <- NULL
  out[["Progress"]] <- NULL

  utils::write.csv(out, path, row.names = FALSE, quote = TRUE)
  invisible(TRUE)
}

.write_batch_working_manifest <- function(df, path, order = character()) {
  if (is.null(path) || is.null(df) || nrow(df) == 0 || "Note" %in% names(df)) {
    return(invisible(FALSE))
  }

  out <- .order_sample_df(df, order, "Batch Queue ID")
  out[["Batch Queue ID"]] <- NULL
  out[["Progress"]] <- NULL

  utils::write.csv(
    out,
    path,
    row.names = FALSE,
    quote = TRUE
  )

  invisible(TRUE)
}

.prepare_merge_inputs <- function(
    mode = c("single", "batch"),
    queued_df = NULL,
    manifest_path = NULL,
    archive_path = NULL,
    extract_dir = NULL
) {
  mode <- match.arg(mode)

  if (mode == "single") {
    if (is.null(queued_df) || nrow(queued_df) == 0) {
      stop("No queued single-sample inputs available.")
    }

    required_cols <- c("Sample Name", "R1 Path", "R2 Path", "Status")
    missing <- setdiff(required_cols, names(queued_df))
    if (length(missing) > 0) {
      stop(
        paste("Single-sample queue is missing required columns:",
              paste(missing, collapse = ", "))
      )
    }

    if (any(queued_df$Status != "Ready")) {
      bad_samples <- queued_df[queued_df$Status != "Ready", "Sample Name", drop = TRUE]
      stop(
        paste(
          "Some queued samples are not ready:",
          paste(bad_samples, collapse = ", ")
        )
      )
    }

    out <- data.frame(
      Sample_Name = queued_df[["Sample Name"]],
      R1_Path = queued_df[["R1 Path"]],
      R2_Path = queued_df[["R2 Path"]],
      Source = if ("Source" %in% names(queued_df)) queued_df[["Source"]] else "single"
    )

    return(out)
  }

  if (is.null(manifest_path) || !file.exists(manifest_path)) {
    stop("Batch manifest path is missing or does not exist.")
  }

  if (is.null(archive_path) || !file.exists(archive_path)) {
    stop("Batch archive path is missing or does not exist.")
  }

  if (is.null(extract_dir)) {
    stop("Batch extract_dir is required.")
  }

  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

  manifest_df <- data.table::fread(
    manifest_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c("Sample Name", "Left (R1) File", "Right (R2) File")
  missing <- setdiff(required_cols, names(manifest_df))
  if (length(missing) > 0) {
    stop(
      paste("Batch manifest is missing required columns:",
            paste(missing, collapse = ", "))
    )
  }

  utils::unzip(archive_path, exdir = extract_dir)

  extracted_files <- list.files(
    extract_dir,
    recursive = TRUE,
    full.names = TRUE
  )
  extracted_basenames <- basename(extracted_files)

  manifest_df$R1_Path <- NA_character_
  manifest_df$R2_Path <- NA_character_

  for (i in seq_len(nrow(manifest_df))) {
    r1_match <- match(manifest_df[i, "Left (R1) File", drop = TRUE], extracted_basenames)
    r2_match <- match(manifest_df[i, "Right (R2) File", drop = TRUE], extracted_basenames)

    if (!is.na(r1_match)) {
      manifest_df[i, "R1_Path"] <- extracted_files[r1_match]
    }
    if (!is.na(r2_match)) {
      manifest_df[i, "R2_Path"] <- extracted_files[r2_match]
    }
  }

  missing_r1 <- is.na(manifest_df$R1_Path)
  missing_r2 <- is.na(manifest_df$R2_Path)

  if (any(missing_r1 | missing_r2)) {
    missing_samples <- manifest_df[missing_r1 | missing_r2,
                                   "Sample Name", drop = TRUE]
    stop(
      paste(
        "Missing FASTQ files for:",
        paste(missing_samples, collapse = ", ")
      )
    )
  }

  data.frame(
    Sample_Name = manifest_df[["Sample Name"]],
    R1_Path = manifest_df[["R1_Path"]],
    R2_Path = manifest_df[["R2_Path"]],
    Source = "batch"
    )
}
