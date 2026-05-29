.fastq_merged_ok <- function(path) {
  grepl("\\.(fastq|fq)(\\.gz)?$", path, ignore.case = TRUE)
}

.mutation_setting_cols <- function() {
  c(
    "Antibody",
    "Reference Sequence",
    "Forward Primer",
    "Reverse Primer",
    "Amplicon Length",
    "Insert Start Coordinate",
    "Insert End Coordinate",
    "Forward Extension Type",
    "Forward Extension Sequence",
    "Reverse Extension Type",
    "Reverse Extension Sequence",
    "FR1 Start",
    "CDR1 Start",
    "FR2 Start",
    "CDR2 Start",
    "FR3 Start",
    "CDR3 Start",
    "FR4 Start",
    "FR4 End"
  )
}

.mutation_default_settings <- function() {
  out <- as.list(rep("", length(.mutation_setting_cols())))
  names(out) <- .mutation_setting_cols()
  out[["Antibody"]] <- FALSE
  out[["Forward Extension Type"]] <- "None"
  out[["Reverse Extension Type"]] <- "None"
  out
}

.mutation_empty_queue <- function(id_col = "Mutation Queue ID") {
  cols <- c(
    id_col,
    "Sample Name",
    "Merged File",
    "Merged Path",
    "Source",
    .mutation_setting_cols(),
    "Status",
    "Progress",
    "Result Path",
    "Output Path",
    "Download Href"
  )
  out <- as.data.frame(setNames(rep(list(character()), length(cols)), cols), check.names = FALSE)
  out
}

.normalize_mutation_row <- function(row) {
  defaults <- .mutation_default_settings()
  for (nm in names(defaults)) {
    if (!nm %in% names(row)) row[[nm]] <- defaults[[nm]]
    if (nm == "Antibody") {
      row[[nm]] <- as.logical(row[[nm]])
      if (is.na(row[[nm]])) row[[nm]] <- FALSE
    } else {
      row[[nm]] <- as.character(row[[nm]])
      row[[nm]][is.na(row[[nm]])] <- ""
    }
  }
  row
}

.mutation_settings_summary <- function(row) {
  antibody <- if (isTRUE(as.logical(row[["Antibody"]]))) "Ab" else "Non-Ab"
  amplicon <- row[["Amplicon Length"]]
  insert <- paste0(row[["Insert Start Coordinate"]], ",  ", row[["Insert End Coordinate"]])
  ext <- paste0(if(row[["Forward Extension Type"]] == "Barcode") "BC" else row[["Forward Extension Type"]],
                " / ",
                if(row[["Reverse Extension Type"]] == "Barcode") "BC" else row[["Reverse Extension Type"]])
  c(antibody, paste(paste0("Len: ", ifelse(nzchar(amplicon), amplicon, "-")), paste0("\nInsert: ", insert), paste0("\nExt: ", ext)))
}

.mutation_settings_button <- function(row_id, input_id, label) {
  sprintf(
    "<button type='button' class='queue-settings-btn' title='Edit settings' aria-label='Edit settings' onclick=\"Shiny.setInputValue('%s', '%s', {priority: 'event'})\">%s</button>",
    .js_string(input_id),
    .js_string(row_id),
    div(
      div(style = "display: flex; justify-content: space-between; font-weight: 700; color: #1f3b53",
          HTML(htmltools::htmlEscape(label[1])),
          icon("pencil")
        ),
      div(HTML(gsub("\n", "<br/>", htmltools::htmlEscape(label[2]))))
    )
  )
}

.make_mutation_table_data <- function(
    df,
    id_col,
    remove_input_id,
    move_up_input_id = NULL,
    move_down_input_id = NULL,
    settings_input_id = "mutation_open_settings"
) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame())
  }

  progress <- if ("Progress" %in% names(df)) df[["Progress"]] else rep("", nrow(df))
  download_href <- if ("Download Href" %in% names(df)) df[["Download Href"]] else rep("", nrow(df))
  download_filename <- if ("Output Path" %in% names(df)) basename(df[["Output Path"]]) else rep("", nrow(df))
  actions <- vapply(df[[id_col]], function(row_id) {
    .queue_actions_html(row_id, remove_input_id, move_up_input_id, move_down_input_id)
  }, character(1))

  settings <- vapply(seq_len(nrow(df)), function(i) {
    .mutation_settings_button(df[[id_col]][i], settings_input_id, .mutation_settings_summary(df[i, , drop = FALSE]))
  }, character(1))

  data.frame(
    "::" = .queue_drag_handle_html(),
    "Display Order" = seq_len(nrow(df)),
    "Sample Name" = df[["Sample Name"]],
    "Merged File" = df[["Merged File"]],
    "Settings" = settings,
    "Status" = vapply(df[["Status"]], .sample_status_html, character(1)),
    "Actions" = actions,
    "Progress" = .sample_progress_html(progress, download_href, download_filename),
    "Sample ID" = df[[id_col]],
    check.names = FALSE
  )
}

.mutation_status_for_row <- function(row, all_samples = character()) {
  row <- .normalize_mutation_row(row)
  sample_name <- trimws(row[["Sample Name"]])
  merged_file <- trimws(row[["Merged File"]])
  ref <- trimws(row[["Reference Sequence"]])
  antibody <- isTRUE(as.logical(row[["Antibody"]]))

  issues <- c(
    if (!nzchar(sample_name)) "Missing sample name",
    if (nzchar(sample_name) && sum(all_samples == sample_name, na.rm = TRUE) > 1) "Duplicate name",
    if (!.fastq_merged_ok(merged_file)) "Invalid merged FASTQ",
    if (!nzchar(ref)) "Missing reference",
    if (nchar(ref) > 500) "Reference > 500 chars",
    if (!nzchar(trimws(row[["Forward Primer"]]))) "Missing forward primer",
    if (!nzchar(trimws(row[["Reverse Primer"]]))) "Missing reverse primer"
  )

  numeric_cols <- c("Amplicon Length", "Insert Start Coordinate", "Insert End Coordinate")
  if (antibody) {
    numeric_cols <- c(
      numeric_cols,
      "FR1 Start",
      "CDR1 Start",
      "FR2 Start",
      "CDR2 Start",
      "FR3 Start",
      "CDR3 Start",
      "FR4 Start",
      "FR4 End"
    )
  }

  for (nm in numeric_cols) {
    val <- trimws(row[[nm]])
    if (!nzchar(val)) {
      issues <- c(issues, paste("Missing", tolower(nm)))
    } else if (is.na(suppressWarnings(as.numeric(val)))) {
      issues <- c(issues, paste("Invalid", tolower(nm)))
    }
  }

  if (row[["Forward Extension Type"]] %in% c("Barcode", "UMI") &&
      !nzchar(trimws(row[["Forward Extension Sequence"]]))) {
    issues <- c(issues, "Missing forward extension sequence")
  }
  if (identical(row[["Forward Extension Type"]], "None") &&
      nzchar(trimws(row[["Forward Extension Sequence"]]))) {
    issues <- c(issues, "Forward extension sequence without extension type")
  }

  if (row[["Reverse Extension Type"]] %in% c("Barcode", "UMI") &&
      !nzchar(trimws(row[["Reverse Extension Sequence"]]))) {
    issues <- c(issues, "Missing reverse extension sequence")
  }
  if (identical(row[["Reverse Extension Type"]], "None") &&
      nzchar(trimws(row[["Reverse Extension Sequence"]]))) {
    issues <- c(issues, "Reverse extension sequence without extension type")
  }

  if (length(issues) == 0) "Ready" else paste(issues, collapse = "<br>")
}

.mutation_queue_with_status <- function(df) {
  if (is.null(df) || nrow(df) == 0 || "Note" %in% names(df)) return(df)
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  samples <- df[["Sample Name"]]
  df$Status <- vapply(seq_len(nrow(df)), function(i) {
    .mutation_status_for_row(df[i, , drop = FALSE], samples)
  }, character(1))
  df
}

.mutation_analyzed_count <- function(df) {
  if (is.null(df) || nrow(df) == 0 || !"Progress" %in% names(df)) {
    return(0L)
  }

  progress <- tolower(trimws(as.character(df[["Progress"]])))
  sum(progress %in% c("done", "complete", "completed", "success", "true"), na.rm = TRUE)
}

.prep_tas_settings <- function(row) {
  data.frame(Name = row[["Sample Name"]],
             IsAntibody = row[["Antibody"]],
             MergedFASTQPath = row[["Merged Path"]],
             ReferenceSequence = row[["Reference Sequence"]],
             ForwardExtensionType = row[["Forward Extension Type"]],
             ForwardExtension = row[["Forward Extension Sequence"]],
             ForwardPrimer = row[["Forward Primer"]],
             ReverseExtensionType = row[["Reverse Extension Type"]],
             ReverseExtension = row[["Reverse Extension Sequence"]],
             ReversePrimer = row[["Reverse Primer"]],
             AmpliconLength = row[["Amplicon Length"]],
             InsertStart = row[["Insert Start Coordinate"]],
             InsertEnd = row[["Insert End Coordinate"]],
             FR1Start = row[["FR1 Start"]],
             CDR1Start = row[["CDR1 Start"]],
             FR2Start = row[["FR2 Start"]],
             CDR2Start = row[["CDR2 Start"]],
             FR3Start = row[["FR3 Start"]],
             CDR3Start = row[["CDR3 Start"]],
             FR4Start = row[["FR4 Start"]],
             FR4End = row[["FR4 End"]]
  )
}

.mutation_settings_df <- function(rows) {
  rows <- as.data.frame(rows, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(rows) || nrow(rows) == 0L) {
    return(data.frame())
  }

  data.table::rbindlist(
    lapply(seq_len(nrow(rows)), function(i) {
      .prep_tas_settings(rows[i, , drop = FALSE])
    }),
    fill = TRUE
  )
}

.mutation_sample_zip_name <- function(sample_name) {
  safe_name <- gsub("[^A-Za-z0-9_.-]", "_", sample_name)
  paste0(safe_name, ".zip")
}

.zip_directory <- function(zipfile, root_dir, compression_level = 1L) {
  if (!dir.exists(root_dir)) {
    stop("Directory not found: ", root_dir, call. = FALSE)
  }

  files <- list.files(root_dir, recursive = TRUE, all.files = FALSE, full.names = FALSE)
  if (length(files) == 0L) {
    stop("No files are available to zip in: ", root_dir, call. = FALSE)
  }

  zip::zipr(
    zipfile = zipfile,
    files = files,
    recurse = FALSE,
    compression_level = compression_level,
    include_directories = TRUE,
    root = root_dir,
    mode = "mirror"
  )
  invisible(zipfile)
}
