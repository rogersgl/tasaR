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
    "CDR4 Start",
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
    "Progress"
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
  insert <- paste0(row[["Insert Start Coordinate"]], "-", row[["Insert End Coordinate"]])
  ext <- paste0(row[["Forward Extension Type"]], "/", row[["Reverse Extension Type"]])
  paste(antibody, paste0("Amp:", ifelse(nzchar(amplicon), amplicon, "-")), paste0("Insert:", insert), paste0("Ext:", ext))
}

.mutation_settings_button <- function(row_id, input_id, label) {
  sprintf(
    "<button type='button' class='queue-settings-btn' title='Edit settings' aria-label='Edit settings' onclick=\"Shiny.setInputValue('%s', '%s', {priority: 'event'})\">%s</button>",
    .js_string(input_id),
    .js_string(row_id),
    htmltools::htmlEscape(label)
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
    "Source" = if ("Source" %in% names(df)) df[["Source"]] else rep("single", nrow(df)),
    "Settings" = settings,
    "Status" = vapply(df[["Status"]], .sample_status_html, character(1)),
    "Actions" = actions,
    "Progress" = .sample_progress_html(progress),
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
      "CDR4 Start",
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
