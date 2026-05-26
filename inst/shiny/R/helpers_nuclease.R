.nuclease_setting_cols <- function() {
  c(
    "Forward Primer",
    "Reverse Primer",
    "Forward Extension Type",
    "Forward Extension Sequence",
    "Reverse Extension Type",
    "Reverse Extension Sequence"
  )
}

.nuclease_default_settings <- function() {
  out <- as.list(rep("", length(.nuclease_setting_cols())))
  names(out) <- .nuclease_setting_cols()
  out[["Forward Extension Type"]] <- "None"
  out[["Reverse Extension Type"]] <- "None"
  out
}

.nuclease_empty_queue <- function(id_col = "Nuclease Queue ID") {
  cols <- c(
    id_col,
    "Sample Name",
    "Control File",
    "Control Path",
    "Experimental File",
    "Experimental Path",
    "gRNA Sequence",
    "Cut Position",
    "Source",
    .nuclease_setting_cols(),
    "Status",
    "Progress"
  )
  as.data.frame(setNames(rep(list(character()), length(cols)), cols), check.names = FALSE)
}

.normalize_nuclease_row <- function(row) {
  defaults <- .nuclease_default_settings()
  for (nm in names(defaults)) {
    if (!nm %in% names(row)) row[[nm]] <- defaults[[nm]]
    row[[nm]] <- as.character(row[[nm]])
    row[[nm]][is.na(row[[nm]])] <- ""
  }

  for (nm in c("Sample Name", "Control File", "Experimental File", "gRNA Sequence", "Cut Position")) {
    if (!nm %in% names(row)) row[[nm]] <- ""
    row[[nm]] <- as.character(row[[nm]])
    row[[nm]][is.na(row[[nm]])] <- ""
  }

  row
}

.nuclease_settings_summary <- function(row) {
  row <- .normalize_nuclease_row(row)
  ext <- paste0(row[["Forward Extension Type"]], "/", row[["Reverse Extension Type"]])
  paste(
    paste0("Cut:", ifelse(nzchar(row[["Cut Position"]]), row[["Cut Position"]], "-")),
    paste0("gRNA:", ifelse(nzchar(row[["gRNA Sequence"]]), row[["gRNA Sequence"]], "-")),
    paste0("Ext:", ext)
  )
}

.make_nuclease_table_data <- function(
    df,
    id_col,
    remove_input_id,
    move_up_input_id = NULL,
    move_down_input_id = NULL,
    settings_input_id = "nuclease_open_settings"
) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame())
  }

  progress <- if ("Progress" %in% names(df)) df[["Progress"]] else rep("", nrow(df))
  actions <- vapply(df[[id_col]], function(row_id) {
    .queue_actions_html(row_id, remove_input_id, move_up_input_id, move_down_input_id)
  }, character(1))

  settings <- vapply(seq_len(nrow(df)), function(i) {
    .mutation_settings_button(df[[id_col]][i], settings_input_id, .nuclease_settings_summary(df[i, , drop = FALSE]))
  }, character(1))

  data.frame(
    "::" = .queue_drag_handle_html(),
    "Display Order" = seq_len(nrow(df)),
    "Sample Name" = df[["Sample Name"]],
    "Control File" = df[["Control File"]],
    "Experimental File" = df[["Experimental File"]],
    "gRNA Sequence" = df[["gRNA Sequence"]],
    "Cut Position" = df[["Cut Position"]],
    "Source" = if ("Source" %in% names(df)) df[["Source"]] else rep("single", nrow(df)),
    "Settings" = settings,
    "Status" = vapply(df[["Status"]], .sample_status_html, character(1)),
    "Actions" = actions,
    "Progress" = .sample_progress_html(progress),
    "Sample ID" = df[[id_col]],
    check.names = FALSE
  )
}

.nuclease_status_for_row <- function(row, all_samples = character()) {
  row <- .normalize_nuclease_row(row)
  sample_name <- trimws(row[["Sample Name"]])
  control_file <- trimws(row[["Control File"]])
  experimental_file <- trimws(row[["Experimental File"]])
  grna <- trimws(row[["gRNA Sequence"]])
  cut_position <- trimws(row[["Cut Position"]])

  issues <- c(
    if (!nzchar(sample_name)) "Missing sample name",
    if (nzchar(sample_name) && sum(all_samples == sample_name, na.rm = TRUE) > 1) "Duplicate name",
    if (!.fastq_merged_ok(control_file)) "Invalid control FASTQ",
    if (!.fastq_merged_ok(experimental_file)) "Invalid experimental FASTQ",
    if (!nzchar(grna)) "Missing gRNA sequence",
    if (!nzchar(cut_position)) "Missing cut position",
    if (nzchar(cut_position) && is.na(suppressWarnings(as.numeric(cut_position)))) "Invalid cut position",
    if (!nzchar(trimws(row[["Forward Primer"]]))) "Missing forward primer",
    if (!nzchar(trimws(row[["Reverse Primer"]]))) "Missing reverse primer"
  )

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

.nuclease_queue_with_status <- function(df) {
  if (is.null(df) || nrow(df) == 0 || "Note" %in% names(df)) return(df)
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  samples <- df[["Sample Name"]]
  df$Status <- vapply(seq_len(nrow(df)), function(i) {
    .nuclease_status_for_row(df[i, , drop = FALSE], samples)
  }, character(1))
  df
}
