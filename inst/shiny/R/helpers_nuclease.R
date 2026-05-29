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

.parse_nuclease_grna <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  chartr("U", "T", toupper(trimws(x)))
}

.valid_nuclease_grna <- function(x) {
  x <- .parse_nuclease_grna(x)
  nzchar(x) & grepl("^[ATCG]{20}$", x)
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

  row[["gRNA Sequence"]] <- .parse_nuclease_grna(row[["gRNA Sequence"]])

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

.nuclease_target_display <- function(row) {
  row <- .normalize_nuclease_row(row)
  cut_position <- trimws(row[["Cut Position"]])
  grna <- trimws(row[["gRNA Sequence"]])

  if (nzchar(cut_position)) {
    cut_position
  } else if (nzchar(grna)) {
    grna
  } else {
    "auto"
  }
}

.nuclease_target_cell <- function(row, row_id, input_id = "nuclease_open_settings") {
  target <- .nuclease_target_display(row)
  sprintf(
    paste0(
      "<div class='nuclease-target-cell'>",
      "<span class='nuclease-target-value'>%s</span>",
      "<button type='button' class='queue-icon-btn queue-settings-edit-btn' title='Edit settings' aria-label='Edit settings' ",
      "onclick=\"Shiny.setInputValue('%s', '%s', {priority: 'event'})\">%s</button>",
      "</div>"
    ),
    htmltools::htmlEscape(target),
    .js_string(input_id),
    .js_string(row_id),
    div(style = "color: #1f3b53; font-weight: 700;",
        shiny::icon("pencil"))
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

  target <- vapply(seq_len(nrow(df)), function(i) {
    .nuclease_target_cell(df[i, , drop = FALSE], df[[id_col]][i], settings_input_id)
  }, character(1))

  data.frame(
    "::" = .queue_drag_handle_html(),
    "Display Order" = seq_len(nrow(df)),
    "Sample Name" = df[["Sample Name"]],
    "Control File" = df[["Control File"]],
    "Experimental File" = df[["Experimental File"]],
    "Target" = target,
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
  has_cut_position <- nzchar(cut_position)
  valid_cut_position <- has_cut_position && !is.na(suppressWarnings(as.numeric(cut_position)))
  valid_grna <- nzchar(grna) && isTRUE(.valid_nuclease_grna(grna))

  issues <- c(
    if (!nzchar(sample_name)) "Missing sample name",
    if (nzchar(sample_name) && sum(all_samples == sample_name, na.rm = TRUE) > 1) "Duplicate name",
    if (!.fastq_merged_ok(control_file)) "Invalid control FASTQ",
    if (!.fastq_merged_ok(experimental_file)) "Invalid experimental FASTQ",
    if (!nzchar(grna) && !has_cut_position) "Missing gRNA sequence or cut position",
    if (!has_cut_position && nzchar(grna) && !valid_grna) "Invalid gRNA sequence",
    if (has_cut_position && !valid_cut_position) "Invalid cut position",
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
  df <- .normalize_nuclease_row(df)
  samples <- df[["Sample Name"]]
  df$Status <- vapply(seq_len(nrow(df)), function(i) {
    .nuclease_status_for_row(df[i, , drop = FALSE], samples)
  }, character(1))
  df
}
