library(shiny)
library(shinyjs)
library(bslib)
library(later)

# ---------- helpers ----------
.save_managed_upload <- function(upload, role, root) {

  if (is.null(upload) || nrow(upload) == 0) {
    return(NULL)
  }

  managed_name <- paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    role,
    "_",
    basename(upload$name)
  )

  managed_path <- file.path(root, managed_name)

  ok <- file.copy(upload$datapath, managed_path, overwrite = TRUE)

  if (!ok) {
    stop("Failed to copy uploaded file.")
  }

  list(
    role = role,
    original_name = upload$name,
    managed_name = managed_name,
    managed_path = managed_path,
    size = upload$size
  )
}

.remove_managed_upload <- function(meta) {
  if (!is.null(meta) &&
      !is.null(meta$managed_path) &&
      file.exists(meta$managed_path)) {
    unlink(meta$managed_path, force = TRUE)
  }
  invisible(NULL)
}

.remove_dir_safely <- function(path) {
  if (!is.null(path) && dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
  invisible(NULL)
}

.read_batch_manifest <- function(path) {
  df <- data.table::fread(path, check.names = FALSE)
  .validate_batch_manifest(df)
  if (nrow(df) == 0) stop("Batch manifest is empty.")
  df
}

.list_archive_contents <- function(zip_path) {
  utils::unzip(zip_path, list = TRUE)
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

.batch_manifest_with_status <- function(df) {
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

.remove_batch_row <- function(df, row_index) {
  row_index <- as.integer(row_index)

  if (is.na(row_index) || row_index < 1 || row_index > nrow(df)) {
    return(list(df = df, removed = NULL, success = FALSE))
  }

  removed <- df[row_index, , drop = FALSE]
  updated <- df[-row_index, , drop = FALSE]

  list(df = updated, removed = removed, success = TRUE)
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
      Source = "single",
      stringsAsFactors = FALSE
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
    missing_samples <- manifest_df[missing_r1 | missing_r2, "Sample Name", drop = TRUE]
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

.run_merge_pipeline <- function(merge_inputs, settings = merge_settings()){
  # TODO: put logic here
}

# Upload limit 50 MB
options(shiny.maxRequestSize = 50 * 1024^2)

ui <- tagList(
  useShinyjs(),
  tags$head(
    tags$style(HTML("
      .bslib-sidebar-resize-handle {
        display: none !important;
      }

      /* Remove numericInput spinner arrows */

/* Chrome, Safari, Edge, Opera */
input[type=number]::-webkit-outer-spin-button,
input[type=number]::-webkit-inner-spin-button {
  -webkit-appearance: none;
  margin: 0;
}

/* Firefox */
input[type=number] {
  -moz-appearance: textfield;
}

      .sidebar {
        resize: none !important;
      }

      .navbar {
        background: rgba(16, 42, 67, 1);
        backdrop-filter: blur(10px);
        border: none !important;
        box-shadow: 0 2px 10px rgba(0,0,0,0.12);
          border-bottom: 1px solid rgba(219,234,254,0.08);

      }

      .navbar-brand {
        color: white !important;
        font-weight: 700;
        font-size: 22px;
        letter-spacing: 0.4px;
      }



      .navbar-nav .nav-link {
        color: rgba(255,255,255,0.88) !important;
        font-weight: 500;
        transition: color 0.2s ease;
      }

      .navbar-nav .nav-link:hover {
        color: white !important;
      }

      .navbar-nav .nav-link.active {
        color: white !important;
        border-bottom: 2px solid white;
      }

      .hero-banner {
        width: 100%;
        margin: 8px 0 20px 0;
        padding: 28px 40px;
        background: linear-gradient(135deg, #1b4965 0%, #2c7da0 50%, #61a5c2 100%);
        border-radius: 18px;
        box-shadow: 0 6px 20px rgba(0,0,0,0.12), inset 0 1px 0 rgba(255,255,255,0.15);
        border: 1px solid rgba(255,255,255,0.12);
        position: relative;
        overflow: hidden;
      }

      .hero-banner::before {
        content: '';
        position: absolute;
        top: -40%;
        left: -10%;
        width: 50%;
        height: 200%;
        background: rgba(255,255,255,0.08);
        transform: rotate(25deg);
      }

      .hero-content {
        position: relative;
        z-index: 2;
        text-align: center;
      }

      .hero-title {
        color: white;
        font-size: 42px;
        font-weight: 800;
        line-height: 1;
        letter-spacing: -0.5px;
        margin-bottom: 8px;
      }

      .hero-subtitle {
        color: rgba(255,255,255,0.92);
        font-size: 20px;
        font-weight: 500;
        letter-spacing: 0.2px;
      }

      .card-header-custom {
        background: linear-gradient(90deg, #1b4965 0%, #2c7da0 100%) !important;
        color: white !important;
        font-weight: 700;
        font-size: 18px !important;
        border-bottom: none !important;
        padding: 14px 18px;
        border-top-left-radius: 12px;
        border-top-right-radius: 12px;

      }

      .card {
        border-radius: 12px;
        border: none;
        box-shadow: 0 4px 14px rgba(0,0,0,0.08);
        overflow: hidden;
      }

      .card-text-body {
        font-size: 15px;
        line-height: 1.7;
        text-align: justify;
        color: #374151;
      }

      .card-text-body h3 {
        font-size: 18px;
        font-weight: 700;
  margin-top: 0 !important;
  margin-bottom: 0.35rem !important;
        color: #102a43;
      }

      .card-text-body p {
  margin-top: 0 !important;
  margin-bottom: 0.55rem !important;
}

.card-text-body ul {
  margin-top: 0.2rem !important;
  margin-bottom: 0.6rem !important;

  padding-left: 1.4rem;
}
.card-text-body li {
  margin-bottom: 0.2rem !important;
}



      .navbar,
      .navbar.navbar-light,
      .navbar.navbar-dark {
        --bs-navbar-toggler-icon-bg: url(\"data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 30 30'%3e%3cpath stroke='rgba%28255,255,255,1%29' stroke-linecap='round' stroke-miterlimit='10' stroke-width='2.5' d='M4 7h22M4 15h22M4 23h22'/%3e%3c/svg%3e\");
        --bs-navbar-toggler-border-color: rgba(255,255,255,0.35);
      }

      .navbar-toggler {
        border: 1px solid rgba(255,255,255,0.35) !important;
        background-color: rgba(255,255,255,0.06) !important;
        border-radius: 10px;
        padding: 0.35rem 0.55rem;
      }

      .navbar-toggler:hover {
        background-color: rgba(255,255,255,0.12) !important;
        border-color: rgba(255,255,255,0.6) !important;
      }

      .navbar-toggler:focus {
        box-shadow: 0 0 0 0.15rem rgba(255,255,255,0.18) !important;
      }

      html, body {
        height: 100%;
      }


    .app-footer {
      width: 100%;
      background: #102a43;
      color: rgba(255,255,255,0.88);

      padding: 12px 24px;

      border-top: 1px solid rgba(255,255,255,0.88);
#border-top: 1px solid #DBEAFE;

      font-size: 14px;
      font-weight: 500;
    }

    .footer-content {
      display: flex;

      justify-content: space-between;

      align-items: center;

      gap: 20px;

      flex-wrap: wrap;
    }

    .footer-left {
      text-align: left;
    }

    .footer-center {
      text-align: center;
      flex-grow: 1;
    }

    .footer-right {
      text-align: right;
    }

        .navbar,
    .navbar .dropdown,
    .navbar .dropdown-menu {
      z-index: 2000 !important;
    }

    .navbar {
      position: relative;
    }

    .navbar .dropdown-menu {
      border: 1px solid rgba(255,255,255,0.12);
      box-shadow: 0 8px 24px rgba(0,0,0,0.25);
    }

    .navbar .dropdown-item {
      color: rgba(255,255,255,0.9) !important;
    }

    .navbar .dropdown-item:hover,
    .navbar .dropdown-item:focus {
      background: rgba(255,255,255,0.08);
      color: white !important;
    }

    .card-text-body {
      font-size: 15px;
      line-height: 1.7;
      color: #374151;
    }

    .card-text-body .doc-block h3 {
      margin: 0 0 0.25rem 0;
      font-size: 18px;
      font-weight: 700;
      color: #102a43;
    }

    .card-text-body .doc-block p {
      margin: 0 0 0.4rem 0;
    }

    .card-text-body .doc-block ul {
      margin: 0.1rem 0 0.7rem 0;
      padding-left: 1.3rem;
    }

    .card-text-body .doc-block li {
      margin: 0 0 0.2rem 0;
    }

        .footnote-divider {
      width: 100%;
      height: 1px;

      background: rgba(0,0,0,0.12);

      margin: 1rem 0 0.7rem 0;
    }

    .footnotes {
      font-size: 0.82rem;
      color: #6b7280;

      line-height: 1.4;
    }

    .footnotes p {
      margin: 0 0 0.25rem 0;
    }

    .footnotes sup {
      font-weight: 700;
      color: #374151;
    }

        .sidebar {
      background: linear-gradient(180deg, #0f2238 0%, #132b45 100%);
      color: rgba(255,255,255,0.92);
      border-right: 1px solid rgba(255,255,255,0.08);

      box-shadow: inset -1px 0 0 rgba(255,255,255,0.03);
    }

    .sidebar .sidebar-title {
      color: white;
      font-weight: 700;
      letter-spacing: 0.2px;
    }

    .sidebar .form-label,
    .sidebar p,
    .sidebar li {
      color: rgba(255,255,255,0.86);
    }

    .sidebar .text-muted {
      color: rgba(255,255,255,0.65) !important;
    }

    .sidebar hr {
      border-top: 1px solid rgba(255,255,255,0.12);
      opacity: 1;
    }

     .sidebar-checklist {
      background: rgba(255,255,255,0.06);
      border: 1px solid rgba(255,255,255,0.12);
      border-radius: 14px;
      padding: 14px 16px;
      margin-top: 8px;
     }


    .sidebar-checklist h5 {
      margin-bottom: 0.5rem;
      color: white;
      font-weight: 700;
    }

    .sidebar-checklist p,
    .sidebar-checklist li {
      color: rgba(255,255,255,0.88);
      margin-bottom: 0.35rem;
    }

    .sidebar-checklist ul {
      padding-left: 1.2rem;
      margin-bottom: 0;
    }

    .sidebar-checklist hr {
      border-top: 1px solid rgba(255,255,255,0.12);
      opacity: 1;
      margin: 0.7rem 0;
    }

    .sidebar-controls h5 {
  color: white;
  font-weight: 700;
}

.sidebar-controls label,
.sidebar-controls .form-label,
.sidebar-controls p,
.sidebar-controls li,
.sidebar-controls .control-label {
  color: rgba(255,255,255,0.88) !important;
}

.sidebar-controls .text-muted {
  color: rgba(255,255,255,0.65) !important;
}

.sidebar-controls hr {
  border: none !important;
  height: 1px !important;
  background-color: rgba(255,255,255,0.5) !important;
  opacity: 1 !important;
  margin: 14px 0 !important;
}

.compact-hr {
  border: none !important;
  height: 1px !important;
  background-color: rgba(255,255,255,0.28) !important;
  opacity: 1 !important;
  margin: 8px 0 !important;
}

.section-hr {
   border: none !important;
  height: 1px !important;
  background-color: rgba(255,255,255,0.5) !important;
  opacity: 1 !important;
  margin: 14px 0 !important;
}


/* --- Merge sidebar text + widget overrides --- */
.sidebar-controls,
.sidebar-controls label,
.sidebar-controls .control-label,
.sidebar-controls h5,
.sidebar-controls p,
.sidebar-controls li,
.sidebar-controls .text-muted {
  color: rgba(255,255,255,0.92) !important;
}

.sidebar-controls hr {
  border-top: 1px solid rgba(255,255,255,0.12) !important;
  opacity: 1 !important;
}

/* Select input */
.sidebar-controls .selectize-control.single .selectize-input,
.sidebar-controls .selectize-control.single .selectize-input input,
.sidebar-controls .selectize-control.single .selectize-dropdown {
  color: #102a43 !important;
}

.sidebar-controls .selectize-control.single .selectize-input {
  background: #ffffff !important;
  border: 1px solid rgba(255,255,255,0.35) !important;
  box-shadow: none !important;
}

.sidebar-controls .selectize-dropdown {
  background: #ffffff !important;
  color: #102a43 !important;
}

/* Numeric input */
.sidebar-controls .form-control {
  color: #102a43 !important;
  background: #ffffff !important;
  border: 1px solid rgba(255,255,255,0.35) !important;
  box-shadow: none !important;
}

.sidebar-controls .input-group-text,
.sidebar-controls .spinner {
  color: #102a43 !important;
  background: #ffffff !important;
  border-color: rgba(255,255,255,0.35) !important;
}

/* Checkbox labels */
.sidebar-controls .shiny-options-group label,
.sidebar-controls .checkbox label,
.sidebar-controls .checkbox-inline label {
  color: rgba(255,255,255,0.92) !important;
}

/* Slider labels / ticks */
.sidebar-controls .irs-min,
.sidebar-controls .irs-max,
.sidebar-controls .irs-single,
.sidebar-controls .irs-grid-text {
  color: rgba(255,255,255,0.92) !important;
}

.sidebar-controls .irs-line {
  background: rgba(255,255,255,0.18) !important;
}

.sidebar-controls .irs-bar {
  background: #1188d6 !important;
}

.sidebar-controls .irs-handle > i:first-child {
  background: #1188d6 !important;
}

/* Run status box from verbatimTextOutput */
.sidebar-controls pre,
.sidebar-controls .shiny-text-output,
.sidebar-controls .shiny-text-output.form-control {
  color: rgba(255,255,255,0.95) !important;
  background: rgba(255,255,255,0.04) !important;
  border: 1px solid rgba(255,255,255,0.45) !important;
  border-radius: 8px !important;
}

/* Keep buttons readable */
.sidebar-controls .btn-primary {
  width: 100%;
}

    .module-header {
  margin-bottom: 18px;
  padding-bottom: 10px;

  border-bottom: 1px solid rgba(219,234,254,0.18);
}

.module-header h1 {
  font-size: 32px;
  font-weight: 750;
  color: #102a43;

  margin-bottom: 0.2rem;
}

.module-header p {
  font-size: 16px;
  color: #6b7280;

  margin-bottom: 0;
}

.minimal-input-group {
      width: 100%;
      margin-bottom: 18px;
    }

    .minimal-input-group .form-control {

      border: none !important;
      border-bottom: 2px solid rgba(16,42,67,0.25) !important;

      border-radius: 0 !important;

      background: transparent !important;

      box-shadow: none !important;

      padding-left: 0;
      padding-right: 0;

      font-size: 16px;

      color: #102A43;
    }
    .shiny-input-container {
  margin-bottom: 0 !important;
}

    .minimal-input-group .form-control:focus {

      border-bottom: 2px solid #2C7DA0 !important;

      box-shadow: none !important;
    }

    .minimal-input-group .form-control::placeholder {

      color: #9CA3AF;
      opacity: 1;

      font-style: italic;
    }

    .minimal-input-label-text {
      margin-top: 2px;
      font-size: 12px;
      color: #6B7280;
      letter-spacing: 0.2px;
    }

    /* keep only one version of this rule */
.minimal-input-label-fileInput {
  display: block;
  margin-top: -10px;     /* pull label closer to the line */
  margin-bottom: 0;
  font-size: 12px;
  color: #6B7280;
  line-height: 1.1;
}

/* remove default spacing from Shiny's fileInput wrapper */
.minimal-file-group .shiny-input-container,
.minimal-file-group .form-group {
  margin-bottom: 0 !important;
}

.minimal-file-group .input-group {
  margin-bottom: 14 !important;
}
.minimal-file-group .shiny-input-container {
  margin-bottom: 0 !important;
  padding-bottom: 18px;  /* room so the warning does not collide */
}

.minimal-file-group .minimal-input-label {
  margin-top: 4px;       /* use positive spacing, not negative */
}

        .minimal-file-group {
      width: 100%;
      margin-bottom: 0px;
    }

    # /* Remove bulky outer spacing */
    # .minimal-file-group .shiny-input-container {
    #   margin-bottom: 0 !important;
    # }




    /* Main file input area */
    .minimal-file-group .form-control {
      border: none !important;
      border-bottom: 2px solid rgba(16,42,67,0.25) !important;

      border-radius: 0 !important;

      background: transparent !important;

      box-shadow: none !important;

      padding-left: 0;
      padding-right: 0;

      color: #102A43;
    }

    /* File text */
    .minimal-file-group input[type='text'] {
      background: transparent !important;
      border: none !important;
      box-shadow: none !important;

      color: #102A43;
      font-size: 15px;
    }

    /* Browse button */
    .minimal-file-group .btn {
      background: transparent !important;

      color: #2C7DA0 !important;

      border: none !important;

      border-radius: 0 !important;

      font-weight: 600;

      padding-left: 0;
      padding-right: 0;
        margin-right: 12px !important;


      box-shadow: none !important;
    }

    .minimal-file-group .btn:hover {
      color: #1B4965 !important;
      background: transparent !important;
    }

    /* Remove rounded input group */
    .minimal-file-group .input-group {
      border-bottom: 2px solid rgba(16,42,67,0.25);
    }

    .minimal-file-group .input-group:focus-within {
      border-bottom: 2px solid #2C7DA0;
    }

.upload-action-row {
  display: flex;

  gap: 12px;

  justify-content: flex-end;

  margin-top: -43px;

  margin-bottom: 12px;

  max-width: 760px;
}

/* Main action button */
.primary-action-btn {
  background: linear-gradient(90deg, #1B4965 0%, #2C7DA0 100%) !important;

  color: white !important;

  border: none !important;

  border-radius: 10px !important;

  padding: 10px 20px !important;

  font-weight: 650 !important;

  letter-spacing: 0.2px;

  box-shadow: 0 4px 12px rgba(44,125,160,0.20);

  transition: all 0.18s ease;
}

.primary-action-btn:hover {
  //transform: translateY(-1px);

  box-shadow: 0 6px 16px rgba(44,125,160,0.28);

  background: linear-gradient(90deg, #163B52 0%, #256B8A 100%) !important;
}

/* Secondary button */
.secondary-action-btn {
  background: rgba(16,42,67,0.06) !important;

  color: #1B4965 !important;

  border: 1px solid rgba(16,42,67,0.12) !important;

  border-radius: 10px !important;

  padding: 10px 20px !important;

  font-weight: 650 !important;

  transition: all 0.18s ease;
}

.secondary-action-btn:hover {
  background: rgba(16,42,67,0.10) !important;

  color: #102A43 !important;
}

.queue-table-wrapper {
  margin-top: 18px;
  max-width: 760px;
}

.queue-table-wrapper table {
  width: 100%;

  border-collapse: collapse;

  font-size: 14px;
    border-radius: 10px;

  overflow: hidden;
}

.queue-table-wrapper th {
  vertical-align: middle !important;
  text-align: left;

   padding-top: 8px;
  padding-bottom: 8px;

  background: rgba(16,42,67,0.08);

  color: #102A43;

  font-weight: 700;

  border-bottom: 1px solid rgba(16,42,67,0.20);
}

.queue-table-wrapper td {
  vertical-align: middle !important;
   padding-top: 8px;
  padding-bottom: 8px;

  border-bottom: 1px solid rgba(16,42,67,0.2);
}



.queue-alert {
  border-radius: 10px;

  padding: 12px 16px;

  // margin-top: 12px;
  // margin-bottom: 8px;

  margin: 0;
  font-size: 14px;

  font-weight: 700;
  width: auto;
  max-width: fit-content;

  transition: opacity 5s ease, transform 5s ease;

}

/* Error */
.queue-alert-error {
  background: rgba(220,38,38,0.08);

  border: 1px solid rgba(220,38,38,0.16);

  color: #991B1B;
}

/* Success */
.queue-alert-success {
  background: rgba(22,163,74,0.08);

  border: 1px solid rgba(22,163,74,0.16);

  color: #166534;
}

/* Info */
.queue-alert-info {
  background: rgba(37,99,235,0.08);

  border: 1px solid rgba(37,99,235,0.16);

  color: #1D4ED8;
}

.message-slot {
  min-height: 52px;

  display: flex;

  align-items: center;

margin-top: 26px;
}

/* Custom slider tick labels */


.custom-slider-ticks {
  position: relative;
  height: 24px;
  margin-top: 0px;
  margin-left: 0;
  margin-right: 0;
  font-size: 12px;
  user-select: none;
}

.custom-slider-ticks span {
  position: absolute;
  top: 0;
  transform: translateX(-50%);
  color: rgba(255,255,255,0.82);
  white-space: nowrap;
}

.custom-slider-ticks span::before {
  content: '';
  position: absolute;
  top: -12px;
  left: 50%;
  transform: translateX(-50%);
  width: 1px;
  height: 8px;
  background: rgba(255,255,255,0.35);
}

.slider-input-row {
  display: flex;
  align-items: flex-start;
  gap: 14px;
  width: 100%;
  margin-bottom: 14px;
}

.numeric-col {
  width: 35px;
  margin-top: -2px;
}

.slider-col {
  flex: 1;
}


/* Remove extra bottom spacing */
.slider-col .shiny-input-container,
.numeric-col .shiny-input-container {
  margin-bottom: 0 !important;
}

#merge_similarity_num,
#min_overlap_num,
#min_length_num,
#max_length_num,
#merge_prefix {
  height: 30px;
  padding: 2px 6px;
  font-size: 12px;
}
.sidebar-controls .control-label,
.sidebar-controls .form-label {
  font-size: 14px !important;

  font-weight: 600;

  margin-bottom: 8px;
}

.sidebar-controls .checkbox label {
  font-size: 14px !important;
}

.sidebar-controls .selectize-input,
.sidebar-controls .selectize-dropdown,
.sidebar-controls .selectize-dropdown-content {
  font-size: 12px !important;
}

.sidebar-controls .selectize-input {
  min-height: 30px !important;

  padding-top: 4px !important;
  padding-bottom: 4px !important;
}

.sidebar-controls .selectize-control.single .selectize-input {
  height: 32px !important;

  min-height: 32px !important;

  max-height: 32px !important;

  padding-top: 4px !important;
  padding-bottom: 4px !important;

  display: flex !important;

  align-items: center !important;

  overflow: hidden !important;
}

.sidebar-controls .selectize-control {
  height: 32px !important;
}

.sidebar-controls .irs-min,
.sidebar-controls .irs-max {
  display: none !important;
}

.sidebar-controls .irs-from,
.sidebar-controls .irs-to {
  display: none !important;
}

.workflow-tabs {
  margin-top: 4px;
}


.workflow-tabs .nav-tabs {
  border-bottom: 1px solid rgba(16, 42, 67, 0.14);

  margin-bottom: 18px;

  position: relative;
}

.workflow-tabs .nav-link {
  color: #1B4965;
  font-weight: 650;

  border-radius: 12px 12px 0 0 !important;

  padding: 12px 22px;

  transition: all 0.18s ease;
}

.workflow-tabs .nav-link.active {
  color: #102A43;

  background-color: #ffffff;

  border-color:
    rgba(16, 42, 67, 0.14)
    rgba(16, 42, 67, 0.14)
    #ffffff;

  box-shadow: 0 -1px 0 rgba(255,255,255,0.65) inset;

  position: relative;

  top: 1px;
}

.workflow-tabs .nav-link:hover {
  color: #102A43;

  background: rgba(44,125,160,0.04);

  border-color: transparent;
}

.workflow-tabs .tab-content {
  padding-top: 8px;
}

.paired-upload-row {
  display: inline-flex;
  gap: 50px;
  align-items: flex-start;
  flex-wrap: wrap;
  width: fit-content;
  max-width: 100%;
}

.paired-upload-col {
  //flex: 0 0 200px;
  //width: 200px;
  max-width: 500px;
}

.upload-card-shell {
  max-width: 900px;

  width: 100%;
}

.workflow-tabs .tab-content {
  padding-top: 18px;
}

.batch-section-help {
  font-size: 13px;
  color: #6B7280;
  margin-bottom: 10px;
}

.minimal-input-label-row {
  display: flex;

  justify-content: space-between;

  align-items: baseline;

  width: 100%;

  margin-top: -10px;
}



.template-download-link {
  font-size: 12px;

  font-weight: 600;

  color: #2C7DA0;

  text-decoration: none;

  transition: color 0.18s ease;
}

.template-download-link:hover {
  color: #1B4965;

  text-decoration: underline;
}

.queue-table td,
.queue-table th {
  vertical-align: middle !important;
  padding-top: 8px;
  padding-bottom: 8px;
  line-height: 1.2;
}

.queue-status-cell {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 12px;
  width: 100%;
}

.queue-status-text {
  flex: 1 1 auto;
  min-width: 0;
}

.queue-remove-btn {
  flex: 0 0 auto;
  border: none;
  background: transparent;
  color: #b91c1c;
  font-weight: 500;
  padding: 0;
  line-height: 1;
  cursor: pointer;
  white-space: nowrap;
}

.queue-remove-btn:hover {
  color: #7f1d1d;
  text-decoration: underline;
}

#forward_primer,
#reverse_primer {
  font-size: 12px !important;
  height: 30px;
  padding: 2px 6px;
}

.sidebar-controls .checkbox {
  margin-bottom: 0 !important;
}

.sidebar-controls .checkbox input {
  margin-top: 0 !important;
}


.sidebar-inline-label {
  font-size: 14px;
  font-weight: 500;
  color: rgba(255,255,255,0.92);
}

.sidebar-controls .checkbox {
  margin-bottom: 0 !important;
}

.sidebar-controls .checkbox input {
  margin-top: 0 !important;
}

.sidebar-checkbox-row {
  display: grid;
  grid-template-columns: auto 1fr auto;
  align-items: center;
  column-gap: 10px;
  margin-bottom: 14px;
}

.sidebar-checkbox-row .checkbox {
  margin: 0 !important;
  padding: 0 !important;
  line-height: 1 !important;
}

.sidebar-checkbox-label {
  font-size: 14px;
  font-weight: 600;
  line-height: 1;
  margin: 0;
  color: rgba(255,255,255,0.92);
  transform: translateY(-2px);
}

.sidebar-checkbox-row .checkbox input[type='checkbox'] {
  margin-top: 0 !important;
}



    "))
  ),
  page_navbar(title = "tasaR",
              id = "main_nav",
              theme = bs_theme(version = 5),
              # potential sidebar background colors or additional accent colors: #DBEAFE, #89C2D9
              sidebar = sidebar(
                title = "Status",
                open = "always",
                resizable = FALSE,
                width = "325px",
                uiOutput("sidebar_ui")
              ),
             nav_panel(title = "",
                       icon = icon("home"),
                      value = "splash.app",
                      fluidRow(
                        div(
                          class = "hero-banner",

                          div(
                            class = "hero-content",

                            div(
                              class = "hero-title",
                              "tasaR"
                            ),

                            div(
                              class = "hero-subtitle",
                              "Toolbox for Amplicon Sequencing Analysis in R"
                            )
                          )
                        )
                      ),
                      fluidRow(card(
                        card_header(
                          "What is tasaR?",
                          class = "card-header-custom"
                        ),
                        card_body(
                          class = "card-text-body",
                          HTML("
    <div class='doc-block'>
      <h3>Overview</h3>
      <div style='height: 4px;'></div>
      <p>
        tasaR is an R package built to support the analysis of Illumina targeted amplicon sequencing experiments.
        This software performs in-depth analysis of single amplicons across multiple modules. tasaR provides 3 key functionalities:
      </p>
      <div style='height: 8px;'></div>
      <ul>
        <li>
          <strong>Merge Reads:</strong> An integrated workflow to merge paired end reads using the software PANDAseq.<sup>1</sup>
        </li>
        <li>
          <strong>Measure Mutations:</strong> Analyzes target sequences and quantifies mutagenesis across the reference sequence.
          Somatic hypermutation is measured by annotating AID hotspot cytosines and quantifying mutation specifically at those sites.
        </li>
        <li>
          <strong>Nuclease Analysis:</strong> Quantifies the DNA mutations created by a targeted nuclease at a specific site in the genome.
          Uses a matched control sequence to determine the reference sequence; supports heterozygous genomic DNA.
        </li>
      </ul>

      <div style='height: 12px;'></div>

      <h3>Features</h3>
      <div style='height: 4px;'></div>
      <ul>
        <li>Automated quality-aware merging of paired end reads.</li>
        <li>Quantify DNA and protein mutations across the sequence.</li>
        <li>Annotate AID hotspot motifs and quantify somatic hypermutation.</li>
        <li>Analyze indels and base change events around a nuclease cut site.</li>
        <li>Predict the usage of DNA repair pathways after nuclease-induced DNA scarring.<sup>2</sup></li>
        <li>Supports unique molecular identifier (UMI) tagging and barcodes to enhance sequence de-multiplexing.</li>
      </ul>
      <div style='height: 4px;'></div>


          <div class='footnote-divider'></div>

    <div class='footnotes'>

      <p>
        <sup>1</sup>
        Masella AP, et al. (2012). PANDAseq: paired-end assembler for Illumina sequences. BMC Bioinformatics 14:13:31.
      </p>

      <p>
        <sup>2</sup>
        Tatiossian KJ, et al. (2021). Rational Selection of CRISPR-Cas9 Guide RNAs for Homology-Directed Genome Editing. Mol Ther 3;29(3):1057-1069.
      </p>

    </div>

    </div>
  ")

                      )))
                      ),
             nav_panel("Merge Reads",
                       value = "merge.app",
                       div(
                         class = "module-header",
                         h1("Merge Reads"),
                         p("Merge paired-end FASTQ reads using integrated workflows.")
                       ),
                       fluidRow(
                         div(
                           class = "upload-card-shell",
                         card(
                           card_header(
                             "Upload Files",
                             class = "card-header-custom"
                           ),
                           card_body(
                             div(
                               class = "workflow-tabs",
                               navset_tab(
                                 id = "upload_mode",
                                 selected = "single_sample",

                                 nav_panel(
                                   title = "Single Sample",
                                   value = "single_sample",

                                   div(
                                     class = "card-text-body",
                                     div(
                                       class = "minimal-input-group",
                                       style = "width: 300px; max-width: 100%;",
                                       textInput(
                                         inputId = "sample_name",
                                         label = NULL,
                                         placeholder = "Enter sample name"
                                       ),
                                       div(
                                         class = "minimal-input-label-text",
                                         "Sample identifier"
                                       )
                                     ),
                                     div(
                                       class = "paired-upload-row",

                                       div(
                                         class = "paired-upload-col",
                                         style = "width: 200px; max-width: 100%;",

                                         div(
                                           class = "minimal-file-group",
                                           fileInput(
                                             inputId = "fastq_r1",
                                             label = NULL,
                                             buttonLabel = "Browse",
                                             placeholder = "No file selected"
                                           ),
                                           div(
                                             class = "minimal-input-label-fileInput",
                                             "Left (R1) FASTQ file"
                                           )
                                         )
                                       ),

                                       div(
                                         class = "paired-upload-col",
                                         style = "width: 200px; max-width: 100%;",

                                         div(
                                           class = "minimal-file-group",
                                           fileInput(
                                             inputId = "fastq_r2",
                                             label = NULL,
                                             buttonLabel = "Browse",
                                             placeholder = "No file selected"
                                           ),
                                           div(
                                             class = "minimal-input-label-fileInput",
                                             "Right (R2) FASTQ file"
                                           )
                                         )
                                       )
                                     ),
                                     div(
                                       class = "message-slot",
                                       uiOutput("merge_message")
                                     ),
                                     div(
                                       class = "upload-action-row",
                                       actionButton(
                                         "queue_sample",
                                         "Add Sample",
                                         class = "secondary-action-btn"
                                       ),
                                       actionButton(
                                         "run_single_merge",
                                         "Start Merge",
                                         class = "primary-action-btn"
                                       )
                                     ),
                                     br(),
                                     div(
                                       class = "queue-table-wrapper",
                                       uiOutput("queued_samples")
                                     )
                                   )
                                 ),
                                 nav_panel(
                                   title = "Batch Import",
                                   value = "batch_import",

                                   div(
                                     class = "card-text-body",

                                     div(
                                       class = "paired-upload-row",

                                       div(
                                         class = "paired-upload-col",
                                         style = "width: 300px; max-width: 100%;",

                                         div(
                                           class = "minimal-file-group",

                                           fileInput(
                                             inputId = "batch_manifest",
                                             label = NULL,
                                             buttonLabel = "Browse",
                                             placeholder = "No file selected",
                                             accept = c(".csv"),
                                             width = "100%"
                                           ),

                                           div(
                                             class = "minimal-input-label-row",

                                             div(
                                               class = "minimal-input-label-fileInput",
                                               "Sample list (.csv)"
                                             ),

                                             tags$a(
                                               href = "samples_to_merge.csv",
                                               download = NA,
                                               class = "template-download-link",
                                               "Download template"
                                             )
                                           )
                                         )
                                       ),

                                       div(
                                         class = "paired-upload-col",
                                         style = "width: 300px; max-width: 100%;",

                                         div(
                                           class = "minimal-file-group",

                                           fileInput(
                                             inputId = "batch_archive",
                                             label = NULL,
                                             buttonLabel = "Browse",
                                             placeholder = "No file selected",
                                             accept = c(".zip"),
                                             width = "100%"
                                           ),

                                           div(
                                             class = "minimal-input-label-fileInput",
                                             "Compressed FASTQ archive (.zip)"
                                           )
                                         )
                                       )
                                     ),

                                     div(
                                       class = "queue-table-wrapper",
                                       uiOutput("batch_queued_samples")
                                         )
                                       ),
                                     div(style = "height: 55px;"),

                                     div(
                                       class = "upload-action-row",

                                       actionButton(
                                         "clear_batch",
                                         "Clear Files",
                                         class = "secondary-action-btn"
                                       ),

                                       actionButton(
                                         "run_batch_merge",
                                         "Start Merge",
                                         class = "primary-action-btn"
                                       )
                                     )
                                   )
                                 )
                               )
                             )
                           )
                         )
                       )),
             nav_panel("Measure Mutations",
                      value = "shm.app"
                      ),
             nav_panel("Nuclease Analysis",
                      value = "cut.app"
                      ),
             nav_spacer(),
             nav_panel("Help",
                      value = "help.app",
                      class = "ms-auto"
                      ),
             nav_menu(title = "Links",
                      align = "right",
                      value = "links.app",
                      nav_item(
                        tags$a(
                          icon("circle-question"),
                          href = "",
                          "About",
                          target = "_blank"
                        )
                      ),
                      nav_item(
                        tags$a(
                          icon("github"),
                          href = "https://github.com/rogersgl/tasaR",
                          "tasaR GitHub",
                          target = "_blank"
                        )
                        ),
                      nav_item(
                        tags$a(
                          icon("docker"),
                          href = "https://TBD_WHEN_DONE",
                          "tasaR Docker",
                          target = "_blank"
                        )
                      ),
                      nav_item(
                        tags$a(
                          icon("github"),
                          href = "https://github.com/neufeld/pandaseq",
                          "PANDAseq GitHub",
                          target = "_blank"
                        )
                      ),
                      nav_item(
                        tags$a(
                          icon("pen-to-square"),
                          href = "",
                          "GPLv3 License",
                          target = "_blank"
                        )
                      ),
                      nav_item(
                        tags$a(
                          icon("user-pen"),
                          href = "",
                          "About the author",
                          target = "_blank"
                        )
                      )
             )),
  div(
    class = "app-footer",

    div(
      class = "footer-content",

      div(
        class = "footer-left",
        "Built with R • Shiny • Bioconductor"
      ),

      div(
        class = "footer-center",
        "Version 0.2.0"
      ),

      div(
        class = "footer-right",
        "tasaR © 2026 Geoffrey L. Rogers"
      )
    )
  )
)

server <- function(input, output, session) {
  output$sidebar_ui <- renderUI({

    switch(
      input$main_nav,

      "splash.app" = div(
        class = "sidebar-checklist",

        h5("Start here"),
        tags$p("Use the top navigation to choose a workflow."),
        tags$hr(class = "compact-hr"),

        tags$ul(
          tags$li(tags$strong("Review"), " the overview section"),
          tags$li(tags$strong("Choose"), " the module you need"),
          tags$li(tags$strong("Load"), " your input data")
        ),

        tags$hr(class = "compact-hr"),
        tags$p(
          class = "text-muted",
          "Supports targeted amplicon sequencing and genome editing analysis."
        )
      ),

      "merge.app" = div(
        class = "sidebar-controls",

        tags$hr(class = "section-hr"),

        h5("Run Settings"),

        verbatimTextOutput("merge_settings"),

        tags$hr(class = "section-hr"),

        h5("Merge Controls"),

        div(style = "height: 8px;"),

        div(
          class = "d-flex align-items-center justify-content-between",
          tags$label("Output prefix", class = "control-label"),
          tooltip(
            icon("circle-question"),
            "Prefix to append to merged file names."
          )
        ),

        textInput(
          "merge_prefix",
          label = NULL,
          value = "merge"
        ),

       br(),

       div(
         class = "d-flex align-items-center justify-content-between",
         tags$label("Merge engine", class = "control-label"),
         tooltip(
           icon("circle-question"),
           "Select the algorithm used to merge paired-end reads. Currently only supports PANDAseq."
         )
       ),

        selectInput(
          "merge_engine",
          label = NULL,
          choices = c(
            "PANDAseq"
          ),
          selected = "PANDAseq"
        ),

        br(),

       div(
         class = "d-flex align-items-center justify-content-between",
         tags$label("Minimum overlap (bp)", class = "control-label"),
         tooltip(
           icon("circle-question"),
           "Minimum number of base pairs that must overlap in the middle for reads to be merged."
         )
       ),

        div(
          class = "slider-input-row",
          div(
            class = "slider-col",
            sliderInput(
              "min_overlap",
              label = NULL,
              min = 5,
              max = 100,
              value = 20,
              step = 5,
              ticks = FALSE
            ),
            div(
              class = "custom-slider-ticks",
              span(style = "left: calc(10px + ((10 - 5) * (100% - 20px) / 95));",  "10"),
              span(style = "left: calc(10px + ((20 - 5) * (100% - 20px) / 95));",  "20"),
              span(style = "left: calc(10px + ((30 - 5) * (100% - 20px) / 95));",  "30"),
              span(style = "left: calc(10px + ((40 - 5) * (100% - 20px) / 95));",  "40"),
              span(style = "left: calc(10px + ((50 - 5) * (100% - 20px) / 95));",  "50"),
              span(style = "left: calc(10px + ((60 - 5) * (100% - 20px) / 95));",  "60"),
              span(style = "left: calc(10px + ((70 - 5) * (100% - 20px) / 95));",  "70"),
              span(style = "left: calc(10px + ((80 - 5) * (100% - 20px) / 95));",  "80"),
              span(style = "left: calc(10px + ((90 - 5) * (100% - 20px) / 95));",  "90"),
              span(style = "left: calc(10px + ((100 - 5) * (100% - 20px) / 95));", "100")
            )
          ),
          div(
            class = "numeric-col",
            numericInput(inputId = "min_overlap_num",
                         label = "",
                         min = 5,
                         max = 100,
                         value = 20)
          )
        ),

        br(),

       div(
         class = "d-flex align-items-center justify-content-between",
         tags$label("Minimum merged length (bp)", class = "control-label"),
         tooltip(
           icon("circle-question"),
           "Minimum total sequence length for a merged read to be valid. Suggested default amplicon size -25."
         )
       ),

        div(
          class = "slider-input-row",
          div(
            class = "slider-col",
            sliderInput(
              "min_length",
              label = NULL,
              min = 50,
              max = 500,
              value = 250,
              step = 25,
              ticks = FALSE
            ),
            div(
              class = "custom-slider-ticks",
              span(style = "left: calc(10px + ((100 - 50) * (100% - 20px) / 450));", "100"),
              span(style = "left: calc(10px + ((200 - 50) * (100% - 20px) / 450));", "200"),
              span(style = "left: calc(10px + ((300 - 50) * (100% - 20px) / 450));", "300"),
              span(style = "left: calc(10px + ((400 - 50) * (100% - 20px) / 450));", "400"),
              span(style = "left: calc(10px + ((500 - 50) * (100% - 20px) / 450));", "500")
            )
          ),
          div(
            class = "numeric-col",
            numericInput(inputId = "min_length_num",
                         label = "",
                         min = 50,
                         max = 500,
                         value = 250)
          )
        ),

        br(),

       div(
         class = "d-flex align-items-center justify-content-between",
         tags$label("Maximum merged length (bp)", class = "control-label"),
         tooltip(
           icon("circle-question"),
           "Max total sequence length for a merged read to be valid. Suggested default amplicon size +25."
         )
       ),

        div(
          class = "slider-input-row",
          div(
            class = "slider-col",
            sliderInput(
              "max_length",
              label = NULL,
              min = 50,
              max = 500,
              value = 475,
              step = 25,
              ticks = FALSE
            ),
            div(
              class = "custom-slider-ticks",
              span(style = "left: calc(10px + ((100 - 50) * (100% - 20px) / 450));", "100"),
              span(style = "left: calc(10px + ((200 - 50) * (100% - 20px) / 450));", "200"),
              span(style = "left: calc(10px + ((300 - 50) * (100% - 20px) / 450));", "300"),
              span(style = "left: calc(10px + ((400 - 50) * (100% - 20px) / 450));", "400"),
              span(style = "left: calc(10px + ((500 - 50) * (100% - 20px) / 450));", "500")
            )
          ),
          div(
            class = "numeric-col",
            numericInput(inputId = "max_length_num",
                         label = "",
                         min = 50,
                         max = 500,
                         value = 475)
          )
        ),

        br(),

       div(
         class = "d-flex align-items-center justify-content-between",
         tags$label("Overlap threshold", class = "control-label"),
         tooltip(
           icon("circle-question"),
           "Minimum accepted value for calculated overlap probability. Default = 0.6. Do not suggest adjusting this without a good reason."
         )
       ),

        div(
          class = "slider-input-row",
          div(
            class = "slider-col",
            sliderInput(
              "merge_similarity",
              label = NULL,
              min = 0.1,
              max = 1,
              value = 0.6,
              step = 0.1,
              ticks = FALSE
            ),
            div(
              class = "custom-slider-ticks",
              span(style = "left: calc(10px + ((0.2 - 0.1) * (100% - 20px) / 0.9));", "0.2"),
              span(style = "left: calc(10px + ((0.4 - 0.1) * (100% - 20px) / 0.9));", "0.4"),
              span(style = "left: calc(10px + ((0.6 - 0.1) * (100% - 20px) / 0.9));", "0.6"),
              span(style = "left: calc(10px + ((0.8 - 0.1) * (100% - 20px) / 0.9));", "0.8"),
              span(style = "left: calc(10px + ((1 - 0.1) * (100% - 20px) / 0.9));", "1.0")
            )
          ),
          div(
            class = "numeric-col",
            numericInput(inputId = "merge_similarity_num",
                         label = "",
                         min = 0,
                         max = 1,
                         value = 0.6)
          )
        ),



       tags$hr(class = "section-hr"),

       div(
         class = "sidebar-checkbox-row",

         checkboxInput(
           "merge_advanced_settings",
           label = NULL,
           value = FALSE
         ),

         tags$span(
           "Advanced settings",
           class = "sidebar-checkbox-label"
         ),

         tooltip(
           icon("circle-question"),
           "Show additional advanced settings for PANDAseq."
         )
       ),

        conditionalPanel(
          condition = "input.merge_advanced_settings === true",
          br(),
                 div(
         class = "d-flex align-items-center justify-content-between",
         tags$label("Forward primer", class = "control-label"),
         tooltip(
           icon("circle-question"),
           "Sequence or length in bp of forward primer to trim. Will also remove any distal adapters (barcodes, UMIs, etc.)"
         )
       ),
          textInput(
            "forward_primer",
            label = NULL,
            placeholder = "Sequence or trim length"
          ),
          div(style = "height: 10px;"),
       div(
         class = "d-flex align-items-center justify-content-between",
         tags$label("Reverse primer", class = "control-label"),
         tooltip(
           icon("circle-question"),
           "Sequence or length in bp of reverse primer to trim. Will also remove any distal adapters (barcodes, UMIs, etc.)"
         )
       ),
          textInput(
            "reverse_primer",
            label = NULL,
            placeholder = "Sequence or trim length"
          ),
       br(),
       div(
         class = "d-flex align-items-center justify-content-between",
         tags$label("Output type", class = "control-label"),
         tooltip(
           icon("circle-question"),
           "Select the file extension for output merged files."
         )
       ),
       selectInput(
         "merge_filetype",
         label = NULL,
         choices = c(
           "FASTQ",
           "FASTA"
         ),
         selected = "FASTQ"
       ),

       br(),

       div(
         class = "d-flex align-items-center justify-content-between",
         tags$label("Additional flags", class = "control-label")),

       div(
         class = "sidebar-checkbox-row",
         checkboxInput(
           "merge_statistics",
           label = NULL,
           value = TRUE
         ),
         tags$span(
           "Show file parsing errors",
           class = "sidebar-checkbox-label"
         ),
         tooltip(
           icon("circle-question"),
           "Provide error about the file parsing."
         )
       ),

       div(style = "height:4px"),

       div(
         class = "sidebar-checkbox-row",
         checkboxInput(
           "merge_statistics",
           label = NULL,
           value = TRUE
         ),
         tags$span(
           "Optional statistics",
           class = "sidebar-checkbox-label"
         ),
         tooltip(
           icon("circle-question"),
           "Show some optional statistics."
         )
       ),

       div(style = "height:4px"),

       div(
         class = "sidebar-checkbox-row",
         checkboxInput(
           "merge_discard_ambig",
           label = NULL,
           value = FALSE
         ),
         tags$span(
           "Discard N's",
           class = "sidebar-checkbox-label"
         ),
         tooltip(
           icon("circle-question"),
           "Discard reads containing ambigious (N) base calls."
         )
       ),


       div(style = "height:4px"),

       div(
         class = "sidebar-checkbox-row",
         checkboxInput(
           "merge_seq_build_log",
           label = NULL,
           value = FALSE
         ),
         tags$span(
           "Sequence build info",
           class = "sidebar-checkbox-label"
         ),
         tooltip(
           icon("circle-question"),
           "Provide information about the building of a sequence."
         )
       ),

       div(style = "height:4px"),

       div(
         class = "sidebar-checkbox-row",
         checkboxInput(
           "merge_seq_reconstruction_log",
           label = NULL,
           value = FALSE
         ),
         tags$span(
           "Sequence reconstruction info",
           class = "sidebar-checkbox-label"
         ),
         tooltip(
           icon("circle-question"),
           "Show excruciating detail about reconstruction. Not recommended."
         )
       ),



       div(style = "height:4px"),

       div(
         class = "sidebar-checkbox-row",
         checkboxInput(
           "merge_kmer",
           label = NULL,
           value = FALSE
         ),
         tags$span(
           "Show k-mer table",
           class = "sidebar-checkbox-label"
         ),
         tooltip(
           icon("circle-question"),
           "Show information about building the k-mer table. Not recommended."
         )
       ),


       div(style = "height:4px"),

       div(
         class = "sidebar-checkbox-row",
         checkboxInput(
           "merge_kmer",
           label = NULL,
           value = FALSE
         ),
         tags$span(
           "Show every mismatch",
           class = "sidebar-checkbox-label"
         ),
         tooltip(
           icon("circle-question"),
           "Show every mismatch. Not recommended."
         )
       ),

       div(style = "height:4px"),

       div(
         class = "sidebar-checkbox-row",
         checkboxInput(
           "require_illumina_tag",
           label = NULL,
           value = FALSE
         ),
         tags$span(
           "Require Illumina tag",
           class = "sidebar-checkbox-label"
         ),
         tooltip(
           icon("circle-question"),
           "Require the Illumina adapter/tag sequence before trimming. If you're not sure, these were probably removed from your files already."
         )
       ),

       )),

      "shm.app" = div(
        class = "sidebar-controls",

        h5("Mutation Controls"),
        tags$p("Measure Mutations controls go here.")
      ),

      "cut.app" = div(
        class = "sidebar-controls",

        h5("Nuclease Controls"),
        tags$p("Nuclease Analysis controls go here.")
      ),

      "help.app" = div(
        class = "sidebar-controls",

        h5("Help"),
        tags$p("Help and documentation go here.")
      )
    )
  })

  ###########
  # merge reads logic

  observe({
    single_ready <- input$upload_mode == "single_sample" &&
      nrow(merge_rv$queued) > 0 &&
      all(merge_rv$queued$Status == "Ready")

    if (single_ready) {
      shinyjs::enable("run_single_merge")
    } else {
      shinyjs::disable("run_single_merge")
    }
  })

  observe({
    batch_ready <- input$upload_mode == "batch_import" &&
      !is.null(batch_rv$archive_meta) &&
      !is.null(batch_rv$working_manifest_path) &&
      file.exists(batch_rv$working_manifest_path)

    if (batch_ready) {
      shinyjs::enable("run_batch_merge")
    } else {
      shinyjs::disable("run_batch_merge")
    }
  })

    # ---------- session-scoped upload root ----------

    upload_root <- file.path(tempdir(), "tasaR_uploads", session$token)
    dir.create(upload_root, recursive = TRUE, showWarnings = FALSE)

    session$onSessionEnded(function() {
      unlink(upload_root, recursive = TRUE, force = TRUE)
    })

    # ---------- single-sample state ----------
    merge_rv <- reactiveValues(
      queued = data.frame(
        "Queue ID" = character(),
        "Sample Name" = character(),
        "Left (R1) File" = character(),
        "Right (R2) File" = character(),
        "R1 Path" = character(),
        "R2 Path" = character(),
        "Status" = character(),
        check.names = FALSE
      ),
      r1_meta = NULL,
      r2_meta = NULL,
      files_ready = FALSE,
      message = NULL,
      message_type = NULL
    )

    refresh_merge_ready <- function() {
      merge_rv$files_ready <- !is.null(merge_rv$r1_meta) && !is.null(merge_rv$r2_meta)
    }

    observeEvent(input$fastq_r1, {
      .remove_managed_upload(merge_rv$r1_meta)
      merge_rv$r1_meta <- .save_managed_upload(input$fastq_r1, "R1", upload_root)
      refresh_merge_ready()
    }, ignoreInit = TRUE)

    observeEvent(input$fastq_r2, {
      .remove_managed_upload(merge_rv$r2_meta)
      merge_rv$r2_meta <- .save_managed_upload(input$fastq_r2, "R2", upload_root)
      refresh_merge_ready()
    }, ignoreInit = TRUE)

    observeEvent(input$queue_sample, {
      sample_name <- trimws(input$sample_name)

      if (!nzchar(sample_name)) {
        merge_rv$message <- "Enter a sample name."
        merge_rv$message_type <- "error"
        return()
      }

      if (!merge_rv$files_ready) {
        merge_rv$message <- "Please upload both FASTQ files."
        merge_rv$message_type <- "error"
        return()
      }

      other_samples <- if (nrow(merge_rv$queued) > 0) merge_rv$queued[["Sample Name"]] else character(0)
      other_files <- if (nrow(merge_rv$queued) > 0) {
        unlist(merge_rv$queued[, c("Left (R1) File", "Right (R2) File")], use.names = FALSE)
      } else {
        character(0)
      }

      status <- .check_fastq_pair(
        sample_name,
        merge_rv$r1_meta$original_name,
        merge_rv$r2_meta$original_name,
        other_samples = other_samples,
        other_files = other_files
      )

      queue_id <- paste0("q_", nrow(merge_rv$queued) + 1L, "_", format(Sys.time(), "%H%M%OS3"))

      merge_rv$queued <- rbind(
        merge_rv$queued,
        data.frame(
          "Queue ID" = queue_id,
          "Sample Name" = sample_name,
          "Left (R1) File" = merge_rv$r1_meta$original_name,
          "Right (R2) File" = merge_rv$r2_meta$original_name,
          "R1 Path" = merge_rv$r1_meta$managed_path,
          "R2 Path" = merge_rv$r2_meta$managed_path,
          "Status" = status,
          check.names = FALSE
        )
      )

      updateTextInput(session, "sample_name", value = "")
      shinyjs::reset("fastq_r1")
      shinyjs::reset("fastq_r2")

      merge_rv$r1_meta <- NULL
      merge_rv$r2_meta <- NULL
      merge_rv$files_ready <- FALSE

      merge_rv$message <- paste("Queued sample", sample_name, "successfully.")
      merge_rv$message_type <- "success"
    })

    observeEvent(input$clear_single, {
      .remove_managed_upload(merge_rv$r1_meta)
      .remove_managed_upload(merge_rv$r2_meta)

      merge_rv$r1_meta <- NULL
      merge_rv$r2_meta <- NULL
      merge_rv$files_ready <- FALSE
      merge_rv$message <- "Cleared uploaded FASTQ files."
      merge_rv$message_type <- "info"

      shinyjs::reset("fastq_r1")
      shinyjs::reset("fastq_r2")
      updateTextInput(session, "sample_name", value = "")
    }, ignoreInit = TRUE)

    observeEvent(merge_rv$message, {
      req(!is.null(merge_rv$message))

      later::later(function() {
        merge_rv$message <- NULL
        merge_rv$message_type <- NULL
      }, 4)
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    output$merge_message <- renderUI({
      req(!is.null(merge_rv$message), !is.null(merge_rv$message_type))

      div(
        class = paste("queue-alert", paste0("queue-alert-", merge_rv$message_type)),
        merge_rv$message
      )
    })

    output$queued_samples <- renderUI({
      if (nrow(merge_rv$queued) == 0) return(NULL)

      rows <- lapply(seq_len(nrow(merge_rv$queued)), function(i) {
        queue_id <- merge_rv$queued[i, "Queue ID", drop = TRUE]

        tags$tr(
          tags$td(merge_rv$queued[i, "Sample Name", drop = TRUE]),
          tags$td(merge_rv$queued[i, "Left (R1) File", drop = TRUE]),
          tags$td(merge_rv$queued[i, "Right (R2) File", drop = TRUE]),
          tags$td(
            tags$div(
              class = "queue-status-cell",
              tags$div(
                class = "queue-status-text",
                HTML(merge_rv$queued[i, "Status", drop = TRUE])
              ),
              tags$button(
                type = "button",
                class = "queue-remove-btn",
                onclick = sprintf(
                  "Shiny.setInputValue('remove_queue_id', '%s', {priority: 'event'})",
                  queue_id
                ),
                "\u00d7 Remove"
              )
            )
          )
        )
      })

      tags$table(
        class = "table queue-table-wrapper queue-table",
        tags$thead(
          tags$tr(
            tags$th("Sample Name"),
            tags$th("Left (R1) File"),
            tags$th("Right (R2) File"),
            tags$th("Status")
          )
        ),
        tags$tbody(rows)
      )
    })

    observeEvent(input$remove_queue_id, {

      res <- .remove_queued_sample(
        merge_rv$queued,
        input$remove_queue_id
      )

      if (!res$success) {
        return()
      }

      row <- res$removed

      .remove_managed_upload(list(
        managed_path = row[["R1 Path"]]
      ))

      .remove_managed_upload(list(
        managed_path = row[["R2 Path"]]
      ))

      merge_rv$queued <- res$queued

      merge_rv$message <- paste(
        "Removed sample",
        row[["Sample Name"]]
      )

      merge_rv$message_type <- "info"

    }, ignoreInit = TRUE)

    # ---------- batch upload state ----------
    batch_rv <- reactiveValues(
      manifest_meta = NULL,
      archive_meta = NULL,
      working_manifest_path = NULL,
      preview_df = NULL,
      extract_dir = NULL,
      message = NULL,
      message_type = NULL
    )


    ensure_batch_extract_dir <- function() {
      if (is.null(batch_rv$archive_meta)) return(NULL)

      if (!is.null(batch_rv$extract_dir) && dir.exists(batch_rv$extract_dir)) {
        return(batch_rv$extract_dir)
      }

      batch_rv$extract_dir <- file.path(
        upload_root,
        "batch_extract"
      )
      dir.create(batch_rv$extract_dir, recursive = TRUE, showWarnings = FALSE)
      batch_rv$extract_dir
    }

    refresh_batch_preview <- function() {
      if (is.null(batch_rv$manifest_meta) || is.null(batch_rv$archive_meta)) {
        batch_rv$preview_df <- data.frame(
          Note = "Upload the manifest CSV and ZIP archive to preview detected samples.",
          check.names = FALSE
        )
        batch_rv$working_manifest_path <- NULL
        return()
      }

      manifest_df <- .read_batch_manifest(batch_rv$manifest_meta$managed_path)
      manifest_df <- .batch_manifest_with_status(manifest_df)

      batch_rv$preview_df <- manifest_df
      batch_rv$working_manifest_path <- file.path(
        upload_root,
        "batch_working_manifest.csv"
      )

      utils::write.csv(
        batch_rv$preview_df,
        batch_rv$working_manifest_path,
        row.names = FALSE,
        quote = TRUE
      )
    }

    observeEvent(input$batch_manifest, {
      .remove_managed_upload(batch_rv$manifest_meta)
      batch_rv$manifest_meta <- .save_managed_upload(
        input$batch_manifest,
        "batch_manifest",
        upload_root
      )

      if (!is.null(batch_rv$manifest_meta)) {
        batch_rv$message <- "Batch manifest uploaded."
        batch_rv$message_type <- "success"
      }

      refresh_batch_preview()
    }, ignoreInit = TRUE)

    observeEvent(input$batch_archive, {
      .remove_managed_upload(batch_rv$archive_meta)
      .remove_dir_safely(batch_rv$extract_dir)

      batch_rv$archive_meta <- .save_managed_upload(
        input$batch_archive,
        "batch_archive",
        upload_root
      )

      batch_rv$extract_dir <- NULL

      if (!is.null(batch_rv$archive_meta)) {
        batch_rv$message <- "Batch archive uploaded."
        batch_rv$message_type <- "success"
      }

      refresh_batch_preview()
    }, ignoreInit = TRUE)

    observeEvent(input$batch_remove_row, {
      req(!is.null(batch_rv$preview_df))

      res <- .remove_batch_row(batch_rv$preview_df, input$batch_remove_row)
      if (!res$success) return()

      batch_rv$preview_df <- if (nrow(res$df) > 0) {
        .batch_manifest_with_status(res$df)
      } else {
        data.frame(
          Note = "Batch manifest is now empty.",
          check.names = FALSE
        )
      }

      if (!is.null(batch_rv$working_manifest_path) && nrow(batch_rv$preview_df) > 0) {
        utils::write.csv(
          batch_rv$preview_df,
          batch_rv$working_manifest_path,
          row.names = FALSE,
          quote = TRUE
        )
      }

      batch_rv$message <- paste("Removed sample", res$removed[["Sample Name"]], "from batch preview.")
      batch_rv$message_type <- "info"
    }, ignoreInit = TRUE)

    observeEvent(input$clear_batch, {

      .remove_managed_upload(batch_rv$manifest_meta)
      .remove_managed_upload(batch_rv$archive_meta)

      .remove_dir_safely(batch_rv$extract_dir)

      if (!is.null(batch_rv$working_manifest_path) &&
          file.exists(batch_rv$working_manifest_path)) {
        unlink(batch_rv$working_manifest_path, force = TRUE)
      }

      batch_rv$manifest_meta <- NULL
      batch_rv$archive_meta <- NULL
      batch_rv$working_manifest_path <- NULL
      batch_rv$preview_df <- NULL
      batch_rv$extract_dir <- NULL

      batch_rv$message <- "Batch inputs cleared."
      batch_rv$message_type <- "info"

      shinyjs::reset("batch_manifest")
      shinyjs::reset("batch_archive")

    }, ignoreInit = TRUE)

    observeEvent(input$run_batch_merge, {

      if (is.null(batch_rv$archive_meta)) {
        batch_rv$message <- "Upload the FASTQ archive first."
        batch_rv$message_type <- "error"
        return()
      }

      if (is.null(batch_rv$working_manifest_path) ||
          !file.exists(batch_rv$working_manifest_path)) {

        batch_rv$message <- "No valid batch manifest available."
        batch_rv$message_type <- "error"
        return()
      }

      extract_dir <- ensure_batch_extract_dir()

      .remove_dir_safely(extract_dir)

      dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

      utils::unzip(
        batch_rv$archive_meta$managed_path,
        exdir = extract_dir
      )

      manifest_df <- data.table::fread(
        batch_rv$working_manifest_path,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      .validate_batch_manifest(manifest_df)

      extracted_files <- list.files(
        extract_dir,
        recursive = TRUE,
        full.names = TRUE
      )

      extracted_basenames <- basename(extracted_files)

      manifest_df$R1_Path <- NA_character_
      manifest_df$R2_Path <- NA_character_

      for (i in seq_len(nrow(manifest_df))) {

        r1_match <- match(
          manifest_df[i, "Left (R1) File", drop = TRUE],
          extracted_basenames
        )

        r2_match <- match(
          manifest_df[i, "Right (R2) File", drop = TRUE],
          extracted_basenames
        )

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

        missing_samples <- manifest_df[
          missing_r1 | missing_r2,
          "Sample Name",
          drop = TRUE
        ]

        batch_rv$message <- paste(
          "Missing FASTQ files for:",
          paste(missing_samples, collapse = ", ")
        )

        batch_rv$message_type <- "error"

        return()
      }

      batch_rv$message <- paste(
        "Batch archive extracted and",
        nrow(manifest_df),
        "samples validated successfully."
      )

      batch_rv$message_type <- "success"

      # -------------------------------------------------
      # YOUR MERGE PIPELINE STARTS HERE
      # -------------------------------------------------
      #
      # Use:
      #
      # manifest_df$R1_Path
      # manifest_df$R2_Path
      #
      # Each row now has fully resolved FASTQ paths.
      #
      # Example:
      #
      # for (i in seq_len(nrow(manifest_df))) {
      #
      #   sample_name <- manifest_df[i, "Sample Name"]
      #   r1 <- manifest_df[i, "R1_Path"]
      #   r2 <- manifest_df[i, "R2_Path"]
      #
      #   ...
      # }
      #
      # -------------------------------------------------

    }, ignoreInit = TRUE)

    output$batch_queued_samples <- renderUI({
      if (is.null(batch_rv$preview_df)) return(NULL)

      if ("Note" %in% names(batch_rv$preview_df)) {
        return(
          div(
            class = "queue-alert queue-alert-info",
            batch_rv$preview_df[["Note"]][1]
          )
        )
      }

      if (nrow(batch_rv$preview_df) == 0) return(NULL)

      rows <- lapply(seq_len(nrow(batch_rv$preview_df)), function(i) {
        tags$tr(
          tags$td(batch_rv$preview_df[i, "Sample Name", drop = TRUE]),
          tags$td(batch_rv$preview_df[i, "Left (R1) File", drop = TRUE]),
          tags$td(batch_rv$preview_df[i, "Right (R2) File", drop = TRUE]),
          tags$td(
            tags$div(
              class = "queue-status-cell",
              tags$div(
                class = "queue-status-text",
                HTML(batch_rv$preview_df[i, "Status", drop = TRUE])
              ),
              tags$button(
                type = "button",
                class = "queue-remove-btn",
                onclick = sprintf(
                  "Shiny.setInputValue('batch_remove_row', '%d', {priority: 'event'})",
                  i
                ),
                "\u00d7 Remove"
              )
            )
          )
        )
      })

      tags$table(
        class = "table queue-table",
        tags$thead(
          tags$tr(
            tags$th("Sample Name"),
            tags$th("Left (R1) File"),
            tags$th("Right (R2) File"),
            tags$th("Status")
          )
        ),
        tags$tbody(rows)
      )
    })

    output$batch_message <- renderUI({
      req(!is.null(batch_rv$message), !is.null(batch_rv$message_type))

      div(
        class = paste("queue-alert", paste0("queue-alert-", batch_rv$message_type)),
        batch_rv$message
      )
    })

  ##########
  output$merge_settings <- renderText({
      paste(
        "Output prefix:", input$merge_prefix,
        "\nEngine:", input$merge_engine,
        "\nMin overlap:", input$min_overlap,
        "\nMin length:", input$min_length,
        "\nMax length:", input$max_length,
        "\nOverlap threshold:", input$merge_similarity
      )
    })

  observeEvent(input$min_overlap_num, {
    val <- suppressWarnings(as.numeric(input$min_overlap_num))
    if (is.na(val)) return()
    val <- max(5, min(100, val))
    if (!identical(val, input$min_overlap)) {
      updateSliderInput(session, "min_overlap", value = val)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$min_overlap, {
    val <- suppressWarnings(as.numeric(input$min_overlap))
    if (is.na(val)) return()
    if (!identical(val, input$min_overlap_num)) {
      updateNumericInput(session, "min_overlap_num", value = val)
    }
  }, ignoreInit = TRUE)


  observeEvent(input$min_length_num, {
    val <- suppressWarnings(as.numeric(input$min_length_num))
    if (is.na(val)) return()
    val <- max(50, min(500, val))
    if (!identical(val, input$min_length)) {
      updateSliderInput(session, "min_length", value = val)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$min_length, {
    val <- suppressWarnings(as.numeric(input$min_length))
    if (is.na(val)) return()
    if (!identical(val, input$min_length_num)) {
      updateNumericInput(session, "min_length_num", value = val)
    }
  }, ignoreInit = TRUE)


  observeEvent(input$max_length_num, {
    val <- suppressWarnings(as.numeric(input$max_length_num))
    if (is.na(val)) return()
    val <- max(50, min(500, val))
    if (!identical(val, input$max_length)) {
      updateSliderInput(session, "max_length", value = val)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$max_length, {
    val <- suppressWarnings(as.numeric(input$max_length))
    if (is.na(val)) return()
    if (!identical(val, input$max_length_num)) {
      updateNumericInput(session, "max_length_num", value = val)
    }
  }, ignoreInit = TRUE)



  observeEvent(input$merge_similarity_num, {
    val <- suppressWarnings(as.numeric(input$merge_similarity_num))
    if (is.na(val)) return()
    val <- max(0, min(1, val))
    if (!identical(val, input$merge_similarity)) {
      updateSliderInput(session, "merge_similarity", value = val)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$merge_similarity, {
    val <- suppressWarnings(as.numeric(input$merge_similarity))
    if (is.na(val)) return()
    if (!identical(val, input$merge_similarity_num)) {
      updateNumericInput(session, "merge_similarity_num", value = val)
    }
  }, ignoreInit = TRUE)

  ###################################################

  observeEvent(input$run_single_merge, {

    merge_inputs <- .prepare_merge_inputs(
      mode = "single",
      queued_df = merge_rv$queued
    )

    ### this function not written yet
    .run_merge_pipeline(
      merge_inputs,
      settings = merge_settings()
    )
  })

  observeEvent(input$run_batch_merge, {

    batch_inputs <- prepare_merge_inputs(
      mode = "batch",
      manifest_path = batch_rv$working_manifest_path,
      archive_path = batch_rv$archive_meta$managed_path,
      extract_dir = ensure_batch_extract_dir()
    )

    ### this function not written yet
    .run_merge_pipeline(
      batch_inputs,
      settings = merge_settings()
    )
  })





}


shinyApp(ui = ui, server = server)
