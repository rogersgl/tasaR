setup_nuclease_server <- function(input, output, session) {
  upload_root <- file.path(tempdir(), "tasaR_nuclease_uploads", session$token)
  dir.create(upload_root, recursive = TRUE, showWarnings = FALSE)

  session$onSessionEnded(function() {
    unlink(upload_root, recursive = TRUE, force = TRUE)
  })

  nuclease_rv <- reactiveValues(
    queued = .nuclease_empty_queue("Nuclease Queue ID"),
    sample_order = character(),
    control_meta = NULL,
    experimental_meta = NULL,
    message = NULL,
    message_type = NULL,
    active_settings_id = NULL
  )

  batch_rv <- reactiveValues(
    manifest_meta = NULL,
    archive_meta = NULL,
    extract_dir = NULL,
    message = NULL,
    message_type = NULL
  )

  setting_map <- c(
    "Forward Primer" = "forward_primer",
    "Reverse Primer" = "reverse_primer",
    "Forward Extension Type" = "forward_extension_type",
    "Forward Extension Sequence" = "forward_extension_sequence",
    "Reverse Extension Type" = "reverse_extension_type",
    "Reverse Extension Sequence" = "reverse_extension_sequence"
  )

  collect_settings <- function(prefix) {
    defaults <- .nuclease_default_settings()
    out <- defaults
    for (nm in names(setting_map)) {
      value <- input[[paste0(prefix, "_", setting_map[[nm]])]]
      if (is.null(value)) value <- defaults[[nm]]
      out[[nm]] <- value
    }
    out
  }

  collect_sample_settings <- function(prefix) {
    c(
      list(
        "gRNA Sequence" = input[[paste0(prefix, "_grna_sequence")]],
        "Cut Position" = input[[paste0(prefix, "_cut_position")]]
      ),
      collect_settings(prefix)
    )
  }

  numeric_value <- function(value) {
    value <- suppressWarnings(as.numeric(value))
    if (is.na(value)) NA else value
  }

  update_settings_inputs <- function(prefix, settings) {
    settings <- modifyList(.nuclease_default_settings(), as.list(settings))
    updateTextInput(session, paste0(prefix, "_forward_primer"), value = settings[["Forward Primer"]])
    updateTextInput(session, paste0(prefix, "_reverse_primer"), value = settings[["Reverse Primer"]])
    updateSelectInput(session, paste0(prefix, "_forward_extension_type"), selected = settings[["Forward Extension Type"]])
    updateTextInput(session, paste0(prefix, "_forward_extension_sequence"), value = settings[["Forward Extension Sequence"]])
    updateSelectInput(session, paste0(prefix, "_reverse_extension_type"), selected = settings[["Reverse Extension Type"]])
    updateTextInput(session, paste0(prefix, "_reverse_extension_sequence"), value = settings[["Reverse Extension Sequence"]])
  }

  update_sample_settings_inputs <- function(prefix, settings) {
    updateTextInput(session, paste0(prefix, "_grna_sequence"), value = settings[["gRNA Sequence"]])
    updateNumericInput(session, paste0(prefix, "_cut_position"), value = numeric_value(settings[["Cut Position"]]))
    update_settings_inputs(prefix, settings)
  }

  reset_main_settings <- function() {
    updateTextInput(session, "nuclease_grna_sequence", value = "")
    updateNumericInput(session, "nuclease_cut_position", value = NA)
    update_settings_inputs("nuclease", .nuclease_default_settings())
  }

  refresh_nuclease_status <- function() {
    nuclease_rv$queued <- .nuclease_queue_with_status(nuclease_rv$queued)
  }

  reconcile_nuclease_order <- function() {
    nuclease_rv$sample_order <- .reconcile_sample_order(
      nuclease_rv$queued[["Nuclease Queue ID"]],
      nuclease_rv$sample_order
    )
  }

  ensure_batch_extract_dir <- function() {
    if (is.null(batch_rv$archive_meta)) return(NULL)

    batch_rv$extract_dir <- file.path(
      upload_root,
      paste0("nuclease_batch_extract_", format(Sys.time(), "%Y%m%d%H%M%OS3"))
    )
    dir.create(batch_rv$extract_dir, recursive = TRUE, showWarnings = FALSE)
    batch_rv$extract_dir
  }

  read_nuclease_manifest <- function(path) {
    df <- data.table::fread(path, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("Sample Name", "Control File", "Experimental File", "gRNA Sequence", "Cut Position")
    missing <- setdiff(required, names(df))
    if (length(missing) > 0) {
      stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
    }
    if (nrow(df) == 0) stop("Nuclease batch manifest is empty.")

    defaults <- .nuclease_default_settings()
    for (nm in names(defaults)) {
      if (!nm %in% names(df)) df[[nm]] <- defaults[[nm]]
    }

    df
  }

  append_nuclease_batch <- function() {
    if (is.null(batch_rv$manifest_meta) || is.null(batch_rv$archive_meta)) {
      return()
    }

    tryCatch({
      extract_dir <- ensure_batch_extract_dir()
      dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

      utils::unzip(batch_rv$archive_meta$managed_path, exdir = extract_dir)

      manifest_df <- read_nuclease_manifest(batch_rv$manifest_meta$managed_path)
      extracted_files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
      extracted_basenames <- basename(extracted_files)
      control_idx <- match(manifest_df[["Control File"]], extracted_basenames)
      experimental_idx <- match(manifest_df[["Experimental File"]], extracted_basenames)
      missing <- is.na(control_idx) | is.na(experimental_idx)

      if (any(missing)) {
        stop(paste(
          "Missing control or experimental FASTQ files for:",
          paste(manifest_df[missing, "Sample Name", drop = TRUE], collapse = ", ")
        ))
      }

      batch_id <- format(Sys.time(), "%Y%m%d%H%M%OS3")
      rows <- data.frame(
        "Nuclease Queue ID" = paste0("nuc_batch_", batch_id, "_", seq_len(nrow(manifest_df))),
        "Sample Name" = manifest_df[["Sample Name"]],
        "Control File" = manifest_df[["Control File"]],
        "Control Path" = extracted_files[control_idx],
        "Experimental File" = manifest_df[["Experimental File"]],
        "Experimental Path" = extracted_files[experimental_idx],
        "gRNA Sequence" = manifest_df[["gRNA Sequence"]],
        "Cut Position" = manifest_df[["Cut Position"]],
        "Source" = "batch",
        check.names = FALSE
      )

      for (nm in .nuclease_setting_cols()) {
        rows[[nm]] <- manifest_df[[nm]]
      }
      rows$Status <- ""
      rows$Progress <- ""

      nuclease_rv$queued <- rbind(nuclease_rv$queued, rows)
      refresh_nuclease_status()
      reconcile_nuclease_order()

      batch_rv$message <- paste("Queued", nrow(rows), "batch samples successfully.")
      batch_rv$message_type <- "success"
    }, error = function(e) {
      batch_rv$message <- conditionMessage(e)
      batch_rv$message_type <- "error"
    })
  }

  observeEvent(input$nuclease_control_fastq, {
    .remove_managed_upload(nuclease_rv$control_meta)
    nuclease_rv$control_meta <- .save_managed_upload(input$nuclease_control_fastq, "nuclease_control_fastq", upload_root)
  }, ignoreInit = TRUE)

  observeEvent(input$nuclease_experimental_fastq, {
    .remove_managed_upload(nuclease_rv$experimental_meta)
    nuclease_rv$experimental_meta <- .save_managed_upload(input$nuclease_experimental_fastq, "nuclease_experimental_fastq", upload_root)
  }, ignoreInit = TRUE)

  observeEvent(input$nuclease_queue_sample, {
    sample_name <- trimws(input$nuclease_sample_name)
    if (is.null(nuclease_rv$control_meta) || is.null(nuclease_rv$experimental_meta)) {
      nuclease_rv$message <- "Please upload both FASTQ files."
      nuclease_rv$message_type <- "error"
      return()
    }

    settings <- collect_sample_settings("nuclease")
    queue_id <- paste0("nuc_single_", nrow(nuclease_rv$queued) + 1L, "_", format(Sys.time(), "%H%M%OS3"))
    row <- data.frame(
      "Nuclease Queue ID" = queue_id,
      "Sample Name" = sample_name,
      "Control File" = nuclease_rv$control_meta$original_name,
      "Control Path" = nuclease_rv$control_meta$managed_path,
      "Experimental File" = nuclease_rv$experimental_meta$original_name,
      "Experimental Path" = nuclease_rv$experimental_meta$managed_path,
      "Source" = "single",
      check.names = FALSE
    )
    for (nm in names(settings)) row[[nm]] <- settings[[nm]]
    row$Status <- ""
    row$Progress <- ""

    nuclease_rv$queued <- rbind(nuclease_rv$queued, row)
    refresh_nuclease_status()
    reconcile_nuclease_order()

    updateTextInput(session, "nuclease_sample_name", value = "")
    shinyjs::reset("nuclease_control_fastq")
    shinyjs::reset("nuclease_experimental_fastq")
    nuclease_rv$control_meta <- NULL
    nuclease_rv$experimental_meta <- NULL
    nuclease_rv$message <- paste("Queued sample", ifelse(nzchar(sample_name), sample_name, "(unnamed)"), "successfully.")
    nuclease_rv$message_type <- "success"
  })

  observeEvent(input$nuclease_clear_settings, {
    reset_main_settings()
    nuclease_rv$message <- "Nuclease settings reset."
    nuclease_rv$message_type <- "info"
  }, ignoreInit = TRUE)

  observeEvent(input$nuclease_batch_manifest, {
    .remove_managed_upload(batch_rv$manifest_meta)
    batch_rv$manifest_meta <- .save_managed_upload(input$nuclease_batch_manifest, "nuclease_manifest", upload_root)
    batch_rv$message <- "Batch manifest uploaded."
    batch_rv$message_type <- "success"
  }, ignoreInit = TRUE)

  observeEvent(input$nuclease_batch_archive, {
    .remove_managed_upload(batch_rv$archive_meta)
    batch_rv$extract_dir <- NULL
    batch_rv$archive_meta <- .save_managed_upload(input$nuclease_batch_archive, "nuclease_archive", upload_root)
    batch_rv$message <- "Batch archive uploaded."
    batch_rv$message_type <- "success"
  }, ignoreInit = TRUE)

  observeEvent(input$nuclease_queue_batch_samples, {
    if (is.null(batch_rv$manifest_meta) || is.null(batch_rv$archive_meta)) {
      batch_rv$message <- "Upload the batch manifest and archive before adding samples."
      batch_rv$message_type <- "error"
      return()
    }

    append_nuclease_batch()
  }, ignoreInit = TRUE)

  observeEvent(input$nuclease_clear_batch, {
    .remove_managed_upload(batch_rv$manifest_meta)
    .remove_managed_upload(batch_rv$archive_meta)

    batch_rv$manifest_meta <- NULL
    batch_rv$archive_meta <- NULL
    batch_rv$extract_dir <- NULL

    shinyjs::reset("nuclease_batch_manifest")
    shinyjs::reset("nuclease_batch_archive")

    batch_rv$message <- "Batch inputs cleared."
    batch_rv$message_type <- "info"
  }, ignoreInit = TRUE)

  observeEvent(input$nuclease_remove_row, {
    res <- .remove_batch_row_by_id(nuclease_rv$queued, input$nuclease_remove_row, id_col = "Nuclease Queue ID")
    if (!res$success) return()

    .remove_managed_upload(list(managed_path = res$removed[["Control Path"]]))
    .remove_managed_upload(list(managed_path = res$removed[["Experimental Path"]]))
    nuclease_rv$queued <- res$df
    refresh_nuclease_status()
    reconcile_nuclease_order()

    nuclease_rv$message <- paste("Removed sample", res$removed[["Sample Name"]])
    nuclease_rv$message_type <- "info"
  }, ignoreInit = TRUE)

  output$nuclease_message <- renderUI({
    .queue_upload_message(nuclease_rv$message, nuclease_rv$message_type)
  })

  output$nuclease_batch_message <- renderUI({
    .queue_upload_message(batch_rv$message, batch_rv$message_type)
  })

  observeEvent(nuclease_rv$message, {
    req(!is.null(nuclease_rv$message))
    .clear_queue_message_later(nuclease_rv)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  observeEvent(batch_rv$message, {
    req(!is.null(batch_rv$message))
    .clear_queue_message_later(batch_rv)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  .setup_reorderable_sample_table(
    input = input,
    output = output,
    output_id = "nuclease_queued_samples",
    data_fn = function() nuclease_rv$queued,
    id_col = "Nuclease Queue ID",
    get_order = function() nuclease_rv$sample_order,
    set_order = function(order) {
      nuclease_rv$sample_order <- order
      refresh_nuclease_status()
    },
    remove_input_id = "nuclease_remove_row",
    move_up_input_id = "nuclease_queued_samples_move_up",
    move_down_input_id = "nuclease_queued_samples_move_down",
    table_builder = function(...) .make_nuclease_table_data(..., settings_input_id = "nuclease_open_settings")
  )

  open_settings_modal <- function(row_id) {
    idx <- match(row_id, nuclease_rv$queued[["Nuclease Queue ID"]])
    if (is.na(idx)) return()

    row <- nuclease_rv$queued[idx, , drop = FALSE]
    nuclease_rv$active_settings_id <- row_id
    showModal(modalDialog(
      title = paste("Settings:", row[["Sample Name"]]),
      size = "l",
      easyClose = FALSE,
      nuclease_sample_settings_ui("nuclease_modal"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("nuclease_modal_save", "Save", class = "primary-action-btn")
      )
    ))
    later::later(function() update_sample_settings_inputs("nuclease_modal", row), 0.1)
  }

  observeEvent(input$nuclease_open_settings, {
    open_settings_modal(input$nuclease_open_settings)
  }, ignoreInit = TRUE)

  observeEvent(input$nuclease_modal_save, {
    req(!is.null(nuclease_rv$active_settings_id))
    idx <- match(nuclease_rv$active_settings_id, nuclease_rv$queued[["Nuclease Queue ID"]])
    if (!is.na(idx)) {
      settings <- collect_sample_settings("nuclease_modal")
      for (nm in names(settings)) nuclease_rv$queued[idx, nm] <- settings[[nm]]
      refresh_nuclease_status()
    }

    removeModal()
  }, ignoreInit = TRUE)

  output$nuclease_run_summary <- renderText({
    df <- nuclease_rv$queued
    if (is.null(df) || nrow(df) == 0) {
      queued <- ready <- issues <- 0
    } else {
      queued <- nrow(df)
      ready <- sum(df$Status == "Ready", na.rm = TRUE)
      issues <- queued - ready
    }

    paste(
      "Queued:", queued,
      "\nReady:", ready,
      "\nIssues:", issues
    )
  })

  observe({
    ready <- nrow(nuclease_rv$queued) > 0 && all(nuclease_rv$queued$Status == "Ready")
    if (ready) {
      shinyjs::enable("run_single_nuclease")
      shinyjs::enable("run_batch_nuclease")
    } else {
      shinyjs::disable("run_single_nuclease")
      shinyjs::disable("run_batch_nuclease")
    }

    batch_ready <- !is.null(batch_rv$manifest_meta) && !is.null(batch_rv$archive_meta)
    if (batch_ready) {
      shinyjs::enable("nuclease_queue_batch_samples")
    } else {
      shinyjs::disable("nuclease_queue_batch_samples")
    }
  })

  run_unified_nuclease <- function(target = c("single", "batch")) {
    target <- match.arg(target)
    ordered <- .order_sample_df(
      nuclease_rv$queued,
      nuclease_rv$sample_order,
      "Nuclease Queue ID"
    )

    if (target == "single") {
      nuclease_rv$message <- paste("Nuclease analysis pipeline is not implemented yet for", nrow(ordered), "queued samples.")
      nuclease_rv$message_type <- "info"
    } else {
      batch_rv$message <- paste("Nuclease analysis pipeline is not implemented yet for", nrow(ordered), "queued samples.")
      batch_rv$message_type <- "info"
    }
  }

  observeEvent(input$run_single_nuclease, {
    run_unified_nuclease("single")
  }, ignoreInit = TRUE)

  observeEvent(input$run_batch_nuclease, {
    run_unified_nuclease("batch")
  }, ignoreInit = TRUE)
}
