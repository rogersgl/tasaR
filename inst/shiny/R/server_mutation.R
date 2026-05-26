setup_mutation_server <- function(input, output, session) {
  upload_root <- file.path(tempdir(), "tasaR_mutation_uploads", session$token)
  dir.create(upload_root, recursive = TRUE, showWarnings = FALSE)

  session$onSessionEnded(function() {
    unlink(upload_root, recursive = TRUE, force = TRUE)
  })

  mutation_rv <- reactiveValues(
    queued = .mutation_empty_queue("Mutation Queue ID"),
    sample_order = character(),
    merged_meta = NULL,
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
    "Antibody" = "antibody",
    "Reference Sequence" = "reference_sequence",
    "Forward Primer" = "forward_primer",
    "Reverse Primer" = "reverse_primer",
    "Amplicon Length" = "amplicon_length",
    "Insert Start Coordinate" = "insert_start_coordinate",
    "Insert End Coordinate" = "insert_end_coordinate",
    "Forward Extension Type" = "forward_extension_type",
    "Forward Extension Sequence" = "forward_extension_sequence",
    "Reverse Extension Type" = "reverse_extension_type",
    "Reverse Extension Sequence" = "reverse_extension_sequence",
    "FR1 Start" = "fr1_start",
    "CDR1 Start" = "cdr1_start",
    "FR2 Start" = "fr2_start",
    "CDR2 Start" = "cdr2_start",
    "FR3 Start" = "fr3_start",
    "CDR4 Start" = "cdr4_start",
    "FR4 Start" = "fr4_start",
    "FR4 End" = "fr4_end"
  )

  collect_settings <- function(prefix) {
    defaults <- .mutation_default_settings()
    out <- defaults
    for (nm in names(setting_map)) {
      value <- input[[paste0(prefix, "_", setting_map[[nm]])]]
      if (is.null(value)) value <- defaults[[nm]]
      out[[nm]] <- value
    }
    out[["Antibody"]] <- isTRUE(out[["Antibody"]])
    out
  }

  update_settings_inputs <- function(prefix, settings) {
    settings <- modifyList(.mutation_default_settings(), as.list(settings))
    numeric_value <- function(value) {
      value <- suppressWarnings(as.numeric(value))
      if (is.na(value)) NA else value
    }

    updateCheckboxInput(session, paste0(prefix, "_antibody"), value = isTRUE(as.logical(settings[["Antibody"]])))
    updateTextAreaInput(session, paste0(prefix, "_reference_sequence"), value = settings[["Reference Sequence"]])
    updateTextInput(session, paste0(prefix, "_forward_primer"), value = settings[["Forward Primer"]])
    updateTextInput(session, paste0(prefix, "_reverse_primer"), value = settings[["Reverse Primer"]])
    updateNumericInput(session, paste0(prefix, "_amplicon_length"), value = numeric_value(settings[["Amplicon Length"]]))
    updateNumericInput(session, paste0(prefix, "_insert_start_coordinate"), value = numeric_value(settings[["Insert Start Coordinate"]]))
    updateNumericInput(session, paste0(prefix, "_insert_end_coordinate"), value = numeric_value(settings[["Insert End Coordinate"]]))
    updateSelectInput(session, paste0(prefix, "_forward_extension_type"), selected = settings[["Forward Extension Type"]])
    updateTextInput(session, paste0(prefix, "_forward_extension_sequence"), value = settings[["Forward Extension Sequence"]])
    updateSelectInput(session, paste0(prefix, "_reverse_extension_type"), selected = settings[["Reverse Extension Type"]])
    updateTextInput(session, paste0(prefix, "_reverse_extension_sequence"), value = settings[["Reverse Extension Sequence"]])
    updateNumericInput(session, paste0(prefix, "_fr1_start"), value = numeric_value(settings[["FR1 Start"]]))
    updateNumericInput(session, paste0(prefix, "_cdr1_start"), value = numeric_value(settings[["CDR1 Start"]]))
    updateNumericInput(session, paste0(prefix, "_fr2_start"), value = numeric_value(settings[["FR2 Start"]]))
    updateNumericInput(session, paste0(prefix, "_cdr2_start"), value = numeric_value(settings[["CDR2 Start"]]))
    updateNumericInput(session, paste0(prefix, "_fr3_start"), value = numeric_value(settings[["FR3 Start"]]))
    updateNumericInput(session, paste0(prefix, "_cdr4_start"), value = numeric_value(settings[["CDR4 Start"]]))
    updateNumericInput(session, paste0(prefix, "_fr4_start"), value = numeric_value(settings[["FR4 Start"]]))
    updateNumericInput(session, paste0(prefix, "_fr4_end"), value = numeric_value(settings[["FR4 End"]]))
  }

  reset_main_settings <- function() {
    update_settings_inputs("mutation", .mutation_default_settings())
  }

  refresh_mutation_status <- function() {
    mutation_rv$queued <- .mutation_queue_with_status(mutation_rv$queued)
  }

  reconcile_mutation_order <- function() {
    mutation_rv$sample_order <- .reconcile_sample_order(
      mutation_rv$queued[["Mutation Queue ID"]],
      mutation_rv$sample_order
    )
  }

  ensure_batch_extract_dir <- function() {
    if (is.null(batch_rv$archive_meta)) return(NULL)

    batch_rv$extract_dir <- file.path(
      upload_root,
      paste0("mutation_batch_extract_", format(Sys.time(), "%Y%m%d%H%M%OS3"))
    )
    dir.create(batch_rv$extract_dir, recursive = TRUE, showWarnings = FALSE)
    batch_rv$extract_dir
  }

  read_mutation_manifest <- function(path) {
    df <- data.table::fread(path, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("Sample Name", "Merged File")
    missing <- setdiff(required, names(df))
    if (length(missing) > 0) {
      stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
    }
    if (nrow(df) == 0) stop("Mutation batch manifest is empty.")

    for (nm in .mutation_setting_cols()) {
      if (!nm %in% names(df)) df[[nm]] <- .mutation_default_settings()[[nm]]
    }

    df
  }

  append_mutation_batch <- function() {
    if (is.null(batch_rv$manifest_meta) || is.null(batch_rv$archive_meta)) {
      return()
    }

    tryCatch({
      extract_dir <- ensure_batch_extract_dir()
      dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

      utils::unzip(batch_rv$archive_meta$managed_path, exdir = extract_dir)

      manifest_df <- read_mutation_manifest(batch_rv$manifest_meta$managed_path)
      extracted_files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
      extracted_basenames <- basename(extracted_files)
      merged_idx <- match(manifest_df[["Merged File"]], extracted_basenames)
      missing <- is.na(merged_idx)

      if (any(missing)) {
        stop(paste(
          "Missing merged FASTQ files for:",
          paste(manifest_df[missing, "Sample Name", drop = TRUE], collapse = ", ")
        ))
      }

      batch_id <- format(Sys.time(), "%Y%m%d%H%M%OS3")
      rows <- data.frame(
        "Mutation Queue ID" = paste0("mut_batch_", batch_id, "_", seq_len(nrow(manifest_df))),
        "Sample Name" = manifest_df[["Sample Name"]],
        "Merged File" = manifest_df[["Merged File"]],
        "Merged Path" = extracted_files[merged_idx],
        "Source" = "batch",
        check.names = FALSE
      )

      for (nm in .mutation_setting_cols()) {
        rows[[nm]] <- manifest_df[[nm]]
      }
      rows$Status <- ""
      rows$Progress <- ""

      mutation_rv$queued <- rbind(mutation_rv$queued, rows)
      refresh_mutation_status()
      reconcile_mutation_order()

      batch_rv$message <- paste("Queued", nrow(rows), "batch samples successfully.")
      batch_rv$message_type <- "success"
    }, error = function(e) {
      batch_rv$message <- conditionMessage(e)
      batch_rv$message_type <- "error"
    })
  }

  observeEvent(input$mutation_merged_fastq, {
    .remove_managed_upload(mutation_rv$merged_meta)
    mutation_rv$merged_meta <- .save_managed_upload(input$mutation_merged_fastq, "merged_fastq", upload_root)
  }, ignoreInit = TRUE)

  observeEvent(input$mutation_queue_sample, {
    sample_name <- trimws(input$mutation_sample_name)
    if (is.null(mutation_rv$merged_meta)) {
      mutation_rv$message <- "Please upload a merged FASTQ file."
      mutation_rv$message_type <- "error"
      return()
    }

    settings <- collect_settings("mutation")
    queue_id <- paste0("mut_single_", nrow(mutation_rv$queued) + 1L, "_", format(Sys.time(), "%H%M%OS3"))
    row <- data.frame(
      "Mutation Queue ID" = queue_id,
      "Sample Name" = sample_name,
      "Merged File" = mutation_rv$merged_meta$original_name,
      "Merged Path" = mutation_rv$merged_meta$managed_path,
      "Source" = "single",
      check.names = FALSE
    )
    for (nm in names(settings)) row[[nm]] <- settings[[nm]]
    row$Status <- ""
    row$Progress <- ""

    mutation_rv$queued <- rbind(mutation_rv$queued, row)
    refresh_mutation_status()
    reconcile_mutation_order()

    updateTextInput(session, "mutation_sample_name", value = "")
    shinyjs::reset("mutation_merged_fastq")
    mutation_rv$merged_meta <- NULL
    mutation_rv$message <- paste("Queued sample", ifelse(nzchar(sample_name), sample_name, "(unnamed)"), "successfully.")
    mutation_rv$message_type <- "success"
  })

  observeEvent(input$mutation_clear_settings, {
    reset_main_settings()
    mutation_rv$message <- "Mutation settings reset."
    mutation_rv$message_type <- "info"
  }, ignoreInit = TRUE)

  observeEvent(input$mutation_batch_manifest, {
    .remove_managed_upload(batch_rv$manifest_meta)
    batch_rv$manifest_meta <- .save_managed_upload(input$mutation_batch_manifest, "mutation_manifest", upload_root)
    batch_rv$message <- "Batch manifest uploaded."
    batch_rv$message_type <- "success"
  }, ignoreInit = TRUE)

  observeEvent(input$mutation_batch_archive, {
    .remove_managed_upload(batch_rv$archive_meta)
    batch_rv$extract_dir <- NULL
    batch_rv$archive_meta <- .save_managed_upload(input$mutation_batch_archive, "mutation_archive", upload_root)
    batch_rv$message <- "Batch archive uploaded."
    batch_rv$message_type <- "success"
  }, ignoreInit = TRUE)

  observeEvent(input$mutation_queue_batch_samples, {
    if (is.null(batch_rv$manifest_meta) || is.null(batch_rv$archive_meta)) {
      batch_rv$message <- "Upload the batch manifest and archive before adding samples."
      batch_rv$message_type <- "error"
      return()
    }

    append_mutation_batch()
  }, ignoreInit = TRUE)

  observeEvent(input$mutation_clear_batch, {
    .remove_managed_upload(batch_rv$manifest_meta)
    .remove_managed_upload(batch_rv$archive_meta)

    batch_rv$manifest_meta <- NULL
    batch_rv$archive_meta <- NULL
    batch_rv$extract_dir <- NULL

    shinyjs::reset("mutation_batch_manifest")
    shinyjs::reset("mutation_batch_archive")

    batch_rv$message <- "Batch inputs cleared."
    batch_rv$message_type <- "info"
  }, ignoreInit = TRUE)

  observeEvent(input$mutation_remove_row, {
    res <- .remove_batch_row_by_id(mutation_rv$queued, input$mutation_remove_row, id_col = "Mutation Queue ID")
    if (!res$success) return()

    .remove_managed_upload(list(managed_path = res$removed[["Merged Path"]]))
    mutation_rv$queued <- res$df
    refresh_mutation_status()
    reconcile_mutation_order()

    mutation_rv$message <- paste("Removed sample", res$removed[["Sample Name"]])
    mutation_rv$message_type <- "info"
  }, ignoreInit = TRUE)

  output$mutation_message <- renderUI({
    .queue_upload_message(mutation_rv$message, mutation_rv$message_type)
  })

  output$mutation_batch_message <- renderUI({
    .queue_upload_message(batch_rv$message, batch_rv$message_type)
  })

  observeEvent(mutation_rv$message, {
    req(!is.null(mutation_rv$message))
    .clear_queue_message_later(mutation_rv)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  observeEvent(batch_rv$message, {
    req(!is.null(batch_rv$message))
    .clear_queue_message_later(batch_rv)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  .setup_reorderable_sample_table(
    input = input,
    output = output,
    output_id = "mutation_queued_samples",
    data_fn = function() mutation_rv$queued,
    id_col = "Mutation Queue ID",
    get_order = function() mutation_rv$sample_order,
    set_order = function(order) {
      mutation_rv$sample_order <- order
      refresh_mutation_status()
    },
    remove_input_id = "mutation_remove_row",
    move_up_input_id = "mutation_queued_samples_move_up",
    move_down_input_id = "mutation_queued_samples_move_down",
    table_builder = function(...) .make_mutation_table_data(..., settings_input_id = "mutation_open_settings")
  )

  open_settings_modal <- function(row_id) {
    idx <- match(row_id, mutation_rv$queued[["Mutation Queue ID"]])
    if (is.na(idx)) return()

    row <- mutation_rv$queued[idx, , drop = FALSE]
    mutation_rv$active_settings_id <- row_id
    showModal(modalDialog(
      title = paste("Settings:", row[["Sample Name"]]),
      size = "l",
      easyClose = FALSE,
      mutation_settings_ui("mutation_modal"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("mutation_modal_save", "Save", class = "primary-action-btn")
      )
    ))
    later::later(function() update_settings_inputs("mutation_modal", row), 0.1)
  }

  observeEvent(input$mutation_open_settings, {
    open_settings_modal(input$mutation_open_settings)
  }, ignoreInit = TRUE)

  observeEvent(input$mutation_modal_save, {
    req(!is.null(mutation_rv$active_settings_id))
    idx <- match(mutation_rv$active_settings_id, mutation_rv$queued[["Mutation Queue ID"]])
    if (!is.na(idx)) {
      settings <- collect_settings("mutation_modal")
      for (nm in names(settings)) mutation_rv$queued[idx, nm] <- settings[[nm]]
      refresh_mutation_status()
    }

    removeModal()
  }, ignoreInit = TRUE)

  selected_row <- reactive({
    df <- mutation_rv$queued
    if (is.null(df) || nrow(df) == 0) return(NULL)

    df <- .order_sample_df(df, mutation_rv$sample_order, "Mutation Queue ID")
    selected <- input$mutation_queued_samples_rows_selected
    if (length(selected) == 1 && selected >= 1 && selected <= nrow(df)) {
      return(df[selected, , drop = FALSE])
    }
    if (!is.null(mutation_rv$active_settings_id)) {
      idx <- match(mutation_rv$active_settings_id, df[["Mutation Queue ID"]])
      if (!is.na(idx)) return(df[idx, , drop = FALSE])
    }
    df[1, , drop = FALSE]
  })

  output$mutation_run_summary <- renderText({
    df <- mutation_rv$queued
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

  output$mutation_selected_summary <- renderUI({
    row <- selected_row()
    if (is.null(row)) return(tags$p(class = "text-muted", "No sample selected."))
    ref <- row[["Reference Sequence"]]
    ref_preview <- if (nchar(ref) > 40) paste0(substr(ref, 1, 40), "...") else ref
    tags$ul(
      class = "sidebar-summary-list",
      tags$li(tags$strong("Sample:"), row[["Sample Name"]]),
      tags$li(tags$strong("Merged:"), row[["Merged File"]]),
      tags$li(tags$strong("Source:"), row[["Source"]]),
      tags$li(tags$strong("Antibody:"), ifelse(isTRUE(as.logical(row[["Antibody"]])), "Yes", "No")),
      tags$li(tags$strong("Amplicon:"), row[["Amplicon Length"]]),
      tags$li(tags$strong("Insert:"), paste(row[["Insert Start Coordinate"]], row[["Insert End Coordinate"]], sep = "-")),
      tags$li(tags$strong("Extensions:"), paste(row[["Forward Extension Type"]], row[["Reverse Extension Type"]], sep = " / ")),
      tags$li(tags$strong("Reference:"), paste0(ref_preview, " (", nchar(ref), " chars)"))
    )
  })

  observe({
    ready <- nrow(mutation_rv$queued) > 0 && all(mutation_rv$queued$Status == "Ready")
    if (ready) {
      shinyjs::enable("run_single_mutation")
      shinyjs::enable("run_batch_mutation")
    } else {
      shinyjs::disable("run_single_mutation")
      shinyjs::disable("run_batch_mutation")
    }

    batch_ready <- !is.null(batch_rv$manifest_meta) && !is.null(batch_rv$archive_meta)
    if (batch_ready) {
      shinyjs::enable("mutation_queue_batch_samples")
    } else {
      shinyjs::disable("mutation_queue_batch_samples")
    }
  })

  run_unified_mutation <- function(target = c("single", "batch")) {
    target <- match.arg(target)
    ordered <- .order_sample_df(
      mutation_rv$queued,
      mutation_rv$sample_order,
      "Mutation Queue ID"
    )

    if (target == "single") {
      mutation_rv$message <- paste("Mutation analysis pipeline is not implemented yet for", nrow(ordered), "queued samples.")
      mutation_rv$message_type <- "info"
    } else {
      batch_rv$message <- paste("Mutation analysis pipeline is not implemented yet for", nrow(ordered), "queued samples.")
      batch_rv$message_type <- "info"
    }
  }

  observeEvent(input$run_single_mutation, {
    run_unified_mutation("single")
  }, ignoreInit = TRUE)

  observeEvent(input$run_batch_mutation, {
    run_unified_mutation("batch")
  }, ignoreInit = TRUE)
}
