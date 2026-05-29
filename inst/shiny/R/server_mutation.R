setup_mutation_server <- function(input, output, session, operation_rv = NULL) {
  upload_root <- file.path(tempdir(), "tasaR_mutation_uploads", session$token)
  mutation_output_root <- file.path(tempdir(), "export", session$token)
  mutation_progress_dir <- file.path(tempdir(), "tasaR_mutation_progress", session$token)
  mutation_resource_prefix <- paste0("mutation-results-", gsub("[^A-Za-z0-9_-]", "_", session$token))
  dir.create(upload_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(mutation_output_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(mutation_progress_dir, recursive = TRUE, showWarnings = FALSE)
  shiny::addResourcePath(mutation_resource_prefix, mutation_output_root)
  if (is.null(operation_rv)) {
    operation_rv <- reactiveValues(active = NULL)
  }

  session$onSessionEnded(function() {
    unlink(upload_root, recursive = TRUE, force = TRUE)
    unlink(mutation_output_root, recursive = TRUE, force = TRUE)
    unlink(mutation_progress_dir, recursive = TRUE, force = TRUE)
    shiny::removeResourcePath(mutation_resource_prefix)
  })

  mutation_rv <- reactiveValues(
    queued = .mutation_empty_queue("Mutation Queue ID"),
    sample_order = character(),
    merged_meta = NULL,
    running = FALSE,
    batch_summarizing = FALSE,
    download_ready_count = 0L,
    download_phase = "hidden",
    results = NULL,
    message = NULL,
    message_type = NULL,
    active_settings_id = NULL
  )
  processed_progress_files <- character()
  active_mutation_run_button <- NULL
  active_mutation_queue_ids <- character()
  current_mutation_run_id <- NULL

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
    "CDR3 Start" = "cdr3_start",
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
    updateNumericInput(session, paste0(prefix, "_cdr3_start"), value = numeric_value(settings[["CDR3 Start"]]))
    updateNumericInput(session, paste0(prefix, "_fr4_start"), value = numeric_value(settings[["FR4 Start"]]))
    updateNumericInput(session, paste0(prefix, "_fr4_end"), value = numeric_value(settings[["FR4 End"]]))
  }

  reset_main_settings <- function() {
    update_settings_inputs("mutation", .mutation_default_settings())
  }

  refresh_mutation_status <- function() {
    mutation_rv$queued <- .mutation_queue_with_status(mutation_rv$queued)
  }

  operation_active <- function() {
    !is.null(operation_rv$active)
  }

  acquire_operation <- function(name) {
    if (operation_active()) return(FALSE)
    operation_rv$active <- name
    TRUE
  }

  release_operation <- function(name) {
    if (identical(operation_rv$active, name)) {
      operation_rv$active <- NULL
    }
    invisible(TRUE)
  }

  start_mutation_button_pulse <- function(button_id) {
    active_mutation_run_button <<- button_id
    shinyjs::addClass(id = button_id, class = "analysis-starting")
    updateActionButton(session, button_id, label = "Preparing")
    invisible(TRUE)
  }

  mark_mutation_button_in_progress <- function() {
    button_id <- active_mutation_run_button
    if (!is.null(button_id) && nzchar(button_id)) {
      shinyjs::removeClass(id = button_id, class = "analysis-starting")
      updateActionButton(session, button_id, label = "In Progress")
    }
    invisible(TRUE)
  }

  clear_mutation_button_pulse <- function() {
    button_id <- active_mutation_run_button
    if (!is.null(button_id) && nzchar(button_id)) {
      shinyjs::removeClass(id = button_id, class = "analysis-starting")
      updateActionButton(session, button_id, label = "Start Analysis")
    }
    active_mutation_run_button <<- NULL
    invisible(TRUE)
  }

  make_mutation_run_id <- function() {
    stamp <- format(Sys.time(), "%Y%m%d%H%M%OS6")
    suffix <- sample.int(.Machine$integer.max, 1L)
    gsub("[^0-9A-Za-z_-]", "", paste0(stamp, "_", suffix))
  }

  mutation_cancel_path <- function(run_id = current_mutation_run_id) {
    if (is.null(run_id) || !nzchar(run_id)) return(NULL)
    file.path(mutation_progress_dir, paste0("cancel_", run_id, ".flag"))
  }

  mutation_canceled_condition <- function() {
    structure(
      list(message = "Mutation analysis canceled.", call = NULL),
      class = c("tasaR_analysis_canceled", "error", "condition")
    )
  }

  check_mutation_canceled <- function(run_id) {
    cancel_path <- mutation_cancel_path(run_id)
    if (!is.null(cancel_path) && file.exists(cancel_path)) {
      stop(mutation_canceled_condition())
    }
    invisible(FALSE)
  }

  mutation_run_active <- function(run_id) {
    !is.null(run_id) && identical(run_id, current_mutation_run_id)
  }

  set_mutation_progress <- function(queue_id, progress) {
    idx <- match(queue_id, mutation_rv$queued[["Mutation Queue ID"]])
    if (is.na(idx)) return(FALSE)

    mutation_rv$queued[idx, "Progress"] <- progress
    TRUE
  }

  make_mutation_progress_callback <- function(run_id) {
    force(run_id)
    function(event) {
      check_mutation_canceled(run_id)
      if (!dir.exists(mutation_progress_dir)) {
        dir.create(mutation_progress_dir, recursive = TRUE, showWarnings = FALSE)
      }

      event$run_id <- run_id
      stamp <- format(Sys.time(), "%Y%m%d%H%M%OS6")
      sample_name <- if (is.null(event$sample_name) || !nzchar(event$sample_name)) "sample" else event$sample_name
      safe_sample <- gsub("[^A-Za-z0-9_-]", "_", sample_name)
      suffix <- paste0(format(Sys.time(), "%OS6"), "_", sample.int(.Machine$integer.max, 1L))
      suffix <- gsub("[^0-9A-Za-z_-]", "", suffix)
      path <- file.path(
        mutation_progress_dir,
        sprintf("%s_%s_%s_%s.rds", stamp, Sys.getpid(), safe_sample, suffix)
      )
      tmp <- paste0(path, ".tmp")
      saveRDS(event, tmp)
      file.rename(tmp, path)
      invisible(TRUE)
    }
  }

  mutation_progress_value <- function(event) {
    status <- if (is.null(event$status)) "step" else event$status
    label <- if (is.null(event$label)) "Processing" else event$label

    if (status %in% c("analyzed", "exporting", "finalizing")) {
      return(paste(status, label, sep = "|"))
    }

    step <- suppressWarnings(as.integer(event$step))
    total <- suppressWarnings(as.integer(event$total))
    if (is.na(step) || is.na(total) || total <= 0L) {
      return("processing")
    }

    paste("step", step, total, label, sep = "|")
  }

  set_mutation_result <- function(queue_id, result_path, zip_path) {
    idx <- match(queue_id, mutation_rv$queued[["Mutation Queue ID"]])
    if (is.na(idx)) return(FALSE)

    mutation_rv$queued[idx, "Result Path"] <- result_path
    mutation_rv$queued[idx, "Output Path"] <- zip_path
    mutation_rv$queued[idx, "Download Href"] <- paste0(
      mutation_resource_prefix,
      "/",
      utils::URLencode(basename(zip_path), reserved = TRUE)
    )
    TRUE
  }

  reset_mutation_progress <- function(queue_ids = mutation_rv$queued[["Mutation Queue ID"]]) {
    if (nrow(mutation_rv$queued) == 0) return(invisible(FALSE))

    idx <- match(queue_ids, mutation_rv$queued[["Mutation Queue ID"]])
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) return(invisible(FALSE))

    mutation_rv$queued[idx, "Progress"] <- ""
    mutation_rv$queued[idx, "Result Path"] <- ""
    mutation_rv$queued[idx, "Output Path"] <- ""
    mutation_rv$queued[idx, "Download Href"] <- ""
    invisible(TRUE)
  }

  reset_mutation_download_state <- function() {
    mutation_rv$download_ready_count <- 0L
    mutation_rv$download_phase <- "hidden"
    invisible(TRUE)
  }

  update_mutation_download_progress <- function() {
    count <- mutation_rv$download_ready_count + 1L
    mutation_rv$download_ready_count <- count
    invisible(count)
  }

  set_mutation_download_phase <- function(phase) {
    mutation_rv$download_phase <- phase
    invisible(phase)
  }

  downloadable_mutation_rows <- function() {
    df <- mutation_rv$queued
    if (is.null(df) || nrow(df) == 0 || !"Progress" %in% names(df)) {
      return(df[0, , drop = FALSE])
    }

    output_path <- if ("Output Path" %in% names(df)) df[["Output Path"]] else rep("", nrow(df))
    result_path <- if ("Result Path" %in% names(df)) df[["Result Path"]] else rep("", nrow(df))
    progress <- df[["Progress"]]
    complete <- tolower(trimws(progress)) %in% c("done", "complete", "completed", "success", "true")
    has_output <- nzchar(output_path) & file.exists(output_path)
    has_result_dir <- nzchar(result_path) & dir.exists(result_path)
    keep <- complete & (has_output | has_result_dir)
    df[keep, , drop = FALSE]
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
      rows[["Result Path"]] <- ""
      rows[["Output Path"]] <- ""
      rows[["Download Href"]] <- ""

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
    if (isTRUE(mutation_rv$running) || operation_active()) {
      mutation_rv$message <- "An analysis is already running."
      mutation_rv$message_type <- "error"
      return()
    }

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
    row[["Result Path"]] <- ""
    row[["Output Path"]] <- ""
    row[["Download Href"]] <- ""

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
    if (isTRUE(mutation_rv$running) || operation_active()) {
      mutation_rv$message <- "An analysis is already running."
      mutation_rv$message_type <- "error"
      return()
    }

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
    if (isTRUE(mutation_rv$running) || operation_active()) {
      batch_rv$message <- "An analysis is already running."
      batch_rv$message_type <- "error"
      return()
    }

    if (is.null(batch_rv$manifest_meta) || is.null(batch_rv$archive_meta)) {
      batch_rv$message <- "Upload the batch manifest and archive before adding samples."
      batch_rv$message_type <- "error"
      return()
    }

    append_mutation_batch()
  }, ignoreInit = TRUE)

  observeEvent(input$mutation_clear_batch, {
    if (isTRUE(mutation_rv$running) || operation_active()) {
      batch_rv$message <- "An analysis is already running."
      batch_rv$message_type <- "error"
      return()
    }

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
    if (isTRUE(mutation_rv$running) || operation_active()) {
      mutation_rv$message <- "An analysis is already running."
      mutation_rv$message_type <- "error"
      return()
    }

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

  output$mutation_downloads <- renderUI({
    phase <- mutation_rv$download_phase
    if (identical(phase, "hidden")) return(NULL)

    download_ready <- identical(phase, "ready")
    button <- if (download_ready) {
      downloadButton(
        "download_all_mutation_results",
        "Download All",
        class = "download-all-btn download-all-ready"
      )
    } else {
      button_class <- "download-all-btn download-all-preparing download-all-pulsing"
      button_label <- "Download All"
      if (identical(phase, "preparing")) {
        button_label <- "Preparing"
      }

      tags$button(
        type = "button",
        class = button_class,
        disabled = "disabled",
        button_label
      )
    }

    div(
      class = "merge-download-row",
      div(class = "merge-download-count", textOutput("mutation_download_count", inline = TRUE)),
      div(class = "action-button-group", button)
    )
  })

  output$mutation_download_count <- renderText({
    ready_count <- mutation_rv$download_ready_count
    paste(ready_count, if (ready_count == 1L) "folder ready" else "folders ready")
  })

  output$download_all_mutation_results <- downloadHandler(
    filename = function() {
      paste0("mutation_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      files <- list.files(mutation_output_root, recursive = TRUE, full.names = FALSE)
      files <- files[!grepl("\\.zip$", files, ignore.case = TRUE)]
      if (length(files) == 0L) {
        stop("No completed mutation result files are available for download.", call. = FALSE)
      }

      zip::zipr(
        zipfile = file,
        files = files,
        recurse = FALSE,
        compression_level = 9,
        include_directories = TRUE,
        root = mutation_output_root,
        mode = "mirror"
      )
    }
  )

  observeEvent(mutation_rv$message, {
    req(!is.null(mutation_rv$message))
    .clear_queue_message_later(mutation_rv)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  observeEvent(batch_rv$message, {
    req(!is.null(batch_rv$message))
    .clear_queue_message_later(batch_rv)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  observe({
    if (!isTRUE(mutation_rv$running)) return()
    invalidateLater(250, session)
    if (!dir.exists(mutation_progress_dir)) return()

    files <- list.files(mutation_progress_dir, pattern = "\\.rds$", full.names = TRUE)
    files <- sort(setdiff(files, processed_progress_files))
    if (length(files) == 0L) return()

    for (path in files) {
      event <- tryCatch(readRDS(path), error = function(e) NULL)
      processed_progress_files <<- c(processed_progress_files, path)
      if (is.null(event) || is.null(event$sample_name) || !nzchar(event$sample_name)) next
      if (!identical(event$run_id, current_mutation_run_id)) next

      idx <- match(event$sample_name, mutation_rv$queued[["Sample Name"]])
      if (is.na(idx)) next

      mark_mutation_button_in_progress()

      current <- tolower(trimws(as.character(mutation_rv$queued[idx, "Progress", drop = TRUE])))
      if (current %in% c("done", "complete", "completed", "success", "true")) next

      mutation_rv$queued[idx, "Progress"] <- mutation_progress_value(event)
    }
  })

  .setup_reorderable_sample_table(
    input = input,
    output = output,
    output_id = "mutation_queued_samples",
    data_fn = function() mutation_rv$queued,
    id_col = "Mutation Queue ID",
    get_order = function() mutation_rv$sample_order,
    set_order = function(order) {
      if (isTRUE(mutation_rv$running) || operation_active()) return()
      mutation_rv$sample_order <- order
      refresh_mutation_status()
    },
    remove_input_id = "mutation_remove_row",
    move_up_input_id = "mutation_queued_samples_move_up",
    move_down_input_id = "mutation_queued_samples_move_down",
    table_builder = function(...) .make_mutation_table_data(..., settings_input_id = "mutation_open_settings"),
    small_text_columns = c("Merged File")
  )

  open_settings_modal <- function(row_id) {
    if (isTRUE(mutation_rv$running) || operation_active()) {
      mutation_rv$message <- "An analysis is already running."
      mutation_rv$message_type <- "error"
      return()
    }

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
    if (isTRUE(mutation_rv$running) || operation_active()) {
      mutation_rv$message <- "An analysis is already running."
      mutation_rv$message_type <- "error"
      removeModal()
      return()
    }

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
      queued <- ready <- issues <- analyzed <- 0
    } else {
      queued <- nrow(df)
      ready <- sum(df$Status == "Ready", na.rm = TRUE)
      issues <- queued - ready
      analyzed <- .mutation_analyzed_count(df)
    }

    paste(
      "Queued:", queued,
      "\nReady:", ready,
      "\nIssues:", issues,
      "\nAnalyzed:", analyzed
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

  output$mutation_stop_control <- renderUI({
    if (!isTRUE(mutation_rv$running)) return(NULL)
    div(
      class = "mutation-stop-row",
      tags$button(
        id = "mutation_stop_analysis",
        type = "button",
        class = "mutation-stop-btn",
        onclick = paste0(
          "if (window.confirm('Stop mutation analysis? Results of in-progress analysis may be lost.')) {",
          "this.disabled = true;",
          "this.textContent = 'Stopping...';",
          "Shiny.setInputValue('mutation_cancel_requested', Date.now(), {priority: 'event'});",
          "}"
        ),
        icon("stop"), "Stop Analysis"
      )
    )
  })

  cancel_mutation_run <- function() {
    run_id <- current_mutation_run_id
    cancel_path <- mutation_cancel_path(run_id)
    if (!is.null(cancel_path)) {
      dir.create(dirname(cancel_path), recursive = TRUE, showWarnings = FALSE)
      writeLines("cancel", cancel_path)
    }

    clear_mutation_processing(active_mutation_queue_ids)
    mutation_rv$running <- FALSE
    mutation_rv$batch_summarizing <- FALSE
    reset_mutation_download_state()
    current_mutation_run_id <<- NULL
    active_mutation_queue_ids <<- character()
    clear_mutation_button_pulse()
    release_operation("mutation")
    mutation_rv$message <- "Mutation analysis canceled."
    mutation_rv$message_type <- "info"
    invisible(TRUE)
  }

  observeEvent(input$mutation_cancel_requested, {
    if (!isTRUE(mutation_rv$running)) return()
    cancel_mutation_run()
  }, ignoreInit = TRUE)

  observe({
    ready <- !isTRUE(mutation_rv$running) && !operation_active() && nrow(mutation_rv$queued) > 0 && all(mutation_rv$queued$Status == "Ready")
    if (ready) {
      shinyjs::enable("run_single_mutation")
      shinyjs::enable("run_batch_mutation")
    } else {
      shinyjs::disable("run_single_mutation")
      shinyjs::disable("run_batch_mutation")
    }

    batch_ready <- !isTRUE(mutation_rv$running) && !operation_active() && !is.null(batch_rv$manifest_meta) && !is.null(batch_rv$archive_meta)
    if (batch_ready) {
      shinyjs::enable("mutation_queue_batch_samples")
    } else {
      shinyjs::disable("mutation_queue_batch_samples")
    }
  })

  set_mutation_run_message <- function(target, message, type) {
    if (target == "single") {
      mutation_rv$message <- message
      mutation_rv$message_type <- type
    } else {
      batch_rv$message <- message
      batch_rv$message_type <- type
    }
  }

  clear_mutation_processing <- function(queue_ids) {
    idx <- match(queue_ids, mutation_rv$queued[["Mutation Queue ID"]])
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0L) return(invisible(FALSE))

    progress_raw <- trimws(as.character(mutation_rv$queued[idx, "Progress", drop = TRUE]))
    progress <- tolower(progress_raw)
    idx <- idx[
      progress %in% c("processing", "running", "in_progress") |
        startsWith(progress_raw, "step|") |
        startsWith(progress_raw, "exporting|") |
        startsWith(progress_raw, "finalizing|") |
        startsWith(progress_raw, "analyzed|")
    ]
    if (length(idx) == 0L) return(invisible(FALSE))

    mutation_rv$queued[idx, "Progress"] <- ""
    mutation_rv$queued[idx, "Result Path"] <- ""
    mutation_rv$queued[idx, "Output Path"] <- ""
    mutation_rv$queued[idx, "Download Href"] <- ""
    invisible(TRUE)
  }

  fail_mutation_run <- function(target, reason, queue_ids, run_id) {
    if (!mutation_run_active(run_id)) return(invisible(NULL))
    clear_mutation_processing(queue_ids)
    mutation_rv$running <- FALSE
    mutation_rv$batch_summarizing <- FALSE
    reset_mutation_download_state()
    current_mutation_run_id <<- NULL
    active_mutation_queue_ids <<- character()
    clear_mutation_button_pulse()
    release_operation("mutation")
    set_mutation_run_message(target, conditionMessage(reason), "error")
    invisible(NULL)
  }

  complete_mutation_run <- function(job_rows, target, samples, amps, summaries, batch_summary = NULL, run_id) {
    if (!mutation_run_active(run_id)) return(invisible(NULL))
    mutation_rv$results <- list(
      samples = samples,
      amps = amps,
      summaries = summaries,
      batch_summary = batch_summary
    )
    mutation_rv$running <- FALSE
    mutation_rv$batch_summarizing <- FALSE
    set_mutation_download_phase("ready")
    current_mutation_run_id <<- NULL
    active_mutation_queue_ids <<- character()
    clear_mutation_button_pulse()
    release_operation("mutation")
    set_mutation_run_message(
      target,
      paste("Mutation analysis complete for", nrow(job_rows), "queued samples."),
      "success"
    )
    invisible(TRUE)
  }

  finish_mutation_outputs <- function(job_rows, target, samples, amps, summaries, run_id) {
    if (!mutation_run_active(run_id)) return(invisible(NULL))
    if (nrow(job_rows) <= 1L) {
      return(complete_mutation_run(job_rows, target, samples, amps, summaries, run_id = run_id))
    }

    mutation_rv$batch_summarizing <- TRUE
    set_mutation_download_phase("preparing")
    batch_job <- promises::future_promise({
      check_mutation_canceled(run_id)
      tasaR::batchSummarize(
        amps,
        export = TRUE,
        path = mutation_output_root,
        suppressConsoleOutput = TRUE
      )
    }, seed = TRUE)

    promises::then(
      batch_job,
      onFulfilled = function(value) {
        if (!mutation_run_active(run_id)) return(invisible(value))
        complete_mutation_run(job_rows, target, samples, amps, summaries, value, run_id = run_id)
        invisible(value)
      },
      onRejected = function(reason) {
        if (mutation_run_active(run_id)) {
          mutation_rv$batch_summarizing <- FALSE
        }
        fail_mutation_run(target, reason, job_rows[["Mutation Queue ID"]], run_id)
      }
    )
  }

  run_next_mutation_output <- function(job_rows, target, amps, samples, summaries, run_id, index = 1L) {
    if (!mutation_run_active(run_id)) return(invisible(NULL))
    if (index > nrow(job_rows)) {
      finish_mutation_outputs(job_rows, target, samples, amps, summaries, run_id)
      return(invisible(TRUE))
    }

    row <- job_rows[index, , drop = FALSE]
    queue_id <- row[["Mutation Queue ID"]]
    amp <- amps[[index]]
    set_mutation_progress(queue_id, "exporting|Writing summary files")
    progress_callback <- make_mutation_progress_callback(run_id)

    output_job <- promises::future_promise({
      progress_callback(list(
        sample_name = row[["Sample Name"]],
        status = "exporting",
        label = "Writing summary files"
      ))
      summary <- tasaR::Summarize(
        amp,
        export = TRUE,
        path = mutation_output_root,
        suppressConsoleOutput = TRUE
      )
      result_path <- file.path(mutation_output_root, row[["Sample Name"]])
      sample_zip <- file.path(mutation_output_root, .mutation_sample_zip_name(row[["Sample Name"]]))
      progress_callback(list(
        sample_name = row[["Sample Name"]],
        status = "finalizing",
        label = "Compressing results"
      ))
      check_mutation_canceled(run_id)
      .zip_directory(sample_zip, result_path)
      check_mutation_canceled(run_id)

      list(
        amp = amp,
        summary = summary,
        result_path = result_path,
        zip_path = sample_zip
      )
    }, seed = TRUE)

    promises::then(
      output_job,
      onFulfilled = function(value) {
        if (!mutation_run_active(run_id)) return(invisible(value))
        samples[[index]] <- value
        summaries[[index]] <- value$summary
        result_ready <- set_mutation_result(queue_id, value$result_path, value$zip_path)
        set_mutation_progress(queue_id, "done")
        if (isTRUE(result_ready)) {
          update_mutation_download_progress()
        }
        run_next_mutation_output(job_rows, target, amps, samples, summaries, run_id, index + 1L)
        invisible(value)
      },
      onRejected = function(reason) {
        fail_mutation_run(target, reason, job_rows[["Mutation Queue ID"]], run_id)
      }
    )
  }

  start_mutation_outputs <- function(job_rows, target, amps, run_id) {
    if (!mutation_run_active(run_id)) return(invisible(NULL))
    sample_names <- job_rows[["Sample Name"]]
    amps <- as.list(amps)
    if (length(amps) != nrow(job_rows)) {
      fail_mutation_run(
        target,
        simpleError("Mutation analysis returned an unexpected number of results."),
        job_rows[["Mutation Queue ID"]],
        run_id
      )
      return(invisible(FALSE))
    }
    if (is.null(names(amps)) || any(!nzchar(names(amps)))) {
      names(amps) <- sample_names
    }

    samples <- vector("list", nrow(job_rows))
    summaries <- vector("list", nrow(job_rows))
    names(samples) <- sample_names
    names(summaries) <- sample_names

    run_next_mutation_output(job_rows, target, amps, samples, summaries, run_id)
  }

  run_single_mutation_analysis <- function(job_rows, target, run_id) {
    queue_id <- job_rows[["Mutation Queue ID"]][[1]]
    set_mutation_progress(queue_id, "processing")
    progress_callback <- make_mutation_progress_callback(run_id)

    analysis_job <- promises::future_promise({
      tasaR::analyzeAmplicon(
        .prep_tas_settings(job_rows[1, , drop = FALSE]),
        progress_callback = progress_callback
      )
    }, seed = TRUE)

    promises::then(
      analysis_job,
      onFulfilled = function(value) {
        if (!mutation_run_active(run_id)) return(invisible(value))
        amps <- list(value)
        names(amps) <- job_rows[["Sample Name"]]
        start_mutation_outputs(job_rows, target, amps, run_id)
        invisible(value)
      },
      onRejected = function(reason) {
        fail_mutation_run(target, reason, job_rows[["Mutation Queue ID"]], run_id)
      }
    )
  }

  run_batch_mutation_analysis <- function(job_rows, target, run_id) {
    progress_callback <- make_mutation_progress_callback(run_id)
    analysis_job <- promises::future_promise({
      tasaR::batchAnalyzeAmplicon(
        .mutation_settings_df(job_rows),
        progress_callback = progress_callback
      )
    }, seed = TRUE)

    promises::then(
      analysis_job,
      onFulfilled = function(value) {
        if (!mutation_run_active(run_id)) return(invisible(value))
        start_mutation_outputs(job_rows, target, value, run_id)
        invisible(value)
      },
      onRejected = function(reason) {
        fail_mutation_run(target, reason, job_rows[["Mutation Queue ID"]], run_id)
      }
    )
  }

  run_unified_mutation <- function(target = c("single", "batch")) {
    target <- match.arg(target)

    if (isTRUE(mutation_rv$running)) {
      mutation_rv$message <- "A mutation analysis is already running."
      mutation_rv$message_type <- "error"
      return(invisible(FALSE))
    }
    if (operation_active()) {
      mutation_rv$message <- "Another operation is already running."
      mutation_rv$message_type <- "error"
      return(invisible(FALSE))
    }

    ordered <- .order_sample_df(
      mutation_rv$queued,
      mutation_rv$sample_order,
      "Mutation Queue ID"
    )

    if (is.null(ordered) || nrow(ordered) == 0L) {
      mutation_rv$message <- "No mutation samples are queued."
      mutation_rv$message_type <- "error"
      return(invisible(FALSE))
    }

    if (!all(ordered$Status == "Ready")) {
      mutation_rv$message <- "Resolve queue issues before starting mutation analysis."
      mutation_rv$message_type <- "error"
      return(invisible(FALSE))
    }

    if (!acquire_operation("mutation")) {
      mutation_rv$message <- "Another operation is already running."
      mutation_rv$message_type <- "error"
      return(invisible(FALSE))
    }

    reset_mutation_progress(ordered[["Mutation Queue ID"]])
    dir.create(mutation_progress_dir, recursive = TRUE, showWarnings = FALSE)
    old_progress_files <- list.files(mutation_progress_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(old_progress_files) > 0L) {
      unlink(old_progress_files, force = TRUE)
    }
    processed_progress_files <<- character()
    current_mutation_run_id <<- make_mutation_run_id()
    active_mutation_queue_ids <<- ordered[["Mutation Queue ID"]]
    mutation_rv$batch_summarizing <- FALSE
    reset_mutation_download_state()
    mutation_rv$running <- TRUE
    mutation_rv$results <- NULL
    if (nrow(ordered) == 1L) {
      run_single_mutation_analysis(ordered, target, current_mutation_run_id)
    } else {
      run_batch_mutation_analysis(ordered, target, current_mutation_run_id)
    }
    invisible(TRUE)
  }

  observeEvent(input$run_single_mutation, {
    start_mutation_button_pulse("run_single_mutation")
    if (!isTRUE(run_unified_mutation("single"))) {
      clear_mutation_button_pulse()
    }
  }, ignoreInit = TRUE)

  observeEvent(input$run_batch_mutation, {
    start_mutation_button_pulse("run_batch_mutation")
    if (!isTRUE(run_unified_mutation("batch"))) {
      clear_mutation_button_pulse()
    }
  }, ignoreInit = TRUE)
}
