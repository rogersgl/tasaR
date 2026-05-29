setup_merge_server <- function(input, output, session, operation_rv = NULL) {
  upload_root <- file.path(tempdir(), "tasaR_uploads", session$token)
  merge_output_root <- file.path(tempdir(), "merge", session$token)
  merge_resource_prefix <- paste0("merge-results-", gsub("[^A-Za-z0-9_-]", "_", session$token))
  dir.create(upload_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(merge_output_root, recursive = TRUE, showWarnings = FALSE)
  shiny::addResourcePath(merge_resource_prefix, merge_output_root)
  if (is.null(operation_rv)) {
    operation_rv <- reactiveValues(active = NULL)
  }

  session$onSessionEnded(function() {
    unlink(upload_root, recursive = TRUE, force = TRUE)
    unlink(merge_output_root, recursive = TRUE, force = TRUE)
    shiny::removeResourcePath(merge_resource_prefix)
  })

  merge_rv <- reactiveValues(
    queued = data.frame(
      "Queue ID" = character(),
      "Sample Name" = character(),
      "Left (R1) File" = character(),
      "Right (R2) File" = character(),
      "R1 Path" = character(),
      "R2 Path" = character(),
      "Source" = character(),
      "Status" = character(),
      "Progress" = character(),
      "Output Path" = character(),
      "Download Href" = character(),
      check.names = FALSE
    ),
    sample_order = character(),
    r1_meta = NULL,
    r2_meta = NULL,
    files_ready = FALSE,
    running = FALSE,
    message = NULL,
    message_type = NULL
  )

  batch_rv <- reactiveValues(
    manifest_meta = NULL,
    archive_meta = NULL,
    extract_dir = NULL,
    message = NULL,
    message_type = NULL
  )

  refresh_merge_ready <- function() {
    merge_rv$files_ready <- !is.null(merge_rv$r1_meta) && !is.null(merge_rv$r2_meta)
  }

  refresh_merge_status <- function() {
    merge_rv$queued <- .sample_table_with_status(merge_rv$queued)
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

  set_merge_progress <- function(queue_id, progress) {
    idx <- match(queue_id, merge_rv$queued[["Queue ID"]])
    if (is.na(idx)) return(FALSE)

    merge_rv$queued[idx, "Progress"] <- progress
    TRUE
  }

  set_merge_result <- function(queue_id, output_path) {
    idx <- match(queue_id, merge_rv$queued[["Queue ID"]])
    if (is.na(idx)) return(FALSE)

    filename <- basename(output_path)
    merge_rv$queued[idx, "Output Path"] <- output_path
    merge_rv$queued[idx, "Download Href"] <- paste0(
      merge_resource_prefix,
      "/",
      utils::URLencode(filename, reserved = TRUE)
    )
    TRUE
  }

  reset_merge_progress <- function(queue_ids = merge_rv$queued[["Queue ID"]]) {
    if (nrow(merge_rv$queued) == 0) return(invisible(FALSE))

    idx <- match(queue_ids, merge_rv$queued[["Queue ID"]])
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) return(invisible(FALSE))

    merge_rv$queued[idx, "Progress"] <- ""
    merge_rv$queued[idx, "Output Path"] <- ""
    merge_rv$queued[idx, "Download Href"] <- ""
    invisible(TRUE)
  }

  downloadable_merge_rows <- function() {
    df <- merge_rv$queued
    if (is.null(df) || nrow(df) == 0 || !"Output Path" %in% names(df)) {
      return(df[0, , drop = FALSE])
    }

    path <- df[["Output Path"]]
    progress <- if ("Progress" %in% names(df)) df[["Progress"]] else rep("", nrow(df))
    keep <- nzchar(path) & file.exists(path) & tolower(trimws(progress)) %in% c("done", "complete", "completed", "success", "true")
    df[keep, , drop = FALSE]
  }

  reconcile_merge_order <- function() {
    merge_rv$sample_order <- .reconcile_sample_order(
      merge_rv$queued[["Queue ID"]],
      merge_rv$sample_order
    )
  }

  ensure_batch_extract_dir <- function() {
    if (is.null(batch_rv$archive_meta)) return(NULL)

    batch_rv$extract_dir <- file.path(
      upload_root,
      paste0("batch_extract_", format(Sys.time(), "%Y%m%d%H%M%OS3"))
    )
    dir.create(batch_rv$extract_dir, recursive = TRUE, showWarnings = FALSE)
    batch_rv$extract_dir
  }

  append_merge_batch <- function() {
    if (is.null(batch_rv$manifest_meta) || is.null(batch_rv$archive_meta)) {
      return()
    }

    tryCatch({
      extract_dir <- ensure_batch_extract_dir()
      dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

      utils::unzip(batch_rv$archive_meta$managed_path, exdir = extract_dir)

      manifest_df <- .read_batch_manifest(batch_rv$manifest_meta$managed_path)
      extracted_files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
      extracted_basenames <- basename(extracted_files)

      r1_idx <- match(manifest_df[["Left (R1) File"]], extracted_basenames)
      r2_idx <- match(manifest_df[["Right (R2) File"]], extracted_basenames)
      missing <- is.na(r1_idx) | is.na(r2_idx)

      if (any(missing)) {
        stop(paste(
          "Missing FASTQ files for:",
          paste(manifest_df[missing, "Sample Name", drop = TRUE], collapse = ", ")
        ))
      }

      batch_id <- format(Sys.time(), "%Y%m%d%H%M%OS3")
      rows <- data.frame(
        "Queue ID" = paste0("batch_", batch_id, "_", seq_len(nrow(manifest_df))),
        "Sample Name" = manifest_df[["Sample Name"]],
        "Left (R1) File" = manifest_df[["Left (R1) File"]],
        "Right (R2) File" = manifest_df[["Right (R2) File"]],
        "R1 Path" = extracted_files[r1_idx],
        "R2 Path" = extracted_files[r2_idx],
        "Source" = "batch",
        "Status" = "",
        "Progress" = "",
        "Output Path" = "",
        "Download Href" = "",
        check.names = FALSE
      )

      merge_rv$queued <- rbind(merge_rv$queued, rows)
      refresh_merge_status()
      reconcile_merge_order()

      batch_rv$message <- paste("Queued", nrow(rows), "batch samples successfully.")
      batch_rv$message_type <- "success"
    }, error = function(e) {
      batch_rv$message <- conditionMessage(e)
      batch_rv$message_type <- "error"
    })
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
    if (isTRUE(merge_rv$running) || operation_active()) {
      merge_rv$message <- "A merge is already running."
      merge_rv$message_type <- "error"
      return()
    }

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

    queue_id <- paste0("single_", nrow(merge_rv$queued) + 1L, "_", format(Sys.time(), "%H%M%OS3"))
    merge_rv$queued <- rbind(
      merge_rv$queued,
      data.frame(
        "Queue ID" = queue_id,
        "Sample Name" = sample_name,
        "Left (R1) File" = merge_rv$r1_meta$original_name,
        "Right (R2) File" = merge_rv$r2_meta$original_name,
        "R1 Path" = merge_rv$r1_meta$managed_path,
        "R2 Path" = merge_rv$r2_meta$managed_path,
        "Source" = "single",
        "Status" = "",
        "Progress" = "",
        "Output Path" = "",
        "Download Href" = "",
        check.names = FALSE
      )
    )
    refresh_merge_status()
    reconcile_merge_order()

    updateTextInput(session, "sample_name", value = "")
    shinyjs::reset("fastq_r1")
    shinyjs::reset("fastq_r2")

    merge_rv$r1_meta <- NULL
    merge_rv$r2_meta <- NULL
    merge_rv$files_ready <- FALSE

    merge_rv$message <- paste("Queued sample", sample_name, "successfully.")
    merge_rv$message_type <- "success"
  })

  observeEvent(input$batch_manifest, {
    .remove_managed_upload(batch_rv$manifest_meta)
    batch_rv$manifest_meta <- .save_managed_upload(input$batch_manifest, "batch_manifest", upload_root)
    batch_rv$message <- "Batch manifest uploaded."
    batch_rv$message_type <- "success"
  }, ignoreInit = TRUE)

  observeEvent(input$batch_archive, {
    .remove_managed_upload(batch_rv$archive_meta)
    batch_rv$extract_dir <- NULL
    batch_rv$archive_meta <- .save_managed_upload(input$batch_archive, "batch_archive", upload_root)
    batch_rv$message <- "Batch archive uploaded."
    batch_rv$message_type <- "success"
  }, ignoreInit = TRUE)

  observeEvent(input$queue_batch_samples, {
    if (isTRUE(merge_rv$running) || operation_active()) {
      batch_rv$message <- "A merge is already running."
      batch_rv$message_type <- "error"
      return()
    }

    if (is.null(batch_rv$manifest_meta) || is.null(batch_rv$archive_meta)) {
      batch_rv$message <- "Upload the batch manifest and archive before adding samples."
      batch_rv$message_type <- "error"
      return()
    }

    append_merge_batch()
  }, ignoreInit = TRUE)

  observeEvent(input$clear_batch, {
    if (isTRUE(merge_rv$running) || operation_active()) {
      batch_rv$message <- "A merge is already running."
      batch_rv$message_type <- "error"
      return()
    }

    .remove_managed_upload(batch_rv$manifest_meta)
    .remove_managed_upload(batch_rv$archive_meta)

    batch_rv$manifest_meta <- NULL
    batch_rv$archive_meta <- NULL
    batch_rv$extract_dir <- NULL

    batch_rv$message <- "Batch inputs cleared."
    batch_rv$message_type <- "info"

    shinyjs::reset("batch_manifest")
    shinyjs::reset("batch_archive")
  }, ignoreInit = TRUE)

  observeEvent(input$remove_queue_id, {
    if (isTRUE(merge_rv$running) || operation_active()) {
      merge_rv$message <- "A merge is already running."
      merge_rv$message_type <- "error"
      return()
    }

    res <- .remove_queued_sample(merge_rv$queued, input$remove_queue_id)
    if (!res$success) return()

    row <- res$removed
    .remove_managed_upload(list(managed_path = row[["R1 Path"]]))
    .remove_managed_upload(list(managed_path = row[["R2 Path"]]))

    merge_rv$queued <- res$queued
    refresh_merge_status()
    reconcile_merge_order()

    merge_rv$message <- paste("Removed sample", row[["Sample Name"]])
    merge_rv$message_type <- "info"
  }, ignoreInit = TRUE)

  observeEvent(merge_rv$message, {
    req(!is.null(merge_rv$message))
    .clear_queue_message_later(merge_rv)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  observeEvent(batch_rv$message, {
    req(!is.null(batch_rv$message))
    .clear_queue_message_later(batch_rv)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  output$merge_message <- renderUI({
    .queue_upload_message(merge_rv$message, merge_rv$message_type)
  })

  output$batch_message <- renderUI({
    .queue_upload_message(batch_rv$message, batch_rv$message_type)
  })

  output$merge_downloads <- renderUI({
    ready_count <- nrow(downloadable_merge_rows())
    if (ready_count == 0L) return(NULL)

    label <- paste(ready_count, if (ready_count == 1L) "file ready" else "files ready")
    button <- downloadButton("download_all_merged", "Download All", class = "download-all-btn")

    div(
      class = "merge-download-row",
      div(class = "merge-download-count", label),
      div(class = "action-button-group", if (ready_count > 0L && !isTRUE(merge_rv$running)) button else shinyjs::disabled(button))
    )
  })

  output$download_all_merged <- downloadHandler(
    filename = function() {
      paste0("merged_reads_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      rows <- downloadable_merge_rows()
      paths <- rows[["Output Path"]]
      paths <- paths[nzchar(paths) & file.exists(paths)]
      if (length(paths) == 0L) {
        stop("No completed merged files are available for download.", call. = FALSE)
      }

      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(dirname(paths[[1]]))
      utils::zip(zipfile = file, files = basename(paths), flags = "-qr9X")
    }
  )

  .setup_reorderable_sample_table(
    input = input,
    output = output,
    output_id = "merge_queued_samples",
    data_fn = function() merge_rv$queued,
    id_col = "Queue ID",
    get_order = function() merge_rv$sample_order,
    set_order = function(order) {
      if (isTRUE(merge_rv$running) || operation_active()) return()
      merge_rv$sample_order <- order
    },
    remove_input_id = "remove_queue_id",
    move_up_input_id = "merge_queued_samples_move_up",
    move_down_input_id = "merge_queued_samples_move_down",
    small_text_columns = c("Left (R1) File", "Right (R2) File")
  )

  observe({
    ready <- !isTRUE(merge_rv$running) && !operation_active() && nrow(merge_rv$queued) > 0 && all(merge_rv$queued$Status == "Ready")
    if (ready) {
      shinyjs::enable("run_single_merge")
      shinyjs::enable("run_batch_merge")
    } else {
      shinyjs::disable("run_single_merge")
      shinyjs::disable("run_batch_merge")
    }
  })

  observe({
    batch_ready <- !isTRUE(merge_rv$running) && !operation_active() && !is.null(batch_rv$manifest_meta) && !is.null(batch_rv$archive_meta)
    if (batch_ready) {
      shinyjs::enable("queue_batch_samples")
    } else {
      shinyjs::disable("queue_batch_samples")
    }
  })

  run_next_merge <- function(job_rows, settings, index = 1L) {
    if (index > nrow(job_rows)) {
      merge_rv$running <- FALSE
      release_operation("merge")
      merge_rv$message <- paste("Merge complete for", nrow(job_rows), "queued samples.")
      merge_rv$message_type <- "success"
      return(invisible(TRUE))
    }

    row <- job_rows[index, , drop = FALSE]
    queue_id <- row[["Queue ID"]]
    set_merge_progress(queue_id, "processing")

    merge_job <- promises::future_promise({
      .run_single_merge_sample(
        row[, c("Sample_Name", "R1_Path", "R2_Path", "Source", "Output_Path", "Log_Path"), drop = FALSE],
        settings = settings,
        merge_fun = tasaR::pandaseq_merge_files
      )
    }, seed = TRUE)

    promises::then(
      merge_job,
      onFulfilled = function(value) {
        set_merge_progress(queue_id, "done")
        set_merge_result(queue_id, row[["Output_Path"]])
        run_next_merge(job_rows, settings, index + 1L)
        invisible(value)
      },
      onRejected = function(reason) {
        set_merge_progress(queue_id, "")
        merge_rv$running <- FALSE
        release_operation("merge")
        merge_rv$message <- conditionMessage(reason)
        merge_rv$message_type <- "error"
        invisible(NULL)
      }
    )
  }

  run_unified_merge <- function() {
    if (isTRUE(merge_rv$running)) {
      merge_rv$message <- "A merge is already running."
      merge_rv$message_type <- "error"
      return(invisible(FALSE))
    }
    if (operation_active()) {
      merge_rv$message <- "Another operation is already running."
      merge_rv$message_type <- "error"
      return(invisible(FALSE))
    }

    ordered_queue <- .order_sample_df(
      merge_rv$queued,
      merge_rv$sample_order,
      "Queue ID"
    )

    merge_inputs <- tryCatch(
      .prepare_merge_inputs(mode = "single", queued_df = ordered_queue),
      error = function(e) {
        merge_rv$message <- conditionMessage(e)
        merge_rv$message_type <- "error"
        NULL
      }
    )

    if (is.null(merge_inputs)) return(invisible(FALSE))

    settings <- .get_merge_settings(input)
    job_rows <- .merge_job_rows(
      ordered_queue = ordered_queue,
      merge_inputs = merge_inputs,
      settings = settings,
      output_dir = merge_output_root
    )

    if (!acquire_operation("merge")) {
      merge_rv$message <- "Another operation is already running."
      merge_rv$message_type <- "error"
      return(invisible(FALSE))
    }

    reset_merge_progress(job_rows[["Queue ID"]])
    merge_rv$running <- TRUE
    run_next_merge(job_rows, settings)
    invisible(TRUE)
  }

  observeEvent(input$run_single_merge, {
    run_unified_merge()
  }, ignoreInit = TRUE)

  observeEvent(input$run_batch_merge, {
    run_unified_merge()
  }, ignoreInit = TRUE)
}
