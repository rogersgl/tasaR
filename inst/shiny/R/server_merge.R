setup_merge_server <- function(input, output, session) {
  upload_root <- file.path(tempdir(), "tasaR_uploads", session$token)
  dir.create(upload_root, recursive = TRUE, showWarnings = FALSE)

  session$onSessionEnded(function() {
    unlink(upload_root, recursive = TRUE, force = TRUE)
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
      check.names = FALSE
    ),
    sample_order = character(),
    r1_meta = NULL,
    r2_meta = NULL,
    files_ready = FALSE,
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
    if (is.null(batch_rv$manifest_meta) || is.null(batch_rv$archive_meta)) {
      batch_rv$message <- "Upload the batch manifest and archive before adding samples."
      batch_rv$message_type <- "error"
      return()
    }

    append_merge_batch()
  }, ignoreInit = TRUE)

  observeEvent(input$clear_batch, {
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

  .setup_reorderable_sample_table(
    input = input,
    output = output,
    output_id = "merge_queued_samples",
    data_fn = function() merge_rv$queued,
    id_col = "Queue ID",
    get_order = function() merge_rv$sample_order,
    set_order = function(order) {
      merge_rv$sample_order <- order
    },
    remove_input_id = "remove_queue_id",
    move_up_input_id = "merge_queued_samples_move_up",
    move_down_input_id = "merge_queued_samples_move_down"
  )

  observe({
    ready <- nrow(merge_rv$queued) > 0 && all(merge_rv$queued$Status == "Ready")
    if (ready) {
      shinyjs::enable("run_single_merge")
      shinyjs::enable("run_batch_merge")
    } else {
      shinyjs::disable("run_single_merge")
      shinyjs::disable("run_batch_merge")
    }
  })

  observe({
    batch_ready <- !is.null(batch_rv$manifest_meta) && !is.null(batch_rv$archive_meta)
    if (batch_ready) {
      shinyjs::enable("queue_batch_samples")
    } else {
      shinyjs::disable("queue_batch_samples")
    }
  })

  run_unified_merge <- function() {

    # returns a data.frame
        # Sample_Name
        # R1_Path
        # R2_Path
        # Source
    merge_inputs <- .prepare_merge_inputs(
      mode = "single",
      queued_df = .order_sample_df(
        merge_rv$queued,
        merge_rv$sample_order,
        "Queue ID"
      )
    )

    .run_merge_pipeline(
      merge_inputs,
      settings = .get_merge_settings(input)
    )
  }

  observeEvent(input$run_single_merge, {
    run_unified_merge()
  }, ignoreInit = TRUE)

  observeEvent(input$run_batch_merge, {
    run_unified_merge()
  }, ignoreInit = TRUE)
}
