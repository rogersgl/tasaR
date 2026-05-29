.js_string <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("'", "\\\\'", x)
  x
}

.reconcile_sample_order <- function(ids, current_order = character()) {
  ids <- as.character(ids)
  current_order <- as.character(current_order)

  current_order <- current_order[current_order %in% ids]
  current_order <- current_order[!duplicated(current_order)]
  new_ids <- ids[!ids %in% current_order]

  c(current_order, new_ids)
}

.order_sample_df <- function(df, order, id_col) {
  if (is.null(df) || nrow(df) == 0 || !id_col %in% names(df)) {
    return(df)
  }

  order <- .reconcile_sample_order(df[[id_col]], order)
  df[match(order, df[[id_col]]), , drop = FALSE]
}

.move_order_item <- function(order, id, direction = c("up", "down")) {
  direction <- match.arg(direction)
  order <- as.character(order)
  id <- as.character(id)

  idx <- match(id, order)
  if (is.na(idx)) return(order)

  swap_idx <- if (direction == "up") idx - 1L else idx + 1L
  if (swap_idx < 1L || swap_idx > length(order)) return(order)

  order[c(idx, swap_idx)] <- order[c(swap_idx, idx)]
  order
}

.remove_batch_row_by_id <- function(df, row_id, id_col = "Batch Queue ID") {
  if (is.null(df) || !id_col %in% names(df)) {
    return(list(df = df, removed = NULL, success = FALSE))
  }

  idx <- match(row_id, df[[id_col]])

  if (is.na(idx)) {
    return(list(df = df, removed = NULL, success = FALSE))
  }

  removed <- df[idx, , drop = FALSE]
  updated <- df[-idx, , drop = FALSE]

  list(df = updated, removed = removed, success = TRUE)
}

.sample_table_button <- function(row_id, input_id, label, class, icon) {
  sprintf(
    "<button type='button' class='%s' title='%s' aria-label='%s' onclick=\"Shiny.setInputValue('%s', '%s', {priority: 'event'})\">%s</button>",
    class,
    label,
    label,
    .js_string(input_id),
    .js_string(row_id),
    icon
  )
}

.queue_actions_html <- function(row_id, remove_input_id, move_up_input_id = NULL, move_down_input_id = NULL) {
  paste0(
    if (!is.null(move_up_input_id)) {
      .sample_table_button(row_id, move_up_input_id, "Move up", "queue-icon-btn queue-move-btn", "&#x25B2;")
    } else {
      ""
    },
    if (!is.null(move_down_input_id)) {
      .sample_table_button(row_id, move_down_input_id, "Move down", "queue-icon-btn queue-move-btn", "&#x25BC;")
    } else {
      ""
    },
    "  ",
    .sample_table_button(row_id, remove_input_id, "Remove", "queue-icon-btn queue-remove-btn", "&#x2715;")
  )
}

.queue_drag_handle_html <- function() {
  "<span class='queue-drag-handle' title='Drag to reorder'>::</span>"
}

.queue_small_text_column_defs <- function(table_df, small_text_columns) {
  if (is.null(table_df) || length(small_text_columns) == 0) {
    return(list())
  }

  targets <- match(small_text_columns, names(table_df)) - 1L
  targets <- targets[!is.na(targets)]

  if (length(targets) == 0) {
    return(list())
  }

  list(list(className = "queue-small-text", targets = targets))
}

.sample_status_html <- function(status) {
  status_class <- ifelse(
    identical(status, "Ready"),
    "queue-status-ready",
    "queue-status-issue"
  )

  sprintf(
    "<div class='queue-status-text %s'>%s</div>",
    status_class,
    status
  )
}

.sample_progress_html <- function(progress, download_href = NULL, download_filename = NULL) {
  progress <- trimws(as.character(progress))
  if (is.null(download_href)) download_href <- rep("", length(progress))
  if (is.null(download_filename)) download_filename <- rep("", length(progress))
  download_href <- rep_len(as.character(download_href), length(progress))
  download_filename <- rep_len(as.character(download_filename), length(progress))
  download_icon <- as.character(shiny::icon("download"))

  vapply(seq_along(progress), function(i) {
    x <- progress[[i]]
    x_norm <- tolower(x)

    if (startsWith(x, "step|")) {
      parts <- strsplit(x, "|", fixed = TRUE)[[1]]
      step <- if (length(parts) >= 2L) suppressWarnings(as.integer(parts[[2]])) else NA_integer_
      total <- if (length(parts) >= 3L) suppressWarnings(as.integer(parts[[3]])) else NA_integer_
      label <- if (length(parts) >= 4L) paste(parts[4:length(parts)], collapse = "|") else "Processing"
      if (is.na(step) || is.na(total) || total <= 0L) {
        step <- 0L
        total <- 1L
      }
      pct <- max(0, min(100, round(100 * step / total)))
      title <- sprintf("%s (%s/%s)", label, step, total)

      return(sprintf(
        paste0(
          "<span class='queue-progress-step' title='%s' aria-label='%s'>",
          "<span class='queue-progress-mini-spinner'></span>",
          "<span class='queue-progress-step-text'>%s/%s</span>",
          "<span class='queue-progress-step-bar'><span style='width:%s%%'></span></span>",
          "</span>"
        ),
        htmltools::htmlEscape(title, attribute = TRUE),
        htmltools::htmlEscape(title, attribute = TRUE),
        htmltools::htmlEscape(step),
        htmltools::htmlEscape(total),
        pct
      ))
    }

    if (startsWith(x, "exporting|")) {
      parts <- strsplit(x, "|", fixed = TRUE)[[1]]
      label <- if (length(parts) >= 2L) paste(parts[2:length(parts)], collapse = "|") else "Exporting"
      return(sprintf(
        paste0(
          "<span class='queue-progress-exporting' title='%s' aria-label='%s'>",
          "<span class='queue-progress-mini-spinner'></span>",
          "<span>Exporting</span>",
          "</span>"
        ),
        htmltools::htmlEscape(label, attribute = TRUE),
        htmltools::htmlEscape(label, attribute = TRUE)
      ))
    }

    if (startsWith(x, "finalizing|")) {
      parts <- strsplit(x, "|", fixed = TRUE)[[1]]
      label <- if (length(parts) >= 2L) paste(parts[2:length(parts)], collapse = "|") else "Finalizing"
      return(sprintf(
        paste0(
          "<span class='queue-progress-finalizing' title='%s' aria-label='%s'>",
          "<span class='queue-progress-mini-spinner'></span>",
          "<span>Finalizing</span>",
          "</span>"
        ),
        htmltools::htmlEscape(label, attribute = TRUE),
        htmltools::htmlEscape(label, attribute = TRUE)
      ))
    }

    if (startsWith(x, "analyzed|")) {
      parts <- strsplit(x, "|", fixed = TRUE)[[1]]
      label <- if (length(parts) >= 2L) paste(parts[2:length(parts)], collapse = "|") else "Analysis complete"
      return(sprintf(
        paste0(
          "<span class='queue-progress-analyzed' title='%s' aria-label='%s'>",
          "<span class='queue-progress-analyzed-icon'>&#10003;</span>",
          "<span>Analyzed</span>",
          "</span>"
        ),
        htmltools::htmlEscape(label, attribute = TRUE),
        htmltools::htmlEscape(label, attribute = TRUE)
      ))
    }

    if (x_norm %in% c("processing", "running", "in_progress")) {
      return("<span class='queue-progress-spinner' title='Processing' aria-label='Processing'></span>")
    }

    if (x_norm %in% c("done", "complete", "completed", "success", "true")) {
      if (nzchar(download_href[[i]])) {
        return(sprintf(
          paste0(
            "<span class='queue-progress-complete'>",
            "<span class='queue-progress-done' title='Complete' aria-label='Complete'>&#10003;</span>",
            "<a class='queue-progress-download' href='%s' download='%s' title='Download merged file' aria-label='Download merged file'>%s</a>",
            "</span>"
          ),
          htmltools::htmlEscape(download_href[[i]], attribute = TRUE),
          htmltools::htmlEscape(download_filename[[i]], attribute = TRUE),
          download_icon
        ))
      }

      return("<span class='queue-progress-done' title='Complete' aria-label='Complete'>&#10003;</span>")
    }

    "<span class='queue-progress-empty' aria-label='Not started'></span>"
  }, character(1))
}

.make_sample_table_data <- function(
    df,
    id_col,
    remove_input_id,
    move_up_input_id = NULL,
    move_down_input_id = NULL
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

  data.frame(
    "::" = .queue_drag_handle_html(),
    "Display Order" = seq_len(nrow(df)),
    "Sample Name" = df[["Sample Name"]],
    "Left (R1) File" = df[["Left (R1) File"]],
    "Right (R2) File" = df[["Right (R2) File"]],
    "Status" = vapply(df[["Status"]], .sample_status_html, character(1)),
    "Actions" = actions,
    "Progress" = .sample_progress_html(progress, download_href, download_filename),
    "Sample ID" = df[[id_col]],
    check.names = FALSE
  )
}

.setup_reorderable_sample_table <- function(
    input,
    output,
    output_id,
    data_fn,
    id_col,
    get_order,
    set_order,
    remove_input_id,
    move_up_input_id,
    move_down_input_id,
    table_builder = .make_sample_table_data,
    small_text_columns = character()
) {
  output[[output_id]] <- DT::renderDT({
    df <- data_fn()

    if (is.null(df) || nrow(df) == 0) {
      return(DT::datatable(
        data.frame(),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    if ("Note" %in% names(df)) {
      return(DT::datatable(
        data.frame(Note = df[["Note"]][1], check.names = FALSE),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    ordered_df <- .order_sample_df(df, get_order(), id_col)
    table_df <- table_builder(
      ordered_df,
      id_col = id_col,
      remove_input_id = remove_input_id,
      move_up_input_id = paste0(output_id, "_move_row_up"),
      move_down_input_id = paste0(output_id, "_move_row_down")
    )
    id_col_index <- match("Sample ID", names(table_df)) - 1L
    order_col_index <- match("Display Order", names(table_df)) - 1L
    progress_col_index <- match("Progress", names(table_df)) - 1L
    progress_col_def <- if (!is.na(progress_col_index)) {
      list(list(className = "queue-progress-cell", targets = progress_col_index))
    } else {
      list()
    }

    column_defs <- c(
      .queue_small_text_column_defs(table_df, small_text_columns),
      progress_col_def,
      list(
        list(visible = FALSE, targets = id_col_index),
        list(visible = FALSE, targets = order_col_index),
        list(orderable = FALSE, targets = "_all")
      )
    )

    DT::datatable(
      table_df,
      rownames = FALSE,
      escape = FALSE,
      extensions = "RowReorder",
      selection = "single",
      class = "queue-table compact",
      options = list(
        dom = "t",
        paging = FALSE,
        searching = FALSE,
        ordering = TRUE,
        order = list(list(order_col_index, "asc")),
        info = FALSE,
        rowReorder = list(
          dataSrc = order_col_index,
          selector = ".queue-drag-handle"
        ),
        columnDefs = column_defs
      ),
      callback = DT::JS(sprintf(
        "table.on('row-reorder', function(e, diff, edit) {
          table.one('draw', function() {
            var ids = table.rows({ order: 'current' }).data().toArray().map(function(row) {
              return row[%d];
            });
            Shiny.setInputValue('%s_row_order', ids, { priority: 'event' });
          });
        });",
        id_col_index,
        output_id
      ))
    )
  }, server = FALSE)

  observeEvent(input[[paste0(output_id, "_row_order")]], {
    df <- data_fn()
    if (is.null(df) || nrow(df) == 0 || !id_col %in% names(df)) return()

    set_order(.reconcile_sample_order(
      df[[id_col]],
      input[[paste0(output_id, "_row_order")]]
    ))
  }, ignoreInit = TRUE)

  observeEvent(input[[move_up_input_id]], {
    df <- .order_sample_df(data_fn(), get_order(), id_col)
    selected <- input[[paste0(output_id, "_rows_selected")]]
    if (is.null(df) || length(selected) != 1 || selected < 1 || selected > nrow(df)) return()

    set_order(.move_order_item(get_order(), df[selected, id_col, drop = TRUE], "up"))
  }, ignoreInit = TRUE)

  observeEvent(input[[move_down_input_id]], {
    df <- .order_sample_df(data_fn(), get_order(), id_col)
    selected <- input[[paste0(output_id, "_rows_selected")]]
    if (is.null(df) || length(selected) != 1 || selected < 1 || selected > nrow(df)) return()

    set_order(.move_order_item(get_order(), df[selected, id_col, drop = TRUE], "down"))
  }, ignoreInit = TRUE)

  observeEvent(input[[paste0(output_id, "_move_row_up")]], {
    df <- data_fn()
    if (is.null(df) || nrow(df) == 0 || !id_col %in% names(df)) return()

    row_id <- input[[paste0(output_id, "_move_row_up")]]
    set_order(.move_order_item(
      .reconcile_sample_order(df[[id_col]], get_order()),
      row_id,
      "up"
    ))
  }, ignoreInit = TRUE)

  observeEvent(input[[paste0(output_id, "_move_row_down")]], {
    df <- data_fn()
    if (is.null(df) || nrow(df) == 0 || !id_col %in% names(df)) return()

    row_id <- input[[paste0(output_id, "_move_row_down")]]
    set_order(.move_order_item(
      .reconcile_sample_order(df[[id_col]], get_order()),
      row_id,
      "down"
    ))
  }, ignoreInit = TRUE)
}
