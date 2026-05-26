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

.sample_progress_html <- function(progress) {
  progress <- tolower(trimws(as.character(progress)))

  vapply(progress, function(x) {
    if (x %in% c("processing", "running", "in_progress")) {
      return("<span class='queue-progress-spinner' title='Processing' aria-label='Processing'></span>")
    }

    if (x %in% c("done", "complete", "completed", "success", "true")) {
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
  actions <- vapply(df[[id_col]], function(row_id) {
    .queue_actions_html(row_id, remove_input_id, move_up_input_id, move_down_input_id)
  }, character(1))

  data.frame(
    "::" = .queue_drag_handle_html(),
    "Display Order" = seq_len(nrow(df)),
    "Sample Name" = df[["Sample Name"]],
    "Left (R1) File" = df[["Left (R1) File"]],
    "Right (R2) File" = df[["Right (R2) File"]],
    "Source" = if ("Source" %in% names(df)) df[["Source"]] else rep("single", nrow(df)),
    "Status" = vapply(df[["Status"]], .sample_status_html, character(1)),
    "Actions" = actions,
    "Progress" = .sample_progress_html(progress),
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
    table_builder = .make_sample_table_data
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
        columnDefs = list(
          list(visible = FALSE, targets = id_col_index),
          list(visible = FALSE, targets = order_col_index),
          list(orderable = FALSE, targets = "_all")
        )
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
