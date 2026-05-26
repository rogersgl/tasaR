.queue_upload_message <- function(message, message_type) {
  if (is.null(message) || is.null(message_type)) {
    return(div(class = "queue-alert message-placeholder", "Message placeholder"))
  }

  div(class = paste("queue-alert", paste0("queue-alert-", message_type)), message)
}

.clear_queue_message_later <- function(rv, delay = 4) {
  later::later(function() {
    rv$message <- NULL
    rv$message_type <- NULL
  }, delay)
}
