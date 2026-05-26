clear_action_button <- function(input_id, label) {
  actionButton(
    input_id,
    tagList(
      tags$span(class = "clear-action-icon", HTML("&#x21BA;")),
      tags$span(label)
    ),
    class = "clear-action-link"
  )
}

drop_upload_file_input <- function(box_label, file_input) {
  htmltools::tagQuery(file_input)$
    find(".input-group")$
    addClass("labeled-drop-upload-input")$
    prepend(div(class = "drop-upload-box-label", box_label))$
    allTags()
}
