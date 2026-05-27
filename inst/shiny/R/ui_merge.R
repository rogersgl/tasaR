merge_page_ui <- function() {
  nav_panel(
    "Merge Reads",
    value = "merge.app",
    div(
      class = "module-header",
      h1("Merge Reads"),
      p("Merge paired-end FASTQ reads with a visual interface.")
    ),
    fluidRow(
      div(
        class = "upload-card-shell",
        card(
          card_header("Upload Files", class = "card-header-custom"),
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
                      class = "upload-input-row upload-input-row-single",
                      div(
                        class = "sample-name-field",
                        div(
                          class = "minimal-input-group",
                          textInput(inputId = "sample_name", label = NULL, placeholder = "Enter sample name"),
                          div(class = "minimal-input-label-text", "Sample identifier")
                        )
                      ),
                      div(
                        class = "upload-zone-col",
                        div(
                        class = "minimal-file-group drop-upload-group",
                          drop_upload_file_input(
                            "R1 reads",
                            fileInput(inputId = "fastq_r1", label = NULL, buttonLabel = "Browse", placeholder = "No file selected")
                          ),
                          div(class = "minimal-input-label-fileInput", "Left (R1) FASTQ file")
                        )
                      ),
                      div(
                        class = "upload-zone-col",
                        div(
                        class = "minimal-file-group drop-upload-group",
                          drop_upload_file_input(
                            "R2 reads",
                            fileInput(inputId = "fastq_r2", label = NULL, buttonLabel = "Browse", placeholder = "No file selected")
                          ),
                          div(class = "minimal-input-label-fileInput", "Right (R2) FASTQ file")
                        )
                      )
                    ),
                    div(class = "input-action-spacer"),
                    div(
                      class = "upload-action-row message-action-row merge-single-action-row",
                      div(class = "message-slot", uiOutput("merge_message")),
                      div(
                        class = "action-button-group",
                        actionButton("queue_sample", "Add Sample", class = "secondary-action-btn"),
                        actionButton("run_single_merge", "Start Merge", class = "primary-action-btn")
                      )
                    )
                  )
                ),
                nav_panel(
                  title = "Batch Import",
                  value = "batch_import",
                  div(
                    class = "card-text-body",
                    div(
                      class = "upload-input-row upload-input-row-batch",
                      div(
                        class = "upload-zone-col",
                        div(
                        class = "minimal-file-group drop-upload-group",
                          drop_upload_file_input(
                            "Manifest (.csv)",
                            fileInput(
                            inputId = "batch_manifest",
                            label = NULL,
                            buttonLabel = "Browse",
                            placeholder = "No file selected",
                            accept = c(".csv"),
                            width = "100%"
                            )
                          ),
                          div(
                            class = "minimal-input-label-row",
                            div(class = "minimal-input-label-fileInput", "List of samples"),
                            tags$a(
                              href = "samples_to_merge.csv",
                              download = NA,
                              class = "template-download-link",
                              "Template"
                            )
                          )
                        )
                      ),
                      div(
                        class = "upload-zone-col",
                        div(
                        class = "minimal-file-group drop-upload-group",
                          drop_upload_file_input(
                            "Archive (.zip)",
                            fileInput(
                            inputId = "batch_archive",
                            label = NULL,
                            buttonLabel = "Browse",
                            placeholder = "No file selected",
                            accept = c(".zip"),
                            width = "100%"
                            )
                          ),
                          div(class = "minimal-input-label-fileInput", "FASTQ archive (.zip)")
                        )
                      ),
                      div(
                        class = "upload-row-action",
                        clear_action_button("clear_batch", "Clear Files")
                      )
                    ),
                    div(class = "input-action-spacer"),
                    div(
                      class = "upload-action-row message-action-row",
                      div(class = "message-slot", uiOutput("batch_message")),
                      div(
                        class = "action-button-group",
                        actionButton("queue_batch_samples", "Add Samples", class = "secondary-action-btn"),
                        actionButton("run_batch_merge", "Start Merge", class = "primary-action-btn")
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    ),
    fluidRow(
      div(
        class = "upload-table-shell",
        div(
          class = "queue-table-wrapper queue-table-section",
          DT::DTOutput("merge_queued_samples")
        ),
        uiOutput("merge_downloads")
      )
    )
  )
}
