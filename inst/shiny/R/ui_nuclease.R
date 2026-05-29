nuclease_parameters_ui <- function(prefix = "nuclease", show_heading = TRUE, sidebar_breaks = FALSE) {
  tagList(
    if (show_heading) div(class = "mutation-settings-heading-spacer"),
    div(
      class = "mutation-settings-section mutation-analysis-section",
      mutation_text_input(paste0(prefix, "_forward_primer"), "Forward primer", "nuclease_forward_primer"),
      if (sidebar_breaks) br(),
      mutation_text_input(paste0(prefix, "_reverse_primer"), "Reverse primer", "nuclease_reverse_primer"),
      if (sidebar_breaks) br(),
      div(class = "mutation-row-spacer"),
      extension_settings_row(prefix, "forward", "nuclease"),
      br(),
      extension_settings_row(prefix, "reverse", "nuclease")
    )
  )
}

nuclease_sample_settings_ui <- function(prefix = "nuclease_modal") {
  tagList(
    div(
      class = "mutation-settings-section mutation-sequence-section",
      div(
        class = "mutation-primer-row",
        div(
          mutation_text_input(paste0(prefix, "_grna_sequence"), "gRNA sequence", "nuclease_grna_sequence", placeholder = "Enter DNA sequence"),
          div(
            class = "minimal-input-label-text mutation-sequence-sublabel",
            style = "white-space: pre;",
            "DNA   (5'-          -3')"
          )
        ),
        mutation_numeric_inline_input(paste0(prefix, "_cut_position"), "Cut Position", "nuclease_cut_position", placeholder = "auto")
      )
    ),
    nuclease_parameters_ui(prefix, show_heading = TRUE)
  )
}

nuclease_page_ui <- function() {
  nav_panel(
    "Nuclease Analysis",
    value = "cut.app",
    div(
      class = "module-header",
      h1("Nuclease Analysis"),
      p("Quantify mutations with paired control and experimental samples.")
    ),
    fluidRow(
      div(
        class = "upload-card-shell mutation-card-shell",
        card(
          card_header("Upload Files", class = "card-header-custom"),
          card_body(
            div(
              class = "workflow-tabs",
              navset_tab(
                id = "nuclease_upload_mode",
                selected = "nuclease_single_sample",
                nav_panel(
                  title = "Single Sample",
                  value = "nuclease_single_sample",
                  div(
                    class = "card-text-body",
                    div(
                      class = "mutation-upload-row upload-input-row-single",
                      div(
                        class = "sample-name-field",
                        div(
                          class = "minimal-input-group",
                          textInput("nuclease_sample_name", label = NULL, placeholder = "Enter sample name"),
                          div(class = "minimal-input-label-text", "Sample identifier")
                        )
                      ),
                      div(
                        class = "upload-zone-col",
                        div(
                          class = "minimal-file-group drop-upload-group",
                          drop_upload_file_input(
                            "Control file",
                            fileInput(
                            "nuclease_control_fastq",
                            label = NULL,
                            buttonLabel = "Browse",
                            placeholder = "No file selected",
                            accept = c(".fastq", ".fastq.gz", ".fq", ".fq.gz", ".gz"),
                            width = "100%"
                            )
                          ),
                          div(class = "minimal-input-label-fileInput", "Control FASTQ file")
                        )
                      ),
                      div(
                        class = "upload-zone-col",
                        div(
                          class = "minimal-file-group drop-upload-group",
                          drop_upload_file_input(
                            "Experimental file",
                            fileInput(
                            "nuclease_experimental_fastq",
                            label = NULL,
                            buttonLabel = "Browse",
                            placeholder = "No file selected",
                            accept = c(".fastq", ".fastq.gz", ".fq", ".fq.gz", ".gz"),
                            width = "100%"
                            )
                          ),
                          div(class = "minimal-input-label-fileInput", "Experimental FASTQ file")
                        )
                      )
                    ),
                    div(class = "mutation-input-settings-divider"),
                    div(
                      class = "nuclease-target-action-row",
                      div(
                        class = "mutation-settings-section mutation-sequence-section nuclease-target-box",
                        div(
                          class = "nuclease-target-fields",
                          div(
                          class = "nuclease-grna-field",
                            mutation_text_input("nuclease_grna_sequence", "gRNA sequence", "nuclease_grna_sequence", placeholder = "Enter DNA sequence"),
                            div(
                              class = "minimal-input-label-text mutation-sequence-sublabel",
                              style = "white-space: pre;",
                              "DNA   (5'-          -3')"
                            )
                          ),
                          mutation_numeric_inline_input("nuclease_cut_position", "Cut Position", "nuclease_cut_position", placeholder = "auto")
                        ),
                        div(
                          class = "clear-action-box-footer",
                          clear_action_button("nuclease_clear_settings", "Clear Settings")
                        )
                      ),
                      div(
                        class = "upload-row-action mutation-single-add-action nuclease-target-action",
                        actionButton("nuclease_queue_sample", "Add Sample", class = "secondary-action-btn")
                      )
                    ),
                    div(class = "input-action-spacer mutation-sequence-action-spacer"),
                    div(
                      class = "upload-action-row message-action-row mutation-action-row",
                      div(class = "message-slot", uiOutput("nuclease_message")),
                      div(
                        class = "action-button-group",
                        actionButton("run_single_nuclease", "Start Analysis", class = "primary-action-btn")
                      )
                    )
                  )
                ),
                nav_panel(
                  title = "Batch Import",
                  value = "nuclease_batch_import",
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
                            "nuclease_batch_manifest",
                            label = NULL,
                            buttonLabel = "Browse",
                            placeholder = "No file selected",
                            accept = c(".csv"),
                            width = "100%"
                            )
                          ),
                          div(
                            class = "minimal-input-label-row",
                            div(class = "minimal-input-label-fileInput", "Samples and settings"),
                            tags$a(
                              href = "samples_to_nuclease_analysis.csv",
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
                            "nuclease_batch_archive",
                            label = NULL,
                            buttonLabel = "Browse",
                            placeholder = "No file selected",
                            accept = c(".zip"),
                            width = "100%"
                            )
                          ),
                          div(class = "minimal-input-label-fileInput", "Merged FASTQ files")
                        )
                      ),
                      div(
                        class = "upload-row-action",
                        clear_action_button("nuclease_clear_batch", "Clear Files")
                      )
                    ),
                    div(class = "input-action-spacer"),
                    div(
                      class = "upload-action-row message-action-row",
                      div(class = "message-slot", uiOutput("nuclease_batch_message")),
                      div(
                        class = "action-button-group",
                        actionButton("nuclease_queue_batch_samples", "Add Samples", class = "secondary-action-btn"),
                        actionButton("run_batch_nuclease", "Start Analysis", class = "primary-action-btn")
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
          DT::DTOutput("nuclease_queued_samples")
        )
      )
    )
  )
}
NA
NA
