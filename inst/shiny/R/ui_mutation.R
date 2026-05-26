mutation_label <- function(label, tooltip_id) {
  div(
    class = "mutation-field-label",
    tags$label(label, class = "control-label"),
    tooltip(
      icon("circle-question"),
      tooltip_text(tooltip_id)
    )
  )
}

tooltip_text <- function(tooltip_id) {
  texts <- c(
    mutation_reference_sequence = "The complete DNA sequence of the region of the amplicon to be analyzed.",
    mutation_forward_primer = "DNA sequence of the portion of the forward primer that binds to the source DNA. Do not include any extensions or partial Illumina adapters.",
    mutation_reverse_primer = "DNA sequence of the portion of the reverse primer that binds to the source DNA. Do not include any extensions or partial Illumina adapters.",
    mutation_insert_start_coordinate = "Number of the first nt of the reference sequence within the total amplicon (start counting from 1).",
    mutation_insert_end_coordinate = "Number of the last nt of the reference sequence relative to the END of the amplicon. The last nt is 0, and count backwards from there with negative numbers.",
    mutation_amplicon_length = "The total length of the amplicon between any Illumina adapters, including extensions such as barcodes or UMIs added during PCR.",
    mutation_forward_extension_type = "Optional barcode or UMI sequence attached to the forward primer.",
    mutation_forward_extension_sequence = "Sequence of the selected forward barcode or UMI.",
    mutation_reverse_extension_type = "Optional barcode or UMI sequence attached to the reverse primer.",
    mutation_reverse_extension_sequence = "Sequence of the selected reverse barcode or UMI.",
    nuclease_forward_primer = "DNA sequence of the portion of the forward primer that binds to the source DNA.",
    nuclease_reverse_primer = "DNA sequence of the portion of the reverse primer that binds to the source DNA.",
    nuclease_grna_sequence = "DNA sequence of an S. pyogenes Cas9 gRNA, if appropriate. Does not check for PAM, will automatically find the cut site.",
    nuclease_cut_position = "Manually set the expected cut site, for instance if using a non-SpyCas9 nuclease. Will override the gRNA sequence for cut site determination.",
    nuclease_forward_extension_type = "Optional barcode or UMI sequence attached to the forward primer.",
    nuclease_forward_extension_sequence = "Sequence of the selected forward barcode or UMI.",
    nuclease_reverse_extension_type = "Optional barcode or UMI sequence attached to the reverse primer.",
    nuclease_reverse_extension_sequence = "Sequence of the selected reverse barcode or UMI."
  )

  text <- texts[[tooltip_id]]
  if (is.null(text)) paste("TODO:", tooltip_id) else text
}

mutation_text_input <- function(input_id, label, tooltip_id, placeholder = "") {
  div(
    class = "minimal-input-group mutation-setting-field",
    mutation_label(label, tooltip_id),
    textInput(input_id, label = NULL, placeholder = placeholder)
  )
}

mutation_numeric_inline_input <- function(input_id, label, tooltip_id, placeholder = NULL) {
  input <- numericInput(input_id, label = NULL, value = NA, step = 1)
  if (!is.null(placeholder)) {
    input <- htmltools::tagQuery(input)$find("input")$addAttrs(placeholder = placeholder)$allTags()
  }

  div(
    class = "minimal-input-group mutation-setting-field mutation-numeric-field mutation-numeric-inline-field",
    mutation_label(label, tooltip_id),
    input
  )
}

mutation_numeric_input_plain <- function(input_id, label) {
  div(
    class = "minimal-input-group mutation-setting-field mutation-numeric-field",
    div(
      class = "mutation-field-label mutation-field-label-plain",
      tags$label(label, class = "control-label")
    ),
    numericInput(input_id, label = NULL, value = NA, step = 1)
  )
}

extension_settings_row <- function(prefix, direction, tooltip_prefix) {
  direction_title <- paste0(tools::toTitleCase(direction), " extension type")
  direction_id <- paste0(direction, "_extension_type")
  sequence_id <- paste0(direction, "_extension_sequence")

  div(
    class = "mutation-extension-row",
    div(
      class = "mutation-setting-field mutation-extension-type",
      mutation_label(direction_title, paste0(tooltip_prefix, "_", direction_id)),
      selectInput(
        paste0(prefix, "_", direction_id),
        label = NULL,
        choices = c("None", "Barcode", "UMI"),
        selected = "None"
      )
    ),
    div(class = "mutation-extension-row spacer1"),
    conditionalPanel(
      condition = sprintf("input.%s_%s == 'Barcode' || input.%s_%s == 'UMI'", prefix, direction_id, prefix, direction_id),
      mutation_text_input(paste0(prefix, "_", sequence_id), "Sequence", paste0(tooltip_prefix, "_", sequence_id))
    ),
    conditionalPanel(
      condition = sprintf("input.%s_%s != 'Barcode' && input.%s_%s != 'UMI'", prefix, direction_id, prefix, direction_id),
      div(class = "mutation-extension-row type-empty")
    ),
    div(class = "mutation-extension-row spacer2")
  )
}

mutation_sequence_settings_ui <- function(prefix = "mutation", show_heading = TRUE, show_sublabels = FALSE, clear_button = NULL) {
  tagList(
    if (show_heading) div(class = "mutation-settings-heading-spacer"),
    div(
      class = "mutation-settings-section mutation-sequence-section",
      if (show_heading) h5("Sequence Definition", class = "mutation-section-title"),
      if (show_heading) div(class = "mutation-section-rule"),
      div(
        class = "minimal-input-group mutation-setting-field mutation-reference-field",
        mutation_label("Reference Sequence", "mutation_reference_sequence"),
        textAreaInput(
          paste0(prefix, "_reference_sequence"),
          label = NULL,
          rows = 1,
          width = "100%",
          placeholder = "Paste reference sequence"
        ),
        if (show_sublabels) div(class = "minimal-input-label-text mutation-reference-sublabel", "Region of sequence to be analyzed")
      ),
      div(class = "mutation-row-spacer"),
      div(
        class = "mutation-primer-row",
        div(
          mutation_text_input(paste0(prefix, "_forward_primer"), "Forward primer", "mutation_forward_primer"),
          if (show_sublabels) div(class = "minimal-input-label-text mutation-sequence-sublabel",
          style = "white-space: pre;",
          "DNA   (5'-          -3')")
        ),
        div(
          mutation_text_input(paste0(prefix, "_reverse_primer"), "Reverse primer", "mutation_reverse_primer"),
          if (show_sublabels) div(class = "minimal-input-label-text mutation-sequence-sublabel",
          style = "white-space: pre;",
           "DNA   (5'-          -3')")
        )
      ),
      div(class = "mutation-row-spacer"),
      div(
        class = "mutation-coordinate-row",
        mutation_numeric_inline_input(paste0(prefix, "_insert_start_coordinate"), "Insert start", "mutation_insert_start_coordinate", placeholder = "+"),
        mutation_numeric_inline_input(paste0(prefix, "_insert_end_coordinate"), "Insert end", "mutation_insert_end_coordinate", placeholder = "-"),
        mutation_numeric_inline_input(paste0(prefix, "_amplicon_length"), "Amplicon length", "mutation_amplicon_length")
      ),
      if (!is.null(clear_button)) div(class = "clear-action-box-footer", clear_button)
    )
  )
}

mutation_analysis_settings_ui <- function(prefix = "mutation", show_heading = TRUE) {
  tagList(
    if (show_heading) div(class = "mutation-settings-heading-spacer"),
    div(
      class = "mutation-settings-section mutation-analysis-section",
      if (show_heading) h5("Analysis Parameters", class = "mutation-section-title"),
      if (show_heading) div(class = "mutation-section-rule"),
      extension_settings_row(prefix, "forward", "mutation"),
      br(),
      extension_settings_row(prefix, "reverse", "mutation"),
      div(class = "mutation-row-spacer"),
      div(
        class = "mutation-antibody-inline-row",
        div(
          class = "mutation-antibody-row",
          checkboxInput(paste0(prefix, "_antibody"), label = NULL, value = FALSE),
          tags$span("Antibody", class = "mutation-checkbox-label"),
          tags$span(
            class = "mutation-antibody-tooltip",
            tooltip(icon("circle-question"), "Is the sequence an antibody? If yes, enter coordinates of the first nt of each region in the boxes that appear below.")
          )
        ),
        div(
          class = "mutation-antibody-coordinate-slot",
          conditionalPanel(
            condition = sprintf("input.%s_antibody === true", prefix),
            div(
              class = "mutation-antibody-coordinate-row",
              mutation_numeric_input_plain(paste0(prefix, "_fr1_start"), "FR1"),
              mutation_numeric_input_plain(paste0(prefix, "_cdr1_start"), "CDR1"),
              mutation_numeric_input_plain(paste0(prefix, "_fr2_start"), "FR2"),
              mutation_numeric_input_plain(paste0(prefix, "_cdr2_start"), "CDR2"),
              mutation_numeric_input_plain(paste0(prefix, "_fr3_start"), "FR3"),
              mutation_numeric_input_plain(paste0(prefix, "_cdr4_start"), "CDR3"),
              mutation_numeric_input_plain(paste0(prefix, "_fr4_start"), "FR4"),
              mutation_numeric_input_plain(paste0(prefix, "_fr4_end"), "End")
            )
          )
        )
      )
    )
  )
}

mutation_settings_ui <- function(prefix = "mutation", show_headings = TRUE) {
  tagList(
    mutation_sequence_settings_ui(prefix, show_heading = show_headings),
    mutation_analysis_settings_ui(prefix, show_heading = show_headings)
  )
}

mutation_page_ui <- function() {
  nav_panel(
    "Measure Mutations",
    value = "shm.app",
    div(
      class = "module-header",
      h1("Measure Mutations"),
      p("Quantify mutations in antibody and non-antibody sequences.")
    ),
    fluidRow(
      div(
        class = "upload-card-shell mutation-card-shell nuclease-card-shell",
        card(
          card_header("Upload Files", class = "card-header-custom"),
          card_body(
            div(
              class = "workflow-tabs",
              navset_tab(
                id = "mutation_upload_mode",
                selected = "mutation_single_sample",
                nav_panel(
                  title = "Single Sample",
                  value = "mutation_single_sample",
                  div(
                    class = "card-text-body",
                    div(
                      class = "mutation-upload-row upload-input-row-single",
                      div(
                        class = "sample-name-field",
                        div(
                          class = "minimal-input-group",
                          textInput("mutation_sample_name", label = NULL, placeholder = "Enter sample name"),
                          div(class = "minimal-input-label-text", "Sample identifier")
                        )
                      ),
                      div(
                        class = "upload-zone-col",
                        div(
                          class = "minimal-file-group drop-upload-group",
                          drop_upload_file_input(
                            "Merged file",
                            fileInput(
                            "mutation_merged_fastq",
                            label = NULL,
                            buttonLabel = "Browse",
                            placeholder = "No file selected",
                            accept = c(".fastq", ".fastq.gz", ".fq", ".fq.gz", ".gz"),
                            width = "100%"
                            )
                          ),
                          div(class = "minimal-input-label-fileInput", "Merged FASTQ file")
                        )
                      ),
                      div(
                        class = "upload-row-action mutation-single-add-action",
                        actionButton("mutation_queue_sample", "Add Sample", class = "secondary-action-btn")
                      )
                    ),
                    div(class = "mutation-input-settings-divider"),
                    mutation_sequence_settings_ui(
                      "mutation",
                      show_heading = FALSE,
                      show_sublabels = TRUE,
                      clear_button = clear_action_button("mutation_clear_settings", "Clear Settings")
                    ),
                    div(class = "input-action-spacer mutation-sequence-action-spacer"),
                    div(
                      class = "upload-action-row message-action-row mutation-action-row",
                      div(class = "message-slot", uiOutput("mutation_message")),
                      div(
                        class = "action-button-group",
                        actionButton("run_single_mutation", "Start Analysis", class = "primary-action-btn")
                      )
                    )
                  )
                ),
                nav_panel(
                  title = "Batch Import",
                  value = "mutation_batch_import",
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
                            "mutation_batch_manifest",
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
                              href = "samples_to_measure_mutations.csv",
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
                            "mutation_batch_archive",
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
                        clear_action_button("mutation_clear_batch", "Clear Files")
                      )
                    ),
                    div(class = "input-action-spacer"),
                    div(
                      class = "upload-action-row message-action-row",
                      div(class = "message-slot", uiOutput("mutation_batch_message")),
                      div(
                        class = "action-button-group",
                        actionButton("mutation_queue_batch_samples", "Add Samples", class = "secondary-action-btn"),
                        actionButton("run_batch_mutation", "Start Analysis", class = "primary-action-btn")
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
          DT::DTOutput("mutation_queued_samples")
        )
      )
    )
  )
}
