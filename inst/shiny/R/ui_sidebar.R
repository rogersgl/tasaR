setup_sidebar_output <- function(input, output) {
  output$sidebar_ui <- renderUI({
    sidebar_ui_for_nav(input$main_nav)
  })
}

sidebar_ui_for_nav <- function(nav) {
  switch(
    nav,

          "splash.app" = div(
            class = "sidebar-checklist",

            h5("Start here"),
            tags$p("Use the top navigation to choose a workflow."),
            tags$hr(class = "compact-hr"),

            tags$ul(
              tags$li(tags$strong("Review"), " the overview section"),
              tags$li(tags$strong("Choose"), " the module you need"),
              tags$li(tags$strong("Load"), " your input data")
            ),

            tags$hr(class = "compact-hr"),
            tags$p(
              class = "text-muted",
              "Supports targeted amplicon sequencing and genome editing analysis."
            )
          ),

          "merge.app" = div(
            class = "sidebar-controls",

            tags$hr(class = "section-hr"),

            h5("Merge Settings"),

            verbatimTextOutput("merge_settings"),

            tags$hr(class = "section-hr"),

            h5("Merge Controls"),

            div(style = "height: 8px;"),

            div(
              class = "d-flex align-items-center justify-content-between",
              tags$label("Output prefix", class = "control-label"),
              tooltip(
                icon("circle-question"),
                "Prefix to append to merged file names."
              )
            ),

            textInput(
              "merge_prefix",
              label = NULL,
              value = "merge"
            ),

           br(),

           div(
             class = "d-flex align-items-center justify-content-between",
             tags$label("Merge algorithm", class = "control-label"),
             tooltip(
               icon("circle-question"),
               "Select the algorithm used to merge paired-end reads. Currently only supports PANDAseq."
             )
           ),

            selectInput(
              "merge_engine",
              label = NULL,
              choices = c(
                "PANDAseq"
              ),
              selected = "PANDAseq"
            ),

            br(),

           div(
             class = "d-flex align-items-center justify-content-between",
             tags$label("Minimum overlap (bp)", class = "control-label"),
             tooltip(
               icon("circle-question"),
               "Minimum number of base pairs that must overlap in the middle for reads to be merged."
             )
           ),

            div(
              class = "slider-input-row",
              div(
                class = "slider-col",
                sliderInput(
                  "min_overlap",
                  label = NULL,
                  min = 5,
                  max = 100,
                  value = 20,
                  step = 5,
                  ticks = FALSE
                ),
                div(
                  class = "custom-slider-ticks",
                  span(style = "left: calc(10px + ((10 - 5) * (100% - 20px) / 95));",  "10"),
                  span(style = "left: calc(10px + ((20 - 5) * (100% - 20px) / 95));",  "20"),
                  span(style = "left: calc(10px + ((30 - 5) * (100% - 20px) / 95));",  "30"),
                  span(style = "left: calc(10px + ((40 - 5) * (100% - 20px) / 95));",  "40"),
                  span(style = "left: calc(10px + ((50 - 5) * (100% - 20px) / 95));",  "50"),
                  span(style = "left: calc(10px + ((60 - 5) * (100% - 20px) / 95));",  "60"),
                  span(style = "left: calc(10px + ((70 - 5) * (100% - 20px) / 95));",  "70"),
                  span(style = "left: calc(10px + ((80 - 5) * (100% - 20px) / 95));",  "80"),
                  span(style = "left: calc(10px + ((90 - 5) * (100% - 20px) / 95));",  "90"),
                  span(style = "left: calc(10px + ((100 - 5) * (100% - 20px) / 95));", "100")
                )
              ),
              div(
                class = "numeric-col",
                numericInput(inputId = "min_overlap_num",
                             label = "",
                             min = 5,
                             max = 100,
                             value = 20)
              )
            ),

            br(),

           div(
             class = "d-flex align-items-center justify-content-between",
             tags$label("Minimum merged length (bp)", class = "control-label"),
             tooltip(
               icon("circle-question"),
               "Minimum total sequence length for a merged read to be valid. Suggested default amplicon size -25."
             )
           ),

            div(
              class = "slider-input-row",
              div(
                class = "slider-col",
                sliderInput(
                  "min_length",
                  label = NULL,
                  min = 50,
                  max = 500,
                  value = 250,
                  step = 25,
                  ticks = FALSE
                ),
                div(
                  class = "custom-slider-ticks",
                  span(style = "left: calc(10px + ((100 - 50) * (100% - 20px) / 450));", "100"),
                  span(style = "left: calc(10px + ((200 - 50) * (100% - 20px) / 450));", "200"),
                  span(style = "left: calc(10px + ((300 - 50) * (100% - 20px) / 450));", "300"),
                  span(style = "left: calc(10px + ((400 - 50) * (100% - 20px) / 450));", "400"),
                  span(style = "left: calc(10px + ((500 - 50) * (100% - 20px) / 450));", "500")
                )
              ),
              div(
                class = "numeric-col",
                numericInput(inputId = "min_length_num",
                             label = "",
                             min = 50,
                             max = 500,
                             value = 250)
              )
            ),

            br(),

           div(
             class = "d-flex align-items-center justify-content-between",
             tags$label("Maximum merged length (bp)", class = "control-label"),
             tooltip(
               icon("circle-question"),
               "Max total sequence length for a merged read to be valid. Suggested default amplicon size +25."
             )
           ),

            div(
              class = "slider-input-row",
              div(
                class = "slider-col",
                sliderInput(
                  "max_length",
                  label = NULL,
                  min = 50,
                  max = 500,
                  value = 475,
                  step = 25,
                  ticks = FALSE
                ),
                div(
                  class = "custom-slider-ticks",
                  span(style = "left: calc(10px + ((100 - 50) * (100% - 20px) / 450));", "100"),
                  span(style = "left: calc(10px + ((200 - 50) * (100% - 20px) / 450));", "200"),
                  span(style = "left: calc(10px + ((300 - 50) * (100% - 20px) / 450));", "300"),
                  span(style = "left: calc(10px + ((400 - 50) * (100% - 20px) / 450));", "400"),
                  span(style = "left: calc(10px + ((500 - 50) * (100% - 20px) / 450));", "500")
                )
              ),
              div(
                class = "numeric-col",
                numericInput(inputId = "max_length_num",
                             label = "",
                             min = 50,
                             max = 500,
                             value = 475)
              )
            ),

            br(),

           div(
             class = "d-flex align-items-center justify-content-between",
             tags$label("Overlap threshold", class = "control-label"),
             tooltip(
               icon("circle-question"),
               "Minimum accepted value for calculated overlap probability. Default = 0.6. Do not suggest adjusting this without a good reason."
             )
           ),

            div(
              class = "slider-input-row",
              div(
                class = "slider-col",
                sliderInput(
                  "merge_similarity",
                  label = NULL,
                  min = 0.1,
                  max = 1,
                  value = 0.6,
                  step = 0.1,
                  ticks = FALSE
                ),
                div(
                  class = "custom-slider-ticks",
                  span(style = "left: calc(10px + ((0.2 - 0.1) * (100% - 20px) / 0.9));", "0.2"),
                  span(style = "left: calc(10px + ((0.4 - 0.1) * (100% - 20px) / 0.9));", "0.4"),
                  span(style = "left: calc(10px + ((0.6 - 0.1) * (100% - 20px) / 0.9));", "0.6"),
                  span(style = "left: calc(10px + ((0.8 - 0.1) * (100% - 20px) / 0.9));", "0.8"),
                  span(style = "left: calc(10px + ((1 - 0.1) * (100% - 20px) / 0.9));", "1.0")
                )
              ),
              div(
                class = "numeric-col",
                numericInput(inputId = "merge_similarity_num",
                             label = "",
                             min = 0,
                             max = 1,
                             value = 0.6)
              )
            ),



           tags$hr(class = "section-hr"),

           div(
             class = "sidebar-checkbox-row",

             checkboxInput(
               "merge_advanced_settings",
               label = NULL,
               value = FALSE
             ),

             tags$span(
               "Advanced settings",
               class = "sidebar-checkbox-label"
             ),

             tooltip(
               icon("circle-question"),
               "Show additional advanced settings for PANDAseq."
             )
           ),

            conditionalPanel(
              condition = "input.merge_advanced_settings === true",
              br(),
                     div(
             class = "d-flex align-items-center justify-content-between",
             tags$label("Forward primer", class = "control-label"),
             tooltip(
               icon("circle-question"),
               "Sequence or length in bp of forward primer to trim. Will also remove any distal adapters (barcodes, UMIs, etc.)"
             )
           ),
              textInput(
                "forward_primer",
                label = NULL,
                placeholder = "Sequence or trim length"
              ),
              div(style = "height: 10px;"),
           div(
             class = "d-flex align-items-center justify-content-between",
             tags$label("Reverse primer", class = "control-label"),
             tooltip(
               icon("circle-question"),
               "Sequence or length in bp of reverse primer to trim. Will also remove any distal adapters (barcodes, UMIs, etc.)"
             )
           ),
              textInput(
                "reverse_primer",
                label = NULL,
                placeholder = "Sequence or trim length"
              ),
           br(),
           div(
             class = "d-flex align-items-center justify-content-between",
             tags$label("Output type", class = "control-label"),
             tooltip(
               icon("circle-question"),
               "Select the file extension for output merged files."
             )
           ),
           selectInput(
             "merge_filetype",
             label = NULL,
             choices = c(
               "FASTQ",
               "FASTA"
             ),
             selected = "FASTQ"
           ),

           br(),

           div(
             class = "d-flex align-items-center justify-content-between",
             tags$label("Additional flags", class = "control-label")),

           div(
             class = "sidebar-checkbox-row",
             checkboxInput(
               "merge_file_errors",
               label = NULL,
               value = TRUE
             ),
             tags$span(
               "Show file parsing errors",
               class = "sidebar-checkbox-label"
             ),
             tooltip(
               icon("circle-question"),
               "Provide error about the file parsing."
             )
           ),

           div(style = "height:4px"),

           div(
             class = "sidebar-checkbox-row",
             checkboxInput(
               "merge_statistics",
               label = NULL,
               value = TRUE
             ),
             tags$span(
               "Optional statistics",
               class = "sidebar-checkbox-label"
             ),
             tooltip(
               icon("circle-question"),
               "Show some optional statistics."
             )
           ),

           div(style = "height:4px"),

           div(
             class = "sidebar-checkbox-row",
             checkboxInput(
               "merge_discard_ambig",
               label = NULL,
               value = FALSE
             ),
             tags$span(
               "Discard N's",
               class = "sidebar-checkbox-label"
             ),
             tooltip(
               icon("circle-question"),
               "Discard reads containing ambiguous (N) base calls."
             )
           ),


           div(style = "height:4px"),

           div(
             class = "sidebar-checkbox-row",
             checkboxInput(
               "merge_seq_build_log",
               label = NULL,
               value = TRUE
             ),
             tags$span(
               "Sequence build info",
               class = "sidebar-checkbox-label"
             ),
             tooltip(
               icon("circle-question"),
               "Provide information about the building of a sequence."
             )
           ),

           div(style = "height:4px"),

           div(
             class = "sidebar-checkbox-row",
             checkboxInput(
               "merge_seq_reconstruction_log",
               label = NULL,
               value = FALSE
             ),
             tags$span(
               "Sequence reconstruction info",
               class = "sidebar-checkbox-label"
             ),
             tooltip(
               icon("circle-question"),
               "Show excruciating detail about reconstruction. Not recommended."
             )
           ),

           div(style = "height:4px"),

           div(
             class = "sidebar-checkbox-row",
             checkboxInput(
               "merge_kmer",
               label = NULL,
               value = FALSE
             ),
             tags$span(
               "Show k-mer table",
               class = "sidebar-checkbox-label"
             ),
             tooltip(
               icon("circle-question"),
               "Show information about building the k-mer table. Not recommended."
             )
           ),
           )),

          "shm.app" = div(
            class = "sidebar-controls",

            tags$hr(class = "section-hr"),

            h5("Sample Summary"),
            verbatimTextOutput("mutation_run_summary"),



            conditionalPanel(
              condition = "input.mutation_upload_mode == 'mutation_single_sample'",
              tags$hr(class = "section-hr"),
              h5("Additional Parameters"),
              div(
                class = "mutation-sidebar-analysis",
                div(
                  class = "d-flex align-items-center justify-content-between",
                  tags$label("Output prefix", class = "control-label"),
                  tooltip(
                    icon("circle-question"),
                    "Prefix reserved for future mutation result files."
                  )
                ),
                textInput("mutation_output_prefix", label = NULL, placeholder = "mutations"),
                br(),
                mutation_analysis_settings_ui("mutation", show_heading = FALSE)
              )
            )
          ),

          "cut.app" = div(
            class = "sidebar-controls",

            tags$hr(class = "section-hr"),

            h5("Sample Summary"),
            verbatimTextOutput("nuclease_run_summary"),

            conditionalPanel(
              condition = "input.nuclease_upload_mode == 'nuclease_single_sample'",
              tags$hr(class = "section-hr"),
              h5("Additional Parameters"),
              div(
                class = "mutation-sidebar-analysis",
                nuclease_parameters_ui("nuclease", show_heading = FALSE, sidebar_breaks = TRUE)
              )
            )
          ),

          "help.app" = div(
            class = "sidebar-controls",

            h5("Help"),
            tags$p("Help and documentation go here.")
          )
  )
}
