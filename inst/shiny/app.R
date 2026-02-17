#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

if (!require("shiny", quietly = TRUE)){
  install("shiny", quietly = TRUE)
  library("shiny", quietly = TRUE)}

if (!require("bslib", quietly = TRUE)){
  install("bslib", quietly = TRUE)
  library("bslib", quietly = TRUE)}

if (!require("shinyjs", quietly = TRUE)){
  install("shinyjs", quietly = TRUE)
  library("shinyjs", quietly = TRUE)}

if (!require("spsComps", quietly = TRUE)){
  install("spsComps", quietly = TRUE)
  library("spsComps", quietly = TRUE)}

if (!require("markdown", quietly = TRUE)){
  install("markdown", quietly = TRUE)
  library("markdown", quietly = TRUE)}

suppressMessages(library(tasAnalyzer))


options(shiny.maxRequestSize = 100*1024^2)

# Define UI for application that draws a histogram
ui <- page_fluid(
    useShinyjs(),
    # Application title
    titlePanel(h1("tasAnalyzer",
                  style = "font-size: 90px;",
                  style = "text-align:center;"
                  ),
               windowTitle = "tasAnalyzer"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            h4("Options"),
            radioButtons(inputId = "merge.reads",
                        "Merge paired-end reads?",
                        choiceNames = c("Yes","No"),
                        choiceValues = c(TRUE,FALSE),
                        selected = TRUE,
                        inline = TRUE
                        ),
            radioButtons(inputId = "measure.shm",
                         "Enable somatic hypermutation module?",
                         choiceNames = c("Yes","No"),
                         choiceValues = c(TRUE,FALSE),
                         selected = TRUE,
                         inline = TRUE
            ),
            radioButtons(inputId = "protein.mutations",
                         "Measure amino acid mutations?",
                         choiceNames = c("Yes","No"),
                         choiceValues = c(TRUE,FALSE),
                         selected = TRUE,
                         inline = TRUE
            ),
            numericInput(inputId = "sequence.alignment.count",
                         label = "Number of sequences to align?",
                         value = 10,
                         min = 1,
                         step = 1
            ),
            numericInput(inputId = "read.frequency.limit",
                         label = "Minimum read frequency for analysis?",
                         value = 0.001,
                         min = 0.0001,
                         max = 5,
                         step = 0.0001
            ),
            h6("(only for non-UMI samples)"
               ),
            radioButtons(inputId = "dna.repair.pathways",
                         "Predict DNA repair pathway use?",
                         choiceNames = c("Yes","No"),
                         choiceValues = c(TRUE,FALSE),
                         selected = FALSE,
                         inline = TRUE
            ),
            radioButtons(inputId = "PhyloTree",
                         "Plot phylogenetic trees of multiple sequence alignments?",
                         choiceNames = c("Yes","No"),
                         choiceValues = c(TRUE,FALSE),
                         selected = FALSE,
                         inline = TRUE
            ),
            radioButtons(inputId = "multicore",
                         "Enable multicore processing?",
                         choiceNames = c("Yes","No"),
                         choiceValues = c(TRUE,FALSE),
                         selected = TRUE,
                         inline = TRUE
            ),
            # actionButton(inputId = "apply",
            #              label = "Apply Settings"
            #              ),
            actionButton(inputId = "run",
                         label = "Begin Analysis")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          navset_card_tab(
            nav_panel("Upload",
                      h4("Upload files"),
                      card(layout_columns(
                          fileInput(inputId = "sheet",
                                    label = "Configuration spreadsheet (.csv)",
                                    accept = ".xlsx"
                          ),
                          fileInput(inputId = "unpaired.reads",
                                    label = "Unmerged fastq files",
                                    multiple = TRUE,
                          ),
                          downloadLink("csv", label = "CSV Configuation Template"),
                          layout_columns(selectInput(inputId = "os", label = "Operating system",
                                      choices = c("Windows" = "bat",
                                                  "MacOS" = "command",
                                                  "Linux" = "sh"),
                                      width = '200px'),
                          downloadLink("script", label = "PANDAseq Automation Script"), col_widths = 12),
                          col_widths = 6)),
              tableOutput("files")
            ),
            nav_panel("Results",
                      uiOutput("graphs")
                      ),
            nav_panel("Download",
                      uiOutput("download.data.ui")
                      ),
            nav_spacer(),
            nav_panel("Help",
                      htmltools::includeMarkdown("help/README.md")),
            nav_menu(
              title = "Links"
            ),
            id = "tabs"),
       hidden(card(textOutput("status"),
                   id = "status.card"))
      )
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {

   shiny.env <- TRUE

   hideTab("tabs","Results")
   hideTab("tabs","Download")

   output$csv <- downloadHandler(filename = function(){paste0("Input.csv")},
                                 content = function(file){
                                   file.copy("templates/Input.csv", file)
                                 }, contentType = "text/csv")

   output$script <- downloadHandler(filename = function(){paste0("run_pandaseq.",input$os)},
                                 content = function(file){
                                   file.copy(paste0("scripts/run_pandaseq.",input$os), file)
                                 })




   # tas_load_dependencies()

   observe({
     if (is.null(input$sheet) || is.null(input$unpaired.reads)) {
       # If files are not uploaded, disable the button
       shinyjs::disable("run")
     } else if (!is.null(input$sheet) && !is.null(input$unpaired.reads)) {
       # If files are uploaded, enable the button
       shinyjs::enable("run")
     }
   })

   fileTable <- reactive({
     req(isTruthy(input$sheet) || isTruthy(input$unpaired.reads))
     x <- rbind(input$sheet, input$unpaired.reads)
     return(x)
   })

   fileTable.display <- reactive({
     req(isTruthy(input$sheet) || isTruthy(input$unpaired.reads))
     x <- rbind(input$sheet, input$unpaired.reads)
     x <- x[,1:2]
     x$size <- x$size/(1024^2)
     colnames(x) <- c("File Name", "Size (MB)")
     return(x)
   })

   output$files <- renderTable({
     fileTable.display()
     })

   observeEvent(input$run,{
     shinyjs::disable("run")
     showTab("tabs","Results")
     showTab("tabs","Download")
     })


   settings <- eventReactive(input$run,{
     list(WorkingDirectory = tempdir(),
          merge.reads = input$merge.reads,
          measure.shm = input$measure.shm,
          protein.mutations = input$protein.mutations,
          sequence.alignment.count = input$sequence.alignment.count,
          read.frequency.limit = input$read.frequency.limit,
          dna.repair.pathways = input$dna.repair.pathways,
          PhyloTree = input$PhyloTree,
          multicore = input$multicore,
          OperatingSystem = if (Sys.info()["sysname"]=="Darwin"){"MacOS"}else{Sys.info()["sysname"]},
          nCores=availableCores())
   })

   results <- eventReactive(input$run,{
     withProgress(tasAnalyzer::tas_analyze(input$sheet$datapath,
                                           config = 'shiny',
                                           shiny.settings = settings(),
                                           shiny.env = TRUE,
                                           shiny.fileTable = fileTable()
       ),
       value = 0, message = "Loading dependencies (step 1 of 8)")
     })


  # #instant results to check events without waiting for the full script to execute
  #  results <- eventReactive(input$run,list(Sequences=list(All=c(1,2))))

  calc <- eventReactive(input$run,{
    req(results(), cancelOutput = TRUE)
    analysis.done <- isolate(results()$Sequences)
    req(analysis.done)
    "Analysis finished."
  })

  observeEvent(calc(),{
    req(results(), cancelOutput = TRUE)
    analysis.done <- isolate(results()$Sequences)
    req(analysis.done)
    shinyjs::toggle("status.card")
    # Sample.Names <- names(results()$Sequences$All)
    # MSA.args <- paste0(lapply(Sample.Names(),function(x){
    #   str_c("renderText('",x,"'), ",
    #         "renderPlot(results()$MSA[[",x,"]]$DNA)")
    # }),collapse=",")
  })

   output$status <- renderText(calc())

   # Server-side logic for the download action
   WD <- tempdir()
   if (substr(WD,nchar(WD),nchar(WD)) != "/"){
     WD <- paste0(WD,"/")
   }
   output$zip <- downloadHandler(filename = function(){paste0("tasAnalyzer results - ",Sys.time(),".zip")},
                                          content = function(file){
                                            file.copy(str_c(WD,"analyzed.zip"), file)},
                                            contentType = "application/zip")

   output$download.data.ui <- renderUI({
     req(results(), cancelOutput = TRUE)
     analysis.done <- isolate(results()$Sequences)
     req(analysis.done)
     downloadButton("zip", "Download Results")
   })






   output$graphs <- renderUI({
     navset_underline(nav_panel(title = "Summary",
                                renderPlot(results()$Mutations$SummaryGraphs$MutAll),
                                if (settings()$measure.shm){
                                renderPlot(results()$Mutations$SummaryGraphs$MutC_All)},
                                if (settings()$dna.repair.pathways){
                                  renderPlot(results()$Mutations$SummaryGraphs$MutTypes)}
                                ))
                      # nav_panel(title = "MSA", paste0(MSA.args))
                      # )

                      # paste0(str_c("title = ",results()$Sequences$All,", renderPlot(", results()$), collapse = ","))
   })
}

# Run the application
shinyApp(ui = ui, server = server)
