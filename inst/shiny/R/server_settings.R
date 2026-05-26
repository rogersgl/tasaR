setup_merge_settings_server <- function(input, output, session) {
    ##########
    output$merge_settings <- renderText({
        paste(
          "Output prefix:", input$merge_prefix,
          "\nAlgorithm:", input$merge_engine,
          "\nMin overlap:", input$min_overlap,
          "\nMin length:", input$min_length,
          "\nMax length:", input$max_length,
          "\nOverlap threshold:", input$merge_similarity
        )
      })
  
    observeEvent(input$min_overlap_num, {
      val <- suppressWarnings(as.numeric(input$min_overlap_num))
      if (is.na(val)) return()
      val <- max(5, min(100, val))
      if (!identical(val, input$min_overlap)) {
        updateSliderInput(session, "min_overlap", value = val)
      }
    }, ignoreInit = TRUE)
  
    observeEvent(input$min_overlap, {
      val <- suppressWarnings(as.numeric(input$min_overlap))
      if (is.na(val)) return()
      if (!identical(val, input$min_overlap_num)) {
        updateNumericInput(session, "min_overlap_num", value = val)
      }
    }, ignoreInit = TRUE)
  
  
    observeEvent(input$min_length_num, {
      val <- suppressWarnings(as.numeric(input$min_length_num))
      if (is.na(val)) return()
      val <- max(50, min(500, val))
      if (!identical(val, input$min_length)) {
        updateSliderInput(session, "min_length", value = val)
      }
    }, ignoreInit = TRUE)
  
    observeEvent(input$min_length, {
      val <- suppressWarnings(as.numeric(input$min_length))
      if (is.na(val)) return()
      if (!identical(val, input$min_length_num)) {
        updateNumericInput(session, "min_length_num", value = val)
      }
    }, ignoreInit = TRUE)
  
  
    observeEvent(input$max_length_num, {
      val <- suppressWarnings(as.numeric(input$max_length_num))
      if (is.na(val)) return()
      val <- max(50, min(500, val))
      if (!identical(val, input$max_length)) {
        updateSliderInput(session, "max_length", value = val)
      }
    }, ignoreInit = TRUE)
  
    observeEvent(input$max_length, {
      val <- suppressWarnings(as.numeric(input$max_length))
      if (is.na(val)) return()
      if (!identical(val, input$max_length_num)) {
        updateNumericInput(session, "max_length_num", value = val)
      }
    }, ignoreInit = TRUE)
  
  
  
    observeEvent(input$merge_similarity_num, {
      val <- suppressWarnings(as.numeric(input$merge_similarity_num))
      if (is.na(val)) return()
      val <- max(0, min(1, val))
      if (!identical(val, input$merge_similarity)) {
        updateSliderInput(session, "merge_similarity", value = val)
      }
    }, ignoreInit = TRUE)
  
    observeEvent(input$merge_similarity, {
      val <- suppressWarnings(as.numeric(input$merge_similarity))
      if (is.na(val)) return()
      if (!identical(val, input$merge_similarity_num)) {
        updateNumericInput(session, "merge_similarity_num", value = val)
      }
    }, ignoreInit = TRUE)
}
