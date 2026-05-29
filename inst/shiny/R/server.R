server <- function(input, output, session) {
  operation_rv <- reactiveValues(active = NULL)
  setup_sidebar_output(input, output)
  setup_merge_server(input, output, session, operation_rv)
  setup_merge_settings_server(input, output, session)
  setup_mutation_server(input, output, session, operation_rv)
  setup_nuclease_server(input, output, session, operation_rv)
}
