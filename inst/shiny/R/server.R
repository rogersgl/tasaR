server <- function(input, output, session) {
  setup_sidebar_output(input, output)
  setup_merge_server(input, output, session)
  setup_merge_settings_server(input, output, session)
  setup_mutation_server(input, output, session)
  setup_nuclease_server(input, output, session)
}
