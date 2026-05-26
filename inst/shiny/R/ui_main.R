app_footer_ui <- function() {
    div(
      class = "app-footer",
  
      div(
        class = "footer-content",
  
        div(
          class = "footer-left",
          "Built with R • Shiny • Bioconductor"
        ),
  
        div(
          class = "footer-center",
          "Version 0.2.0"
        ),
  
        div(
          class = "footer-right",
          "tasaR © 2026 Geoffrey L. Rogers"
        )
      )
    )
}

build_ui <- function() {
tagList(
  useShinyjs(),
  tags$head(
    includeCSS("www/css/base.css"),
    includeCSS("www/css/sidebar.css"),
    includeCSS("www/css/merge.css"),
    tags$script(HTML("
      function resizeMutationReferenceTextareas() {
        document.querySelectorAll('textarea[id$=\"_reference_sequence\"]').forEach(function(el) {
          el.style.height = 'auto';
          var nextHeight = Math.min(160, Math.max(38, el.scrollHeight));
          el.style.height = nextHeight + 'px';
          el.style.overflowY = el.scrollHeight > 160 ? 'auto' : 'hidden';
        });
      }
      document.addEventListener('input', function(e) {
        if (e.target && e.target.matches('textarea[id$=\"_reference_sequence\"]')) {
          resizeMutationReferenceTextareas();
        }
      });
      document.addEventListener('wheel', function(e) {
        if (e.target && e.target.matches('input[type=\"number\"]')) {
          e.target.blur();
        }
      }, { passive: true });
      document.addEventListener('shiny:value', resizeMutationReferenceTextareas);
      document.addEventListener('shown.bs.modal', resizeMutationReferenceTextareas);
      document.addEventListener('DOMContentLoaded', resizeMutationReferenceTextareas);
    "))
  ),
  page_navbar(title = "tasaR",
              id = "main_nav",
              theme = bs_theme(version = 5),
              # potential sidebar background colors or additional accent colors: #DBEAFE, #89C2D9
              sidebar = sidebar(
                title = "Status",
                open = "always",
                resizable = FALSE,
                width = "325px",
                uiOutput("sidebar_ui")
              ),
             home_page_ui(),
             merge_page_ui(),
             mutation_page_ui(),
             nuclease_page_ui(),
             nav_spacer(),
             nav_panel("Help",
                      value = "help.app",
                      class = "ms-auto"
                      ),
             nav_menu(title = "Links",
                      align = "right",
                      value = "links.app",
                      nav_item(
                        tags$a(
                          icon("circle-question"),
                          href = "",
                          "About",
                          target = "_blank"
                        )
                      ),
                      nav_item(
                        tags$a(
                          icon("github"),
                          href = "https://github.com/rogersgl/tasaR",
                          "tasaR GitHub",
                          target = "_blank"
                        )
                        ),
                      nav_item(
                        tags$a(
                          icon("docker"),
                          href = "https://TBD_WHEN_DONE",
                          "tasaR Docker",
                          target = "_blank"
                        )
                      ),
                      nav_item(
                        tags$a(
                          icon("github"),
                          href = "https://github.com/neufeld/pandaseq",
                          "PANDAseq GitHub",
                          target = "_blank"
                        )
                      ),
                      nav_item(
                        tags$a(
                          icon("pen-to-square"),
                          href = "",
                          "GPLv3 License",
                          target = "_blank"
                        )
                      ),
                      nav_item(
                        tags$a(
                          icon("user-pen"),
                          href = "",
                          "About the author",
                          target = "_blank"
                        )
                      )
             )),
  app_footer_ui()
)
}
