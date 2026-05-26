library(shiny)
library(shinyjs)
library(bslib)
library(later)
library(DT)

# Upload limit 50 MB
options(shiny.maxRequestSize = 50 * 1024^2)
