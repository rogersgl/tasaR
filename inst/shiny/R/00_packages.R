library(shiny)
library(shinyjs)
library(bslib)
library(later)
library(DT)
library(future)
library(promises)
library(tasaR)

future::plan(future::multisession)

# Upload limit 50 MB
options(shiny.maxRequestSize = 50 * 1024^2)
