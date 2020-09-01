library(shiny)
library(datasets)
# library(ggpairs)
require(GGally)
require(ggplot2)
library(tidyverse)
library(gridExtra)
library(plotly)
library(DT)
library(reshape2)
library(table1)
library(pals)
library(ggthemes)

source("data_ui.R")
source("qc_ui.R")
source("manifest_ui.R")

source("data_server.R")
source("qc_server.R")
source("manifest_server.R")

ui <- shinyUI(fluidPage(
  titlePanel("CAAPA2 QC Metrics"),
  tabsetPanel( id = "top_tabs", data_ui, qc_ui, manifest_ui )
  
  # tabsetPanel( id = "top_tabs", qc_ui, manifest_ui )
  # tabsetPanel( id = "top_tabs", data_ui, manifest_ui )
  # tabsetPanel( id = "top_tabs", data_ui, qc_ui )
))

server <- shinyServer(function(input, output, session) {
  data_server(input, output, session)
  # qc_server(input, output, session)
  manifest_server(input, output, session)
})

shinyApp(ui, server)