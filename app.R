library(shiny)
# library(datasets)
# library(ggpairs)
# require(GGally)
require(ggplot2)
library(tidyverse)
library(gridExtra)
library(plotly)
library(DT)
library(reshape2)
library(table1)
# library(pals)
# library(ggthemes)
library(xlsx)
library(readxl)


source("data_ui.R", local = TRUE)
source("qc_ui.R", local = TRUE)
source("manifest_ui.R", local = TRUE)

ui <- shinyUI(fluidPage(
  titlePanel("CAAPA2 QC Metrics"),
  tabsetPanel( id = "top_tabs", data_ui, qc_ui, manifest_ui )
  
  # tabsetPanel( id = "top_tabs", qc_ui, manifest_ui )
  # tabsetPanel( id = "top_tabs", data_ui, manifest_ui )
  # tabsetPanel( id = "top_tabs", data_ui, qc_ui )
))

source("data_server.R", local = TRUE)
source("qc_server.R", local = TRUE)
source("manifest_server.R", local = TRUE)

server <- shinyServer(function(input, output, session) {
  data_server(input, output, session)
  qc_server(input, output, session)
  manifest_server(input, output, session)
})

shinyApp(ui, server)