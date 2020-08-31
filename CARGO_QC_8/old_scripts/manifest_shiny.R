library(xlsx)
library(tidyverse)
library(tibble)
library(shiny)

source("manifest_functions.R")

manifest_ui <- function() { tabPanel( titlePanel("Plate Samples"), mainPanel( manifest_main() ) ) }
manifest_main <- function() { sidebarLayout( manifest_sidebar(), manifest_tables() ) }

manifest_sidebar <- function() {
  sidebarPanel(
    tags$head(
      tags$style(type="text/css", "select { min-width: 220px; max-width: 220px; }"),
      tags$style(type="text/css", ".span4 { min-width: 220px; max-width: 220px; }"),
      tags$style(type="text/css", ".well { min-width: 220px; max-width: 220px; }")
    ),
    actionButton("getPassedSamples", "Load Passed Samples"),
    uiOutput("ui"),
    downloadButton("downloadManifest", "Download Manifest")
  )
}

manifest_tables <- function() {
  mainPanel(
    tabsetPanel( id = "tabs", type = "pills",
      tabPanel( "Passed QC", DT::dataTableOutput("passedQC") ),
      tabPanel( "Manifest", DT::dataTableOutput("manifestTable")),
      tabPanel( "Plate Layout", )
    ),
    textOutput("debug")
  )
}

manifest_server <- function(input, output, session) {
  observeEvent(input$getPassedSamples, {
    field_names <- colnames(get_data())
    updateActionButton(session, "getPassedSamples", label = "Reload Passed Samples")
    output$ui <- renderUI({ tagList(
      selectInput("id_col", "Sample ID Column", choices = field_names),
      checkboxGroupInput("by_cols", "Balance by Columns", choices = set_names(field_names),
                         select = c("site", "Age_category", "Asthma")),
      checkboxGroupInput("add_cols", "Add to Manifest", choices = set_names(field_names))
    )})
  })
  
  get_data <<- eventReactive(input$getPassedSamples, {
    load("savePassed.RData")
    forCorey
  })
  
  get_manifest <- reactive({
    req(input$id_col)
    controls <- c("Hypo-methylated Control", "Hyper-methylated control")
    col_vals = c(Species = 'Homo sapiens', `Tissue Source` = 'Whole Blood') #SHINY get default values from user.
    create_manifest_2(get_data(), controls, input$id_col, input$by_cols, col_vals, input$add_cols)
  })
  
  output$passedQC <- DT::renderDataTable({ get_data() }, options = list(pageLength = 96))
  output$manifestTable <- DT::renderDataTable({ get_manifest() }, options = list(pageLength = 96))
  
  output$downloadManifest <- downloadHandler(
    filename = function() { paste('Manifest-', Sys.Date(), '.xlsx', sep='') },
    content = function(con) { write.xlsx(get_manifest(), file = con, showNA = FALSE) })
  
  output$debug <- renderText({ # req(input$fileUploaded)
    input$by_cols
  })
}