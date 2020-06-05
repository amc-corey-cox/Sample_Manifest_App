#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(tidyverse)
library(readxl)
library(xlsx)
library(rlang)

# source("helpers.R")
source("manifest.R")

# make_names <- function(p, s) {
#   paste0(p, str_remove_all(s, "[^[:alpha:]]"))
# }

upload_ctrl <- function() {
  # fileInput("files", "Upload Phenotype File", multiple = TRUE,
  fileInput("files", "Upload Phenotype File", multiple = FALSE,
    accept = c(".csv", ".xls", ".xlsx", "text/csv",
      "text/comma-separated-values,text/plain",
      "application/vnd.ms-excel",
      "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    )
}

file_ctrl <- function() {
  conditionalPanel( condition = "output.fileUploaded",
    checkboxInput("col_names", "Column names in header", TRUE),
    numericInput("skip", "Skip lines:", 0, min = 0)
  )
}

excel_ctrl <- function() {
  conditionalPanel( condition = "output.fileExcel",
    selectInput("sheet", "Sheet", choices = NULL)
  )
}

txt_ctrl <- function() {
  conditionalPanel( condition = "output.fileUploaded && ! output.fileExcel",
    radioButtons("delim", "Separator", selected = "\t",
                 choices = c(Comma = ",", Semicolon = ";", Tab = "\t")),
    radioButtons("quote", "Quote", selected = '"',
                 choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'")),
  )
}

id_ctrl <- function() {
  conditionalPanel( condition = "output.fileUploaded",
    # radioButtons("id_radio", "Sample ID Column", choices = c(1))
    selectInput("id_col", "Sample ID Column", choices = NULL)
  )
}

ui <- fluidPage(
  splitLayout(cellWidths = c("60%", "20%", "20%"),
    titlePanel("CARGO Manifest Generator", windowTitle = "CARGO Manifest Generator"),
    radioButtons( "disp", "Display", selected = "head", inline = TRUE,
                  choices = c(Head = "head", All = "all")),
    downloadButton("downloadManifest", "Download"),
    tags$style(type='text/css', "#disp { margin-top: 15px; }"),
    tags$style(type='text/css', "#downloadManifest { margin-top: 20px; }")
  ),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel( condition = "input.tabs == 'Phenotype'",
        upload_ctrl(), file_ctrl(), excel_ctrl(), txt_ctrl()
      ),
      conditionalPanel( condition = "input.tabs == 'Phenotype' || input.tabs == 'Cleanup'", id_ctrl() ),
      conditionalPanel( condition = "input.tabs == 'Cleanup'",
        checkboxGroupInput("clean_cols", "Convert Columns to Number", choices = NULL),
        p("Select Columns to remove all non-digits and convert to number.")
      ),
      conditionalPanel(
        condition = "input.tabs == 'Manifest'",
        checkboxGroupInput("by_cols", "Balance by Columns", choices = NULL)
      )
    ),
    mainPanel(
      tabsetPanel( id = "tabs", type = "pills",
        # tabPanel( "Phenotype", uiOutput('uploadUI') ),
        tabPanel( "Phenotype", tableOutput("pheno") ),
        tabPanel("Cleanup", tableOutput("cleaned_pheno")),
        tabPanel("Manifest", tableOutput("manifest"))
      )
    )
  ),
  hr(),
  textOutput("debug")
)

server <- function(input, output, session) {
  output$fileUploaded <- reactive({ return(!is.null(input$files)) })
  outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)

  output$fileExcel <- reactive({
    if (!is.null(input$files) &&
        !is.na(excel_format(input$files$datapath))) {
      updateSelectInput(session, "sheet", choices = excel_sheets(input$files$datapath))
      return(TRUE)
    }
    updateSelectInput(session, "sheet", choices = character(0))
    return(FALSE)
  })
  outputOptions(output, 'fileExcel', suspendWhenHidden = FALSE)
  
  # output$uploadUI <- renderUI({ req(input$files)
  #   nTabs = nrow(input$files)
  #   myTabs <- map(input$files$name, tabPanel)
  #   # myTabs = lapply(paste('Tab', 1: nTabs), tabPanel)
  #   do.call(tabsetPanel, myTabs)
  # })

  read_pheno <- reactive({
    # return()
    if (is.na(excel_format(input$files$datapath))) {
      pheno <- read_delim(input$files$datapath, col_names = input$col_names,
                 delim = input$delim, quote = input$quote, skip = input$skip)
    } else {
      pheno <- read_excel( input$files$datapath, skip = input$skip, col_names = input$col_names)
    }
    updateSelectInput(session, "id_col", choices = colnames(pheno))
    map(list("align_cols", "clean_cols", "filter_cols", "by_cols"),
        ~ updateCheckboxGroupInput(session, .x, choices = set_names(colnames(pheno))))
    return(pheno)
  })
  
  get_pheno <- reactive({ req(input$files)
    pheno <- read_pheno()
    req(input$id_col %in% colnames(pheno))
    pheno %>% rename("Sample ID" = input$id_col) %>%
      mutate(`Sample ID` = as.character(`Sample ID`))
  })
  
  clean_pheno <- reactive({
    # clean_cols <- c("Avg [ng/ul]", "% CV")
    clean <- function(s) { str_remove_all(s, "[^[:digit:].]") %>% as.double() }
    get_pheno() %>%
      mutate_at(input$clean_cols, clean)
  })

  filter_pheno <- reactive({
    filters <- c("") #SHINY use paste to create this from user input
    clean_pheno() %>% filter(!!! parse_exprs(filters))
  })

  get_manifest <- reactive({ req(input$by_cols)
    #RSHINY Get this from UI
    controls <- c("Hypo-methylated Control", "Hyper-methylated control")
    by_cols <- input$by_cols
    col_vals = c(Species = 'Homo sapiens', `Tissue Source` = 'Whole Blood') #SHINY get default values from user.
    create_manifest(filter_pheno(), controls, by_cols, col_vals)
  })

  createTableOutput <- function(df) {
    if (input$disp == "head") { return(head(df)) }
    return(df)
  }

  output$pheno <- renderTable({ createTableOutput(get_pheno()) })
  output$cleaned_pheno <- renderTable({ createTableOutput(clean_pheno()) })
  output$filtered_pheno <- renderTable({ createTableOutput(filter_pheno()) })
  output$manifest <- renderTable({ createTableOutput(get_manifest()) })

  output$downloadManifest <- downloadHandler(
    filename = function() {
      paste('data-', Sys.Date(), '.xlsx', sep='')
    },
    content = function(con) {
      write.xlsx(get_manifest(), file = con, showNA = FALSE)
    }
  )

  output$debug <- renderText({ # req(input$fileUploaded)
    input$id_col
  })
}
# Run the app ----
shinyApp(ui, server)