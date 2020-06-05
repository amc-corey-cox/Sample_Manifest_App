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

source("helpers.R")
make_names <- function(p, s) {
  paste0(p, str_remove_all(s, "[^[:alpha:]]"))
}

source("manifest.R")

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

tab_upload <- function() {
  tabPanel( "Read",
    fileInput("pheno_file", "Choose CSV File", multiple = FALSE,
      accept = c(".csv", ".xls", ".xlsx", "text/csv",
       "text/comma-separated-values,text/plain",
       "application/vnd.ms-excel",
       "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")),
    file_ctrl(),
    excel_ctrl(),
    txt_ctrl()
  )
}

tab_align <- function() {
  tabPanel( "Align",
    conditionalPanel( condition = "output.fileUploaded",
      checkboxGroupInput("align_cols", "Columns to Align", choices = NULL),
      uiOutput("align_ui")
    )
  )
}

tab_cleanup <- function() {
  tabPanel( "Cleanup",
    conditionalPanel( condition = "output.fileUploaded",
      checkboxGroupInput("clean_cols", "Columns to Cleanup", choices = NULL) #,
      # uiOutput("cleanup_ui")
    )
  )
}

tab_filter <- function() {
  tabPanel( "Filter",
    conditionalPanel( condition = "output.fileUploaded",
      checkboxGroupInput("filter_cols", "Columns to Filter", choices = NULL)
    )
  )
}

tab_balance <- function() {
    tabPanel( "Balance",
      conditionalPanel( condition = "output.fileUploaded",
        checkboxGroupInput("by_cols", "Columns to Balance", choices = NULL)
      )
    )
}

main_ctrl <- function() {
  splitLayout( 
    radioButtons( "disp", "Display", selected = "head", inline = TRUE,
                  choices = c(Head = "head", All = "all")),
    downloadButton("downloadManifest", "Download", class = 'rightAlign')
  )
}

ui <- fluidPage(
  titlePanel("CARGO Manifest Generator"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel( type = "pills", tab_upload(), tab_align(), tab_cleanup(), tab_filter(), tab_balance() ),
      textOutput("debug")
    ),
    mainPanel( main_ctrl(),
               tabsetPanel( type = "pills",
                 tabPanel("Phenotype File", tableOutput("pheno")),
                 tabPanel("Cleaned Phenotype", tableOutput("cleaned_pheno")),
                 tabPanel("Filtered Phenotype", tableOutput("filtered_pheno")),
                 tabPanel("Manifest", tableOutput("manifest"))
               )
    )
  )
)

server <- function(input, output, session) {
  output$fileUploaded <- reactive({ return(!is.null(input$pheno_file)) })
  outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)
  
  output$fileExcel <- reactive({
    if (!is.null(input$pheno_file) &&
        !is.na(excel_format(input$pheno_file$datapath))) {
      updateSelectInput(session, "sheet", choices = excel_sheets(input$pheno_file$datapath))
      return(TRUE)
    }
    updateSelectInput(session, "sheet", choices = character(0))
    return(FALSE)
  })
  outputOptions(output, 'fileExcel', suspendWhenHidden = FALSE)
  
  get_pheno <- reactive({
    if (is.na(excel_format(input$pheno_file$datapath))) {
      pheno <- read_delim(input$pheno_file$datapath, col_names = input$col_names,
                 delim = input$delim, quote = input$quote, skip = input$skip)
    } else {
      pheno <- read_excel( input$pheno_file$datapath, skip = input$skip, col_names = input$col_names)
    }
    map(list("align_cols", "clean_cols", "filter_cols", "by_cols"),
        ~ updateCheckboxGroupInput(session, .x, choices = set_names(colnames(pheno))))
    return(pheno)
  })
  
  # output$cleanup_ui <- renderUI({
  #   map(input$clean_cols, function(i) {
  #     verticalLayout(
  #      h4(i),
  #       splitLayout(
  #         textInput(paste0("clean_txt_", i), "Remove"),
  #         selectInput(paste0("clean_type_", i), label = "Type", multiple = FALSE,
  #                     choices = c("logical", "integer", "double", "character"), selected = "double")
  #       )
  #     )
  #   })
  # })
  
  output$align_ui <- renderUI({
    map(input$align_cols, function(n) {
      splitLayout(
        selectInput(make_names("align_input_", n), label = n, multiple = FALSE,
                    selectize = FALSE, choices = col_names)
      )
    })
  })
  
  get_align_cols <- reactive({
    req(input$align_cols, cancelOutput = TRUE)
    # %>% set_names()
    input$align_cols %>% set_names() %>%
      map_chr(function(n) {
        pluck(input, make_names("align_input_", n))
      }) %>% split(names(.), .)
  })
  
  get_clean_cols <- reactive({
    req(input$clean_cols)
    input$clean_cols %>% set_names()
  })
  
  clean_pheno <- reactive({
    # align_cols <- c("Sample ID" = "Sample Number")
    align_cols <- get_align_cols()
    clean_cols <- c("Avg [ng/ul]", "% CV")
    clean <- function(s) { str_remove_all(s, "[^[:digit:].]") %>% as.double() }
    get_pheno() %>%
      # rename(`Sample ID` = `Sample Number`) %>%
      rename(!!! get_align_cols()) %>%
      mutate_at(input$clean_cols, clean) %>%
      # mutate(`Avg [ng/ul]` = str_remove_all(`Avg [ng/ul]`, ">") %>% as.double(),
      #        `% CV` = str_remove_all(`% CV`, "!") %>% as.double()) %>%
      mutate(`Sample ID` = as.character(`Sample ID`))
  })
  
  filter_pheno <- reactive({
    filters <- c("") #SHINY use paste to create this from user input 
    clean_pheno() %>% filter(!!! parse_exprs(filters))
  })
  
  get_manifest <- reactive({
    #RSHINY Get this from UI
    controls <- c("Hypo-methylated Control", "Hyper-methylated control")
    by_cols <- c("Covid +/-")
    col_vals = c(Species = 'Homo sapiens', `Tissue Source` = 'Whole Blood') #SHINY get default values from user.
    create_manifest(filter_pheno(), controls, by_cols, col_vals)
  })
  
  createTableOutput <- function(df) {
    if (input$disp == "head") { return(head(df)) }
    return(df)
  }
  
  output$pheno <- renderTable({ req(input$pheno_file)
    createTableOutput(get_pheno())
  })
  
  output$cleaned_pheno <- renderTable({
    cleaned_pheno <- clean_pheno()
    if(input$disp == "head") { return(head(cleaned_pheno) ) }
    else { return(cleaned_pheno) }
  })
  
  output$filtered_pheno <- renderTable({
    filtered_pheno <- filter_pheno()
    if(input$disp == "head") { return(head(filtered_pheno) ) }
    else { return(filtered_pheno) }
  })
  
  output$manifest <- renderTable({ req(input$pheno_file)
    manifest <- get_manifest()
    if(input$disp == "head") { return(head(manifest) ) }
    else { return(manifest) }
  })
  
  output$downloadManifest <- downloadHandler(
    filename = function() {
      paste('data-', Sys.Date(), '.xlsx', sep='')
    },
    content = function(con) {
      write.xlsx(get_manifest(), file = con, showNA = FALSE)
    }
  )
  
  output$debug <- renderPrint({
    req(input$clean_cols)
    get_clean_cols()
  })
}
# Run the app ----
shinyApp(ui, server)