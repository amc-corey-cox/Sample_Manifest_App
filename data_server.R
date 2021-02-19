library(tidyverse)
library(shiny)

data_server <- function(input, output, session) {
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
    updateSelectInput(session, "d_id_col", choices = colnames(pheno))
    map(list("align_cols", "clean_cols", "filter_cols", "by_cols"),
        ~ updateCheckboxGroupInput(session, .x, choices = set_names(colnames(pheno))))
    return(pheno)
  })
  
  get_pheno <- reactive({ req(input$files)
    pheno <- read_pheno()
    req(input$d_id_col %in% colnames(pheno))
    pheno %>% rename("Sample ID" = input$d_id_col) %>%
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
  
  get_manifest_d <- reactive({ req(input$by_cols)
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
  output$manifest <- renderTable({ createTableOutput(get_manifest_d()) })
  
  output$d_downloadManifest <- downloadHandler(
    filename = function() {
      paste('data-', Sys.Date(), '.xlsx', sep='')
    },
    content = function(con) {
      write.xlsx(get_manifest_d(), file = con, showNA = FALSE)
    }
  )
  
  output$debug <- renderText({ # req(input$fileUploaded)
    input$d_id_col
  })
}