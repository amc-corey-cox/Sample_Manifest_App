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

#SHINY Get this from UI.
col_names <- c("Row", "Group ID", "Sample ID", "Plate Barcode or Number", "Plate", "Well", "Sample Well", "Species",
               "Gender (M/F/U)", "Volume (ul)", "Concentration (ng/ul)", "OD 260/280", "Tissue Source",
               "Extraction Method", "Ethnicity", "Parent 1 ID", "Parent 2 ID", "Replicate(s) ID", "Cancer Sample (Y/N)")

get_n_plates <- function (samples, controls, n_plate_wells, n_chip_wells) {
  n_samples <- nrow(samples)
  n_controls <- length(controls)
  
  n_plates <- ceiling( n_samples / ( n_plate_wells - n_controls ) )
  total_controls <- n_plates * n_controls
  total_samples <- n_samples + total_controls
  
  n_chips <- ceiling( total_samples / n_chip_wells )
  empty_wells <- n_chips * n_chip_wells - total_samples
  
  ( total_samples + empty_wells ) / n_plate_wells
}

add_column_na <- function(d, col_names) {
  add_cols <- col_names[!col_names %in% colnames(d)]
  
  if(length(add_cols) != 0) d[add_cols] <- NA
  d
}

get_plates <- function(n_plates, n_wells = 96) {
  rep(1:ceiling(n_plates), each = n_wells, length.out = n_wells * n_plates)
}

get_plate_wells <- function(n = 1, r = 12, c = 8) {
  if (n == 0) { return(NULL); }
  else if (n < 1) { return( get_plate_wells(n = 1, r = r * n, c) ) }
  else { return( c(
    paste0(rep(LETTERS[1:8], each = r), formatC(rep(1:r, times = c), width = 2, flag = "0")),
    get_plate_wells(n-1, r, c)
  ) ) }
}

get_plate_chips <- function(wells) {
  as.numeric(substring(wells, 2))
}

get_controls <- function(control_list, whole_plates, fill = 0) {
  n_controls <- whole_plates * length(control_list) + fill
  
  tibble("Sample ID" = rep(control_list, length.out = n_controls)) %>%
    mutate(r = sample(rep(sample(LETTERS[1:8]), length.out = n_controls)),
           c = sample(1:12, size = n_controls),
           n = rep(1:whole_plates, each = n_controls),
           Well = paste0(r, formatC(c, width = 2, flag = "0"))
    ) %>%
    select (-r, -c) %>% group_split(n, keep = FALSE) %>%
    imap(~ mutate(.x, Plate = .y, Chip = as.numeric(substring(Well, 2))))
}

get_open_wells <- function(all_controls, n_plates) {
  tibble(Plate = get_plates(n_plates), Well = get_plate_wells(n_plates), Chip = get_plate_chips(Well)) %>%
    group_split(Plate) %>%
    map2(all_controls, ~ filter(.x, ! Well %in% .y$Well)) %>%
    reduce(rbind) %>%
    arrange(Well)
}

plate_samples <- function(pheno, open_wells, by_cols, col_vals) {
  pheno %>% slice(sample(1:n(), n())) %>%
    mutate(!!! col_vals) %>%
    arrange(!!! quos(!!! syms(by_cols))) %>%
    cbind(open_wells) %>%
    group_split(Plate, Chip) %>%
    map(~ slice(., sample(1:n(), n()))) %>%
    reduce(rbind)
}

create_manifest <- function(pheno, controls, by_cols, col_vals) {
  set.seed(42)  # RSHINY Add control to allow changing seed
  
  wells_per_plate <- 96
  samples_per_chip <- 8
  
  n_plates <- get_n_plates(pheno, controls, 96, 8)
  whole_plates <- ceiling (n_plates)
  fill_wells <- (96 - length(controls)) * whole_plates - nrow(pheno)
  
  all_controls <- get_controls(controls, whole_plates, fill_wells)
  open_wells <- get_open_wells(all_controls, n_plates)
  plated_samples <- plate_samples(pheno, open_wells, by_cols, col_vals)
  formatted_plates <- plated_samples %>% add_column_na(col_names) %>% select(union(col_names, by_cols))
  return(formatted_plates)
}

upload_ctrl <- function() {
  fileInput("pheno_file", "Upload Phenotype File", multiple = FALSE,
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
    radioButtons("id_radio", "Select Sample ID Column", choices = c(1))
  )
}

ui <- fluidPage(
  titlePanel("CARGO Manifest Generator"),
  sidebarLayout(
    sidebarPanel(
      radioButtons( "disp", "Display", selected = "head", inline = TRUE,
                    choices = c(Head = "head", All = "all")),
      conditionalPanel( condition = "input.tabs == 'Phenotype'",
        upload_ctrl(), file_ctrl(), excel_ctrl(), txt_ctrl()
      ),
      conditionalPanel( condition = "input.tabs == 'Phenotype' || input.tabs == 'Cleanup'", id_ctrl() ),
      conditionalPanel( condition = "input.tabs == 'Cleanup'",
        checkboxGroupInput("clean_cols", "Columns to Cleanup", choices = NULL)
      ),
      conditionalPanel(
        condition = "input.tabs == 'Manifest'",
        checkboxGroupInput("by_cols", "Columns to Balance", choices = NULL)
      ),
      downloadButton("downloadManifest", "Download", class = 'rightAlign')
    ),
    mainPanel(
      tabsetPanel( id = "tabs", type = "pills",
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

  read_pheno <- reactive({
    if (is.na(excel_format(input$pheno_file$datapath))) {
      pheno <- read_delim(input$pheno_file$datapath, col_names = input$col_names,
                 delim = input$delim, quote = input$quote, skip = input$skip)
    } else {
      pheno <- read_excel( input$pheno_file$datapath, skip = input$skip, col_names = input$col_names)
    }
    updateRadioButtons(session, "id_radio", choices = colnames(pheno))
    map(list("align_cols", "clean_cols", "filter_cols", "by_cols"),
        ~ updateCheckboxGroupInput(session, .x, choices = set_names(colnames(pheno))))
    return(pheno)
  })
  
  get_pheno <- reactive({ req(input$pheno_file)
    pheno <- read_pheno()
    req(input$id_radio %in% colnames(pheno))
    pheno %>% rename("Sample ID" = input$id_radio) %>%
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

  output$debug <- renderPrint({ req(input$by_cols)
    input$by_cols
  })
}
# Run the app ----
shinyApp(ui, server)