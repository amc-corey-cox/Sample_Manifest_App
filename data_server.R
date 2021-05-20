data_server <- function(input, output, session) {
  output$fileUploaded <- reactive({ return(!is.null(input$files)) })
  outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)
  
  output$fileExcel <- reactive({
    if (!is.null(input$files) &&
        !is.na(excel_format(input$files$datapath))) {
      updateSelectInput(session, "sheet", choices = excel_sheets(input$files$datapath))
      # updateSelectInput(session, "clean_col_1", choices = c("None", colnames(get_pheno())))
      return(TRUE)
    }
    updateSelectInput(session, "sheet", choices = character(0))
    return(FALSE)
  })
  outputOptions(output, 'fileExcel', suspendWhenHidden = FALSE)
  vals <- reactiveValues(data = NULL)
  
  observe({ req(input$files)
    output$ui_clean_cols <- renderUI({ req(input$files)
      tagList(
        varSelectInput("clean_cols", "Column(s) to Number", multiple = TRUE, data = vals$data),
        p("Select Columns to remove all non-digits and convert to number."),
        varSelectInput("na_cols", "NA to Missing", multiple = TRUE, data = vals$data),
        p('Select Columns to convert "NA" to "Missing"'))
    })
  })
  
  # Move this to a utilities file?
  import_file <<- function(file, col_names, delim, quote, skip) {
    if (is.na(excel_format(file))) {
      file <- read_delim(file, col_names = col_names,
                          delim = delim, quote = quote, skip = skip)
    } else {
      file <- read_excel(file, skip = skip, col_names = col_names)
    }
    file
  }
  
  # Have to put these in global environment for now. Rewrite using moduleServer and nested servers.
  get_pheno <<- reactive({ req(input$files)
    # pheno <- read_pheno()
    pheno <- import_file(input$files$datapath, input$col_names, input$delim, input$quote, input$skip)
    updateSelectInput(session, "d_id_col", choices = colnames(pheno))
    updateSelectInput(session, "m_id_col_2", choices = colnames(pheno))
    map(list("align_cols", "filter_cols", "by_cols"),
        ~ updateCheckboxGroupInput(session, .x, choices = set_names(colnames(pheno))))
    vals$data <- pheno
    req(input$d_id_col %in% colnames(pheno))
    pheno %>% rename("Sample ID" = input$d_id_col) %>%
      mutate(`Sample ID` = as.character(`Sample ID`))
  })
  
  clean_pheno <- reactive({
    make_numeric <- function(x) { str_remove_all(x, "[^[:digit:].]") %>% as.double() }
    na_to_missing <- function(x) { ifelse(is.na(x), "Missing", x) }
    
    get_pheno() %>% mutate(
      across(as.character(input$clean_cols), make_numeric),
      across(as.character(input$na_cols), na_to_missing))
  })
  
  filter_pheno <- reactive({
    filters <- c("") #SHINY use paste to create this from user input
    clean_pheno() # %>% filter(!!! parse_exprs(filters))
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