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
  import_file <- function(file, col_names, delim, quote, skip) {
    if (is.na(excel_format(file))) {
      file <- read_delim(file, col_names = col_names,
                          delim = delim, quote = quote, skip = skip)
    } else {
      file <- read_excel(file, skip = skip, col_names = col_names)
    }
    file
  }
  
  read_pheno <- reactive({ req(input$files)
    pheno <- import_file(input$files$datapath, input$col_names, input$delim, input$quote, input$skip)
    updateSelectInput(session, "d_id_col", choices = colnames(pheno))
    updateSelectInput(session, "m_id_col_2", choices = colnames(pheno))
    map(list("align_cols", "filter_cols", "by_cols"),
        ~ updateCheckboxGroupInput(session, .x, choices = set_names(colnames(pheno))))
    vals$data <- pheno
  })
  
  get_pheno <- reactive({ req(input$files)
    # pheno <- read_pheno()
    pheno <- read_pheno()
    req(input$d_id_col %in% colnames(pheno))
    pheno %>%
      possibly(rename, otherwise = ., )("Sample ID_old" = "Sample ID") %>%
      rename("Sample ID" = input$d_id_col) %>%
      mutate(`Sample ID` = as.character(`Sample ID`))
  })
  
  clean_pheno <- reactive({
    make_numeric <- function(x) { str_remove_all(x, "[^[:digit:].]") %>% as.double() }
    na_to_missing <- function(x) { ifelse(is.na(x), "Missing", x) }
    
    get_pheno() %>%
      when(length(input$clean_cols) < 1 ~ .,
           mutate(., across(as.character(input$clean_cols), make_numeric))) %>%
      when(length(input$na_cols) < 1 ~ .,
           mutate(., across(as.character(input$na_cols), na_to_missing)))
  })

  observeEvent(input$getFilterColumns, {
    # filter_choices <- colnames(clean_pheno())
    showModal(modalDialog(
      checkboxGroupInput("filter_cols", "Select Columns for Filtering", choices = colnames(clean_pheno())),
      actionButton("filterSelectedColumns", "Filter Selected"),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  filter_names <- reactive({ req(input$filter_cols)
    input$filter_cols %>%
      map(~ str_c("filt_", make.names(.))) %>%
      set_names(input$filter_cols, .)
    })
  
  has_filter <- reactive({ req(input$filter_cols)
    filter_names() %>% names() %>%
      map( ~ pluck(input, .) %>% length() > 0) %>%
      reduce(all)
  })
  
  observeEvent(input$filterSelectedColumns , {
    output$ui_filter_cols <- renderUI({
      tagList(
        filter_names() %>%
          imap(~ checkboxGroupInput( inputId = .y, label = .x,
            choices = pull(clean_pheno(), .x) %>% unique(),
            selected = pull(clean_pheno(), .x) %>% unique())),
      )
    })
    removeModal()
  })
  
  filter_pheno <- reactive({ req(filter_names)
    if (! has_filter()) { clean_pheno() }
    else {
    data_filters <- filter_names() %>%
      imap(~ str_c(.x, " %in% c(", str_c('"', pluck(input, .y), '"', collapse = ', '), ")")) %>%
      unlist()
    clean_pheno() %>%
      filter_(data_filters)
    }
  })
  
  # Have to put this in global environment for now. Rewrite using moduleServer and nested servers.
  pheno <<- reactive({
    if(isTruthy(input$filter_cols)) { filter_pheno() } else { clean_pheno() }
  })
  
  createTableOutput <- function(df) {
    if (input$d_disp == "head") { return(head(df)) }
    return(df)
  }
  
  output$pheno <- renderTable({ createTableOutput(get_pheno()) })
  output$cleaned_pheno <- renderTable({ createTableOutput(clean_pheno()) })
  output$filtered_pheno <- renderTable({ createTableOutput(filter_pheno()) })
  
  output$downloadCleanData <- downloadHandler(
    filename = function() { str_c('cleaned_and_filtered_data-', Sys.Date(), '.tsv') },
    content = function(con) { write_tsv(filter_pheno(), file = con) })
  
  output$debug <- renderText({ # req(input$fileUploaded)
    input$d_id_col
  })
}