source("manifest_functions.R", local = TRUE)

my_colors <- c(stepped2(), stepped3(4), "#DDDDDD")
unique_n <- function(x) { unique(x) %>% length() }

manifest_server <- function(input, output, session) {
  output$templateUploaded <- reactive({ return(!is.null(input$template)) })
  outputOptions(output, 'templateUploaded', suspendWhenHidden = FALSE)
  
  output$templateExcel <- reactive({
    if (!is.null(input$template) &&
        !is.na(excel_format(input$template$datapath))) {
      updateSelectInput(session, "tmp_sheet", choices = excel_sheets(input$template$datapath))
      return(TRUE)
    }
    updateSelectInput(session, "tmp_sheet", choices = character(0))
    return(FALSE)
  })
  outputOptions(output, 'fileExcel', suspendWhenHidden = FALSE)
  
  observe({ req(input$files)
    field_names <- colnames(get_data())
    default_cols <- get_default_cols()
    updateActionButton(session, "getPassedSamples", label = "Reload Passed Samples")
    output$manifest_controls <- renderUI({
      tagList(
        radioButtons("control_type", "Control Type", choices = c("None", "Epic", "MEGA")),
        radioButtons("empty_wells", "Empty Wells", choices = c("Use Controls", "Leave Empty")),
        conditionalPanel( "input.mtabs == 'Plate Layout' || input.mtabs == 'Layout Facets'",
          checkboxInput("show_ids", "Show IDs", value = TRUE),
          numericInput("layout_plate", "Select Plate", value = 1, min = 1, max = 5, width = "100px")),
        conditionalPanel("input.mtabs != 'Layout Facets'",
          radioButtons("bal_type", "Balance Type", choices = c("Grouped Disperse", "Simple Disperse", "Randomize")),
          selectInput("id_col", "Sample ID Column", choices = field_names),
          checkboxGroupInput("m_by_cols", "Balance by Columns", choices = set_names(field_names),
                             select = default_cols),
          checkboxGroupInput("add_cols", "Add to Manifest", choices = set_names(field_names)),
          # selectInput("col_vals", "Select Values", choices = c(Species = 'Homo sapiens', `Tissue Source` = 'Whole Blood'),
          #             selected = c("Species", "Tissue Source")),
          numericInput("seed", "Set Random Seed", value = 44),
          numericInput("num_plates", "Number of Plates", value = NA_integer_),
        ),
        conditionalPanel("input.mtabs == 'Layout Facets'",
          checkboxGroupInput("layout_cols", "Facet Columns", choices = set_names(field_names),
            select = default_cols))
      )
    })
  })
  
  output$facet_UI <- renderUI ({
    input$layout_cols %>% str_c("layout_", .) %>%
      map(~ div(DT::dataTableOutput(str_c(., "_key")), DT::dataTableOutput(.)))
  })
  
  observe({ req(input$layout_cols)
    input$layout_cols %>% set_names(str_c("layout_", ., "_key")) %>%
      iwalk( function(by_col, out_name) { 
        output[[out_name]] <- DT::renderDataTable( get_layout_key(get_plates(), by_col))
        })
    input$layout_cols %>% set_names(str_c("layout_", .)) %>%
      iwalk(function(by_col, out_name) {
        output[[out_name]] <- DT::renderDataTable( get_layout(get_plates(), by_col, input$layout_plate, input$show_ids) )
        })
  })
  
  # Use for returning from qc or data portion
  # get_data <- eventReactive(input$getPassedSamples, {
  #   if (input$dataSource == "Data") { return (clean_pheno()) }
  #   load("savePassed.RData")
  #   forCorey
  # })
  
  # Get data only from data portion of app
  get_data <- reactive({ filter_pheno() })
  
  get_controls <- reactive({
    if(input$control_type == "MEGA") {controls <- c("HapMap Control", "HapMap Control", "HapMap Control", "Duplicate", "Duplicate") }
    else if(input$control_type == "Epic") { controls <- c("Hypo-Methylated Control", "Hyper-Methylated Control") }
    else { controls <- character() }
    controls
  })
  
  get_plates <- reactive({ req(input$id_col)
    controls <- get_controls()
    
    if (input$bal_type == "Grouped Disperse") {
      plates <- grouped_disperse(get_data(), controls, input$seed, input$id_col, input$m_by_cols, input$empty_wells)
    } else if (input$bal_type == "Simple Disperse") {
       plates <- simple_disperse(get_data(), controls, input$seed, input$id_col, input$m_by_cols, input$empty_wells)
    } else {
      plates <- plate_randomize(get_data(), controls, input$seed, input$id_col, input$m_by_cols, input$empty_wells)
    }
    if (is.na(input$num_plates)) {
      return (plates)
    } else {
      plates %>% filter(Plate <= input$num_plates) %>% return
    }
  })
  
  get_default_cols <- reactive({
    get_data() %>%
      select(where(~ 1 < unique_n(.x) && unique_n(.x) < 5)) %>%
      select(matches(c("Site", "Age", "Asthma", "Gender", "Status", "Race", "Ethnicity", "Symptom")), everything()) %>%
      select(1:4) %>% colnames()
  })
  
  get_manifest_m <- reactive({
    ### TODO: Get this from UI
    col_vals <- c(Species = 'Homo sapiens', `Tissue Source` = 'Whole Blood') #SHINY get default values from user.
    # col_vals for WGS, this is a quick hack for now... try to do better in the future!
    # col_vals <- c('Specimen type' = 'Saliva', Instrument = 'NovaSeq') #SHINY get default values from user.
    if (! is.null(input$template)) {
      col_names <- import_file(input$template$datapath, col_names = TRUE, input$tmp_delim, input$tmp_quote, input$tmp_skip) %>%
        colnames
    } else {
      col_names <- c("Row", "Group ID", "Sample ID", "Plate Barcode or Number", "Plate", "Well", "Sample Well", "Species",
                     "Gender (M/F/U)", "Volume (ul)", "Concentration (ng/ul)", "OD 260/280", "Tissue Source",
                     "Extraction Method", "Ethnicity", "Parent 1 ID", "Parent 2 ID", "Replicate(s) ID", "Cancer Sample (Y/N)")
    }
    updateSelectInput(session, "m_id_col_1", choices = col_names)
    format_manifest(get_plates(), input$m_by_cols, input$add_cols, col_vals, col_names) 
  })
  
  get_layout_types <- function(plates, m_by_cols) {
    m_by_cols %>% set_names %>%
      map( ~ pull(plates, .) %>% na.omit() %>% unique() ) %>%
      # map( ~ pull(plates, .) %>% replace_na(as.list(rep("Missing", length(m_by_cols))) %>% set_names(m_by_cols)) %>% unique() ) %>%
      expand.grid() %>% unite(Sample_Type, all_of(m_by_cols)) %>%
      arrange(Sample_Type) %>% pull(Sample_Type) %>% c("NA_NA_NA")
  }
  
  get_layout_colors <- function(num_types) { n = num_types - 1
    if (n < 9) { return (c(get_layout_colors(30)[c(1, 5, 9, 13, 17, 21, 25, 29)][1:n], "#DDDDDD")) }
    c(c(stepped2(), stepped3(), stepped(20))[1:(n)], "#DDDDDD")
  }
  
  get_layout_ids <- function(plates, plate_num, rm_string) {
    plates %>%
      mutate(`Sample ID` = str_replace(`Sample ID`, rm_string, "")) %>%
      group_split(Plate) %>%
      map(~ pull(., `Sample ID`) ) %>%
      map(~ matrix(c(., rep("", 96 - length(.))), nrow = 8, dimnames = list(LETTERS[1:8], 1:12))) %>%
      pluck(plate_num) # Move to where we use the data?
  }
  
  get_id_types <- function(plates, m_by_cols, plate_num) {
    plates %>%
      group_split(Plate) %>%
      map(~ unite(., Sample_Type, all_of(m_by_cols))) %>% map(~ pull(., Sample_Type)) %>%
      map(~ matrix(c(., rep("", 96 - length(.))), nrow = 8, dimnames = list(LETTERS[1:8], 1:12))) %>%
      pluck(plate_num) # Move to where we use the data?
  }
  
  get_layout <- function (plates, m_by_cols, plate_num, show_ids = TRUE) {
    types <- get_layout_types(plates, m_by_cols)
    rm_str <- "-methyl.*"
    
    ids <- get_layout_ids(plates, plate_num, rm_str)
    if(show_ids) { display_ids <- ids } else {
      display_ids <- ids %>% str_replace(".*", " ") %>% matrix(nrow = 8, dimnames = list(LETTERS[1:8], 1:12))
    }
    id_types <- get_id_types(plates, m_by_cols, plate_num)
    layout_colors <- get_layout_colors(length(types))
    
    cbind(display_ids, id_types) %>%
      datatable(class = 'cell-border stripe', height = "100%", width = "100%",
        options = list(dom = 't', pageLength = 8, columnDefs = list(list(visible = FALSE, targets = 13:24)))) %>%
      formatStyle(columns = 1:12, valueColumns = 13:24,
                  backgroundColor = styleEqual(levels = types, values = layout_colors))
  }
  
  get_layout_key <- function(plates, m_by_cols) {
    all_types <- get_layout_types(plates, m_by_cols)
    
    types <-  plates %>% unite(Sample_Type, m_by_cols) %>%
      arrange(Sample_Type) %>% pull(Sample_Type) %>% unique() %>% intersect(all_types, .)
    
    type_names <- types # %>% str_replace("NA.*", "Control") %>% str_replace_all("_", " ")
    col_nums <- 1:length(types)
    
    dt_opts <- list(dom = 't', columnDefs = list(list(visible = FALSE, targets = col_nums - 1)))
    bg_color <- styleEqual(levels = all_types, values = get_layout_colors(length(all_types)))
    
    matrix(c(types, type_names), nrow = 1) %>%
      datatable(colnames = rep("", ncol(.)), options = dt_opts, height = "100%", width = "100%") %>%
      formatStyle(columns = col_nums + length(types), valueColumns = col_nums,
                  backgroundColor = bg_color, color = "#FFFFFF")
  }
  
  createTableOutput_m <- function(df) {
    if (input$m_disp == "head") { return(head(df)) }
    return(df)
  }
  
  output$passedQC <- DT::renderDataTable({ get_data() %>% createTableOutput_m() }, options = list(pageLength = 96))
  output$manifestTable <- DT::renderDataTable({ get_manifest_m() %>% createTableOutput_m() }, options = list(pageLength = 96))
  output$downloadManifest <- downloadHandler(
    filename = function() { str_c('Manifest-', Sys.Date(), '.xlsx') },
    content = function(con) { write.xlsx(get_manifest_m(), file = con, showNA = FALSE) })
  
  output$manifestReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = function() { str_c('Manifest_Report-', Sys.Date(), '.html') },
    content = function(con) {
      params <- lst(input = input,
                    info = get_info(get_data(), get_controls(), 96, 8),
                    num_plates = ifelse(is.na(input$num_plates), info$total_plates, input$num_plates),
                    plate_layouts = 1:num_plates %>% set_names %>%
                      map(~ get_layout(get_plates(), input$m_by_cols, ., input$show_ids)),
                    layout_key = get_layout_key(get_plates(), input$m_by_cols),
                    types = get_layout_types(get_plates(), input$m_by_cols),
                    layout_colors = layout_colors <- get_layout_colors(length(types)),
                    facet_plates = input$layout_cols %>% set_names %>%
                      map(~ map2(., 1:num_plates, ~ get_layout(get_plates(), .x, .y, input$show_ids) )),
                    facet_keys = input$layout_cols %>% set_names(str_c(., "_key")) %>%
                      map( ~ get_layout_key(get_plates(), .)))
      
      render("Manifest_Report.Rmd", output_file = con,
             params = params,
             # Can I just use globalenv() and not pass params? 
             envir = new.env(parent = globalenv()))
    })
  
  output$plateLayout <- DT::renderDataTable(get_layout(get_plates(), input$m_by_cols, input$layout_plate, input$show_ids))
  output$layoutKey <- DT::renderDataTable(get_layout_key(get_plates(), input$m_by_cols))
}
