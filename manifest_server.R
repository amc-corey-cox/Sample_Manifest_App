source("manifest_functions.R", local = TRUE)

my_colors <- c(stepped2(), stepped3(4), "#DDDDDD")

manifest_server <- function(input, output, session) {
  observeEvent(input$getPassedSamples, {
    field_names <- colnames(get_data())
    updateActionButton(session, "getPassedSamples", label = "Reload Passed Samples")
    output$controls <- renderUI({
      tagList(
        conditionalPanel( "input.mtabs == 'Plate Layout' || input.mtabs == 'Layout Facets'",
          checkboxInput("show_ids", "Show IDs", value = TRUE),
          numericInput("layout_plate", "Select Plate", value = 1, min = 1, max = 5, width = "100px")),
        conditionalPanel("input.mtabs != 'Layout Facets'",
          radioButtons("bal_type", "Balance Type", choices = c("Disperse", "Randomize")),
          selectInput("id_col", "Sample ID Column", choices = field_names),
          checkboxGroupInput("m_by_cols", "Balance by Columns", choices = set_names(field_names),
                         select = c("site", "Age_category", "Asthma")),
          # checkboxGroupInput("m_by_cols", "Balance by Columns", choices = set_names(field_names)),
          checkboxGroupInput("add_cols", "Add to Manifest", choices = set_names(field_names)),
          # selectInput("col_vals", "Select Values", choices = c(Species = 'Homo sapiens', `Tissue Source` = 'Whole Blood'),
          #             selected = c("Species", "Tissue Source")),
          numericInput("seed", "Set Random Seed", value = 44),
          downloadButton("downloadManifest", "Download Manifest")
        ),
        # conditionalPanel("input.mtabs == 'Layout Facets'",
        #     checkboxGroupInput("layout_cols", "Balance by Columns", choices = set_names(field_names)))
        conditionalPanel("input.mtabs == 'Layout Facets'",
          checkboxGroupInput("layout_cols", "Balance by Columns", choices = set_names(field_names),
            select = c("site", "Age_category", "Asthma", "Gender"))
        )
    )})
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
  
  get_data <- eventReactive(input$getPassedSamples, {
    if (input$dataSource == "Data") { return (get_pheno()) }
    load("savePassed.RData")
    forCorey
  })
  
  get_plates <- reactive({ req(input$id_col)
    ### TODO: Get this from UI
    # controls <- c("Hypo-Methylated Control", "Hyper-Methylated Control")
    controls <- c("HapMap Control", "HapMap Control", "HapMap Control", "Duplicate", "Duplicate")
    
    set.seed(input$seed)
    if (input$bal_type == "Disperse") { manifest <- plate_disperse(input, get_data(), controls) }
    else { manifest <- plate_randomize(input, get_data(), controls) }
    manifest
  })
  
  get_manifest_m <- reactive({
    ### TODO: Get this from UI
    col_vals <- c(Species = 'Homo sapiens', `Tissue Source` = 'Whole Blood') #SHINY get default values from user.
    format_manifest(get_plates(), input$m_by_cols, input$add_cols, col_vals)
  })
  
  get_layout_types <- function(plates, m_by_cols) {
    m_by_cols %>% set_names %>%
      map( ~ pull(plates, .) %>% na.omit() %>% unique() ) %>%
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
    m_by_cols <- m_by_cols
    types <- get_layout_types(plates, m_by_cols)
    rm_str <- "-methyl.*"
    
    ids <- get_layout_ids(plates, plate_num, rm_str)
    if(show_ids) { display_ids <- ids } else {
      display_ids <- ids %>% str_replace(".*", " ") %>% matrix(nrow = 8, dimnames = list(LETTERS[1:8], 1:12))
    }
    id_types <- get_id_types(plates, m_by_cols, plate_num)
    layout_colors <- get_layout_colors(length(types))
    
    cbind(display_ids, id_types) %>%
      datatable(class = 'cell-border stripe', 
        options = list(dom = 't', pageLength = 8, columnDefs = list(list(visible = FALSE, targets = 13:24)))) %>%
      formatStyle(columns = 1:12, valueColumns = 13:24,
                  backgroundColor = styleEqual(levels = types, values = layout_colors))
  }
  
  get_layout_key <- function(plates, m_by_cols) {
    all_types <- get_layout_types(plates, m_by_cols)
    
    types <-  plates %>% unite(Sample_Type, m_by_cols) %>%
      arrange(Sample_Type) %>% pull(Sample_Type) %>% unique() %>% intersect(all_types, .)
    
    type_names <- types %>% str_replace("NA.*", "Control") %>% str_replace_all("_", " ")
    col_nums <- 1:length(types)
    
    dt_opts <- list(dom = 't', columnDefs = list(list(visible = FALSE, targets = col_nums - 1)))
    bg_color <- styleEqual(levels = all_types, values = get_layout_colors(length(all_types)))
    
    matrix(c(types, type_names), nrow = 1) %>%
      datatable(colnames = rep("", ncol(.)), options = dt_opts) %>%
      formatStyle(columns = col_nums + length(types), valueColumns = col_nums,
                  backgroundColor = bg_color, color = "#FFFFFF")
  }
  
  output$passedQC <- DT::renderDataTable({ get_data() }, options = list(pageLength = 96))
  output$manifestTable <- DT::renderDataTable({ get_manifest_m() }, options = list(pageLength = 96))
  output$downloadManifest <- downloadHandler(
    filename = function() { paste('Manifest-', Sys.Date(), '.xlsx', sep='') },
    content = function(con) { write.xlsx(get_manifest_m(), file = con, showNA = FALSE) })
  
  output$plateLayout <- DT::renderDataTable(get_layout(get_plates(), input$m_by_cols, input$layout_plate, input$show_ids))
  output$layoutKey <- DT::renderDataTable(get_layout_key(get_plates(), input$m_by_cols))

  # output$debug <- renderText({ # req(input$fileUploaded)
  #   input$mtabs
  # })
}