library(pals)
library(ggthemes)
library(xlsx)
library(tidyverse)
library(tibble)
library(shiny)

source("manifest_functions.R")

my_colors <- c(stepped2(), stepped3(4), "#DDDDDD")

manifest_ui <- function(input, output, session) { 
  tabPanel( titlePanel("Plate Samples"), mainPanel( manifest_main(input, output, session) ) ) }
manifest_main <- function(input, output, session) { flowLayout( manifest_sidebar, manifest_tabs ) }

manifest_sidebar <- div(
    style = "min-width: 220px; max-width: 220px;",
    wellPanel( # sidebar
      actionButton("getPassedSamples", "Load Passed Samples"),
      # textOutput("debug"),
      uiOutput("controls")
) )

manifest_tabs <- div( style="display:inline-block; min-width: 400px; padding-left:10px; padding-top:10px",
  tabsetPanel( 
  id = "mtabs", type = "pills",
  tabPanel( "Passed QC", br(), DT::dataTableOutput("passedQC") ),
  tabPanel( "Manifest", br(), DT::dataTableOutput("manifestTable")),
  tabPanel( "Plate Layout", br(),
    div(DT::dataTableOutput("layoutKey"), style = "font-size:80%"),
    DT::dataTableOutput("plateLayout"),
    uiOutput("layoutFacets")),
  tabPanel( "Layout Facets", uiOutput("facet_UI"))
) )

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
          checkboxGroupInput("by_cols", "Balance by Columns", choices = set_names(field_names),
                         select = c("site", "Age_category", "Asthma")),
          checkboxGroupInput("add_cols", "Add to Manifest", choices = set_names(field_names)),
          # selectInput("col_vals", "Select Values", choices = c(Species = 'Homo sapiens', `Tissue Source` = 'Whole Blood'),
          #             selected = c("Species", "Tissue Source")),
          numericInput("seed", "Set Random Seed", value = 44),
          downloadButton("downloadManifest", "Download Manifest")
        ),
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
  
  observe({ req(input$by_cols)
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
    load("savePassed.RData")
    # passedData <- forCorey
    forCorey
  })
  
  get_plates <- reactive({ req(input$id_col)
    ### TODO: Get this from UI
    controls <- c("Hypo-methylated Control", "Hyper-methylated Control")
    
    set.seed(input$seed)
    if (input$bal_type == "Disperse") { manifest <- plate_disperse(input, get_data(), controls) }
    else { manifest <- plate_randomize(input, get_data(), controls) }
    manifest
  })
  
  get_manifest <- reactive({
    ### TODO: Get this from UI
    col_vals <- c(Species = 'Homo sapiens', `Tissue Source` = 'Whole Blood') #SHINY get default values from user.
    format_manifest(get_plates(), input$by_cols, input$add_cols, col_vals)
  })
  
  get_layout_types <- function(plates, by_cols) {
    by_cols %>% set_names %>%
      map( ~ pull(plates, .) %>% na.omit() %>% unique() ) %>%
      expand.grid() %>% unite(Sample_Type, all_of(by_cols)) %>%
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
  
  get_id_types <- function(plates, by_cols, plate_num) {
    plates %>%
      group_split(Plate) %>%
      map(~ unite(., Sample_Type, by_cols)) %>% map(~ pull(., Sample_Type)) %>%
      map(~ matrix(c(., rep("", 96 - length(.))), nrow = 8, dimnames = list(LETTERS[1:8], 1:12))) %>%
      pluck(plate_num) # Move to where we use the data?
  }
  
  get_layout <- function (plates, by_cols, plate_num, show_ids = TRUE) {
    by_cols <- by_cols
    types <- get_layout_types(plates, by_cols)
    rm_str <- "-methyl.*"
    
    ids <- get_layout_ids(plates, plate_num, rm_str)
    if(show_ids) { display_ids <- ids } else {
      display_ids <- ids %>% str_replace(".*", " ") %>% matrix(nrow = 8, dimnames = list(LETTERS[1:8], 1:12))
    }
    id_types <- get_id_types(plates, by_cols, plate_num)
    layout_colors <- get_layout_colors(length(types))
    
    cbind(display_ids, id_types) %>%
      datatable(class = 'cell-border stripe', 
        options = list(dom = 't', pageLength = 8, columnDefs = list(list(visible = FALSE, targets = 13:24)))) %>%
      formatStyle(columns = 1:12, valueColumns = 13:24,
                  backgroundColor = styleEqual(levels = types, values = layout_colors))
  }
  
  get_layout_key <- function(plates, by_cols) {
    all_types <- get_layout_types(plates, by_cols)
    
    types <-  plates %>% unite(Sample_Type, by_cols) %>%
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
  output$manifestTable <- DT::renderDataTable({ get_manifest() }, options = list(pageLength = 96))
  output$downloadManifest <- downloadHandler(
    filename = function() { paste('Manifest-', Sys.Date(), '.xlsx', sep='') },
    content = function(con) { write.xlsx(get_manifest(), file = con, showNA = FALSE) })
  
  output$plateLayout <- DT::renderDataTable(get_layout(get_plates(), input$by_cols, input$layout_plate, input$show_ids))
  output$layoutKey <- DT::renderDataTable(get_layout_key(get_plates(), input$by_cols))

  # output$debug <- renderText({ # req(input$fileUploaded)
  #   input$mtabs
  # })
}