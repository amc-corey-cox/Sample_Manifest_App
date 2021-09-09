manifest_panel <- div(
  style = "min-width: 240px; max-width: 240px;",
  wellPanel( # sidebar
    # conditionalPanel(condition = "input.mtabs == 'Passed QC'",
    #   radioButtons("dataSource", "Data Source", choices = c("QC", "Data"), selected = "QC"),
    #   actionButton("getPassedSamples", "Load Passed Samples")),
    # textOutput("debug"),
    conditionalPanel(condition = "input.mtabs == 'Manifest'",
      fileInput("template", "Upload Manifest Template", multiple = FALSE,
        accept = c(".csv", ".xls", ".xlsx", "text/csv", "text/comma-separated-values,text/plain", "application/vnd.ms-excel",
          "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")),
        # conditionalPanel(condition = "output.templateUploaded",
        #   numericInput("tmp_skip", "Skip lines:", 0, min = 0),
        #   p("Align Data and Manifest Columns"),
        #   splitLayout(
        #     selectInput("m_id_col_1", "Manifest Col", choices = NULL),
        #     selectInput("m_id_col_2", "Data Col", choices = NULL)
        #   )),
      conditionalPanel(condition = "output.templateUploaded",
        numericInput("tmp_skip", "Skip lines:", 0, min = 0)),
    ),
    uiOutput("manifest_controls")
))

manifest_tabs <- div(
  style="display:inline-block; min-width: 400px; padding-left:25px; padding-top:10px",
  tabsetPanel(
    id = "mtabs", type = "pills",
    tabPanel( "Imported Data", br(), DT::dataTableOutput("passedQC") ),
    tabPanel( "Manifest", br(), DT::dataTableOutput("manifestTable")),
    tabPanel( "Plate Layout", br(),
      div(DT::dataTableOutput("layoutKey"), style = "font-size:80%"),
      DT::dataTableOutput("plateLayout")),
    tabPanel( "Layout Facets", uiOutput("facet_UI"))
))

manifest_ui <- tabPanel(
  titlePanel("Plate Samples"),
  splitLayout(
    cellWidths = c("40%", "20%", "20%", "20%"),
    titlePanel("Plating & Manifest"),
    radioButtons("m_disp", "Display", selected = "head", inline = TRUE, choices = c(Head = "head", All = "all")),
    downloadButton("downloadManifest", "Get Manifest"),
    downloadButton("manifestReport", "Get Manifest Report"),
    tags$style(type = 'text/css', "#m_disp { margin-top: 15px; }"),
    tags$style(type = 'text/css', "#downloadManifest { margin-top: 20px; }"),
    tags$style(type = 'text/css', "#manifestReport { margin-top: 20px; }")
  ),
  flowLayout(manifest_panel, manifest_tabs))