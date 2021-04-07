manifest_sidebar <- div(
  style = "min-width: 240px; max-width: 240px;",
  wellPanel( # sidebar
    radioButtons("dataSource", "Data Source", choices = c("QC", "Data"), selected = "QC"),
    actionButton("getPassedSamples", "Load Passed Samples"),
    # textOutput("debug"),
    uiOutput("manifest_controls")
) )

manifest_tabs <- div( style="display:inline-block; min-width: 400px; padding-left:25px; padding-top:10px",
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

manifest_main <- flowLayout( manifest_sidebar, manifest_tabs )

manifest_ui <- tabPanel( titlePanel("Plate Samples"), mainPanel( manifest_main ) )